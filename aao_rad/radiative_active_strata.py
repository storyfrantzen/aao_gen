#!/usr/bin/env python3
"""Build an evidence-audited active-stratum mask for radiative mode 4.

The classifier deliberately separates physical support, analysis relevance,
and generator readiness.  A zero Monte Carlo count is never sufficient to
declare a stratum structurally empty or negligible.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable

import radiative_guards
import radiative_mode4


EVIDENCE_SCHEMA = "aao-rad-stratum-relevance-evidence-v1"
MASK_SCHEMA = "aao-rad-active-stratum-mask-v1"
PHYSICAL_STATUSES = {
    "unknown",
    "nonempty",
    "partially_accessible",
    "structurally_empty",
}
ACTIVE_STATUSES = {
    "active_ready",
    "active_needs_pilot_validation",
    "active_needs_optimization",
}
ALL_STATUSES = ACTIVE_STATUSES | {
    "closure_only",
    "structurally_empty",
    "needs_relevance_assessment",
}


class ActiveStratumError(RuntimeError):
    """Raised when classification evidence is incomplete or inconsistent."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path: Path, payload: object) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _bounds(stratum: radiative_guards.Stratum) -> dict[str, list[float]]:
    return {
        "Q2": list(stratum.q2),
        "xB": list(stratum.xb),
        "minus_t": list(stratum.minus_t),
        "phi_deg": list(stratum.phi_deg),
    }


def _indices(stratum: radiative_guards.Stratum) -> dict[str, int]:
    return {
        "Q2": stratum.iq2,
        "xB": stratum.ixb,
        "minus_t": stratum.it,
        "phi_deg": stratum.iphi,
    }


def _normalized_indices(value: object) -> tuple[int, int, int, int]:
    if not isinstance(value, dict):
        raise ValueError("indices must be an object")
    formats = (
        ("Q2", "xB", "minus_t", "phi_deg"),
        ("iq2", "ixb", "it", "iphi"),
    )
    for names in formats:
        if set(value) != set(names):
            continue
        entries = tuple(value[name] for name in names)
        if any(
            isinstance(entry, bool) or not isinstance(entry, int)
            for entry in entries
        ):
            raise ValueError("indices must be integers")
        return entries
    raise ValueError(
        "indices must use Q2/xB/minus_t/phi_deg or iq2/ixb/it/iphi"
    )


def _selected_strata(
    config: dict, flat_indices: Iterable[int] | None
) -> list[radiative_guards.Stratum]:
    catalog = radiative_guards.enumerate_strata(config)
    requested = list(flat_indices or [])
    if not requested:
        return catalog
    if len(set(requested)) != len(requested):
        raise ValueError("duplicate --flat-index selection")
    lookup = {item.flat_index: item for item in catalog}
    unknown = sorted(set(requested) - set(lookup))
    if unknown:
        raise ValueError(
            "flat indices lie outside the analysis catalog: "
            + ", ".join(str(value) for value in unknown)
        )
    return [lookup[value] for value in requested]


def _parse_source(specification: str) -> tuple[str, Path]:
    if "=" not in specification:
        raise ValueError(
            "--evidence-source must have the form IDENTIFIER=PATH"
        )
    identifier, raw_path = specification.split("=", 1)
    identifier = identifier.strip()
    if not identifier or not raw_path.strip():
        raise ValueError(
            "--evidence-source must have the form IDENTIFIER=PATH"
        )
    return identifier, Path(raw_path).expanduser().resolve()


def create_template(args: argparse.Namespace) -> Path:
    config_path = args.config.expanduser().resolve()
    config, config_sha256 = radiative_guards.load_analysis_config(config_path)
    output = args.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"{output} already exists")
    output.parent.mkdir(parents=True, exist_ok=True)

    sources: dict[str, dict[str, object]] = {}
    for specification in args.evidence_sources or []:
        identifier, path = _parse_source(specification)
        if identifier in sources:
            raise ValueError(f"duplicate evidence source {identifier}")
        if not path.is_file():
            raise FileNotFoundError(path)
        sources[identifier] = {
            "identifier": identifier,
            "path": str(path),
            "sha256": _sha256(path),
        }

    strata = _selected_strata(config, args.flat_indices)
    payload = {
        "schema": EVIDENCE_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "analysis_config": str(config_path),
        "analysis_config_sha256": config_sha256,
        "instructions": {
            "physical_status_values": sorted(PHYSICAL_STATUSES),
            "fractions_are_relative_to_the_full_analysis_prediction": True,
            "zero_mc_targets_do_not_prove_structural_emptiness": True,
            "nondefault_claims_require_rationale_and_source_ids": True,
        },
        "sources": list(sources.values()),
        "strata": [
            {
                "stratum_id": stratum.identifier,
                "flat_index": stratum.flat_index,
                "indices": _indices(stratum),
                "bounds": _bounds(stratum),
                "physical_status": "unknown",
                "analysis_included": None,
                "data_events": None,
                "model_cross_section_fraction": None,
                "maximum_feed_in_fraction": None,
                "global_closure_impact_fraction": None,
                "force_active": False,
                "rationale": None,
                "source_ids": [],
            }
            for stratum in strata
        ],
    }
    _write_json(output, payload)
    return output


def _load_relevance(
    path: Path, config_sha256: str
) -> tuple[dict[str, dict[str, object]], list[dict[str, object]]]:
    resolved = path.expanduser().resolve()
    payload = json.loads(resolved.read_text(encoding="utf-8"))
    if payload.get("schema") != EVIDENCE_SCHEMA:
        raise ActiveStratumError(
            f"{resolved}: expected schema {EVIDENCE_SCHEMA}"
        )
    if payload.get("analysis_config_sha256") != config_sha256:
        raise ActiveStratumError(
            f"{resolved}: analysis configuration hash differs"
        )

    sources: dict[str, dict[str, object]] = {}
    for source in payload.get("sources", []):
        identifier = str(source.get("identifier", ""))
        if not identifier or identifier in sources:
            raise ActiveStratumError(
                f"{resolved}: malformed or duplicate evidence source"
            )
        source_path = Path(str(source.get("path", ""))).expanduser().resolve()
        expected_hash = str(source.get("sha256", ""))
        if not source_path.is_file():
            raise FileNotFoundError(source_path)
        if _sha256(source_path) != expected_hash:
            raise ActiveStratumError(
                f"{resolved}: evidence source {identifier} changed"
            )
        sources[identifier] = {
            "identifier": identifier,
            "path": str(source_path),
            "sha256": expected_hash,
        }

    records: dict[str, dict[str, object]] = {}
    fraction_fields = (
        "model_cross_section_fraction",
        "maximum_feed_in_fraction",
        "global_closure_impact_fraction",
    )
    numeric_fields = ("data_events",) + fraction_fields
    for raw in payload.get("strata", []):
        record = dict(raw)
        stratum_id = str(record.get("stratum_id", ""))
        if not stratum_id or stratum_id in records:
            raise ActiveStratumError(
                f"{resolved}: malformed or duplicate stratum evidence"
            )
        physical_status = str(record.get("physical_status", "unknown"))
        if physical_status not in PHYSICAL_STATUSES:
            raise ActiveStratumError(
                f"{resolved}: {stratum_id} has invalid physical_status"
            )
        analysis_included = record.get("analysis_included")
        if analysis_included not in (None, True, False):
            raise ActiveStratumError(
                f"{resolved}: {stratum_id} analysis_included must be "
                "true, false, or null"
            )
        for field in numeric_fields:
            value = record.get(field)
            if value is None:
                continue
            if isinstance(value, bool):
                raise ActiveStratumError(
                    f"{resolved}: {stratum_id} {field} is not numeric"
                )
            value = float(value)
            if not math.isfinite(value) or value < 0.0:
                raise ActiveStratumError(
                    f"{resolved}: {stratum_id} {field} must be finite and "
                    "nonnegative"
                )
            record[field] = value
        for field in fraction_fields:
            value = record.get(field)
            if value is not None and float(value) > 1.0:
                raise ActiveStratumError(
                    f"{resolved}: {stratum_id} {field} must not exceed one"
                )
        if record.get("data_events") is not None:
            value = float(record["data_events"])
            if not value.is_integer():
                raise ActiveStratumError(
                    f"{resolved}: {stratum_id} data_events must be integral"
                )
            record["data_events"] = int(value)

        source_ids = [str(value) for value in record.get("source_ids", [])]
        if len(set(source_ids)) != len(source_ids):
            raise ActiveStratumError(
                f"{resolved}: {stratum_id} repeats a source identifier"
            )
        unknown_sources = sorted(set(source_ids) - set(sources))
        if unknown_sources:
            raise ActiveStratumError(
                f"{resolved}: {stratum_id} cites unknown sources "
                + ", ".join(unknown_sources)
            )
        force_active = record.get("force_active", False)
        if force_active not in (True, False):
            raise ActiveStratumError(
                f"{resolved}: {stratum_id} force_active must be boolean"
            )
        has_claim = (
            physical_status != "unknown"
            or analysis_included is not None
            or any(record.get(field) is not None for field in numeric_fields)
            or force_active
        )
        if has_claim and (
            not str(record.get("rationale") or "").strip() or not source_ids
        ):
            raise ActiveStratumError(
                f"{resolved}: {stratum_id} nondefault claims require a "
                "rationale and at least one source_id"
            )
        record["source_ids"] = source_ids
        record["force_active"] = force_active
        records[stratum_id] = record
    return records, list(sources.values())


def _load_calibrations(
    paths: Iterable[Path], config_sha256: str
) -> tuple[dict[str, dict[str, object]], list[dict[str, object]]]:
    records: dict[str, dict[str, object]] = {}
    sources: list[dict[str, object]] = []
    for raw_path in paths:
        path = raw_path.expanduser().resolve()
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema") != radiative_mode4.CALIBRATION_SCHEMA:
            raise ActiveStratumError(
                f"{path}: expected {radiative_mode4.CALIBRATION_SCHEMA}"
            )
        if payload.get("analysis_config_sha256") != config_sha256:
            raise ActiveStratumError(
                f"{path}: analysis configuration hash differs"
            )
        digest = _sha256(path)
        sources.append({"path": str(path), "sha256": digest})
        for raw in payload.get("strata", []):
            record = dict(raw)
            stratum_id = str(record.get("stratum_id", ""))
            if stratum_id in records:
                raise ActiveStratumError(
                    f"multiple canonical calibration reports contain "
                    f"{stratum_id}; pass only the report to use"
                )
            record["source_path"] = str(path)
            record["source_sha256"] = digest
            records[stratum_id] = record
    return records, sources


def _load_validations(
    paths: Iterable[Path],
    calibration_sources: list[dict[str, object]],
) -> tuple[dict[str, dict[str, object]], list[dict[str, object]]]:
    calibration_hashes = {str(item["sha256"]) for item in calibration_sources}
    records: dict[str, dict[str, object]] = {}
    sources: list[dict[str, object]] = []
    for raw_path in paths:
        path = raw_path.expanduser().resolve()
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema") != radiative_mode4.PILOT_VALIDATION_SCHEMA:
            raise ActiveStratumError(
                f"{path}: expected {radiative_mode4.PILOT_VALIDATION_SCHEMA}"
            )
        calibration_hash = str(payload.get("calibration_report_sha256", ""))
        if calibration_hash not in calibration_hashes:
            raise ActiveStratumError(
                f"{path}: its calibration report was not supplied"
            )
        digest = _sha256(path)
        sources.append({"path": str(path), "sha256": digest})
        for raw in payload.get("strata", []):
            record = dict(raw)
            stratum_id = str(record.get("stratum_id", ""))
            if stratum_id in records:
                raise ActiveStratumError(
                    f"multiple pilot-validation reports contain "
                    f"{stratum_id}; pass only the report to use"
                )
            record["source_path"] = str(path)
            record["source_sha256"] = digest
            record["calibration_report_sha256"] = calibration_hash
            records[stratum_id] = record
    return records, sources


def _calibration_positive(record: dict[str, object] | None) -> bool:
    if record is None:
        return False
    if float(record.get("integrated_cross_section_microbarn") or 0.0) > 0.0:
        return True
    inside = record.get("inside_guard") or {}
    complement = record.get("guard_complement") or {}
    return (
        int(inside.get("target_candidates") or 0) > 0
        or int(complement.get("target_candidates") or 0) > 0
    )


def _validate_catalog_metadata(
    stratum: radiative_guards.Stratum,
    record: dict[str, object] | None,
    *,
    label: str,
    required: bool,
) -> None:
    if record is None:
        return
    expected = {
        "flat_index": stratum.flat_index,
        "indices": _indices(stratum),
        "bounds": _bounds(stratum),
    }
    for name, value in expected.items():
        if name not in record:
            if required:
                raise ActiveStratumError(
                    f"{stratum.identifier}: {label} lacks {name}"
                )
            continue
        recorded_value = record[name]
        try:
            differs = (
                _normalized_indices(recorded_value)
                != _normalized_indices(value)
                if name == "indices"
                else recorded_value != value
            )
        except ValueError as error:
            raise ActiveStratumError(
                f"{stratum.identifier}: {label} has malformed indices: "
                f"{error}"
            ) from error
        if differs:
            raise ActiveStratumError(
                f"{stratum.identifier}: {label} {name} differs from the "
                "analysis catalog"
            )


def _classify_record(
    stratum: radiative_guards.Stratum,
    relevance: dict[str, object] | None,
    calibration: dict[str, object] | None,
    validation: dict[str, object] | None,
    *,
    minimum_data_events: int,
    maximum_model_fraction: float,
    maximum_feed_in_fraction: float,
    maximum_closure_fraction: float,
    reference_cross_section: float | None,
) -> dict[str, object]:
    relevance = dict(relevance or {})
    _validate_catalog_metadata(
        stratum, relevance, label="relevance evidence", required=bool(relevance)
    )
    _validate_catalog_metadata(
        stratum, calibration, label="calibration", required=False
    )
    _validate_catalog_metadata(
        stratum, validation, label="pilot validation", required=False
    )
    if validation is not None and (
        calibration is None
        or validation.get("calibration_report_sha256")
        != calibration.get("source_sha256")
    ):
        raise ActiveStratumError(
            f"{stratum.identifier}: pilot validation does not reference its "
            "supplied canonical calibration"
        )
    explicit_physical = str(relevance.get("physical_status", "unknown"))
    positive_calibration = _calibration_positive(calibration)
    if explicit_physical == "structurally_empty" and positive_calibration:
        raise ActiveStratumError(
            f"{stratum.identifier}: structural-empty claim contradicts a "
            "positive calibration target or cross section"
        )
    physical_status = explicit_physical
    physical_basis = "external_relevance_evidence"
    if physical_status == "unknown" and positive_calibration:
        physical_status = "nonempty_observed"
        physical_basis = "positive_mode4_calibration"

    model_fraction = relevance.get("model_cross_section_fraction")
    if (
        model_fraction is None
        and reference_cross_section is not None
        and calibration is not None
    ):
        model_fraction = (
            float(calibration.get("integrated_cross_section_microbarn") or 0.0)
            / reference_cross_section
        )
    data_events = relevance.get("data_events")
    analysis_included = relevance.get("analysis_included")
    feed_in = relevance.get("maximum_feed_in_fraction")
    closure = relevance.get("global_closure_impact_fraction")
    if explicit_physical == "structurally_empty" and (
        analysis_included is True
        or (data_events is not None and int(data_events) > 0)
        or bool(relevance.get("force_active", False))
        or any(
            value is not None and float(value) > 0.0
            for value in (model_fraction, feed_in, closure)
        )
    ):
        raise ActiveStratumError(
            f"{stratum.identifier}: structural-empty claim contradicts "
            "positive relevance evidence"
        )
    reasons: list[str] = []
    if bool(relevance.get("force_active", False)):
        reasons.append("explicit_force_active")
    if analysis_included is True:
        reasons.append("included_analysis_bin")
    if data_events is not None and int(data_events) >= minimum_data_events:
        reasons.append("data_occupancy")
    if (
        model_fraction is not None
        and float(model_fraction) > maximum_model_fraction
    ):
        reasons.append("material_model_cross_section")
    if feed_in is not None and float(feed_in) > maximum_feed_in_fraction:
        reasons.append("material_feed_in")
    if closure is not None and float(closure) > maximum_closure_fraction:
        reasons.append("material_global_closure_impact")

    complete_for_exclusion = (
        physical_status != "unknown"
        and analysis_included is False
        and data_events is not None
        and model_fraction is not None
        and feed_in is not None
        and closure is not None
    )
    calibration_ready = bool(
        calibration is not None
        and calibration.get("pilot_readiness")
        in ("ready", "ready_provisional_zero_complement")
        and calibration.get("recommended_envelope") is not None
    )
    pilot_passed = bool(validation is not None and validation.get("passed"))

    if explicit_physical == "structurally_empty":
        status = "structurally_empty"
        reasons = ["externally_proven_structural_emptiness"]
    elif reasons:
        if calibration_ready and pilot_passed:
            status = "active_ready"
        elif calibration_ready:
            status = "active_needs_pilot_validation"
        else:
            status = "active_needs_optimization"
    elif complete_for_exclusion:
        status = "closure_only"
        reasons = ["complete_relevance_evidence_below_all_thresholds"]
    else:
        status = "needs_relevance_assessment"
        reasons = ["relevance_evidence_incomplete"]

    if status not in ALL_STATUSES:
        raise AssertionError(status)
    return {
        "stratum_id": stratum.identifier,
        "flat_index": stratum.flat_index,
        "indices": _indices(stratum),
        "bounds": _bounds(stratum),
        "status": status,
        "active": status in ACTIVE_STATUSES,
        "production_ready": status == "active_ready",
        "reasons": reasons,
        "physical_status": physical_status,
        "physical_status_basis": physical_basis,
        "relevance_complete_for_exclusion": complete_for_exclusion,
        "relevance": {
            "analysis_included": analysis_included,
            "data_events": data_events,
            "model_cross_section_fraction": model_fraction,
            "maximum_feed_in_fraction": feed_in,
            "global_closure_impact_fraction": closure,
            "force_active": bool(relevance.get("force_active", False)),
            "rationale": relevance.get("rationale"),
            "source_ids": relevance.get("source_ids", []),
        },
        "calibration": (
            {
                "source_path": calibration.get("source_path"),
                "source_sha256": calibration.get("source_sha256"),
                "integrated_cross_section_microbarn": calibration.get(
                    "integrated_cross_section_microbarn"
                ),
                "integrated_cross_section_sem_microbarn": calibration.get(
                    "integrated_cross_section_sem_microbarn"
                ),
                "recommendation_status": calibration.get(
                    "recommendation_status"
                ),
                "pilot_readiness": calibration.get("pilot_readiness"),
                "recommended_sigr_max": (
                    calibration.get("recommended_envelope") or {}
                ).get("sigr_max"),
            }
            if calibration is not None
            else None
        ),
        "pilot_validation": (
            {
                "source_path": validation.get("source_path"),
                "source_sha256": validation.get("source_sha256"),
                "passed": validation.get("passed"),
                "recommendation": validation.get("recommendation"),
                "relative_cross_section_difference": validation.get(
                    "relative_cross_section_difference"
                ),
                "cross_section_difference_z_score": validation.get(
                    "cross_section_difference_z_score"
                ),
            }
            if validation is not None
            else None
        ),
    }


def classify(args: argparse.Namespace) -> Path:
    config_path = args.config.expanduser().resolve()
    config, config_sha256 = radiative_guards.load_analysis_config(config_path)
    output = args.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"{output} already exists")
    for name, value in (
        ("maximum model fraction", args.maximum_model_fraction),
        ("maximum feed-in fraction", args.maximum_feed_in_fraction),
        ("maximum closure fraction", args.maximum_closure_fraction),
    ):
        if not 0.0 <= value < 1.0:
            raise ValueError(f"{name} must lie in [0,1)")
    if args.minimum_data_events < 1:
        raise ValueError("--minimum-data-events must be positive")
    if (
        args.reference_cross_section_microbarn is not None
        and args.reference_cross_section_microbarn <= 0.0
    ):
        raise ValueError(
            "--reference-cross-section-microbarn must be positive"
        )

    relevance_path = args.relevance.expanduser().resolve()
    relevance, relevance_sources = _load_relevance(
        relevance_path, config_sha256
    )
    calibrations, calibration_sources = _load_calibrations(
        args.calibrations or [], config_sha256
    )
    validations, validation_sources = _load_validations(
        args.pilot_validations or [], calibration_sources
    )
    strata = _selected_strata(config, args.flat_indices)
    catalog_ids = {item.identifier for item in strata}
    unknown_evidence = sorted(set(relevance) - catalog_ids)
    if unknown_evidence:
        raise ActiveStratumError(
            "relevance evidence contains unselected or unknown strata: "
            + ", ".join(unknown_evidence)
        )
    records = [
        _classify_record(
            stratum,
            relevance.get(stratum.identifier),
            calibrations.get(stratum.identifier),
            validations.get(stratum.identifier),
            minimum_data_events=args.minimum_data_events,
            maximum_model_fraction=args.maximum_model_fraction,
            maximum_feed_in_fraction=args.maximum_feed_in_fraction,
            maximum_closure_fraction=args.maximum_closure_fraction,
            reference_cross_section=args.reference_cross_section_microbarn,
        )
        for stratum in strata
    ]
    counts = {
        status: sum(record["status"] == status for record in records)
        for status in sorted(ALL_STATUSES)
    }
    payload = {
        "schema": MASK_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "classifier_revision": radiative_mode4._source_revision(),
        "classifier_source_sha256": _sha256(Path(__file__).resolve()),
        "analysis_config": str(config_path),
        "analysis_config_sha256": config_sha256,
        "relevance_evidence": str(relevance_path),
        "relevance_evidence_sha256": _sha256(relevance_path),
        "evidence_sources": relevance_sources,
        "calibration_sources": calibration_sources,
        "pilot_validation_sources": validation_sources,
        "thresholds": {
            "minimum_data_events": args.minimum_data_events,
            "maximum_model_cross_section_fraction_for_closure_only": (
                args.maximum_model_fraction
            ),
            "maximum_feed_in_fraction_for_closure_only": (
                args.maximum_feed_in_fraction
            ),
            "maximum_global_closure_impact_fraction_for_closure_only": (
                args.maximum_closure_fraction
            ),
            "reference_cross_section_microbarn": (
                args.reference_cross_section_microbarn
            ),
        },
        "zero_mc_targets_never_prove_structural_emptiness": True,
        "closure_only_requires_complete_external_relevance_evidence": True,
        "status_counts": counts,
        "strata": records,
    }
    output.mkdir(parents=True)
    json_path = output / "active_stratum_mask.json"
    _write_json(json_path, payload)
    fields = (
        "flat_index",
        "stratum_id",
        "status",
        "active",
        "production_ready",
        "physical_status",
        "analysis_included",
        "data_events",
        "model_cross_section_fraction",
        "maximum_feed_in_fraction",
        "global_closure_impact_fraction",
        "calibration_cross_section_microbarn",
        "calibration_pilot_readiness",
        "pilot_validation_passed",
        "reasons",
    )
    with (output / "active_stratum_mask.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        writer = csv.DictWriter(destination, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for record in records:
            calibration = record["calibration"] or {}
            validation = record["pilot_validation"] or {}
            rel = record["relevance"]
            writer.writerow(
                {
                    "flat_index": record["flat_index"],
                    "stratum_id": record["stratum_id"],
                    "status": record["status"],
                    "active": record["active"],
                    "production_ready": record["production_ready"],
                    "physical_status": record["physical_status"],
                    "analysis_included": rel["analysis_included"],
                    "data_events": rel["data_events"],
                    "model_cross_section_fraction": rel[
                        "model_cross_section_fraction"
                    ],
                    "maximum_feed_in_fraction": rel[
                        "maximum_feed_in_fraction"
                    ],
                    "global_closure_impact_fraction": rel[
                        "global_closure_impact_fraction"
                    ],
                    "calibration_cross_section_microbarn": calibration.get(
                        "integrated_cross_section_microbarn"
                    ),
                    "calibration_pilot_readiness": calibration.get(
                        "pilot_readiness"
                    ),
                    "pilot_validation_passed": validation.get("passed"),
                    "reasons": ",".join(record["reasons"]),
                }
            )
    for status in sorted(ALL_STATUSES):
        values = [
            str(record["flat_index"])
            for record in records
            if record["status"] == status
        ]
        (output / f"{status}_flat_indices.txt").write_text(
            "".join(f"{value}\n" for value in values), encoding="utf-8"
        )
    active = [
        str(record["flat_index"]) for record in records if record["active"]
    ]
    ready = [
        str(record["flat_index"])
        for record in records
        if record["production_ready"]
    ]
    (output / "active_flat_indices.txt").write_text(
        "".join(f"{value}\n" for value in active), encoding="utf-8"
    )
    (output / "production_ready_flat_indices.txt").write_text(
        "".join(f"{value}\n" for value in ready), encoding="utf-8"
    )
    return json_path


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    template = subparsers.add_parser(
        "create-template",
        help="create an evidence-hashed relevance-assessment template",
    )
    template.add_argument("--config", type=Path, required=True)
    template.add_argument("--output", type=Path, required=True)
    template.add_argument(
        "--flat-index", dest="flat_indices", type=int, action="append"
    )
    template.add_argument(
        "--evidence-source",
        dest="evidence_sources",
        action="append",
        help="external evidence artifact as IDENTIFIER=PATH",
    )

    classifier = subparsers.add_parser(
        "classify",
        help="classify physical support, relevance, and generator readiness",
    )
    classifier.add_argument("--config", type=Path, required=True)
    classifier.add_argument("--relevance", type=Path, required=True)
    classifier.add_argument(
        "--calibration",
        dest="calibrations",
        type=Path,
        action="append",
    )
    classifier.add_argument(
        "--pilot-validation",
        dest="pilot_validations",
        type=Path,
        action="append",
    )
    classifier.add_argument("--output", type=Path, required=True)
    classifier.add_argument(
        "--flat-index", dest="flat_indices", type=int, action="append"
    )
    classifier.add_argument("--minimum-data-events", type=int, default=1)
    classifier.add_argument(
        "--maximum-model-fraction", type=float, required=True
    )
    classifier.add_argument(
        "--maximum-feed-in-fraction", type=float, required=True
    )
    classifier.add_argument(
        "--maximum-closure-fraction", type=float, required=True
    )
    classifier.add_argument(
        "--reference-cross-section-microbarn", type=float
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    result = (
        create_template(args)
        if args.command == "create-template"
        else classify(args)
    )
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
