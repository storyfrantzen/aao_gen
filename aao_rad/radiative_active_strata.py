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

import numpy as np

import radiative_guards
import radiative_migrations
import radiative_mode4


DATA_OCCUPANCY_SCHEMA = "aao-rad-data-occupancy-v1"
SURVEY_MODEL_SCHEMA = "aao-rad-survey-model-evidence-v1"
EVIDENCE_SCHEMA = "aao-rad-stratum-relevance-evidence-v1"
MASK_SCHEMA = "aao-rad-active-stratum-mask-v1"
CUMULATIVE_QUEUE_SCHEMA = "aao-rad-cumulative-stratum-queue-v1"
STRATIFIED_BATCH_SCHEMA = "aao-rad-stratified-calibration-batch-v1"
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


def _one_dimensional_array(
    sample: np.lib.npyio.NpzFile,
    name: str,
    expected_rows: int | None = None,
) -> np.ndarray:
    if name not in sample.files:
        raise ActiveStratumError(f"data event sample lacks {name}")
    values = np.asarray(sample[name])
    if values.ndim != 1:
        raise ActiveStratumError(f"data event sample {name} is not one-dimensional")
    if expected_rows is not None and values.size != expected_rows:
        raise ActiveStratumError(
            f"data event sample {name} has {values.size} rows; "
            f"expected {expected_rows}"
        )
    return values


def _array_bin_indices(values: np.ndarray, edges: Iterable[float]) -> np.ndarray:
    edge_array = np.asarray(tuple(edges), dtype=float)
    indices = np.searchsorted(edge_array, values, side="right") - 1
    valid = (
        np.isfinite(values)
        & (indices >= 0)
        & (indices < edge_array.size - 1)
    )
    return np.where(valid, indices, -1).astype(np.int64, copy=False)


def _event_key_audit(
    run: np.ndarray | None,
    event: np.ndarray | None,
    selected: np.ndarray,
) -> dict[str, object]:
    if run is None or event is None:
        return {
            "available": False,
            "selected_rows": int(selected.sum()),
            "unique_keys": None,
            "duplicate_rows": None,
        }
    if run.shape != selected.shape or event.shape != selected.shape:
        raise ActiveStratumError(
            "data event run/event arrays do not match the coordinate rows"
        )
    try:
        run_values = np.asarray(run[selected], dtype=np.int64)
        event_values = np.asarray(event[selected], dtype=np.int64)
    except (TypeError, ValueError, OverflowError) as error:
        raise ActiveStratumError(
            "data event run/event arrays cannot be interpreted as integers"
        ) from error
    keys = np.column_stack((run_values, event_values))
    unique = int(np.unique(keys, axis=0).shape[0])
    rows = int(keys.shape[0])
    return {
        "available": True,
        "selected_rows": rows,
        "unique_keys": unique,
        "duplicate_rows": rows - unique,
    }


def build_data_evidence(args: argparse.Namespace) -> Path:
    config_path = args.config.expanduser().resolve()
    config, config_sha256 = radiative_guards.load_analysis_config(config_path)
    data_path = args.data_events.expanduser().resolve()
    mask_path = args.selection_mask.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"{output} already exists")
    if not data_path.is_file():
        raise FileNotFoundError(data_path)
    if not mask_path.is_file():
        raise FileNotFoundError(mask_path)

    if "target_mass" not in config:
        raise ActiveStratumError("analysis config is missing target_mass")
    target_mass = float(config["target_mass"])
    beam_energy = float(config["beam_energy"])
    if not math.isfinite(target_mass) or target_mass <= 0.0:
        raise ActiveStratumError("analysis target_mass must be positive")
    if not math.isfinite(beam_energy) or beam_energy <= 0.0:
        raise ActiveStratumError("analysis beam_energy must be positive")

    with np.load(data_path, allow_pickle=False) as sample:
        q2 = np.asarray(_one_dimensional_array(sample, "rec_Q2"), dtype=float)
        rows = q2.size
        xb = np.asarray(
            _one_dimensional_array(sample, "rec_xB", rows), dtype=float
        )
        minus_t = np.asarray(
            _one_dimensional_array(sample, "rec_minus_t", rows), dtype=float
        )
        phi_rad = np.asarray(
            _one_dimensional_array(sample, "rec_trento_phi", rows), dtype=float
        )
        run = (
            _one_dimensional_array(sample, "run", rows)
            if "run" in sample.files
            else None
        )
        event = (
            _one_dimensional_array(sample, "event", rows)
            if "event" in sample.files
            else None
        )

    raw_mask = np.asarray(np.load(mask_path, allow_pickle=False))
    if raw_mask.ndim != 1 or raw_mask.size != rows:
        raise ActiveStratumError(
            f"selection mask has shape {raw_mask.shape}; expected ({rows},)"
        )
    if raw_mask.dtype != np.bool_:
        raise ActiveStratumError("selection mask must have boolean dtype")
    selected = np.ones(rows, dtype=bool)

    cut_flow: list[dict[str, object]] = [
        {
            "stage": "input_rows",
            "applied": True,
            "rows_before": rows,
            "rows_removed": 0,
            "rows_remaining": rows,
        }
    ]

    def apply_cut(
        name: str,
        condition: np.ndarray,
        *,
        applied: bool = True,
    ) -> None:
        nonlocal selected
        before = int(selected.sum())
        if applied:
            selected &= condition
        remaining = int(selected.sum())
        cut_flow.append(
            {
                "stage": name,
                "applied": applied,
                "rows_before": before,
                "rows_removed": before - remaining,
                "rows_remaining": remaining,
            }
        )

    apply_cut("exclusivity_selection_mask", raw_mask)
    finite_coordinates = (
        np.isfinite(q2)
        & np.isfinite(xb)
        & np.isfinite(minus_t)
        & np.isfinite(phi_rad)
    )
    apply_cut("finite_analysis_coordinates", finite_coordinates)
    apply_cut("positive_xB", xb > 0.0)

    phase_space = config.get("phase_space", {})
    q2_minimum = phase_space.get("Q2_min")
    if q2_minimum is None:
        apply_cut("Q2_min", np.ones(rows, dtype=bool), applied=False)
    else:
        q2_minimum = float(q2_minimum)
        if not math.isfinite(q2_minimum) or q2_minimum < 0.0:
            raise ActiveStratumError(
                "phase_space.Q2_min must be finite and nonnegative"
            )
        apply_cut("Q2_min", q2 >= q2_minimum)

    with np.errstate(divide="ignore", invalid="ignore"):
        y = q2 / (2.0 * target_mass * beam_energy * xb)
        w_squared = target_mass * target_mass + q2 * (1.0 / xb - 1.0)
    w_minimum = float(phase_space["W_min"])
    if not math.isfinite(w_minimum) or w_minimum <= 0.0:
        raise ActiveStratumError("phase_space.W_min must be positive")
    apply_cut(
        "W_min",
        np.isfinite(w_squared) & (w_squared >= w_minimum * w_minimum),
    )

    y_maximum = phase_space.get("y_max")
    if y_maximum is None:
        apply_cut("y_max", np.ones(rows, dtype=bool), applied=False)
    else:
        y_maximum = float(y_maximum)
        if not 0.0 < y_maximum <= 1.0:
            raise ActiveStratumError("phase_space.y_max must lie in (0,1]")
        apply_cut("y_max", np.isfinite(y) & (y <= y_maximum))

    binning = config["binning"]
    iq2 = _array_bin_indices(q2, binning["Q2"])
    ixb = _array_bin_indices(xb, binning["xB"])
    it = _array_bin_indices(minus_t, binning["minus_t"])
    phi_edges = tuple(float(value) for value in binning["phi_deg"])
    phi_origin = phi_edges[0]
    phi_deg = phi_origin + np.mod(
        np.degrees(phi_rad) - phi_origin, 360.0
    )
    iphi = _array_bin_indices(phi_deg, phi_edges)
    inside_binning = (iq2 >= 0) & (ixb >= 0) & (it >= 0) & (iphi >= 0)
    apply_cut("analysis_binning", inside_binning)

    nq2 = len(binning["Q2"]) - 1
    nt = len(binning["minus_t"]) - 1
    nphi = len(binning["phi_deg"]) - 1
    flat = ixb * nq2 * nphi * nt + iq2 * nphi * nt + iphi * nt + it
    flat = np.where(inside_binning, flat, -1).astype(np.int64, copy=False)
    catalog = radiative_guards.enumerate_strata(config)
    counts = np.bincount(flat[selected], minlength=len(catalog)).astype(
        np.int64, copy=False
    )
    if counts.size != len(catalog):
        raise AssertionError("data occupancy does not match the catalog size")
    selected_rows = int(selected.sum())
    if int(counts.sum()) != selected_rows:
        raise AssertionError("data occupancy does not conserve selected rows")

    key_audit = _event_key_audit(run, event, selected)
    if (
        not args.allow_duplicate_event_keys
        and int(key_audit.get("duplicate_rows") or 0) > 0
    ):
        raise ActiveStratumError(
            "selected data contain duplicate (run,event) keys; pass "
            "--allow-duplicate-event-keys only after auditing the overlap"
        )

    strata = []
    # Python 3.9 on the JLab farm predates zip(..., strict=True).  The size
    # equality is already asserted above, so ordinary zip has the same
    # fail-closed behavior here while retaining farm compatibility.
    for stratum, count_value in zip(catalog, counts):
        count = int(count_value)
        strata.append(
            {
                "stratum_id": stratum.identifier,
                "flat_index": stratum.flat_index,
                "indices": _indices(stratum),
                "bounds": _bounds(stratum),
                "data_events": count,
                "data_event_fraction": (
                    count / selected_rows if selected_rows else 0.0
                ),
            }
        )

    occupancy_payload = {
        "schema": DATA_OCCUPANCY_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "builder_revision": radiative_mode4._source_revision(),
        "builder_source_sha256": _sha256(Path(__file__).resolve()),
        "analysis_config": str(config_path),
        "analysis_config_sha256": config_sha256,
        "data_events": {
            "path": str(data_path),
            "sha256": _sha256(data_path),
        },
        "selection_mask": {
            "path": str(mask_path),
            "sha256": _sha256(mask_path),
        },
        "coordinate_definition": {
            "source": "reconstructed selected EPPI0 candidate",
            "Q2": "rec_Q2",
            "xB": "rec_xB",
            "minus_t": "rec_minus_t",
            "phi": "rec_trento_phi in radians, wrapped periodically",
            "y_formula": "Q2/(2*target_mass*beam_energy*xB)",
            "W2_formula": "target_mass^2+Q2*(1/xB-1)",
        },
        "selection_policy": {
            "Q2_min": q2_minimum,
            "W_min": w_minimum,
            "y_max": y_maximum,
            "y_max_applied": y_maximum is not None,
            "analysis_bin_edges_are_half_open": True,
        },
        "cut_flow": cut_flow,
        "event_key_audit": key_audit,
        "catalog_summary": {
            "strata_total": len(catalog),
            "selected_data_events": selected_rows,
            "occupied_strata": int(np.count_nonzero(counts)),
            "strata_with_at_least_5_events": int(np.count_nonzero(counts >= 5)),
            "strata_with_at_least_10_events": int(
                np.count_nonzero(counts >= 10)
            ),
            "strata_with_at_least_50_events": int(
                np.count_nonzero(counts >= 50)
            ),
        },
        "zero_data_events_do_not_prove_structural_emptiness": True,
        "strata": strata,
    }

    output.mkdir(parents=True)
    occupancy_path = output / "data_occupancy.json"
    _write_json(occupancy_path, occupancy_payload)
    occupancy_sha256 = _sha256(occupancy_path)
    source_identifier = "selected_data_occupancy"
    relevance_payload = {
        "schema": EVIDENCE_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "analysis_config": str(config_path),
        "analysis_config_sha256": config_sha256,
        "instructions": {
            "physical_status_values": sorted(PHYSICAL_STATUSES),
            "fractions_are_relative_to_the_full_analysis_prediction": True,
            "zero_mc_targets_do_not_prove_structural_emptiness": True,
            "zero_data_events_do_not_prove_structural_emptiness": True,
            "data_occupancy_alone_does_not_define_analysis_inclusion": True,
            "nondefault_claims_require_rationale_and_source_ids": True,
        },
        "sources": [
            {
                "identifier": source_identifier,
                "path": str(occupancy_path),
                "sha256": occupancy_sha256,
            }
        ],
        "strata": [
            {
                "stratum_id": record["stratum_id"],
                "flat_index": record["flat_index"],
                "indices": record["indices"],
                "bounds": record["bounds"],
                "physical_status": "unknown",
                "analysis_included": None,
                "data_events": record["data_events"],
                "model_cross_section_fraction": None,
                "maximum_feed_in_fraction": None,
                "global_closure_impact_fraction": None,
                "force_active": False,
                "rationale": (
                    "Observed exclusivity-selected reconstructed data "
                    "occupancy; a zero count is finite-sample evidence only "
                    "and does not establish structural emptiness."
                ),
                "source_ids": [source_identifier],
            }
            for record in strata
        ],
    }
    relevance_path = output / "relevance_evidence.json"
    _write_json(relevance_path, relevance_payload)

    fields = (
        "flat_index",
        "stratum_id",
        "Q2",
        "xB",
        "minus_t",
        "phi_deg",
        "data_events",
        "data_event_fraction",
    )
    with (output / "data_occupancy.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        writer = csv.DictWriter(destination, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for record in strata:
            bounds = record["bounds"]
            writer.writerow(
                {
                    "flat_index": record["flat_index"],
                    "stratum_id": record["stratum_id"],
                    "Q2": ":".join(str(value) for value in bounds["Q2"]),
                    "xB": ":".join(str(value) for value in bounds["xB"]),
                    "minus_t": ":".join(
                        str(value) for value in bounds["minus_t"]
                    ),
                    "phi_deg": ":".join(
                        str(value) for value in bounds["phi_deg"]
                    ),
                    "data_events": record["data_events"],
                    "data_event_fraction": record["data_event_fraction"],
                }
            )

    for threshold in (1, 5, 10, 50):
        label = "occupied" if threshold == 1 else f"at_least_{threshold}_events"
        values = [
            str(record["flat_index"])
            for record in strata
            if int(record["data_events"]) >= threshold
        ]
        (output / f"{label}_flat_indices.txt").write_text(
            "".join(f"{value}\n" for value in values), encoding="utf-8"
        )
    return relevance_path


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


def _moment_from_metrics(
    metrics: dict[str, object] | None,
    *,
    label: str,
) -> radiative_guards.Moment:
    if metrics is None:
        return radiative_guards.Moment()
    fields = {
        "count": "contributing_rows",
        "total": "sum_trial_contributions_microbarn",
        "square_total": "sum_squared_trial_contributions_microbarn2",
        "maximum": "largest_trial_contribution_microbarn",
    }
    values: dict[str, float | int] = {}
    try:
        values["count"] = int(metrics[fields["count"]])
        for name in ("total", "square_total", "maximum"):
            values[name] = float(metrics[fields[name]])
    except (KeyError, TypeError, ValueError, OverflowError) as error:
        raise ActiveStratumError(f"{label}: malformed fixed-trial metrics") from error
    count = int(values["count"])
    numeric = [float(values[name]) for name in ("total", "square_total", "maximum")]
    if count < 0 or any(not math.isfinite(value) or value < 0.0 for value in numeric):
        raise ActiveStratumError(
            f"{label}: fixed-trial sufficient statistics must be nonnegative"
        )
    total, square_total, maximum = numeric
    if count == 0 and any(value != 0.0 for value in numeric):
        raise ActiveStratumError(
            f"{label}: nonzero sufficient statistics have zero contributing rows"
        )
    return radiative_guards.Moment(
        count=count,
        total=total,
        square_total=square_total,
        maximum=maximum,
    )


def _add_moments(
    first: radiative_guards.Moment,
    second: radiative_guards.Moment,
) -> radiative_guards.Moment:
    return radiative_guards.Moment(
        count=first.count + second.count,
        total=first.total + second.total,
        square_total=first.square_total + second.square_total,
        maximum=max(first.maximum, second.maximum),
    )


def _selection_matches_config(selection: dict, config: dict) -> None:
    phase_space = config["phase_space"]
    expected_q2 = float(
        phase_space.get("Q2_min", config["binning"]["Q2"][0])
    )
    expected_w = float(phase_space["W_min"])
    expected_y = phase_space.get("y_max")
    expected_apply_y = expected_y is not None
    checks = (
        ("q2_minimum", expected_q2),
        ("w_minimum", expected_w),
    )
    for name, expected in checks:
        try:
            recorded = float(selection[name])
        except (KeyError, TypeError, ValueError) as error:
            raise ActiveStratumError(
                f"migration manifest has malformed analysis_selection.{name}"
            ) from error
        if not math.isclose(recorded, expected, rel_tol=0.0, abs_tol=1.0e-12):
            raise ActiveStratumError(
                f"migration analysis_selection.{name} differs from the "
                "analysis configuration"
            )
    if bool(selection.get("apply_y_max")) != expected_apply_y:
        raise ActiveStratumError(
            "migration y_max policy differs from the analysis configuration"
        )
    recorded_y = selection.get("y_maximum")
    if expected_apply_y and (
        recorded_y is None
        or not math.isclose(
            float(recorded_y), float(expected_y), rel_tol=0.0, abs_tol=1.0e-12
        )
    ):
        raise ActiveStratumError(
            "migration analysis_selection.y_maximum differs from the "
            "analysis configuration"
        )
    if not expected_apply_y and recorded_y is not None:
        raise ActiveStratumError(
            "migration unexpectedly records an active y_maximum"
        )


def augment_survey_evidence(args: argparse.Namespace) -> Path:
    config_path = args.config.expanduser().resolve()
    config, config_sha256 = radiative_guards.load_analysis_config(config_path)
    base_path = args.base_relevance.expanduser().resolve()
    manifest_path = args.migration_manifest.expanduser().resolve()
    validation_path = args.migration_validation.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"{output} already exists")
    for path in (base_path, manifest_path, validation_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    base_records, base_sources = _load_relevance(base_path, config_sha256)
    base_payload = json.loads(base_path.read_text(encoding="utf-8"))
    if any(
        record.get("model_cross_section_fraction") is not None
        for record in base_records.values()
    ):
        raise ActiveStratumError(
            "base relevance already contains model_cross_section_fraction; "
            "refusing to overwrite existing model evidence"
        )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != radiative_migrations.MANIFEST_SCHEMA:
        raise ActiveStratumError(
            f"{manifest_path}: expected {radiative_migrations.MANIFEST_SCHEMA}"
        )
    if validation.get("schema") != radiative_migrations.VALIDATION_SCHEMA:
        raise ActiveStratumError(
            f"{validation_path}: expected {radiative_migrations.VALIDATION_SCHEMA}"
        )
    if manifest.get("analysis_config_sha256") != config_sha256:
        raise ActiveStratumError(
            "migration manifest analysis configuration hash differs"
        )
    manifest_sha256 = _sha256(manifest_path)
    if validation.get("manifest_sha256") != manifest_sha256:
        raise ActiveStratumError(
            "migration validation does not reference the supplied manifest"
        )
    _selection_matches_config(manifest.get("analysis_selection") or {}, config)

    catalog = radiative_guards.enumerate_strata(config)
    catalog_ids = {stratum.identifier for stratum in catalog}
    for label, identifiers in (
        ("base relevance", set(base_records)),
        ("migration manifest", set(manifest.get("strata", {}))),
    ):
        if identifiers != catalog_ids:
            missing = len(catalog_ids - identifiers)
            extra = len(identifiers - catalog_ids)
            raise ActiveStratumError(
                f"{label} is not the full analysis catalog "
                f"(missing={missing}, extra={extra})"
            )
    validation_records = validation.get("strata", {})
    if not isinstance(validation_records, dict) or not set(
        validation_records
    ).issubset(catalog_ids):
        raise ActiveStratumError(
            "migration validation contains unknown or malformed strata"
        )

    try:
        training_proposals = int(manifest["training"]["total_proposals"])
        holdout_proposals = int(validation["validation"]["total_proposals"])
    except (KeyError, TypeError, ValueError) as error:
        raise ActiveStratumError(
            "migration artifacts lack fixed-trial proposal totals"
        ) from error
    if training_proposals <= 0 or holdout_proposals <= 0:
        raise ActiveStratumError(
            "migration fixed-trial proposal totals must be positive"
        )
    pooled_proposals = training_proposals + holdout_proposals
    training_inside = _moment_from_metrics(
        manifest["training"].get("inside_analysis_partition"),
        label="training inside-analysis partition",
    )
    holdout_inside = _moment_from_metrics(
        validation["validation"].get("inside_analysis_partition"),
        label="holdout inside-analysis partition",
    )
    pooled_inside = _add_moments(training_inside, holdout_inside)
    if pooled_inside.total <= 0.0:
        raise ActiveStratumError(
            "pooled surveys have zero inside-analysis cross section"
        )

    support_counts = {
        "independent_support": 0,
        "training_only": 0,
        "holdout_only": 0,
        "no_survey_contribution": 0,
    }
    coverage_counts = {
        "passed": 0,
        "failed": 0,
        "no_holdout_contribution": 0,
        "not_assessed": 0,
    }
    survey_records: list[dict[str, object]] = []
    pooled_strata_moment = radiative_guards.Moment()
    for stratum in catalog:
        training_record = manifest["strata"][stratum.identifier]
        _validate_catalog_metadata(
            stratum,
            training_record,
            label="migration training",
            required=True,
        )
        validation_record = validation_records.get(stratum.identifier)
        training_moment = _moment_from_metrics(
            training_record.get("training_total"),
            label=f"{stratum.identifier} training",
        )
        holdout_moment = _moment_from_metrics(
            (validation_record or {}).get("holdout_total"),
            label=f"{stratum.identifier} holdout",
        )
        pooled_moment = _add_moments(training_moment, holdout_moment)
        pooled_strata_moment = _add_moments(
            pooled_strata_moment, pooled_moment
        )
        training_nonzero = training_moment.total > 0.0
        holdout_nonzero = holdout_moment.total > 0.0
        if training_nonzero and holdout_nonzero:
            support = "independent_support"
        elif training_nonzero:
            support = "training_only"
        elif holdout_nonzero:
            support = "holdout_only"
        else:
            support = "no_survey_contribution"
        support_counts[support] += 1

        coverage_value = (
            validation_record.get("coverage_passed")
            if validation_record is not None
            else None
        )
        if coverage_value is True:
            coverage = "passed"
        elif coverage_value is False:
            coverage = "failed"
        elif validation_record is not None:
            coverage = "no_holdout_contribution"
        else:
            coverage = "not_assessed"
        coverage_counts[coverage] += 1

        pooled_metrics = radiative_guards._metrics(
            pooled_moment, pooled_proposals
        )
        model_fraction = pooled_moment.total / pooled_inside.total
        survey_records.append(
            {
                "stratum_id": stratum.identifier,
                "flat_index": stratum.flat_index,
                "indices": _indices(stratum),
                "bounds": _bounds(stratum),
                "data_events": base_records[stratum.identifier].get(
                    "data_events"
                ),
                "training_status": training_record.get("status"),
                "training": radiative_guards._metrics(
                    training_moment, training_proposals
                ),
                "holdout": radiative_guards._metrics(
                    holdout_moment, holdout_proposals
                ),
                "pooled": pooled_metrics,
                "model_cross_section_fraction": model_fraction,
                "support_status": support,
                "parent_coverage_status": coverage,
                "training_holdout_difference_z_score": (
                    validation_record.get(
                        "training_holdout_difference_z_score"
                    )
                    if validation_record is not None
                    else None
                ),
            }
        )

    if not math.isclose(
        pooled_strata_moment.total,
        pooled_inside.total,
        rel_tol=5.0e-12,
        abs_tol=1.0e-12,
    ):
        raise ActiveStratumError(
            "pooled per-stratum contributions do not close to the "
            "inside-analysis partition"
        )
    fraction_sum = sum(
        float(record["model_cross_section_fraction"])
        for record in survey_records
    )
    if not math.isclose(fraction_sum, 1.0, rel_tol=5.0e-12, abs_tol=1.0e-12):
        raise ActiveStratumError(
            "pooled model cross-section fractions do not sum to one"
        )

    source_identifier = "pooled_radiative_survey_model"
    if source_identifier in {source["identifier"] for source in base_sources}:
        raise ActiveStratumError(
            f"base relevance already defines source {source_identifier}"
        )

    survey_payload = {
        "schema": SURVEY_MODEL_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "builder_revision": radiative_mode4._source_revision(),
        "builder_source_sha256": _sha256(Path(__file__).resolve()),
        "analysis_config": str(config_path),
        "analysis_config_sha256": config_sha256,
        "base_relevance": {
            "path": str(base_path),
            "sha256": _sha256(base_path),
        },
        "migration_manifest": {
            "path": str(manifest_path),
            "sha256": manifest_sha256,
            "generator_revision": manifest.get("generator_revision"),
            "generator_revision_source": manifest.get(
                "generator_revision_source"
            ),
        },
        "migration_validation": {
            "path": str(validation_path),
            "sha256": _sha256(validation_path),
            "passed": validation.get("passed"),
            "global_training_holdout_difference_z_score": validation.get(
                "global_training_holdout_difference_z_score"
            ),
        },
        "analysis_selection": manifest["analysis_selection"],
        "pooling": {
            "estimator": (
                "sum of fixed-trial contribution sums divided by the sum "
                "of fixed-trial proposal counts"
            ),
            "training_proposals": training_proposals,
            "holdout_proposals": holdout_proposals,
            "pooled_proposals": pooled_proposals,
            "training_inside_analysis": radiative_guards._metrics(
                training_inside, training_proposals
            ),
            "holdout_inside_analysis": radiative_guards._metrics(
                holdout_inside, holdout_proposals
            ),
            "pooled_inside_analysis": radiative_guards._metrics(
                pooled_inside, pooled_proposals
            ),
            "model_fraction_sum": fraction_sum,
        },
        "support_counts": support_counts,
        "parent_coverage_counts": coverage_counts,
        "zero_survey_contribution_does_not_prove_structural_emptiness": True,
        "parent_coverage_failure_does_not_make_a_stratum_irrelevant": True,
        "detector_feed_in_is_not_measured_by_this_artifact": True,
        "strata": survey_records,
    }

    output.mkdir(parents=True)
    survey_path = output / "survey_model_evidence.json"
    _write_json(survey_path, survey_payload)
    survey_sha256 = _sha256(survey_path)

    survey_lookup = {
        str(record["stratum_id"]): record for record in survey_records
    }
    augmented_records = []
    survey_rationale = (
        "Pooled independent fixed-trial radiative surveys provide the model "
        "cross-section fraction; zero survey contribution does not establish "
        "structural emptiness, and parent-coverage failure does not remove "
        "analysis relevance."
    )
    for stratum in catalog:
        base = dict(base_records[stratum.identifier])
        survey = survey_lookup[stratum.identifier]
        rationale = str(base.get("rationale") or "").strip()
        base["rationale"] = (
            f"{rationale} {survey_rationale}".strip()
        )
        base["source_ids"] = list(base.get("source_ids", [])) + [
            source_identifier
        ]
        base["model_cross_section_fraction"] = survey[
            "model_cross_section_fraction"
        ]
        base["survey_model_evidence"] = {
            "support_status": survey["support_status"],
            "parent_coverage_status": survey["parent_coverage_status"],
            "pooled_cross_section_microbarn": survey["pooled"][
                "cross_section_microbarn"
            ],
            "pooled_sem_microbarn": survey["pooled"][
                "cross_section_sem_microbarn"
            ],
            "pooled_ess": survey["pooled"][
                "importance_effective_sample_size"
            ],
            "training_holdout_difference_z_score": survey[
                "training_holdout_difference_z_score"
            ],
        }
        augmented_records.append(base)

    instructions = dict(base_payload.get("instructions") or {})
    instructions.update(
        {
            "zero_survey_contribution_does_not_prove_structural_emptiness": True,
            "parent_coverage_failure_does_not_make_a_stratum_irrelevant": True,
            "detector_feed_in_requires_separate_GEMC_evidence": True,
        }
    )
    relevance_payload = {
        **base_payload,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "instructions": instructions,
        "sources": base_sources
        + [
            {
                "identifier": source_identifier,
                "path": str(survey_path),
                "sha256": survey_sha256,
            }
        ],
        "augmentation": {
            "base_relevance": str(base_path),
            "base_relevance_sha256": _sha256(base_path),
            "survey_model_evidence": str(survey_path),
            "survey_model_evidence_sha256": survey_sha256,
        },
        "strata": augmented_records,
    }
    relevance_path = output / "relevance_evidence.json"
    _write_json(relevance_path, relevance_payload)

    fields = (
        "flat_index",
        "stratum_id",
        "data_events",
        "training_status",
        "training_cross_section_microbarn",
        "training_sem_microbarn",
        "training_ess",
        "holdout_cross_section_microbarn",
        "holdout_sem_microbarn",
        "holdout_ess",
        "pooled_cross_section_microbarn",
        "pooled_sem_microbarn",
        "pooled_ess",
        "model_cross_section_fraction",
        "support_status",
        "parent_coverage_status",
        "training_holdout_difference_z_score",
    )
    with (output / "survey_model_evidence.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        writer = csv.DictWriter(destination, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for record in survey_records:
            training = record["training"]
            holdout = record["holdout"]
            pooled = record["pooled"]
            writer.writerow(
                {
                    "flat_index": record["flat_index"],
                    "stratum_id": record["stratum_id"],
                    "data_events": record["data_events"],
                    "training_status": record["training_status"],
                    "training_cross_section_microbarn": training[
                        "cross_section_microbarn"
                    ],
                    "training_sem_microbarn": training[
                        "cross_section_sem_microbarn"
                    ],
                    "training_ess": training[
                        "importance_effective_sample_size"
                    ],
                    "holdout_cross_section_microbarn": holdout[
                        "cross_section_microbarn"
                    ],
                    "holdout_sem_microbarn": holdout[
                        "cross_section_sem_microbarn"
                    ],
                    "holdout_ess": holdout[
                        "importance_effective_sample_size"
                    ],
                    "pooled_cross_section_microbarn": pooled[
                        "cross_section_microbarn"
                    ],
                    "pooled_sem_microbarn": pooled[
                        "cross_section_sem_microbarn"
                    ],
                    "pooled_ess": pooled[
                        "importance_effective_sample_size"
                    ],
                    "model_cross_section_fraction": record[
                        "model_cross_section_fraction"
                    ],
                    "support_status": record["support_status"],
                    "parent_coverage_status": record[
                        "parent_coverage_status"
                    ],
                    "training_holdout_difference_z_score": record[
                        "training_holdout_difference_z_score"
                    ],
                }
            )

    list_specs = {
        "survey_nonzero_flat_indices.txt": lambda record: float(
            record["model_cross_section_fraction"]
        )
        > 0.0,
        "independent_survey_support_flat_indices.txt": lambda record: record[
            "support_status"
        ]
        == "independent_support",
        "parent_coverage_failed_flat_indices.txt": lambda record: record[
            "parent_coverage_status"
        ]
        == "failed",
    }
    for filename, predicate in list_specs.items():
        values = [
            str(record["flat_index"])
            for record in survey_records
            if predicate(record)
        ]
        (output / filename).write_text(
            "".join(f"{value}\n" for value in values), encoding="utf-8"
        )
    return relevance_path


def _queue_work_category(record: dict[str, object]) -> str:
    survey = record["survey_model_evidence"]
    support = str(survey["support_status"])
    coverage = str(survey["parent_coverage_status"])
    if support == "no_survey_contribution":
        return "targeted_discovery"
    if support == "independent_support" and coverage == "passed":
        return "supported_calibration"
    return "guard_refinement"


def build_cumulative_queue(args: argparse.Namespace) -> Path:
    """Select data bins plus the model tail needed for global closure."""

    config_path = args.config.expanduser().resolve()
    config, config_sha256 = radiative_guards.load_analysis_config(config_path)
    relevance_path = args.relevance.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"{output} already exists")
    if args.minimum_data_events < 1:
        raise ValueError("--minimum-data-events must be positive")
    if not 0.0 <= args.maximum_global_model_residual_fraction < 1.0:
        raise ValueError(
            "--maximum-global-model-residual-fraction must lie in [0,1)"
        )

    relevance, relevance_sources = _load_relevance(
        relevance_path, config_sha256
    )
    catalog = radiative_guards.enumerate_strata(config)
    catalog_ids = {stratum.identifier for stratum in catalog}
    if set(relevance) != catalog_ids:
        missing = len(catalog_ids - set(relevance))
        extra = len(set(relevance) - catalog_ids)
        raise ActiveStratumError(
            "cumulative selection requires full-catalog relevance evidence "
            f"(missing={missing}, extra={extra})"
        )

    records: list[dict[str, object]] = []
    model_fraction_sum = 0.0
    allowed_support = {
        "independent_support",
        "training_only",
        "holdout_only",
        "no_survey_contribution",
    }
    allowed_coverage = {
        "passed",
        "failed",
        "no_holdout_contribution",
        "not_assessed",
    }
    for stratum in catalog:
        evidence = relevance[stratum.identifier]
        _validate_catalog_metadata(
            stratum,
            evidence,
            label="relevance evidence",
            required=True,
        )
        data_events = evidence.get("data_events")
        model_fraction = evidence.get("model_cross_section_fraction")
        survey = evidence.get("survey_model_evidence")
        if data_events is None or model_fraction is None:
            raise ActiveStratumError(
                f"{stratum.identifier}: cumulative selection requires "
                "data_events and model_cross_section_fraction"
            )
        if not isinstance(survey, dict):
            raise ActiveStratumError(
                f"{stratum.identifier}: cumulative selection requires "
                "survey_model_evidence"
            )
        support = str(survey.get("support_status", ""))
        coverage = str(survey.get("parent_coverage_status", ""))
        if support not in allowed_support or coverage not in allowed_coverage:
            raise ActiveStratumError(
                f"{stratum.identifier}: malformed survey support metadata"
            )
        fraction = float(model_fraction)
        model_fraction_sum += fraction
        records.append(
            {
                "stratum_id": stratum.identifier,
                "flat_index": stratum.flat_index,
                "indices": _indices(stratum),
                "bounds": _bounds(stratum),
                "data_events": int(data_events),
                "model_cross_section_fraction": fraction,
                "survey_model_evidence": dict(survey),
            }
        )
    if not math.isclose(
        model_fraction_sum, 1.0, rel_tol=5.0e-12, abs_tol=1.0e-12
    ):
        raise ActiveStratumError(
            "full-catalog model_cross_section_fraction values do not sum "
            f"to one (sum={model_fraction_sum})"
        )

    data_selected = [
        record
        for record in records
        if int(record["data_events"]) >= args.minimum_data_events
    ]
    zero_data_ranked = sorted(
        (
            record
            for record in records
            if int(record["data_events"]) < args.minimum_data_events
        ),
        key=lambda record: (
            -float(record["model_cross_section_fraction"]),
            int(record["flat_index"]),
        ),
    )
    initial_residual = math.fsum(
        float(record["model_cross_section_fraction"])
        for record in zero_data_ranked
    )
    residual = initial_residual
    model_selected: list[dict[str, object]] = []
    model_rank: dict[int, int] = {}
    cumulative_selected = 0.0
    cumulative_by_flat_index: dict[int, float] = {}
    tolerance = 5.0e-15
    for rank, record in enumerate(zero_data_ranked, start=1):
        flat_index = int(record["flat_index"])
        model_rank[flat_index] = rank
        if residual <= args.maximum_global_model_residual_fraction + tolerance:
            continue
        model_selected.append(record)
        fraction = float(record["model_cross_section_fraction"])
        cumulative_selected = math.fsum((cumulative_selected, fraction))
        cumulative_by_flat_index[flat_index] = cumulative_selected
        residual = math.fsum(
            float(item["model_cross_section_fraction"])
            for item in zero_data_ranked[rank:]
        )

    data_selected = sorted(
        data_selected,
        key=lambda record: (
            -int(record["data_events"]),
            -float(record["model_cross_section_fraction"]),
            int(record["flat_index"]),
        ),
    )
    selected = data_selected + model_selected
    selected_flat_indices = {
        int(record["flat_index"]) for record in selected
    }
    omitted = [
        record
        for record in zero_data_ranked
        if int(record["flat_index"]) not in selected_flat_indices
    ]
    actual_residual = math.fsum(
        float(record["model_cross_section_fraction"])
        for record in omitted
    )
    if actual_residual > (
        args.maximum_global_model_residual_fraction + tolerance
    ):
        raise ActiveStratumError(
            "cumulative selection failed to meet the requested residual"
        )

    queue_records: list[dict[str, object]] = []
    priority_rank = {
        int(record["flat_index"]): rank
        for rank, record in enumerate(selected, start=1)
    }
    model_selected_indices = {
        int(record["flat_index"]) for record in model_selected
    }
    for record in records:
        flat_index = int(record["flat_index"])
        is_selected = flat_index in selected_flat_indices
        if int(record["data_events"]) >= args.minimum_data_events:
            basis = "data_occupancy"
        elif flat_index in model_selected_indices:
            basis = "cumulative_model_tail"
        else:
            basis = "omitted_global_residual"
        queue_records.append(
            {
                **record,
                "selected": is_selected,
                "selection_basis": basis,
                "calibration_priority_rank": priority_rank.get(flat_index),
                "zero_data_model_rank": model_rank.get(flat_index),
                "zero_data_cumulative_selected_model_fraction": (
                    cumulative_by_flat_index.get(flat_index)
                ),
                "work_category": (
                    _queue_work_category(record) if is_selected else "omitted"
                ),
            }
        )

    category_counts = {
        category: sum(
            record["work_category"] == category for record in queue_records
        )
        for category in (
            "supported_calibration",
            "guard_refinement",
            "targeted_discovery",
            "omitted",
        )
    }
    selected_model_fraction = math.fsum(
        float(record["model_cross_section_fraction"])
        for record in selected
    )
    payload = {
        "schema": CUMULATIVE_QUEUE_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "builder_revision": radiative_mode4._source_revision(),
        "builder_source_sha256": _sha256(Path(__file__).resolve()),
        "analysis_config": str(config_path),
        "analysis_config_sha256": config_sha256,
        "relevance_evidence": str(relevance_path),
        "relevance_evidence_sha256": _sha256(relevance_path),
        "evidence_sources": relevance_sources,
        "selection_policy": {
            "minimum_data_events": args.minimum_data_events,
            "maximum_global_model_residual_fraction": (
                args.maximum_global_model_residual_fraction
            ),
            "data_occupied_strata_are_always_selected": True,
            "zero_data_order": (
                "descending pooled model cross-section fraction, then "
                "ascending flat index"
            ),
            "detector_feed_in_not_yet_included": True,
            "omitted_does_not_mean_structurally_empty": True,
        },
        "summary": {
            "catalog_strata": len(records),
            "data_occupied_strata": len(data_selected),
            "model_required_zero_data_strata": len(model_selected),
            "selected_strata": len(selected),
            "omitted_strata": len(omitted),
            "data_occupied_model_fraction": math.fsum(
                float(record["model_cross_section_fraction"])
                for record in data_selected
            ),
            "initial_zero_data_model_fraction": initial_residual,
            "selected_zero_data_model_fraction": math.fsum(
                float(record["model_cross_section_fraction"])
                for record in model_selected
            ),
            "selected_total_model_fraction": selected_model_fraction,
            "actual_global_model_residual_fraction": actual_residual,
            "work_category_counts": category_counts,
        },
        "strata": queue_records,
    }

    output.mkdir(parents=True)
    json_path = output / "cumulative_stratum_queue.json"
    _write_json(json_path, payload)
    fields = (
        "calibration_priority_rank",
        "flat_index",
        "stratum_id",
        "selected",
        "selection_basis",
        "work_category",
        "data_events",
        "model_cross_section_fraction",
        "zero_data_model_rank",
        "zero_data_cumulative_selected_model_fraction",
        "support_status",
        "parent_coverage_status",
        "pooled_cross_section_microbarn",
        "pooled_sem_microbarn",
        "pooled_ess",
        "Q2",
        "xB",
        "minus_t",
        "phi_deg",
    )
    with (output / "cumulative_stratum_queue.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        writer = csv.DictWriter(destination, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        ordered_records = sorted(
            queue_records,
            key=lambda record: (
                not bool(record["selected"]),
                record["calibration_priority_rank"]
                if record["calibration_priority_rank"] is not None
                else int(record["zero_data_model_rank"] or len(records) + 1),
                int(record["flat_index"]),
            ),
        )
        for record in ordered_records:
            survey = record["survey_model_evidence"]
            bounds = record["bounds"]
            writer.writerow(
                {
                    "calibration_priority_rank": record[
                        "calibration_priority_rank"
                    ],
                    "flat_index": record["flat_index"],
                    "stratum_id": record["stratum_id"],
                    "selected": record["selected"],
                    "selection_basis": record["selection_basis"],
                    "work_category": record["work_category"],
                    "data_events": record["data_events"],
                    "model_cross_section_fraction": record[
                        "model_cross_section_fraction"
                    ],
                    "zero_data_model_rank": record["zero_data_model_rank"],
                    "zero_data_cumulative_selected_model_fraction": record[
                        "zero_data_cumulative_selected_model_fraction"
                    ],
                    "support_status": survey["support_status"],
                    "parent_coverage_status": survey[
                        "parent_coverage_status"
                    ],
                    "pooled_cross_section_microbarn": survey.get(
                        "pooled_cross_section_microbarn"
                    ),
                    "pooled_sem_microbarn": survey.get(
                        "pooled_sem_microbarn"
                    ),
                    "pooled_ess": survey.get("pooled_ess"),
                    "Q2": ":".join(str(value) for value in bounds["Q2"]),
                    "xB": ":".join(str(value) for value in bounds["xB"]),
                    "minus_t": ":".join(
                        str(value) for value in bounds["minus_t"]
                    ),
                    "phi_deg": ":".join(
                        str(value) for value in bounds["phi_deg"]
                    ),
                }
            )

    list_specs = {
        "selected_flat_indices.txt": selected,
        "data_occupied_flat_indices.txt": data_selected,
        "model_required_zero_data_flat_indices.txt": model_selected,
        "omitted_flat_indices.txt": omitted,
    }
    for category in (
        "supported_calibration",
        "guard_refinement",
        "targeted_discovery",
    ):
        list_specs[f"{category}_flat_indices.txt"] = [
            record
            for record in selected
            if _queue_work_category(record) == category
        ]
    for filename, selected_records in list_specs.items():
        (output / filename).write_text(
            "".join(
                f"{int(record['flat_index'])}\n"
                for record in selected_records
            ),
            encoding="utf-8",
        )
    return json_path


def _normalized_stratum_distance_squared(
    first: tuple[int, int, int, int],
    second: tuple[int, int, int, int],
    axis_bins: tuple[int, int, int, int],
) -> float:
    """Distance in bin-index space, with periodic phi."""

    _, first_xb, first_t, first_phi = first
    _, second_xb, second_t, second_phi = second
    _, xb_bins, t_bins, phi_bins = axis_bins
    xb_scale = max(1, xb_bins - 1)
    t_scale = max(1, t_bins - 1)
    phi_scale = max(1.0, phi_bins / 2.0)
    phi_difference = abs(first_phi - second_phi)
    phi_difference = min(phi_difference, phi_bins - phi_difference)
    return (
        ((first_xb - second_xb) / xb_scale) ** 2
        + ((first_t - second_t) / t_scale) ** 2
        + (phi_difference / phi_scale) ** 2
    )


def select_stratified_batch(args: argparse.Namespace) -> Path:
    """Select a deterministic maximin calibration-scale batch."""

    queue_path = args.queue.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"{output} already exists")
    if args.representatives_per_q2_basis < 1:
        raise ValueError("--representatives-per-q2-basis must be positive")
    if not queue_path.is_file():
        raise FileNotFoundError(queue_path)
    queue_sha256 = _sha256(queue_path)
    queue = json.loads(queue_path.read_text(encoding="utf-8"))
    if queue.get("schema") != CUMULATIVE_QUEUE_SCHEMA:
        raise ActiveStratumError(
            f"{queue_path}: expected schema {CUMULATIVE_QUEUE_SCHEMA}"
        )

    config_path = Path(str(queue.get("analysis_config", ""))).expanduser().resolve()
    if not config_path.is_file():
        raise FileNotFoundError(config_path)
    config, config_sha256 = radiative_guards.load_analysis_config(config_path)
    if queue.get("analysis_config_sha256") != config_sha256:
        raise ActiveStratumError(
            "cumulative queue analysis configuration hash differs"
        )
    catalog = radiative_guards.enumerate_strata(config)
    catalog_lookup = {stratum.identifier: stratum for stratum in catalog}
    raw_records = queue.get("strata")
    if not isinstance(raw_records, list) or len(raw_records) != len(catalog):
        raise ActiveStratumError(
            "cumulative queue does not contain the full analysis catalog"
        )
    records: list[dict[str, object]] = []
    identifiers: set[str] = set()
    for raw in raw_records:
        if not isinstance(raw, dict):
            raise ActiveStratumError("cumulative queue has a malformed record")
        record = dict(raw)
        identifier = str(record.get("stratum_id", ""))
        if identifier not in catalog_lookup or identifier in identifiers:
            raise ActiveStratumError(
                "cumulative queue has an unknown or duplicate stratum"
            )
        identifiers.add(identifier)
        _validate_catalog_metadata(
            catalog_lookup[identifier],
            record,
            label="cumulative queue",
            required=True,
        )
        if record.get("selected") not in (True, False):
            raise ActiveStratumError(
                f"{identifier}: cumulative queue selected must be boolean"
            )
        fraction = record.get("model_cross_section_fraction")
        data_events = record.get("data_events")
        if fraction is None or data_events is None:
            raise ActiveStratumError(
                f"{identifier}: cumulative queue lacks ranking fields"
            )
        fraction = float(fraction)
        events = float(data_events)
        if (
            not math.isfinite(fraction)
            or fraction < 0.0
            or not math.isfinite(events)
            or events < 0.0
            or not events.is_integer()
        ):
            raise ActiveStratumError(
                f"{identifier}: cumulative queue has invalid ranking fields"
            )
        record["model_cross_section_fraction"] = fraction
        record["data_events"] = int(events)
        records.append(record)
    if identifiers != set(catalog_lookup):
        raise ActiveStratumError(
            "cumulative queue stratum identifiers differ from the catalog"
        )
    recorded_selected = int((queue.get("summary") or {}).get("selected_strata", -1))
    actual_selected = sum(record["selected"] is True for record in records)
    if recorded_selected != actual_selected:
        raise ActiveStratumError(
            "cumulative queue selected count differs from its summary"
        )

    default_bases = ("data_occupancy", "cumulative_model_tail")
    bases = tuple(args.selection_bases or default_bases)
    if len(set(bases)) != len(bases):
        raise ValueError("duplicate --selection-basis")
    unknown_bases = sorted(set(bases) - set(default_bases))
    if unknown_bases:
        raise ValueError(
            "unsupported --selection-basis values: "
            + ", ".join(unknown_bases)
        )
    axis_bins = tuple(
        len(config["binning"][name]) - 1
        for name in ("Q2", "xB", "minus_t", "phi_deg")
    )
    q2_bins = axis_bins[0]
    eligible = [
        record
        for record in records
        if record["selected"] is True
        and record.get("work_category") == args.work_category
        and record.get("selection_basis") in bases
    ]
    groups: dict[tuple[str, int], list[dict[str, object]]] = {
        (basis, iq2): [] for basis in bases for iq2 in range(q2_bins)
    }
    normalized_indices: dict[int, tuple[int, int, int, int]] = {}
    for record in eligible:
        indices = _normalized_indices(record["indices"])
        flat_index = int(record["flat_index"])
        normalized_indices[flat_index] = indices
        key = (str(record["selection_basis"]), indices[0])
        groups[key].append(record)

    allow_basis_fallback = bool(
        getattr(args, "allow_basis_fallback", False)
    )

    def choose_maximin(
        candidates: list[dict[str, object]],
        count: int,
        context: list[dict[str, object]],
        *,
        anchor_role: str,
        additional_role: str,
        basis_fallback: bool,
    ) -> list[dict[str, object]]:
        remaining = list(candidates)
        local_selected: list[dict[str, object]] = []
        details: list[dict[str, object]] = []
        while len(local_selected) < count:
            references = context + local_selected
            if not references:
                next_record = max(
                    remaining,
                    key=lambda record: (
                        float(record["model_cross_section_fraction"]),
                        int(record["data_events"]),
                        -int(record["flat_index"]),
                    ),
                )
                role = anchor_role
                distance = None
            else:
                distance_by_flat_index = {
                    int(candidate["flat_index"]): min(
                        _normalized_stratum_distance_squared(
                            normalized_indices[int(candidate["flat_index"])],
                            normalized_indices[int(reference["flat_index"])],
                            axis_bins,
                        )
                        for reference in references
                    )
                    for candidate in remaining
                }
                next_record = max(
                    remaining,
                    key=lambda record: (
                        distance_by_flat_index[int(record["flat_index"])],
                        float(record["model_cross_section_fraction"]),
                        int(record["data_events"]),
                        -int(record["flat_index"]),
                    ),
                )
                role = additional_role
                distance = distance_by_flat_index[int(next_record["flat_index"])]
            remaining.remove(next_record)
            local_selected.append(next_record)
            details.append(
                {
                    "record": next_record,
                    "selection_role": role,
                    "minimum_normalized_distance_squared": distance,
                    "basis_fallback": basis_fallback,
                }
            )
        return details

    chosen: list[dict[str, object]] = []
    group_summaries: list[dict[str, object]] = []
    required_per_basis = args.representatives_per_q2_basis
    required_per_q2 = required_per_basis * len(bases)
    for iq2 in range(q2_bins):
        q2_details: list[dict[str, object]] = []
        primary_details: dict[str, list[dict[str, object]]] = {}
        for basis in bases:
            candidates = groups[(basis, iq2)]
            if len(candidates) < required_per_basis and not allow_basis_fallback:
                raise ActiveStratumError(
                    f"{basis}, Q2 index {iq2}: only {len(candidates)} "
                    f"eligible {args.work_category} strata; "
                    f"{required_per_basis} required"
                )
            count = min(len(candidates), required_per_basis)
            details = choose_maximin(
                candidates,
                count,
                [],
                anchor_role="largest_model_contribution_anchor",
                additional_role="normalized_index_maximin",
                basis_fallback=False,
            )
            primary_details[basis] = details
            q2_details.extend(details)

        fallback_needed = required_per_q2 - len(q2_details)
        if fallback_needed:
            selected_indices = {
                int(detail["record"]["flat_index"]) for detail in q2_details
            }
            fallback_candidates = [
                record
                for basis in bases
                for record in groups[(basis, iq2)]
                if int(record["flat_index"]) not in selected_indices
            ]
            if len(fallback_candidates) < fallback_needed:
                raise ActiveStratumError(
                    f"Q2 index {iq2}: basis fallback needs {fallback_needed} "
                    f"additional {args.work_category} strata but only "
                    f"{len(fallback_candidates)} remain"
                )
            q2_details.extend(
                choose_maximin(
                    fallback_candidates,
                    fallback_needed,
                    [detail["record"] for detail in q2_details],
                    anchor_role="basis_fallback_model_anchor",
                    additional_role="basis_fallback_normalized_index_maximin",
                    basis_fallback=True,
                )
            )

        within_basis_counts = {basis: 0 for basis in bases}
        for detail in q2_details:
            source = detail["record"]
            actual_basis = str(source["selection_basis"])
            within_basis_counts[actual_basis] += 1
            chosen.append(
                {
                    **source,
                    "batch_rank": len(chosen) + 1,
                    "within_group_rank": within_basis_counts[actual_basis],
                    "selection_role": detail["selection_role"],
                    "minimum_normalized_distance_squared": detail[
                        "minimum_normalized_distance_squared"
                    ],
                    "basis_fallback": detail["basis_fallback"],
                }
            )

        for basis in bases:
            primary = primary_details[basis]
            fallback = [
                detail
                for detail in q2_details
                if detail["basis_fallback"]
                and detail["record"]["selection_basis"] == basis
            ]
            group_summaries.append(
                {
                    "selection_basis": basis,
                    "q2_index": iq2,
                    "q2_bounds": list(config["binning"]["Q2"][iq2 : iq2 + 2]),
                    "eligible_strata": len(groups[(basis, iq2)]),
                    "requested_representatives": required_per_basis,
                    "unfilled_primary_quota": max(
                        0, required_per_basis - len(primary)
                    ),
                    "primary_selected_flat_indices": [
                        int(detail["record"]["flat_index"])
                        for detail in primary
                    ],
                    "basis_fallback_selected_flat_indices": [
                        int(detail["record"]["flat_index"])
                        for detail in fallback
                    ],
                    "selected_flat_indices": [
                        int(detail["record"]["flat_index"])
                        for detail in primary + fallback
                    ],
                }
            )

    if len({int(record["flat_index"]) for record in chosen}) != len(chosen):
        raise ActiveStratumError("stratified selection produced duplicates")
    basis_counts = {
        basis: sum(record["selection_basis"] == basis for record in chosen)
        for basis in bases
    }
    q2_counts = {
        str(iq2): sum(
            _normalized_indices(record["indices"])[0] == iq2
            for record in chosen
        )
        for iq2 in range(q2_bins)
    }
    payload = {
        "schema": STRATIFIED_BATCH_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "builder_revision": radiative_mode4._source_revision(),
        "builder_source_sha256": _sha256(Path(__file__).resolve()),
        "analysis_config": str(config_path),
        "analysis_config_sha256": config_sha256,
        "cumulative_queue": str(queue_path),
        "cumulative_queue_sha256": queue_sha256,
        "cumulative_queue_summary": queue.get("summary"),
        "selection_policy": {
            "work_category": args.work_category,
            "selection_bases": list(bases),
            "representatives_per_q2_basis": (
                args.representatives_per_q2_basis
            ),
            "allow_basis_fallback": allow_basis_fallback,
            "basis_fallback": (
                "preserve the requested total per Q2 by maximin selection "
                "from another requested basis in the same Q2 interval"
                if allow_basis_fallback
                else None
            ),
            "anchor": (
                "largest model fraction, then data events, then lowest "
                "flat index"
            ),
            "additional_representatives": (
                "maximin squared distance in normalized xB/-t/periodic-phi "
                "bin-index space; ties use model fraction, data events, "
                "then lowest flat index"
            ),
            "no_random_selection": True,
        },
        "summary": {
            "selected_strata": len(chosen),
            "basis_fallback_strata": sum(
                bool(record["basis_fallback"]) for record in chosen
            ),
            "basis_counts": basis_counts,
            "q2_index_counts": q2_counts,
            "selected_model_cross_section_fraction": math.fsum(
                float(record["model_cross_section_fraction"])
                for record in chosen
            ),
        },
        "groups": group_summaries,
        "strata": chosen,
    }

    output.mkdir(parents=True)
    json_path = output / "stratified_batch.json"
    _write_json(json_path, payload)
    fields = (
        "batch_rank",
        "flat_index",
        "stratum_id",
        "selection_basis",
        "work_category",
        "within_group_rank",
        "selection_role",
        "basis_fallback",
        "minimum_normalized_distance_squared",
        "data_events",
        "model_cross_section_fraction",
        "pooled_ess",
        "Q2",
        "xB",
        "minus_t",
        "phi_deg",
    )
    with (output / "stratified_batch.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        writer = csv.DictWriter(destination, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for record in chosen:
            bounds = record["bounds"]
            survey = record.get("survey_model_evidence") or {}
            writer.writerow(
                {
                    "batch_rank": record["batch_rank"],
                    "flat_index": record["flat_index"],
                    "stratum_id": record["stratum_id"],
                    "selection_basis": record["selection_basis"],
                    "work_category": record["work_category"],
                    "within_group_rank": record["within_group_rank"],
                    "selection_role": record["selection_role"],
                    "basis_fallback": record["basis_fallback"],
                    "minimum_normalized_distance_squared": record[
                        "minimum_normalized_distance_squared"
                    ],
                    "data_events": record["data_events"],
                    "model_cross_section_fraction": record[
                        "model_cross_section_fraction"
                    ],
                    "pooled_ess": survey.get("pooled_ess"),
                    "Q2": ":".join(str(value) for value in bounds["Q2"]),
                    "xB": ":".join(str(value) for value in bounds["xB"]),
                    "minus_t": ":".join(
                        str(value) for value in bounds["minus_t"]
                    ),
                    "phi_deg": ":".join(
                        str(value) for value in bounds["phi_deg"]
                    ),
                }
            )
    (output / "selected_flat_indices.txt").write_text(
        "".join(f"{int(record['flat_index'])}\n" for record in chosen),
        encoding="utf-8",
    )
    for basis in bases:
        (output / f"{basis}_flat_indices.txt").write_text(
            "".join(
                f"{int(record['flat_index'])}\n"
                for record in chosen
                if record["selection_basis"] == basis
            ),
            encoding="utf-8",
        )
    return json_path


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

    data = subparsers.add_parser(
        "build-data-evidence",
        help=(
            "build full-catalog reconstructed-data occupancy and compatible "
            "relevance evidence"
        ),
    )
    data.add_argument("--config", type=Path, required=True)
    data.add_argument("--data-events", type=Path, required=True)
    data.add_argument("--selection-mask", type=Path, required=True)
    data.add_argument("--output", type=Path, required=True)
    data.add_argument(
        "--allow-duplicate-event-keys",
        action="store_true",
        help=(
            "permit repeated (run,event) keys after selection; disabled by "
            "default because overlaps would double-count data"
        ),
    )

    survey = subparsers.add_parser(
        "augment-survey-evidence",
        help=(
            "pool migration training/holdout surveys and add full-catalog "
            "model relevance to existing evidence"
        ),
    )
    survey.add_argument("--config", type=Path, required=True)
    survey.add_argument("--base-relevance", type=Path, required=True)
    survey.add_argument("--migration-manifest", type=Path, required=True)
    survey.add_argument("--migration-validation", type=Path, required=True)
    survey.add_argument("--output", type=Path, required=True)

    queue = subparsers.add_parser(
        "build-cumulative-queue",
        help=(
            "select all data-occupied strata plus the ranked zero-data "
            "model tail required by a global residual budget"
        ),
    )
    queue.add_argument("--config", type=Path, required=True)
    queue.add_argument("--relevance", type=Path, required=True)
    queue.add_argument("--output", type=Path, required=True)
    queue.add_argument("--minimum-data-events", type=int, default=1)
    queue.add_argument(
        "--maximum-global-model-residual-fraction",
        type=float,
        required=True,
    )

    batch = subparsers.add_parser(
        "select-stratified-batch",
        help=(
            "select a deterministic maximin scale batch from a cumulative "
            "work queue"
        ),
    )
    batch.add_argument("--queue", type=Path, required=True)
    batch.add_argument("--output", type=Path, required=True)
    batch.add_argument(
        "--work-category", default="supported_calibration"
    )
    batch.add_argument(
        "--selection-basis", dest="selection_bases", action="append"
    )
    batch.add_argument(
        "--representatives-per-q2-basis", type=int, default=2
    )
    batch.add_argument(
        "--allow-basis-fallback",
        action="store_true",
        help=(
            "when one requested basis is sparse, preserve the total per Q2 "
            "with maximin representatives from another requested basis"
        ),
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
    if args.command == "create-template":
        result = create_template(args)
    elif args.command == "build-data-evidence":
        result = build_data_evidence(args)
    elif args.command == "augment-survey-evidence":
        result = augment_survey_evidence(args)
    elif args.command == "build-cumulative-queue":
        result = build_cumulative_queue(args)
    elif args.command == "select-stratified-batch":
        result = select_stratified_batch(args)
    else:
        result = classify(args)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
