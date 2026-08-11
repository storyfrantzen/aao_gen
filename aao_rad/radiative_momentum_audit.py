#!/usr/bin/env python3
"""Audit final-electron momentum retention in radiative analysis strata.

The audit consumes the same fixed-trial survey artifacts used to learn the
radiative guards.  It estimates, with survey cross-section contributions,
what fraction of each observed analysis stratum also satisfies a stricter
final-LUND-electron momentum threshold.  Raw survey-row fractions are retained
only as diagnostics; they are not used as physics estimates.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import radiative_guards
import radiative_mode4
import radiative_survey


AUDIT_SCHEMA = "aao-rad-electron-momentum-audit-v1"
ELECTRON_MASS_GEV = 0.00051099895
REPORT_FILENAME = "electron_momentum_audit.json"
TABLE_FILENAME = "electron_momentum_audit.tsv"


class MomentumAuditError(RuntimeError):
    """Raised when survey artifacts cannot define a trustworthy audit."""


@dataclass
class GroupAudit:
    proposals: int
    replicas: list[dict[str, object]]
    denominator: dict[str, radiative_guards.Moment] = field(
        default_factory=dict
    )
    nominal: dict[str, radiative_guards.Moment] = field(default_factory=dict)
    below_threshold: dict[str, radiative_guards.Moment] = field(
        default_factory=dict
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, payload: object) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _copy_moment(moment: radiative_guards.Moment | None) -> radiative_guards.Moment:
    if moment is None:
        return radiative_guards.Moment()
    return radiative_guards.Moment(
        count=moment.count,
        total=moment.total,
        square_total=moment.square_total,
        maximum=moment.maximum,
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


def _sum_moments(
    moments: Iterable[radiative_guards.Moment],
) -> radiative_guards.Moment:
    result = radiative_guards.Moment()
    for moment in moments:
        result = _add_moments(result, moment)
    return result


def _electron_momentum(row: dict[str, str]) -> float:
    try:
        components = tuple(
            float(row[name])
            for name in ("final_e_px", "final_e_py", "final_e_pz")
        )
    except (KeyError, TypeError, ValueError) as error:
        raise MomentumAuditError(
            "survey row has malformed final-electron momentum components"
        ) from error
    momentum_squared = sum(value * value for value in components)
    if not math.isfinite(momentum_squared) or momentum_squared < 0.0:
        raise MomentumAuditError(
            "survey row has nonfinite final-electron momentum"
        )
    return math.sqrt(momentum_squared)


def _close_moment(
    observed: radiative_guards.Moment,
    expected: radiative_guards.Moment,
    *,
    label: str,
) -> None:
    if observed.count != expected.count:
        raise MomentumAuditError(
            f"{label}: row count does not close ({observed.count} versus "
            f"{expected.count})"
        )
    for name in ("total", "square_total", "maximum"):
        left = float(getattr(observed, name))
        right = float(getattr(expected, name))
        if not math.isclose(left, right, rel_tol=5.0e-13, abs_tol=1.0e-15):
            raise MomentumAuditError(
                f"{label}: {name} does not close ({left} versus {right})"
            )


def _validate_threshold_support(
    directories: list[Path], electron_p_min: float
) -> None:
    threshold_energy = math.hypot(electron_p_min, ELECTRON_MASS_GEV)
    for directory in directories:
        norm = radiative_guards._validated_norm(directory.resolve())
        try:
            generator_energy_min = float(norm["ep_min"])
        except (KeyError, ValueError) as error:
            raise MomentumAuditError(
                f"{directory}: normalization metadata lacks ep_min"
            ) from error
        if threshold_energy + 1.0e-10 < generator_energy_min:
            raise MomentumAuditError(
                f"{directory}: requested p_e >= {electron_p_min:.12g} GeV "
                f"corresponds to E_e >= {threshold_energy:.12g} GeV, below "
                f"the survey generator cutoff {generator_energy_min:.12g} "
                "GeV; the missing support cannot be audited"
            )


def _validate_cross_group_compatibility(
    requested_groups: dict[str, list[Path]], config: dict
) -> None:
    reference_norm: dict[str, str] | None = None
    reference_signature: str | None = None
    replica_ids: set[int] = set()
    expected_binning = {
        name: [float(value) for value in config["binning"][name]]
        for name in ("Q2", "xB", "minus_t", "phi_deg")
    }
    for directories in requested_groups.values():
        for directory in directories:
            norm = radiative_guards._validated_norm(directory)
            signature, _, settings = radiative_guards._legacy_input_metadata(
                directory / "survey_input.inp"
            )
            if reference_norm is None:
                reference_norm = norm
                reference_signature = signature
            else:
                radiative_guards._compatible_settings(reference_norm, norm)
                if signature != reference_signature:
                    raise MomentumAuditError(
                        "training and validation surveys have different "
                        "physical or proposal settings"
                    )
                if norm.get("survey_schema") != reference_norm.get(
                    "survey_schema"
                ):
                    raise MomentumAuditError(
                        "training and validation survey schemas differ"
                    )
            replica = int(norm["survey_replica"])
            if replica in replica_ids:
                raise MomentumAuditError(f"duplicate survey replica ID {replica}")
            replica_ids.add(replica)
            proposal = settings.get("survey_proposal", {})
            if proposal.get("mode") == 1 and proposal.get("binning") != expected_binning:
                raise MomentumAuditError(
                    f"{directory}: balanced survey binning differs from the "
                    "analysis config"
                )


def _audit_group(
    directories: list[Path],
    config: dict,
    *,
    apply_y_max: bool,
    electron_p_min: float,
) -> GroupAudit:
    if not directories:
        raise MomentumAuditError("an audit group cannot be empty")
    _validate_threshold_support(directories, electron_p_min)

    # This performs the complete existing survey validation, compatibility,
    # fixed-trial closure, and observed-stratum assignment before the momentum
    # split is evaluated in a second streaming pass.
    validated = radiative_guards.aggregate_surveys(
        directories,
        config,
        radiative_guards.parse_partition(None),
        apply_y_max=apply_y_max,
    )
    nominal: dict[str, radiative_guards.Moment] = {}
    below: dict[str, radiative_guards.Moment] = {}
    tolerance = 8.0 * math.ulp(max(1.0, electron_p_min))

    for requested_directory in directories:
        directory = requested_directory.resolve()
        survey_path = directory / radiative_survey.SURVEY_FILENAME
        for row in radiative_guards._survey_rows(survey_path):
            stratum_id, _ = radiative_guards.assign_stratum(
                row, config, apply_y_max=apply_y_max
            )
            if stratum_id is None:
                continue
            try:
                weight = float(row["trial_xsec_observed_microbarn"])
            except (KeyError, TypeError, ValueError) as error:
                raise MomentumAuditError(
                    f"{survey_path}: malformed observed contribution"
                ) from error
            destination = (
                nominal
                if _electron_momentum(row) + tolerance >= electron_p_min
                else below
            )
            destination.setdefault(
                stratum_id, radiative_guards.Moment()
            ).add(weight)

    for stratum_id, denominator in validated.strata.items():
        split = _add_moments(
            nominal.get(stratum_id, radiative_guards.Moment()),
            below.get(stratum_id, radiative_guards.Moment()),
        )
        _close_moment(
            split,
            denominator,
            label=f"{stratum_id} nominal/below-threshold partition",
        )

    return GroupAudit(
        proposals=validated.proposals,
        replicas=validated.replicas,
        denominator={
            key: _copy_moment(value) for key, value in validated.strata.items()
        },
        nominal=nominal,
        below_threshold=below,
    )


def _combine_groups(groups: Iterable[GroupAudit]) -> GroupAudit:
    result = GroupAudit(proposals=0, replicas=[])
    for group in groups:
        result.proposals += group.proposals
        result.replicas.extend(group.replicas)
        for attribute in ("denominator", "nominal", "below_threshold"):
            destination = getattr(result, attribute)
            for stratum_id, moment in getattr(group, attribute).items():
                destination[stratum_id] = _add_moments(
                    destination.get(stratum_id, radiative_guards.Moment()),
                    moment,
                )
    result.replicas.sort(key=lambda item: int(item["replica"]))
    return result


def _fraction_metrics(
    nominal: radiative_guards.Moment,
    below: radiative_guards.Moment,
    proposals: int,
) -> dict[str, object]:
    denominator = _add_moments(nominal, below)
    total = denominator.total
    fraction = nominal.total / total if total > 0.0 else None
    fraction_sem: float | None = None
    if fraction is not None and proposals > 1:
        # Delta-method variance of sum(a_i)/sum(a_i+b_i).  A proposal can
        # contribute to only one side, so the pass/fail cross-product sum is
        # exactly zero.  The fixed-trial zero contributions are represented by
        # `proposals` even though they are absent from the survey CSV.
        residual_square_sum = (
            (1.0 - fraction) ** 2 * nominal.square_total
            + fraction**2 * below.square_total
        )
        fraction_sem = math.sqrt(
            proposals
            * residual_square_sum
            / ((proposals - 1) * total * total)
        )

    if total <= 0.0:
        status = "no_survey_contribution"
    elif nominal.total <= 0.0:
        status = "no_nominal_contribution_observed"
    elif below.total <= 0.0:
        status = "all_observed_contribution_passes"
    elif min(
        radiative_guards._ess(nominal), radiative_guards._ess(below)
    ) < 5.0:
        status = "mixed_low_effective_support"
    else:
        status = "mixed_quantified"

    return {
        "status": status,
        "denominator": radiative_guards._metrics(denominator, proposals),
        "nominal": radiative_guards._metrics(nominal, proposals),
        "below_threshold": radiative_guards._metrics(below, proposals),
        "nominal_cross_section_fraction": fraction,
        "nominal_fraction_sem_delta_method": fraction_sem,
        "raw_nominal_row_fraction_diagnostic_only": (
            nominal.count / denominator.count if denominator.count else None
        ),
        "loose_events_per_expected_nominal_event": (
            1.0 / fraction if fraction is not None and fraction > 0.0 else None
        ),
        "finite_survey_boundary": (
            total > 0.0 and (nominal.total <= 0.0 or below.total <= 0.0)
        ),
    }


def _difference_z_score(
    first: dict[str, object], second: dict[str, object]
) -> float | None:
    first_fraction = first["nominal_cross_section_fraction"]
    second_fraction = second["nominal_cross_section_fraction"]
    first_sem = first["nominal_fraction_sem_delta_method"]
    second_sem = second["nominal_fraction_sem_delta_method"]
    if None in (first_fraction, second_fraction, first_sem, second_sem):
        return None
    uncertainty = math.hypot(float(first_sem), float(second_sem))
    if uncertainty <= 0.0:
        return None
    return (float(second_fraction) - float(first_fraction)) / uncertainty


def _stratum_metadata(stratum: radiative_guards.Stratum) -> dict[str, object]:
    return {
        "stratum_id": stratum.identifier,
        "flat_index": stratum.flat_index,
        "indices": {
            "Q2": stratum.iq2,
            "xB": stratum.ixb,
            "minus_t": stratum.it,
            "phi_deg": stratum.iphi,
        },
        "bounds": {
            "Q2": list(stratum.q2),
            "xB": list(stratum.xb),
            "minus_t": list(stratum.minus_t),
            "phi_deg": list(stratum.phi_deg),
        },
    }


def _load_selection(
    path: Path | None, catalog: list[radiative_guards.Stratum]
) -> tuple[set[int], dict[str, object] | None]:
    catalog_indices = {stratum.flat_index for stratum in catalog}
    if path is None:
        return catalog_indices, None
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    try:
        values = [int(token) for token in resolved.read_text().split()]
    except ValueError as error:
        raise MomentumAuditError(
            f"{resolved}: selection must contain integer flat indices"
        ) from error
    if len(values) != len(set(values)):
        raise MomentumAuditError(f"{resolved}: duplicate flat indices")
    unknown = sorted(set(values) - catalog_indices)
    if unknown:
        raise MomentumAuditError(
            f"{resolved}: unknown flat indices, beginning with {unknown[:5]}"
        )
    return set(values), {
        "path": str(resolved),
        "sha256": _sha256(resolved),
        "selected_strata": len(values),
    }


def _group_summary(
    group: GroupAudit,
    identifiers: set[str],
) -> dict[str, object]:
    denominator = _sum_moments(
        group.denominator.get(identifier, radiative_guards.Moment())
        for identifier in identifiers
    )
    nominal = _sum_moments(
        group.nominal.get(identifier, radiative_guards.Moment())
        for identifier in identifiers
    )
    below = _sum_moments(
        group.below_threshold.get(identifier, radiative_guards.Moment())
        for identifier in identifiers
    )
    return _fraction_metrics(nominal, below, group.proposals)


def _group_sources(
    directories: list[Path], group: GroupAudit
) -> dict[str, object]:
    by_replica = {int(item["replica"]): item for item in group.replicas}
    sources: list[dict[str, object]] = []
    for requested in directories:
        directory = requested.resolve()
        norm = radiative_guards._validated_norm(directory)
        replica = int(norm["survey_replica"])
        sources.append(
            {
                "directory": str(directory),
                "replica": replica,
                "proposals": int(norm["ntries"]),
                "survey_csv_sha256": _sha256(
                    directory / radiative_survey.SURVEY_FILENAME
                ),
                "normalization_sha256": _sha256(
                    directory / radiative_survey.NORM_FILENAME
                ),
                "survey_input_sha256": _sha256(
                    directory / "survey_input.inp"
                ),
                "validated_summary": by_replica[replica],
            }
        )
    return {"proposals": group.proposals, "surveys": sources}


def _analysis_selection(
    config: dict, *, apply_y_max: bool, electron_p_min: float
) -> dict[str, object]:
    phase_space = config["phase_space"]
    result: dict[str, object] = {
        "coordinate_definition": "final_lund_analysis",
        "electron_momentum_definition": (
            "sqrt(final_e_px^2 + final_e_py^2 + final_e_pz^2)"
        ),
        "q2_minimum": float(
            phase_space.get("Q2_min", config["binning"]["Q2"][0])
        ),
        "w_minimum": float(phase_space["W_min"]),
        "apply_y_max_to_denominator": apply_y_max,
        "denominator_y_maximum": (
            float(phase_space["y_max"]) if apply_y_max else None
        ),
        "no_implicit_y_minimum": True,
        "nominal_electron_momentum_minimum_GeV": electron_p_min,
        "threshold_edge_convention": "lower_inclusive",
    }
    return result


def audit(args: argparse.Namespace) -> Path:
    config_path = args.config.expanduser().resolve()
    config, config_sha256 = radiative_guards.load_analysis_config(config_path)
    phase_space = config["phase_space"]
    configured_threshold = phase_space.get("electron_p_min")
    if args.electron_p_min is None:
        if configured_threshold is None:
            raise MomentumAuditError(
                "phase_space.electron_p_min is absent; provide "
                "--electron-p-min explicitly"
            )
        electron_p_min = float(configured_threshold)
        threshold_source = "analysis_config"
    else:
        electron_p_min = float(args.electron_p_min)
        threshold_source = "command_line"
    if not math.isfinite(electron_p_min) or electron_p_min < 0.0:
        raise MomentumAuditError("electron momentum threshold must be finite and nonnegative")

    if args.apply_y_max and phase_space.get("y_max") is None:
        raise MomentumAuditError(
            "--apply-y-max requires a finite phase_space.y_max"
        )
    if args.apply_y_max and not math.isfinite(float(phase_space["y_max"])):
        raise MomentumAuditError("phase_space.y_max must be finite")

    if args.surveys:
        if args.training_surveys or args.validation_surveys:
            raise MomentumAuditError(
                "use either --survey or the training/validation survey split"
            )
        requested_groups = {"pooled": args.surveys}
        split_groups: dict[str, GroupAudit] = {}
    else:
        if not args.training_surveys or not args.validation_surveys:
            raise MomentumAuditError(
                "provide --survey, or provide both --training-survey and "
                "--validation-survey"
            )
        requested_groups = {
            "training": args.training_surveys,
            "validation": args.validation_surveys,
        }
        split_groups = {}

    resolved_seen: set[Path] = set()
    for label, directories in requested_groups.items():
        resolved = [path.expanduser().resolve() for path in directories]
        if len(resolved) != len(set(resolved)):
            raise MomentumAuditError(f"{label}: duplicate survey directory")
        overlap = resolved_seen.intersection(resolved)
        if overlap:
            raise MomentumAuditError(
                f"survey directory appears in multiple groups: {sorted(overlap)[0]}"
            )
        resolved_seen.update(resolved)
        requested_groups[label] = resolved

    _validate_cross_group_compatibility(requested_groups, config)
    for label, resolved in requested_groups.items():
        split_groups[label] = _audit_group(
            resolved,
            config,
            apply_y_max=args.apply_y_max,
            electron_p_min=electron_p_min,
        )

    pooled = (
        split_groups["pooled"]
        if "pooled" in split_groups
        else _combine_groups(
            (split_groups["training"], split_groups["validation"])
        )
    )

    catalog = radiative_guards.enumerate_strata(config)
    selected_indices, selection_source = _load_selection(
        args.flat_index_file, catalog
    )
    selected_ids = {
        stratum.identifier
        for stratum in catalog
        if stratum.flat_index in selected_indices
    }

    records: list[dict[str, object]] = []
    status_counts: dict[str, int] = {}
    selected_status_counts: dict[str, int] = {}
    fraction_bands = {
        "zero": 0,
        "greater_than_zero_below_0p1": 0,
        "0p1_to_below_0p5": 0,
        "0p5_to_below_0p8": 0,
        "0p8_to_below_0p95": 0,
        "0p95_to_one": 0,
        "no_survey_contribution": 0,
    }
    for stratum in catalog:
        identifier = stratum.identifier
        pooled_metrics = _fraction_metrics(
            pooled.nominal.get(identifier, radiative_guards.Moment()),
            pooled.below_threshold.get(identifier, radiative_guards.Moment()),
            pooled.proposals,
        )
        record = {
            **_stratum_metadata(stratum),
            "selected": stratum.flat_index in selected_indices,
            "pooled": pooled_metrics,
        }
        if "training" in split_groups:
            training_metrics = _fraction_metrics(
                split_groups["training"].nominal.get(
                    identifier, radiative_guards.Moment()
                ),
                split_groups["training"].below_threshold.get(
                    identifier, radiative_guards.Moment()
                ),
                split_groups["training"].proposals,
            )
            validation_metrics = _fraction_metrics(
                split_groups["validation"].nominal.get(
                    identifier, radiative_guards.Moment()
                ),
                split_groups["validation"].below_threshold.get(
                    identifier, radiative_guards.Moment()
                ),
                split_groups["validation"].proposals,
            )
            record["training"] = training_metrics
            record["validation"] = validation_metrics
            record["training_validation_fraction_difference_z_score"] = (
                _difference_z_score(training_metrics, validation_metrics)
            )
        records.append(record)

        status = str(pooled_metrics["status"])
        status_counts[status] = status_counts.get(status, 0) + 1
        if not record["selected"]:
            continue
        selected_status_counts[status] = selected_status_counts.get(status, 0) + 1
        fraction = pooled_metrics["nominal_cross_section_fraction"]
        if fraction is None:
            band = "no_survey_contribution"
        elif float(fraction) <= 0.0:
            band = "zero"
        elif float(fraction) < 0.1:
            band = "greater_than_zero_below_0p1"
        elif float(fraction) < 0.5:
            band = "0p1_to_below_0p5"
        elif float(fraction) < 0.8:
            band = "0p5_to_below_0p8"
        elif float(fraction) < 0.95:
            band = "0p8_to_below_0p95"
        else:
            band = "0p95_to_one"
        fraction_bands[band] += 1

    threshold_energy = math.hypot(electron_p_min, ELECTRON_MASS_GEV)
    beam_energy = float(config["beam_energy"])
    equivalent_y_maximum = 1.0 - threshold_energy / beam_energy
    payload: dict[str, object] = {
        "schema": AUDIT_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "builder_revision": radiative_mode4._source_revision(),
        "builder_source_sha256": _sha256(Path(__file__).resolve()),
        "analysis_config": str(config_path),
        "analysis_config_sha256": config_sha256,
        "analysis_selection": _analysis_selection(
            config,
            apply_y_max=args.apply_y_max,
            electron_p_min=electron_p_min,
        ),
        "threshold": {
            "source": threshold_source,
            "electron_momentum_minimum_GeV": electron_p_min,
            "electron_energy_minimum_GeV": threshold_energy,
            "equivalent_observed_y_maximum": equivalent_y_maximum,
            "electron_mass_GeV": ELECTRON_MASS_GEV,
        },
        "selection": (
            selection_source
            if selection_source is not None
            else {
                "path": None,
                "sha256": None,
                "selected_strata": len(catalog),
                "basis": "full_analysis_catalog",
            }
        ),
        "survey_groups": {
            label: _group_sources(directories, split_groups[label])
            for label, directories in requested_groups.items()
        },
        "pooling": {
            "estimator": (
                "sum of fixed-trial survey contributions divided by total "
                "fixed proposals; nominal fraction is the ratio of nominal "
                "and denominator contribution sums"
            ),
            "proposals": pooled.proposals,
            "raw_row_counts_are_not_used_as_cross_section_weights": True,
            "fraction_uncertainty": (
                "delta method with the exact disjoint pass/fail covariance "
                "and the full fixed-trial proposal denominator"
            ),
        },
        "summary": {
            "catalog_strata": len(catalog),
            "selected_strata": len(selected_indices),
            "catalog_status_counts": status_counts,
            "selected_status_counts": selected_status_counts,
            "selected_fraction_bands": fraction_bands,
            "selected_cross_section": _group_summary(pooled, selected_ids),
        },
        "interpretation": {
            "nominal_fraction": (
                "expected retained fraction of a loose, bin-conditional "
                "unweighted sample after applying the final-electron "
                "momentum threshold"
            ),
            "loose_events_per_expected_nominal_event": (
                "inverse nominal fraction; a planning diagnostic, not an "
                "event-by-event weight"
            ),
            "finite_survey_boundary_warning": (
                "an observed fraction of zero or one is finite-survey "
                "evidence, not proof of exact zero or one"
            ),
            "constant_per_stratum_weights_remain_valid": True,
        },
        "strata": records,
    }

    output = args.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"{output} already exists")
    output.mkdir(parents=True)
    report_path = output / REPORT_FILENAME
    _write_json(report_path, payload)

    fields = [
        "flat_index",
        "stratum_id",
        "selected",
        "Q2",
        "xB",
        "minus_t",
        "phi_deg",
        "status",
        "denominator_cross_section_microbarn",
        "denominator_sem_microbarn",
        "denominator_ess",
        "denominator_rows",
        "nominal_cross_section_microbarn",
        "nominal_sem_microbarn",
        "nominal_ess",
        "nominal_rows",
        "below_threshold_cross_section_microbarn",
        "below_threshold_sem_microbarn",
        "below_threshold_ess",
        "below_threshold_rows",
        "nominal_cross_section_fraction",
        "nominal_fraction_sem_delta_method",
        "raw_nominal_row_fraction_diagnostic_only",
        "loose_events_per_expected_nominal_event",
        "finite_survey_boundary",
        "training_nominal_fraction",
        "validation_nominal_fraction",
        "training_validation_fraction_difference_z_score",
    ]
    with (output / TABLE_FILENAME).open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        writer = csv.DictWriter(
            destination, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for record in records:
            pooled_metrics = record["pooled"]
            denominator = pooled_metrics["denominator"]
            nominal = pooled_metrics["nominal"]
            below = pooled_metrics["below_threshold"]
            bounds = record["bounds"]
            writer.writerow(
                {
                    "flat_index": record["flat_index"],
                    "stratum_id": record["stratum_id"],
                    "selected": record["selected"],
                    "Q2": ":".join(str(value) for value in bounds["Q2"]),
                    "xB": ":".join(str(value) for value in bounds["xB"]),
                    "minus_t": ":".join(
                        str(value) for value in bounds["minus_t"]
                    ),
                    "phi_deg": ":".join(
                        str(value) for value in bounds["phi_deg"]
                    ),
                    "status": pooled_metrics["status"],
                    "denominator_cross_section_microbarn": denominator[
                        "cross_section_microbarn"
                    ],
                    "denominator_sem_microbarn": denominator[
                        "cross_section_sem_microbarn"
                    ],
                    "denominator_ess": denominator[
                        "importance_effective_sample_size"
                    ],
                    "denominator_rows": denominator["contributing_rows"],
                    "nominal_cross_section_microbarn": nominal[
                        "cross_section_microbarn"
                    ],
                    "nominal_sem_microbarn": nominal[
                        "cross_section_sem_microbarn"
                    ],
                    "nominal_ess": nominal[
                        "importance_effective_sample_size"
                    ],
                    "nominal_rows": nominal["contributing_rows"],
                    "below_threshold_cross_section_microbarn": below[
                        "cross_section_microbarn"
                    ],
                    "below_threshold_sem_microbarn": below[
                        "cross_section_sem_microbarn"
                    ],
                    "below_threshold_ess": below[
                        "importance_effective_sample_size"
                    ],
                    "below_threshold_rows": below["contributing_rows"],
                    "nominal_cross_section_fraction": pooled_metrics[
                        "nominal_cross_section_fraction"
                    ],
                    "nominal_fraction_sem_delta_method": pooled_metrics[
                        "nominal_fraction_sem_delta_method"
                    ],
                    "raw_nominal_row_fraction_diagnostic_only": pooled_metrics[
                        "raw_nominal_row_fraction_diagnostic_only"
                    ],
                    "loose_events_per_expected_nominal_event": pooled_metrics[
                        "loose_events_per_expected_nominal_event"
                    ],
                    "finite_survey_boundary": pooled_metrics[
                        "finite_survey_boundary"
                    ],
                    "training_nominal_fraction": (
                        record.get("training", {}).get(
                            "nominal_cross_section_fraction"
                        )
                    ),
                    "validation_nominal_fraction": (
                        record.get("validation", {}).get(
                            "nominal_cross_section_fraction"
                        )
                    ),
                    "training_validation_fraction_difference_z_score": record.get(
                        "training_validation_fraction_difference_z_score"
                    ),
                }
            )

    low_threshold = float(args.low_nominal_fraction)
    low = []
    unresolved = []
    no_nominal = []
    for record in records:
        if not record["selected"]:
            continue
        metrics = record["pooled"]
        fraction = metrics["nominal_cross_section_fraction"]
        if fraction is None:
            unresolved.append(int(record["flat_index"]))
        elif float(fraction) <= 0.0:
            no_nominal.append(int(record["flat_index"]))
        elif float(fraction) < low_threshold:
            low.append(int(record["flat_index"]))
    for name, values in (
        ("low_nominal_fraction_flat_indices.txt", low),
        ("no_nominal_contribution_observed_flat_indices.txt", no_nominal),
        ("unresolved_no_survey_contribution_flat_indices.txt", unresolved),
    ):
        (output / name).write_text(
            "".join(f"{value}\n" for value in values), encoding="utf-8"
        )
    if args.flat_index_file is not None:
        (output / "selected_flat_indices.txt").write_text(
            "".join(f"{value}\n" for value in sorted(selected_indices)),
            encoding="utf-8",
        )
    return report_path


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    command = subparsers.add_parser(
        "audit",
        help=(
            "estimate per-stratum retention under a final-electron momentum cut"
        ),
    )
    command.add_argument("--config", type=Path, required=True)
    command.add_argument(
        "--survey",
        dest="surveys",
        type=Path,
        nargs="+",
        help="survey directories to pool without an independent split",
    )
    command.add_argument(
        "--training-survey",
        dest="training_surveys",
        type=Path,
        nargs="+",
        help="training survey directories",
    )
    command.add_argument(
        "--validation-survey",
        dest="validation_surveys",
        type=Path,
        nargs="+",
        help="independent validation survey directories",
    )
    command.add_argument(
        "--electron-p-min",
        type=float,
        help="GeV; defaults to phase_space.electron_p_min in the config",
    )
    command.add_argument(
        "--apply-y-max",
        action="store_true",
        help="apply phase_space.y_max to the loose denominator",
    )
    command.add_argument(
        "--flat-index-file",
        type=Path,
        help="optional selected-stratum list for campaign summaries and lists",
    )
    command.add_argument(
        "--low-nominal-fraction",
        type=float,
        default=0.5,
        help="threshold for the low-retention planning list (default: 0.5)",
    )
    command.add_argument("--output", type=Path, required=True)
    return parser


def main() -> int:
    args = _parser().parse_args()
    if not 0.0 < args.low_nominal_fraction <= 1.0:
        raise SystemExit("--low-nominal-fraction must lie in (0,1]")
    result = audit(args)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
