#!/usr/bin/env python3
"""Learn and validate hard-parent radiative migration footprints.

This milestone-2b diagnostic relates the Born-like hard-vertex coordinates of
an unrestricted fixed-trial survey to the final-LUND analysis stratum.  It is
diagnostic only: it does not change AAO sampling, perform unweighting, or emit
LUND events.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import datetime as dt
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import radiative_guards as guards
import radiative_survey


MANIFEST_SCHEMA = "aao-rad-migration-v1"
VALIDATION_SCHEMA = "aao-rad-migration-validation-v1"
PLOT_SCHEMA = "aao-rad-migration-plots-v1"
HARD_COORDINATE_DEFINITION = "hard_vertex_born_parent"
OBSERVED_COORDINATE_DEFINITION = "final_lund_analysis"
PARENT_AXES = ("Q2", "xB", "minus_t", "phi_deg")
HARD_COLUMNS = {
    "Q2": "q2_hard",
    "xB": "xb_hard",
    "minus_t": "minus_t_hard",
    "phi_deg": "phi_cm_deg",
}


class MigrationError(RuntimeError):
    """Raised when inputs cannot define a reproducible migration campaign."""


@dataclass(frozen=True, order=True)
class ParentIndex:
    iq2: int
    ixb: int
    it: int
    iphi: int


@dataclass(frozen=True)
class ParentGrid:
    q2_edges: tuple[float, ...]
    xb_edges: tuple[float, ...]
    minus_t_edges: tuple[float, ...]
    phi_edges: tuple[float, ...]

    @classmethod
    def from_config(cls, config: dict) -> "ParentGrid":
        return cls(
            q2_edges=guards._strict_edges(config["binning"]["Q2"], "Q2"),
            xb_edges=guards._strict_edges(config["binning"]["xB"], "xB"),
            minus_t_edges=guards._strict_edges(
                config["binning"]["minus_t"], "minus_t"
            ),
            phi_edges=guards._strict_edges(
                config["binning"]["phi_deg"], "phi_deg"
            ),
        )

    @property
    def edges(self) -> dict[str, tuple[float, ...]]:
        return {
            "Q2": self.q2_edges,
            "xB": self.xb_edges,
            "minus_t": self.minus_t_edges,
            "phi_deg": self.phi_edges,
        }

    @property
    def shape(self) -> tuple[int, int, int, int]:
        return (
            len(self.q2_edges) + 1,
            len(self.xb_edges) + 1,
            len(self.minus_t_edges) + 1,
            len(self.phi_edges) - 1,
        )


@dataclass
class MigrationCampaign:
    proposals: int
    rows: int
    final_valid_rows: int
    replicas: list[dict[str, object]]
    global_observed: guards.Moment
    outside: dict[str, guards.Moment]
    strata: dict[str, guards.Moment]
    migrations: dict[str, dict[tuple[str, int], guards.Moment]]
    parent_totals: dict[tuple[str, int], guards.Moment]
    relationships: dict[str, guards.Moment]
    channels: dict[int, guards.Moment]
    norm_reference: dict[str, str]
    input_signature: str
    legacy_input_settings: dict[str, object]
    channel_probabilities: dict[str, float]


def _now() -> str:
    return (
        dt.datetime.now(dt.timezone.utc)
        .replace(microsecond=0)
        .isoformat()
    )


def _parent_grid(config_or_grid: dict | ParentGrid) -> ParentGrid:
    if isinstance(config_or_grid, ParentGrid):
        return config_or_grid
    return ParentGrid.from_config(config_or_grid)


def _config_edges(
    config_or_grid: dict | ParentGrid,
) -> dict[str, tuple[float, ...]]:
    return _parent_grid(config_or_grid).edges


def _parent_shape(
    config_or_grid: dict | ParentGrid,
) -> tuple[int, int, int, int]:
    return _parent_grid(config_or_grid).shape


def _extended_index(value: float, edges: tuple[float, ...]) -> int:
    """Return -1/regular-bin/n for underflow/interior/overflow."""
    if not math.isfinite(value):
        raise MigrationError("hard-parent coordinates must be finite")
    if value < edges[0]:
        return -1
    if value >= edges[-1]:
        return len(edges) - 1
    return bisect.bisect_right(edges, value) - 1


def parent_from_row(
    row: dict[str, str], config_or_grid: dict | ParentGrid
) -> ParentIndex:
    edges = _config_edges(config_or_grid)
    try:
        q2 = float(row[HARD_COLUMNS["Q2"]])
        xb = float(row[HARD_COLUMNS["xB"]])
        minus_t = float(row[HARD_COLUMNS["minus_t"]])
        phi = float(row[HARD_COLUMNS["phi_deg"]])
        phi_base = float(row["hadron_phi_base"])
    except (KeyError, TypeError, ValueError) as error:
        raise MigrationError("survey row has malformed hard coordinates") from error
    if not all(math.isfinite(value) for value in (q2, xb, minus_t, phi)):
        raise MigrationError("survey row has nonfinite hard coordinates")
    expected_phi = 360.0 * phi_base
    phi_difference = abs((phi - expected_phi + 180.0) % 360.0 - 180.0)
    if phi_difference > 3.0e-5:
        raise MigrationError(
            "phi_cm_deg is inconsistent with hadron_phi_base"
        )
    phi_edges = edges["phi_deg"]
    wrapped_phi = (phi - phi_edges[0]) % 360.0 + phi_edges[0]
    iphi = guards._bin_index(wrapped_phi, phi_edges)
    if iphi is None:
        raise MigrationError("periodic hard phi could not be assigned")
    return ParentIndex(
        _extended_index(q2, edges["Q2"]),
        _extended_index(xb, edges["xB"]),
        _extended_index(minus_t, edges["minus_t"]),
        iphi,
    )


def parent_identifier(
    parent: ParentIndex, config_or_grid: dict | ParentGrid
) -> str:
    shape = _parent_shape(config_or_grid)
    encoded = (
        parent.iq2 + 1,
        parent.ixb + 1,
        parent.it + 1,
        parent.iphi,
    )
    for index, size in zip(encoded, shape):
        if not 0 <= index < size:
            raise MigrationError(f"hard-parent index {parent} is out of range")
    flat = 0
    for index, size in zip(encoded, shape):
        flat = flat * size + index
    width = max(5, len(str(math.prod(shape) - 1)))
    return f"h{flat:0{width}d}"


def parent_from_identifier(
    identifier: str, config_or_grid: dict | ParentGrid
) -> ParentIndex:
    if not identifier.startswith("h"):
        raise MigrationError(f"invalid hard-parent identifier {identifier!r}")
    try:
        flat = int(identifier[1:])
    except ValueError as error:
        raise MigrationError(
            f"invalid hard-parent identifier {identifier!r}"
        ) from error
    shape = _parent_shape(config_or_grid)
    if not 0 <= flat < math.prod(shape):
        raise MigrationError(f"hard-parent identifier {identifier!r} is out of range")
    encoded = [0] * len(shape)
    for position in range(len(shape) - 1, -1, -1):
        encoded[position] = flat % shape[position]
        flat //= shape[position]
    return ParentIndex(
        encoded[0] - 1,
        encoded[1] - 1,
        encoded[2] - 1,
        encoded[3],
    )


def component_identifier(parent_id: str, channel: int) -> str:
    if channel not in range(1, 7):
        raise MigrationError(f"invalid radiative channel {channel}")
    return f"{parent_id}:c{channel}"


def component_from_identifier(
    identifier: str, config_or_grid: dict | ParentGrid
) -> tuple[str, int]:
    try:
        parent_id, channel_text = identifier.rsplit(":c", 1)
        channel = int(channel_text)
    except (ValueError, TypeError) as error:
        raise MigrationError(
            f"invalid parent-component identifier {identifier!r}"
        ) from error
    parent_from_identifier(parent_id, config_or_grid)
    if channel not in range(1, 7):
        raise MigrationError(
            f"invalid parent-component identifier {identifier!r}"
        )
    return parent_id, channel


def _axis_bin_metadata(
    index: int, edges: tuple[float, ...]
) -> dict[str, object]:
    regular_bins = len(edges) - 1
    if index == -1:
        return {
            "index": index,
            "region": "underflow",
            "bounds": [None, edges[0]],
        }
    if index == regular_bins:
        return {
            "index": index,
            "region": "overflow",
            "bounds": [edges[-1], None],
        }
    if not 0 <= index < regular_bins:
        raise MigrationError(f"parent index {index} is invalid")
    return {
        "index": index,
        "region": "analysis_bin",
        "bounds": [edges[index], edges[index + 1]],
    }


def parent_metadata(
    parent: ParentIndex, config_or_grid: dict | ParentGrid
) -> dict[str, object]:
    edges = _config_edges(config_or_grid)
    return {
        "parent_id": parent_identifier(parent, config_or_grid),
        "coordinates": {
            "Q2": _axis_bin_metadata(parent.iq2, edges["Q2"]),
            "xB": _axis_bin_metadata(parent.ixb, edges["xB"]),
            "minus_t": _axis_bin_metadata(parent.it, edges["minus_t"]),
            "phi_deg": _axis_bin_metadata(parent.iphi, edges["phi_deg"]),
        },
        "has_underflow_or_overflow": any(
            item["region"] != "analysis_bin"
            for item in (
                _axis_bin_metadata(parent.iq2, edges["Q2"]),
                _axis_bin_metadata(parent.ixb, edges["xB"]),
                _axis_bin_metadata(parent.it, edges["minus_t"]),
            )
        ),
    }


def _component_neighbors(
    component: tuple[str, int], config_or_grid: dict | ParentGrid
) -> Iterator[tuple[str, int]]:
    parent_id, channel = component
    parent = parent_from_identifier(parent_id, config_or_grid)
    indices = [parent.iq2 + 1, parent.ixb + 1, parent.it + 1, parent.iphi]
    shape = _parent_shape(config_or_grid)
    for dimension in range(4):
        for shift in (-1, 1):
            neighbor = indices.copy()
            candidate = neighbor[dimension] + shift
            if dimension == 3:
                candidate %= shape[dimension]
            elif not 0 <= candidate < shape[dimension]:
                continue
            neighbor[dimension] = candidate
            neighbor_parent = ParentIndex(
                neighbor[0] - 1,
                neighbor[1] - 1,
                neighbor[2] - 1,
                neighbor[3],
            )
            yield parent_identifier(neighbor_parent, config_or_grid), channel


def dilate_components(
    seeds: set[tuple[str, int]],
    config_or_grid: dict | ParentGrid,
    radius: int,
) -> set[tuple[str, int]]:
    if radius < 0:
        raise MigrationError("parent dilation must be nonnegative")
    result = set(seeds)
    frontier = set(seeds)
    for _ in range(radius):
        additions = {
            neighbor
            for component in frontier
            for neighbor in _component_neighbors(
                component, config_or_grid
            )
        }
        additions -= result
        result.update(additions)
        frontier = additions
        if not frontier:
            break
    return result


def _target_lookup(config: dict) -> dict[str, guards.Stratum]:
    return {
        stratum.identifier: stratum
        for stratum in guards.enumerate_strata(config)
    }


def _periodic_index_distance(left: int, right: int, bins: int) -> int:
    difference = abs(left - right)
    return min(difference, bins - difference)


def parent_relationship(
    parent: ParentIndex,
    target: guards.Stratum,
    config_or_grid: dict | ParentGrid,
) -> tuple[str, dict[str, int | None]]:
    edges = _config_edges(config_or_grid)
    outside = (
        parent.iq2 not in range(len(edges["Q2"]) - 1)
        or parent.ixb not in range(len(edges["xB"]) - 1)
        or parent.it not in range(len(edges["minus_t"]) - 1)
    )
    if outside:
        return "hard_parent_underflow_or_overflow", {
            "Q2": None,
            "xB": None,
            "minus_t": None,
            "phi_deg": _periodic_index_distance(
                parent.iphi, target.iphi, len(edges["phi_deg"]) - 1
            ),
        }
    deltas = {
        "Q2": parent.iq2 - target.iq2,
        "xB": parent.ixb - target.ixb,
        "minus_t": parent.it - target.it,
        "phi_deg": (
            (parent.iphi - target.iphi + (len(edges["phi_deg"]) - 1) // 2)
            % (len(edges["phi_deg"]) - 1)
            - (len(edges["phi_deg"]) - 1) // 2
        ),
    }
    if all(value == 0 for value in deltas.values()):
        return "same_hard_stratum", deltas
    if all(abs(value) <= 1 for value in deltas.values()):
        return "neighboring_hard_stratum", deltas
    return "nonlocal_hard_stratum", deltas


def _combined(moments: Iterable[guards.Moment]) -> guards.Moment:
    return guards._combined_moment(moments)


def _subset(
    moments: dict[tuple[str, int], guards.Moment],
    selected: set[tuple[str, int]],
) -> guards.Moment:
    return _combined(
        moment for key, moment in moments.items() if key in selected
    )


def _select_seed_components(
    moments: dict[tuple[str, int], guards.Moment], target_fraction: float
) -> set[tuple[str, int]]:
    total = sum(moment.total for moment in moments.values())
    if total <= 0.0:
        return set()
    selected: set[tuple[str, int]] = set()
    retained = 0.0
    for key, moment in sorted(
        moments.items(), key=lambda item: (-item[1].total, item[0])
    ):
        if moment.total <= 0.0:
            continue
        selected.add(key)
        retained += moment.total
        if retained >= target_fraction * total:
            break
    return selected


def _component_ids(components: set[tuple[str, int]]) -> list[str]:
    return [
        component_identifier(parent_id, channel)
        for parent_id, channel in sorted(components)
    ]


def _components_from_record(
    record: dict,
    config_or_grid: dict | ParentGrid,
    *,
    extra_steps: int = 0,
) -> set[tuple[str, int]]:
    definition = record["parent_footprint"]
    seeds = {
        component_from_identifier(identifier, config_or_grid)
        for identifier in definition["seed_component_ids"]
    }
    radius = int(definition["dilation_axis_steps"]) + extra_steps
    return dilate_components(seeds, config_or_grid, radius)


def _selected_purity(
    target_moments: dict[tuple[str, int], guards.Moment],
    parent_totals: dict[tuple[str, int], guards.Moment],
    selected: set[tuple[str, int]],
) -> float | None:
    numerator = _subset(target_moments, selected).total
    denominator = sum(
        parent_totals.get(component, guards.Moment()).total
        for component in selected
    )
    return numerator / denominator if denominator > 0.0 else None


def aggregate_surveys(
    directories: list[Path],
    config: dict,
    *,
    apply_y_max: bool = False,
) -> MigrationCampaign:
    if not directories:
        raise MigrationError("at least one survey directory is required")
    proposals = rows = final_valid_rows = 0
    replicas: list[dict[str, object]] = []
    replica_ids: set[int] = set()
    global_observed = guards.Moment()
    outside: dict[str, guards.Moment] = {}
    strata: dict[str, guards.Moment] = {}
    migrations: dict[
        str, dict[tuple[str, int], guards.Moment]
    ] = {}
    parent_totals: dict[tuple[str, int], guards.Moment] = {}
    relationships: dict[str, guards.Moment] = {}
    channels: dict[int, guards.Moment] = {}
    norm_reference: dict[str, str] | None = None
    input_signature: str | None = None
    legacy_input_settings: dict[str, object] | None = None
    channel_probabilities: dict[str, float] | None = None
    targets = _target_lookup(config)
    parent_grid = ParentGrid.from_config(config)

    for requested_directory in directories:
        directory = requested_directory.resolve()
        try:
            norm = guards._validated_norm(directory)
        except guards.GuardLearningError as error:
            raise MigrationError(str(error)) from error
        if norm_reference is None:
            norm_reference = norm
        else:
            try:
                guards._compatible_settings(norm_reference, norm)
            except guards.GuardLearningError as error:
                raise MigrationError(str(error)) from error
        try:
            signature, probabilities, settings = guards._legacy_input_metadata(
                directory / "survey_input.inp"
            )
        except guards.GuardLearningError as error:
            raise MigrationError(str(error)) from error
        if input_signature is None:
            input_signature = signature
            channel_probabilities = probabilities
            legacy_input_settings = settings
        elif signature != input_signature:
            raise MigrationError(
                "survey legacy inputs differ; migration replicas must have "
                "identical physical and proposal settings"
            )

        try:
            ntrials = int(norm["ntries"])
            requested = int(norm["survey_ntrials_requested"])
            expected_rows = int(norm["survey_rows"])
            expected_valid = int(norm["survey_final_valid"])
            replica = int(norm["survey_replica"])
        except (KeyError, ValueError) as error:
            raise MigrationError(
                f"{directory}: incomplete survey normalization metadata"
            ) from error
        if ntrials != requested or ntrials <= 0:
            raise MigrationError(f"{directory}: incomplete fixed-trial survey")
        if replica in replica_ids:
            raise MigrationError(f"duplicate survey replica ID {replica}")
        replica_ids.add(replica)

        local = guards.Moment()
        local_rows = local_valid = 0
        previous_trial = 0
        phase_volume = float(norm["survey_phase_volume"])
        beam_energy = float(norm["ebeam"])
        try:
            survey_rows = guards._survey_rows(
                directory / radiative_survey.SURVEY_FILENAME
            )
            for row in survey_rows:
                local_rows += 1
                guards._validate_guard_row(
                    row,
                    directory=directory,
                    phase_volume=phase_volume,
                    beam_energy=beam_energy,
                )
                trial = int(row["trial"])
                row_replica = int(row["replica"])
                valid = int(row["final_valid"])
                channel = int(row["intreg"])
                weight = float(row["trial_xsec_observed_microbarn"])
                if trial <= previous_trial or trial > ntrials:
                    raise MigrationError(
                        f"{directory}: trial IDs must increase within [1,ntries]"
                    )
                previous_trial = trial
                if row_replica != replica:
                    raise MigrationError(
                        f"{directory}: CSV replica ID mismatch"
                    )
                if valid not in (0, 1) or channel not in range(1, 7):
                    raise MigrationError(
                        f"{directory}: invalid status or channel"
                    )
                if not math.isfinite(weight) or weight < 0.0:
                    raise MigrationError(
                        f"{directory}: invalid observed contribution"
                    )

                rows += 1
                component: tuple[str, int] | None = None
                if valid == 1:
                    final_valid_rows += 1
                    local_valid += 1
                    local.add(weight)
                    global_observed.add(weight)
                    parent = parent_from_row(row, parent_grid)
                    component = (
                        parent_identifier(parent, parent_grid),
                        channel,
                    )
                    parent_totals.setdefault(
                        component, guards.Moment()
                    ).add(weight)

                stratum_id, outside_reason = guards.assign_stratum(
                    row, config, apply_y_max=apply_y_max
                )
                if stratum_id is None:
                    outside.setdefault(
                        outside_reason or "unknown", guards.Moment()
                    ).add(weight)
                    continue
                if component is None:
                    raise MigrationError(
                        f"{directory}: assigned stratum lacks hard parent"
                    )
                strata.setdefault(stratum_id, guards.Moment()).add(weight)
                migrations.setdefault(stratum_id, {}).setdefault(
                    component, guards.Moment()
                ).add(weight)
                channel_moment = channels.setdefault(
                    channel, guards.Moment()
                )
                channel_moment.add(weight)
                parent = parent_from_identifier(
                    component[0], parent_grid
                )
                relationship, _ = parent_relationship(
                    parent, targets[stratum_id], parent_grid
                )
                relationships.setdefault(
                    relationship, guards.Moment()
                ).add(weight)
        except (guards.GuardLearningError, ValueError, TypeError) as error:
            raise MigrationError(f"{directory}: {error}") from error

        if local_rows != expected_rows or local_valid != expected_valid:
            raise MigrationError(
                f"{directory}: CSV row counts disagree with normalization metadata"
            )
        observed_estimate = local.total / ntrials
        expected_estimate = float(norm["survey_observed_sig_sum"])
        if not math.isclose(
            observed_estimate,
            expected_estimate,
            rel_tol=5.0e-6,
            abs_tol=1.0e-15,
        ):
            raise MigrationError(
                f"{directory}: CSV observed integral "
                f"{observed_estimate:.12g} does not match norm "
                f"{expected_estimate:.12g}"
            )
        proposals += ntrials
        replicas.append(
            {
                "directory": str(directory),
                "replica": replica,
                "seed": int(norm["survey_seed"]),
                "proposals": ntrials,
                "rows": local_rows,
                "final_valid_rows": local_valid,
                "observed_cross_section_microbarn": observed_estimate,
                "observed_cross_section_sem_microbarn": guards._sem(
                    local, ntrials
                ),
            }
        )

    return MigrationCampaign(
        proposals=proposals,
        rows=rows,
        final_valid_rows=final_valid_rows,
        replicas=sorted(replicas, key=lambda item: int(item["replica"])),
        global_observed=global_observed,
        outside=outside,
        strata=strata,
        migrations=migrations,
        parent_totals=parent_totals,
        relationships=relationships,
        channels=channels,
        norm_reference=norm_reference or {},
        input_signature=input_signature or "",
        legacy_input_settings=legacy_input_settings or {},
        channel_probabilities=channel_probabilities or {},
    )


def _fraction_metrics(
    moments: dict[object, guards.Moment],
    proposals: int,
) -> dict[str, dict[str, object]]:
    total = sum(moment.total for moment in moments.values())
    result: dict[str, dict[str, object]] = {}
    for name, moment in sorted(moments.items(), key=lambda item: str(item[0])):
        result[str(name)] = {
            **guards._metrics(moment, proposals),
            "fraction_of_inside_analysis_cross_section": (
                moment.total / total if total > 0.0 else None
            ),
        }
    return result


def _leading_components(
    moments: dict[tuple[str, int], guards.Moment],
    config_or_grid: dict | ParentGrid,
    total: float,
    *,
    limit: int = 10,
) -> list[dict[str, object]]:
    leading: list[dict[str, object]] = []
    for (parent_id, channel), moment in sorted(
        moments.items(), key=lambda item: (-item[1].total, item[0])
    )[:limit]:
        parent = parent_from_identifier(parent_id, config_or_grid)
        leading.append(
            {
                "component_id": component_identifier(parent_id, channel),
                "parent_id": parent_id,
                "channel": channel,
                "hard_indices": {
                    "Q2": parent.iq2,
                    "xB": parent.ixb,
                    "minus_t": parent.it,
                    "phi_deg": parent.iphi,
                },
                "fraction_of_stratum_cross_section": (
                    moment.total / total if total > 0.0 else None
                ),
                "contributing_rows": moment.count,
                "sum_trial_contributions_microbarn": moment.total,
                "sum_squared_trial_contributions_microbarn2": (
                    moment.square_total
                ),
                "largest_trial_contribution_microbarn": moment.maximum,
                "importance_effective_sample_size": guards._ess(moment),
            }
        )
    return leading


def build_manifest(
    *,
    config: dict,
    config_path: Path,
    config_sha256: str,
    campaign: MigrationCampaign,
    target_parent_fraction: float,
    parent_dilation: int,
    iteration: int,
    generator_revision: str,
    generator_revision_source: str,
    learner_revision: str,
    minimum_training_rows: int,
    minimum_training_ess: float,
    apply_y_max: bool,
) -> dict:
    catalog = guards.enumerate_strata(config)
    parent_grid = ParentGrid.from_config(config)
    strata_output: dict[str, dict] = {}
    learned = low_support = empty = 0
    for stratum in catalog:
        total = campaign.strata.get(stratum.identifier, guards.Moment())
        moments = campaign.migrations.get(stratum.identifier, {})
        record = guards._stratum_metadata(stratum)
        record["training_total"] = guards._metrics(
            total, campaign.proposals
        )
        if total.total <= 0.0:
            record.update(
                {
                    "status": "no_training_contribution",
                    "parent_footprint": {
                        "representation": (
                            "axis_dilation_of_weighted_hard_parent_components"
                        ),
                        "seed_component_ids": [],
                        "dilation_axis_steps": parent_dilation,
                        "derived_component_count": 0,
                    },
                    "estimated_selected_parent_fraction": 0.0,
                    "estimated_tail_fraction": 0.0,
                    "selected_parent_cross_section_purity_within_global_observed": (
                        None
                    ),
                    "leading_training_parent_components": [],
                }
            )
            empty += 1
            strata_output[stratum.identifier] = record
            continue

        seeds = _select_seed_components(moments, target_parent_fraction)
        selected = dilate_components(
            seeds, parent_grid, parent_dilation
        )
        expanded = dilate_components(
            seeds, parent_grid, parent_dilation + 1
        )
        selected_moment = _subset(moments, selected)
        tail_moment = _combined(
            moment for key, moment in moments.items() if key not in selected
        )
        expanded_moment = _subset(moments, expanded)
        sufficient = (
            total.count >= minimum_training_rows
            and guards._ess(total) >= minimum_training_ess
        )
        status = "learned" if sufficient else "learned_low_support"
        learned += int(sufficient)
        low_support += int(not sufficient)
        record.update(
            {
                "status": status,
                "parent_footprint": {
                    "representation": (
                        "axis_dilation_of_weighted_hard_parent_components"
                    ),
                    "seed_component_ids": _component_ids(seeds),
                    "dilation_axis_steps": parent_dilation,
                    "derived_component_count": len(selected),
                },
                "training_selected_parents": guards._metrics(
                    selected_moment, campaign.proposals
                ),
                "training_tail": guards._metrics(
                    tail_moment, campaign.proposals
                ),
                "training_one_more_dilation": guards._metrics(
                    expanded_moment, campaign.proposals
                ),
                "estimated_selected_parent_fraction": (
                    selected_moment.total / total.total
                ),
                "estimated_tail_fraction": tail_moment.total / total.total,
                "one_more_dilation_fraction": (
                    expanded_moment.total / total.total
                ),
                "selected_parent_cross_section_purity_within_global_observed": (
                    _selected_purity(
                        moments, campaign.parent_totals, selected
                    )
                ),
                "one_more_dilation_cross_section_purity_within_global_observed": (
                    _selected_purity(
                        moments, campaign.parent_totals, expanded
                    )
                ),
                "leading_training_parent_components": _leading_components(
                    moments, parent_grid, total.total
                ),
            }
        )
        strata_output[stratum.identifier] = record

    inside = _combined(campaign.strata.values())
    outside_total = _combined(campaign.outside.values())
    return {
        "schema": MANIFEST_SCHEMA,
        "created_utc": _now(),
        "production_ready": False,
        "production_readiness_note": (
            "Milestone 2b diagnostic artifact; it does not define or activate "
            "the mode-4 proposal."
        ),
        "coordinate_definitions": {
            "hard_parent": HARD_COORDINATE_DEFINITION,
            "hard_columns": HARD_COLUMNS,
            "observed_target": OBSERVED_COORDINATE_DEFINITION,
        },
        "artifacts": {
            "training_parent_coverage": "training_parent_coverage.csv",
            "training_migrations": "training_migrations.csv",
            "hash_sidecar_suffix": ".sha256",
        },
        "migration_iteration": iteration,
        "generator_revision": generator_revision,
        "generator_revision_source": generator_revision_source,
        "migration_learner_revision": learner_revision,
        "survey_schema": radiative_survey.SURVEY_SCHEMA,
        "analysis_config": config,
        "analysis_config_source": str(config_path.resolve()),
        "analysis_config_sha256": config_sha256,
        "analysis_selection": {
            "edge_convention": "lower_inclusive_upper_exclusive",
            "phi_periodic": True,
            "q2_minimum": float(
                config["phase_space"].get(
                    "Q2_min", config["binning"]["Q2"][0]
                )
            ),
            "w_minimum": float(config["phase_space"]["W_min"]),
            "apply_y_max": apply_y_max,
            "y_maximum": (
                float(config["phase_space"]["y_max"])
                if apply_y_max
                else None
            ),
            "no_implicit_y_minimum": True,
        },
        "hard_parent_definition": {
            "axes": [
                {
                    "name": name,
                    "survey_column": HARD_COLUMNS[name],
                    "analysis_edges": list(parent_grid.edges[name]),
                    "underflow_parent": name != "phi_deg",
                    "overflow_parent": name != "phi_deg",
                    "periodic": name == "phi_deg",
                }
                for name in PARENT_AXES
            ],
            "grid_shape_including_underflow_overflow": list(
                parent_grid.shape
            ),
            "radiative_channel_is_separate": True,
            "underflow_overflow_note": (
                "Hard Q2, xB, and -t values outside the analysis edges are "
                "retained as explicit parents so radiative feed-in is not lost."
            ),
        },
        "legacy_generator_settings": campaign.legacy_input_settings,
        "learning": {
            "target_parent_fraction_before_dilation": target_parent_fraction,
            "parent_dilation_axis_steps": parent_dilation,
            "minimum_training_rows_for_supported_label": minimum_training_rows,
            "minimum_training_ess_for_supported_label": minimum_training_ess,
            "learned_strata": learned,
            "learned_low_support_strata": low_support,
            "no_training_contribution_strata": empty,
        },
        "training": {
            "replica_ids": [item["replica"] for item in campaign.replicas],
            "replicas": campaign.replicas,
            "total_proposals": campaign.proposals,
            "recorded_rows": campaign.rows,
            "final_valid_rows": campaign.final_valid_rows,
            "legacy_input_sha256": campaign.input_signature,
            "global_observed": guards._metrics(
                campaign.global_observed, campaign.proposals
            ),
            "inside_analysis_partition": guards._metrics(
                inside, campaign.proposals
            ),
            "outside_analysis_partition": {
                name: guards._metrics(moment, campaign.proposals)
                for name, moment in sorted(campaign.outside.items())
            },
            "hard_parent_relationships_inside_analysis": _fraction_metrics(
                campaign.relationships, campaign.proposals
            ),
            "radiative_channels_inside_analysis": _fraction_metrics(
                campaign.channels, campaign.proposals
            ),
            "cross_section_partition_closure_microbarn": (
                campaign.global_observed.total
                - inside.total
                - outside_total.total
            )
            / campaign.proposals,
        },
        "strata": strata_output,
    }


def validate_manifest(
    manifest: dict,
    manifest_path: Path,
    campaign: MigrationCampaign,
    minimum_parent_coverage: float,
) -> dict:
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise MigrationError(
            f"manifest schema is {manifest.get('schema')!r}, "
            f"expected {MANIFEST_SCHEMA}"
        )
    if campaign.input_signature != manifest["training"]["legacy_input_sha256"]:
        raise MigrationError(
            "validation surveys do not match the manifest's legacy input settings"
        )
    training_ids = {
        int(value) for value in manifest["training"]["replica_ids"]
    }
    validation_ids = {int(item["replica"]) for item in campaign.replicas}
    overlap = training_ids & validation_ids
    if overlap:
        raise MigrationError(
            f"validation replicas overlap training replicas: {sorted(overlap)}"
        )

    config = manifest["analysis_config"]
    parent_grid = ParentGrid.from_config(config)
    strata_output: dict[str, dict] = {}
    assessed = passed = failed = no_holdout = 0
    for stratum_id, training_record in manifest["strata"].items():
        training_total = float(
            training_record["training_total"][
                "sum_trial_contributions_microbarn"
            ]
        )
        holdout_total = campaign.strata.get(
            stratum_id, guards.Moment()
        )
        if training_total <= 0.0 and holdout_total.total <= 0.0:
            continue
        moments = campaign.migrations.get(stratum_id, {})
        selected = _components_from_record(
            training_record, parent_grid
        )
        expanded = _components_from_record(
            training_record, parent_grid, extra_steps=1
        )
        selected_moment = _subset(moments, selected)
        tail_moment = _combined(
            moment for key, moment in moments.items() if key not in selected
        )
        expanded_moment = _subset(moments, expanded)
        if holdout_total.total > 0.0:
            selected_fraction = selected_moment.total / holdout_total.total
            expanded_fraction = expanded_moment.total / holdout_total.total
            coverage_passed = (
                selected_fraction >= minimum_parent_coverage
            )
            assessed += 1
            passed += int(coverage_passed)
            failed += int(not coverage_passed)
        else:
            selected_fraction = expanded_fraction = None
            coverage_passed = None
            no_holdout += 1
        largest_tail = []
        for (parent_id, channel), moment in sorted(
            (
                (key, moment)
                for key, moment in moments.items()
                if key not in selected and moment.total > 0.0
            ),
            key=lambda item: (-item[1].total, item[0]),
        )[:5]:
            parent = parent_from_identifier(parent_id, parent_grid)
            largest_tail.append(
                {
                    "component_id": component_identifier(parent_id, channel),
                    "parent_id": parent_id,
                    "channel": channel,
                    "hard_indices": {
                        "Q2": parent.iq2,
                        "xB": parent.ixb,
                        "minus_t": parent.it,
                        "phi_deg": parent.iphi,
                    },
                    **guards._metrics(moment, campaign.proposals),
                    "fraction_of_holdout_cross_section": (
                        moment.total / holdout_total.total
                        if holdout_total.total > 0.0
                        else None
                    ),
                }
            )
        training_metrics = training_record["training_total"]
        holdout_metrics = guards._metrics(
            holdout_total, campaign.proposals
        )
        strata_output[stratum_id] = {
            "training_status": training_record["status"],
            "training_selected_parent_fraction": training_record.get(
                "estimated_selected_parent_fraction", 0.0
            ),
            "holdout_total": holdout_metrics,
            "holdout_selected_parents": guards._metrics(
                selected_moment, campaign.proposals
            ),
            "holdout_tail": guards._metrics(
                tail_moment, campaign.proposals
            ),
            "holdout_one_more_dilation": guards._metrics(
                expanded_moment, campaign.proposals
            ),
            "holdout_selected_parent_fraction": selected_fraction,
            "holdout_tail_fraction": (
                1.0 - selected_fraction
                if selected_fraction is not None
                else None
            ),
            "holdout_one_more_dilation_fraction": expanded_fraction,
            "holdout_selected_parent_cross_section_purity_within_global_observed": (
                _selected_purity(
                    moments, campaign.parent_totals, selected
                )
            ),
            "holdout_one_more_dilation_cross_section_purity_within_global_observed": (
                _selected_purity(
                    moments, campaign.parent_totals, expanded
                )
            ),
            "largest_tail_parent_components": largest_tail,
            "training_holdout_difference_z_score": guards._difference_z_score(
                float(training_metrics["cross_section_microbarn"]),
                training_metrics["cross_section_sem_microbarn"],
                float(holdout_metrics["cross_section_microbarn"]),
                holdout_metrics["cross_section_sem_microbarn"],
            ),
            "coverage_passed": coverage_passed,
        }

    inside = _combined(campaign.strata.values())
    outside_total = _combined(campaign.outside.values())
    validation_global = guards._metrics(
        campaign.global_observed, campaign.proposals
    )
    training_global = manifest["training"]["global_observed"]
    return {
        "schema": VALIDATION_SCHEMA,
        "created_utc": _now(),
        "manifest": str(manifest_path.resolve()),
        "manifest_sha256": hashlib.sha256(
            manifest_path.read_bytes()
        ).hexdigest(),
        "coordinate_definitions": manifest["coordinate_definitions"],
        "artifacts": {
            "validation_parent_coverage": "validation_parent_coverage.csv",
            "validation_migrations": "validation_migrations.csv",
            "hash_sidecar_suffix": ".sha256",
        },
        "minimum_parent_coverage": minimum_parent_coverage,
        "passed": failed == 0 and assessed > 0,
        "global_training_holdout_difference_z_score": (
            guards._difference_z_score(
                float(training_global["cross_section_microbarn"]),
                training_global["cross_section_sem_microbarn"],
                float(validation_global["cross_section_microbarn"]),
                validation_global["cross_section_sem_microbarn"],
            )
        ),
        "coverage_summary": {
            "assessed_strata": assessed,
            "passed_strata": passed,
            "failed_strata": failed,
            "training_strata_without_holdout_contribution": no_holdout,
        },
        "validation": {
            "replica_ids": [item["replica"] for item in campaign.replicas],
            "replicas": campaign.replicas,
            "total_proposals": campaign.proposals,
            "recorded_rows": campaign.rows,
            "final_valid_rows": campaign.final_valid_rows,
            "global_observed": validation_global,
            "inside_analysis_partition": guards._metrics(
                inside, campaign.proposals
            ),
            "outside_analysis_partition": {
                name: guards._metrics(moment, campaign.proposals)
                for name, moment in sorted(campaign.outside.items())
            },
            "hard_parent_relationships_inside_analysis": _fraction_metrics(
                campaign.relationships, campaign.proposals
            ),
            "radiative_channels_inside_analysis": _fraction_metrics(
                campaign.channels, campaign.proposals
            ),
            "cross_section_partition_closure_microbarn": (
                campaign.global_observed.total
                - inside.total
                - outside_total.total
            )
            / campaign.proposals,
        },
        "strata": strata_output,
    }


def _write_json_artifact(
    output: Path, filename: str, payload: dict
) -> Path:
    if output.exists():
        raise FileExistsError(
            f"{output} already exists; migration artifacts are immutable"
        )
    output.mkdir(parents=True)
    path = output / filename
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    path.write_text(serialized, encoding="utf-8")
    guards._write_sha256(path)
    return path


def _write_coverage_csv(
    path: Path, payload: dict, *, validation: bool
) -> None:
    if validation:
        fields = (
            "stratum_id",
            "training_status",
            "holdout_cross_section_microbarn",
            "holdout_sem_microbarn",
            "holdout_ess",
            "holdout_selected_parent_fraction",
            "holdout_tail_fraction",
            "holdout_one_more_dilation_fraction",
            "holdout_selected_parent_purity",
            "coverage_passed",
            "largest_tail_parent_fraction",
        )
    else:
        fields = (
            "stratum_id",
            "flat_index",
            "status",
            "training_cross_section_microbarn",
            "training_sem_microbarn",
            "training_ess",
            "training_rows",
            "selected_parent_fraction",
            "tail_fraction",
            "one_more_dilation_fraction",
            "selected_parent_purity",
            "seed_parent_components",
            "derived_parent_components",
        )
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        for stratum_id, record in payload["strata"].items():
            if validation:
                total = record["holdout_total"]
                largest = record["largest_tail_parent_components"]
                writer.writerow(
                    {
                        "stratum_id": stratum_id,
                        "training_status": record["training_status"],
                        "holdout_cross_section_microbarn": total[
                            "cross_section_microbarn"
                        ],
                        "holdout_sem_microbarn": total[
                            "cross_section_sem_microbarn"
                        ],
                        "holdout_ess": total[
                            "importance_effective_sample_size"
                        ],
                        "holdout_selected_parent_fraction": record[
                            "holdout_selected_parent_fraction"
                        ],
                        "holdout_tail_fraction": record[
                            "holdout_tail_fraction"
                        ],
                        "holdout_one_more_dilation_fraction": record[
                            "holdout_one_more_dilation_fraction"
                        ],
                        "holdout_selected_parent_purity": record[
                            "holdout_selected_parent_cross_section_purity_within_global_observed"
                        ],
                        "coverage_passed": record["coverage_passed"],
                        "largest_tail_parent_fraction": (
                            largest[0][
                                "fraction_of_holdout_cross_section"
                            ]
                            if largest
                            else ""
                        ),
                    }
                )
            else:
                total = record["training_total"]
                footprint = record["parent_footprint"]
                writer.writerow(
                    {
                        "stratum_id": stratum_id,
                        "flat_index": record["flat_index"],
                        "status": record["status"],
                        "training_cross_section_microbarn": total[
                            "cross_section_microbarn"
                        ],
                        "training_sem_microbarn": total[
                            "cross_section_sem_microbarn"
                        ],
                        "training_ess": total[
                            "importance_effective_sample_size"
                        ],
                        "training_rows": total["contributing_rows"],
                        "selected_parent_fraction": record[
                            "estimated_selected_parent_fraction"
                        ],
                        "tail_fraction": record["estimated_tail_fraction"],
                        "one_more_dilation_fraction": record.get(
                            "one_more_dilation_fraction", ""
                        ),
                        "selected_parent_purity": record[
                            "selected_parent_cross_section_purity_within_global_observed"
                        ],
                        "seed_parent_components": len(
                            footprint["seed_component_ids"]
                        ),
                        "derived_parent_components": footprint[
                            "derived_component_count"
                        ],
                    }
                )


def _write_migration_csv(
    path: Path,
    *,
    manifest: dict,
    campaign: MigrationCampaign,
    validation: bool,
) -> None:
    config = manifest["analysis_config"]
    parent_grid = ParentGrid.from_config(config)
    targets = _target_lookup(config)
    prefix = "holdout" if validation else "training"
    fields = (
        "stratum_id",
        "observed_q2_index",
        "observed_xb_index",
        "observed_minus_t_index",
        "observed_phi_index",
        "component_id",
        "parent_id",
        "channel",
        "hard_q2_index",
        "hard_xb_index",
        "hard_minus_t_index",
        "hard_phi_index",
        "hard_q2_region",
        "hard_xb_region",
        "hard_minus_t_region",
        "delta_q2_index",
        "delta_xb_index",
        "delta_minus_t_index",
        "delta_phi_index",
        "relationship",
        "contributing_rows",
        f"{prefix}_sum_trial_contributions_microbarn",
        f"{prefix}_sum_squared_trial_contributions_microbarn2",
        f"{prefix}_largest_trial_contribution_microbarn",
        f"{prefix}_cross_section_microbarn",
        "fraction_of_stratum_cross_section",
        "selected_seed",
        "in_selected_parents",
        "in_one_more_dilation",
    )
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        for stratum_id, moments in sorted(campaign.migrations.items()):
            target = targets[stratum_id]
            record = manifest["strata"][stratum_id]
            seeds = {
                component_from_identifier(identifier, parent_grid)
                for identifier in record["parent_footprint"][
                    "seed_component_ids"
                ]
            }
            selected = _components_from_record(record, parent_grid)
            expanded = _components_from_record(
                record, parent_grid, extra_steps=1
            )
            total = campaign.strata[stratum_id].total
            for (parent_id, channel), moment in sorted(moments.items()):
                parent = parent_from_identifier(parent_id, parent_grid)
                metadata = parent_metadata(
                    parent, parent_grid
                )["coordinates"]
                relationship, deltas = parent_relationship(
                    parent, target, parent_grid
                )
                component = (parent_id, channel)
                writer.writerow(
                    {
                        "stratum_id": stratum_id,
                        "observed_q2_index": target.iq2,
                        "observed_xb_index": target.ixb,
                        "observed_minus_t_index": target.it,
                        "observed_phi_index": target.iphi,
                        "component_id": component_identifier(
                            parent_id, channel
                        ),
                        "parent_id": parent_id,
                        "channel": channel,
                        "hard_q2_index": parent.iq2,
                        "hard_xb_index": parent.ixb,
                        "hard_minus_t_index": parent.it,
                        "hard_phi_index": parent.iphi,
                        "hard_q2_region": metadata["Q2"]["region"],
                        "hard_xb_region": metadata["xB"]["region"],
                        "hard_minus_t_region": metadata["minus_t"][
                            "region"
                        ],
                        "delta_q2_index": (
                            "" if deltas["Q2"] is None else deltas["Q2"]
                        ),
                        "delta_xb_index": (
                            "" if deltas["xB"] is None else deltas["xB"]
                        ),
                        "delta_minus_t_index": (
                            ""
                            if deltas["minus_t"] is None
                            else deltas["minus_t"]
                        ),
                        "delta_phi_index": deltas["phi_deg"],
                        "relationship": relationship,
                        "contributing_rows": moment.count,
                        f"{prefix}_sum_trial_contributions_microbarn": (
                            moment.total
                        ),
                        f"{prefix}_sum_squared_trial_contributions_microbarn2": (
                            moment.square_total
                        ),
                        f"{prefix}_largest_trial_contribution_microbarn": (
                            moment.maximum
                        ),
                        f"{prefix}_cross_section_microbarn": (
                            moment.total / campaign.proposals
                        ),
                        "fraction_of_stratum_cross_section": (
                            moment.total / total if total > 0.0 else ""
                        ),
                        "selected_seed": int(component in seeds),
                        "in_selected_parents": int(component in selected),
                        "in_one_more_dilation": int(component in expanded),
                    }
                )


def learn(args: argparse.Namespace) -> dict:
    if args.output.resolve().exists():
        raise FileExistsError(
            f"{args.output.resolve()} already exists; "
            "migration artifacts are immutable"
        )
    if not 0.0 < args.target_parent_fraction <= 1.0:
        raise MigrationError("--target-parent-fraction must be in (0,1]")
    if args.parent_dilation < 0:
        raise MigrationError("--parent-dilation must be nonnegative")
    if args.minimum_training_rows < 1 or args.minimum_training_ess < 0.0:
        raise MigrationError("training support thresholds are invalid")
    try:
        config, config_sha = guards.load_analysis_config(
            args.config.resolve()
        )
    except guards.GuardLearningError as error:
        raise MigrationError(str(error)) from error
    campaign = aggregate_surveys(
        args.survey, config, apply_y_max=args.apply_y_max
    )
    if not math.isclose(
        float(config["beam_energy"]),
        float(campaign.norm_reference["ebeam"]),
        rel_tol=2.0e-7,
    ):
        raise MigrationError("analysis and survey beam energies differ")
    learner_revision = guards._current_revision()
    generator_revision = args.generator_revision or learner_revision
    manifest = build_manifest(
        config=config,
        config_path=args.config,
        config_sha256=config_sha,
        campaign=campaign,
        target_parent_fraction=args.target_parent_fraction,
        parent_dilation=args.parent_dilation,
        iteration=args.iteration,
        generator_revision=generator_revision,
        generator_revision_source=(
            "user_supplied"
            if args.generator_revision
            else "assumed_current_checkout"
        ),
        learner_revision=learner_revision,
        minimum_training_rows=args.minimum_training_rows,
        minimum_training_ess=args.minimum_training_ess,
        apply_y_max=args.apply_y_max,
    )
    path = _write_json_artifact(
        args.output.resolve(), "migration_manifest.json", manifest
    )
    coverage_path = path.parent / "training_parent_coverage.csv"
    migration_path = path.parent / "training_migrations.csv"
    _write_coverage_csv(coverage_path, manifest, validation=False)
    _write_migration_csv(
        migration_path,
        manifest=manifest,
        campaign=campaign,
        validation=False,
    )
    guards._write_sha256(coverage_path)
    guards._write_sha256(migration_path)
    return {
        "passed": True,
        "schema": MANIFEST_SCHEMA,
        "manifest": str(path),
        "manifest_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "training_replicas": manifest["training"]["replica_ids"],
        "total_proposals": campaign.proposals,
        **manifest["learning"],
    }


def validate(args: argparse.Namespace) -> dict:
    if args.output.resolve().exists():
        raise FileExistsError(
            f"{args.output.resolve()} already exists; "
            "migration artifacts are immutable"
        )
    if not 0.0 <= args.minimum_parent_coverage <= 1.0:
        raise MigrationError("--minimum-parent-coverage must be in [0,1]")
    manifest_path = args.manifest.resolve()
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise MigrationError(
            f"cannot read manifest {manifest_path}: {error}"
        ) from error
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise MigrationError(
            f"manifest schema is {manifest.get('schema')!r}, "
            f"expected {MANIFEST_SCHEMA}"
        )
    config = manifest["analysis_config"]
    campaign = aggregate_surveys(
        args.survey,
        config,
        apply_y_max=bool(
            manifest.get("analysis_selection", {}).get(
                "apply_y_max", False
            )
        ),
    )
    report = validate_manifest(
        manifest,
        manifest_path,
        campaign,
        args.minimum_parent_coverage,
    )
    path = _write_json_artifact(
        args.output.resolve(), "migration_validation.json", report
    )
    coverage_path = path.parent / "validation_parent_coverage.csv"
    migration_path = path.parent / "validation_migrations.csv"
    _write_coverage_csv(coverage_path, report, validation=True)
    _write_migration_csv(
        migration_path,
        manifest=manifest,
        campaign=campaign,
        validation=True,
    )
    guards._write_sha256(coverage_path)
    guards._write_sha256(migration_path)
    return {
        "passed": report["passed"],
        "schema": VALIDATION_SCHEMA,
        "validation": str(path),
        "validation_replicas": report["validation"]["replica_ids"],
        **report["coverage_summary"],
    }


def _aggregate_serialized_moments(
    records: Iterable[dict], key: str
) -> guards.Moment:
    return guards._aggregate_record_moments(records, key)


def _weighted_parent_coverage(
    records: list[dict],
    *,
    proposals: int,
    total_key: str,
    selected_key: str,
    tail_key: str,
    expanded_key: str,
) -> dict[str, object]:
    total = _aggregate_serialized_moments(records, total_key)
    selected = _aggregate_serialized_moments(records, selected_key)
    tail = _aggregate_serialized_moments(records, tail_key)
    expanded = _aggregate_serialized_moments(records, expanded_key)
    return {
        "contributing_strata": sum(
            1
            for record in records
            if float(
                record.get(total_key, {}).get(
                    "sum_trial_contributions_microbarn", 0.0
                )
            )
            > 0.0
        ),
        "total": guards._compact_metrics(total, proposals),
        "selected_parents": guards._compact_metrics(
            selected, proposals
        ),
        "tail": guards._compact_metrics(tail, proposals),
        "one_more_dilation": guards._compact_metrics(
            expanded, proposals
        ),
        "cross_section_weighted_selected_parent_fraction": (
            selected.total / total.total if total.total > 0.0 else None
        ),
        "cross_section_weighted_tail_fraction": (
            tail.total / total.total if total.total > 0.0 else None
        ),
        "cross_section_weighted_one_more_dilation_fraction": (
            expanded.total / total.total if total.total > 0.0 else None
        ),
        "extra_cross_section_fraction_recovered_by_one_more_dilation": (
            (expanded.total - selected.total) / total.total
            if total.total > 0.0
            else None
        ),
    }


def _material_parent_coverage(
    records: list[dict],
    *,
    proposals: int,
    targets: tuple[float, ...] = (0.5, 0.9, 0.95, 0.99),
) -> dict[str, object]:
    contributing = [
        record
        for record in records
        if float(
            record["holdout_total"][
                "sum_trial_contributions_microbarn"
            ]
        )
        > 0.0
    ]
    contributing.sort(
        key=lambda record: -float(
            record["holdout_total"][
                "sum_trial_contributions_microbarn"
            ]
        )
    )
    total = sum(
        float(
            record["holdout_total"][
                "sum_trial_contributions_microbarn"
            ]
        )
        for record in contributing
    )
    result: dict[str, object] = {}
    for target in targets:
        retained = 0.0
        selected_records: list[dict] = []
        for record in contributing:
            selected_records.append(record)
            retained += float(
                record["holdout_total"][
                    "sum_trial_contributions_microbarn"
                ]
            )
            if retained >= target * total:
                break
        coverage = _weighted_parent_coverage(
            selected_records,
            proposals=proposals,
            total_key="holdout_total",
            selected_key="holdout_selected_parents",
            tail_key="holdout_tail",
            expanded_key="holdout_one_more_dilation",
        )
        result[f"top_{100.0 * target:g}_percent_cross_section"] = {
            "strata": len(selected_records),
            "actual_fraction_of_inside_cross_section": (
                retained / total if total > 0.0 else None
            ),
            "cross_section_weighted_selected_parent_fraction": coverage[
                "cross_section_weighted_selected_parent_fraction"
            ],
            "cross_section_weighted_one_more_dilation_fraction": coverage[
                "cross_section_weighted_one_more_dilation_fraction"
            ],
        }
    return result


def summarize(args: argparse.Namespace) -> dict:
    if args.limit < 0:
        raise MigrationError("--limit must be nonnegative")
    manifest_path = args.manifest.resolve()
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise MigrationError(
            f"cannot read manifest {manifest_path}: {error}"
        ) from error
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise MigrationError(
            f"manifest schema is {manifest.get('schema')!r}, "
            f"expected {MANIFEST_SCHEMA}"
        )
    manifest_hash = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    training_records = list(manifest["strata"].values())
    training_proposals = int(manifest["training"]["total_proposals"])
    summary: dict[str, object] = {
        "manifest_schema": manifest["schema"],
        "manifest_sha256": manifest_hash,
        "analysis_selection": manifest["analysis_selection"],
        "hard_parent_definition": manifest["hard_parent_definition"],
        "training_replicas": manifest["training"]["replica_ids"],
        "learning": manifest["learning"],
        "weighted_parent_coverage": {
            "training_inside_analysis_partition": (
                _weighted_parent_coverage(
                    training_records,
                    proposals=training_proposals,
                    total_key="training_total",
                    selected_key="training_selected_parents",
                    tail_key="training_tail",
                    expanded_key="training_one_more_dilation",
                )
            )
        },
        "training_hard_parent_relationships": manifest["training"][
            "hard_parent_relationships_inside_analysis"
        ],
    }
    if not args.validation:
        return summary

    validation_path = args.validation.resolve()
    try:
        report = json.loads(validation_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise MigrationError(
            f"cannot read validation {validation_path}: {error}"
        ) from error
    if report.get("schema") != VALIDATION_SCHEMA:
        raise MigrationError(
            f"validation schema is {report.get('schema')!r}, "
            f"expected {VALIDATION_SCHEMA}"
        )
    if report.get("manifest_sha256") != manifest_hash:
        raise MigrationError(
            "validation artifact does not reference the supplied manifest"
        )
    validation_records = list(report["strata"].values())
    validation_proposals = int(report["validation"]["total_proposals"])
    weighted_validation = _weighted_parent_coverage(
        validation_records,
        proposals=validation_proposals,
        total_key="holdout_total",
        selected_key="holdout_selected_parents",
        tail_key="holdout_tail",
        expanded_key="holdout_one_more_dilation",
    )
    by_status = {
        status: _weighted_parent_coverage(
            [
                record
                for record in validation_records
                if record["training_status"] == status
            ],
            proposals=validation_proposals,
            total_key="holdout_total",
            selected_key="holdout_selected_parents",
            tail_key="holdout_tail",
            expanded_key="holdout_one_more_dilation",
        )
        for status in (
            "learned",
            "learned_low_support",
            "no_training_contribution",
        )
    }
    assessed = [
        (stratum_id, record)
        for stratum_id, record in report["strata"].items()
        if record["holdout_selected_parent_fraction"] is not None
    ]
    worst = sorted(
        assessed,
        key=lambda item: (
            item[1]["holdout_selected_parent_fraction"],
            item[0],
        ),
    )[: args.limit]
    training_inside = manifest["training"]["inside_analysis_partition"]
    validation_inside = report["validation"][
        "inside_analysis_partition"
    ]
    summary["weighted_parent_coverage"].update(
        {
            "validation_inside_analysis_partition": weighted_validation,
            "inside_training_holdout_relative_difference": (
                (
                    float(validation_inside["cross_section_microbarn"])
                    - float(training_inside["cross_section_microbarn"])
                )
                / float(training_inside["cross_section_microbarn"])
                if float(training_inside["cross_section_microbarn"]) > 0.0
                else None
            ),
            "inside_training_holdout_difference_z_score": (
                guards._difference_z_score(
                    float(training_inside["cross_section_microbarn"]),
                    training_inside["cross_section_sem_microbarn"],
                    float(validation_inside["cross_section_microbarn"]),
                    validation_inside["cross_section_sem_microbarn"],
                )
            ),
            "validation_by_training_status": by_status,
            "validation_material_strata": _material_parent_coverage(
                validation_records, proposals=validation_proposals
            ),
        }
    )
    summary.update(
        {
            "validation_schema": report["schema"],
            "validation_passed": report["passed"],
            "validation_replicas": report["validation"]["replica_ids"],
            "validation_hard_parent_relationships": report["validation"][
                "hard_parent_relationships_inside_analysis"
            ],
            "coverage_summary": report["coverage_summary"],
            "worst_holdout_parent_fractions": [
                {
                    "stratum_id": stratum_id,
                    "training_status": record["training_status"],
                    "holdout_cross_section_microbarn": record[
                        "holdout_total"
                    ]["cross_section_microbarn"],
                    "holdout_selected_parent_fraction": record[
                        "holdout_selected_parent_fraction"
                    ],
                    "holdout_one_more_dilation_fraction": record[
                        "holdout_one_more_dilation_fraction"
                    ],
                    "largest_tail_parent_components": record[
                        "largest_tail_parent_components"
                    ],
                }
                for stratum_id, record in worst
            ],
        }
    )
    return summary


def _load_plot_inputs(
    manifest_path: Path, validation_path: Path
) -> tuple[dict, dict, Path]:
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        report = json.loads(validation_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise MigrationError(f"cannot read plotting inputs: {error}") from error
    manifest_hash = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise MigrationError("plot manifest has the wrong schema")
    if report.get("schema") != VALIDATION_SCHEMA:
        raise MigrationError("plot validation has the wrong schema")
    if report.get("manifest_sha256") != manifest_hash:
        raise MigrationError(
            "validation artifact does not reference the supplied manifest"
        )
    migration_path = validation_path.parent / report["artifacts"][
        "validation_migrations"
    ]
    if not migration_path.is_file():
        raise MigrationError(f"missing validation migration CSV {migration_path}")
    return manifest, report, migration_path


def plot(args: argparse.Namespace) -> dict:
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(
            f"{output} already exists; migration plot artifacts are immutable"
        )
    if args.top_strata < 0:
        raise MigrationError("--top-strata must be nonnegative")
    manifest_path = args.manifest.resolve()
    validation_path = args.validation.resolve()
    manifest, report, migration_path = _load_plot_inputs(
        manifest_path, validation_path
    )
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError as error:
        raise MigrationError(
            "plotting requires matplotlib"
        ) from error

    output.mkdir(parents=True)
    pdf_path = output / "migration_diagnostics.pdf"
    records = [
        (stratum_id, record)
        for stratum_id, record in report["strata"].items()
        if float(record["holdout_total"]["cross_section_microbarn"]) > 0.0
    ]
    records.sort(
        key=lambda item: -float(
            item[1]["holdout_total"]["cross_section_microbarn"]
        )
    )
    selected_strata = {
        stratum_id for stratum_id, _ in records[: args.top_strata]
    }
    representative: dict[
        str, dict[tuple[int, int], float]
    ] = {stratum_id: {} for stratum_id in selected_strata}
    delta_qx: dict[tuple[int, int], float] = {}
    delta_t: dict[int, float] = {}
    delta_phi: dict[int, float] = {}
    with migration_path.open(encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source)
        for row in reader:
            contribution = float(
                row["holdout_sum_trial_contributions_microbarn"]
            )
            if row["delta_q2_index"] and row["delta_xb_index"]:
                key = (
                    int(row["delta_q2_index"]),
                    int(row["delta_xb_index"]),
                )
                delta_qx[key] = delta_qx.get(key, 0.0) + contribution
            if row["delta_minus_t_index"]:
                index = int(row["delta_minus_t_index"])
                delta_t[index] = delta_t.get(index, 0.0) + contribution
            if row["delta_phi_index"]:
                index = int(row["delta_phi_index"])
                delta_phi[index] = delta_phi.get(index, 0.0) + contribution
            if row["stratum_id"] in selected_strata:
                key = (
                    int(row["hard_q2_index"]),
                    int(row["hard_xb_index"]),
                )
                target = representative[row["stratum_id"]]
                target[key] = target.get(key, 0.0) + contribution

    with PdfPages(pdf_path) as pdf:
        figure, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        cross_sections = [
            float(record["holdout_total"]["cross_section_microbarn"])
            for _, record in records
        ]
        coverages = [
            float(record["holdout_selected_parent_fraction"])
            for _, record in records
        ]
        axes[0].scatter(cross_sections, coverages, s=8, alpha=0.5)
        axes[0].set_xscale("log")
        axes[0].set_ylim(-0.02, 1.02)
        axes[0].set_xlabel("Held-out stratum cross section [microbarn]")
        axes[0].set_ylabel("Selected-parent coverage")
        axes[0].grid(alpha=0.25)
        axes[1].hist(
            coverages,
            bins=20,
            range=(0.0, 1.0),
            weights=cross_sections,
        )
        axes[1].set_xlabel("Selected-parent coverage")
        axes[1].set_ylabel("Held-out cross section [microbarn]")
        axes[1].grid(alpha=0.25)
        figure.suptitle("Hard-parent migration coverage")
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        relationships = report["validation"][
            "hard_parent_relationships_inside_analysis"
        ]
        labels = list(relationships)
        values = [
            relationships[label][
                "fraction_of_inside_analysis_cross_section"
            ]
            or 0.0
            for label in labels
        ]
        figure, axis = plt.subplots(figsize=(9, 4.5))
        axis.bar(range(len(labels)), values)
        axis.set_xticks(
            range(len(labels)),
            [label.replace("_", "\n") for label in labels],
        )
        axis.set_ylabel("Fraction of inside-analysis cross section")
        axis.set_ylim(0.0, 1.0)
        axis.grid(axis="y", alpha=0.25)
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        if delta_qx:
            q_values = [key[0] for key in delta_qx]
            x_values = [key[1] for key in delta_qx]
            q_min, q_max = min(q_values), max(q_values)
            x_min, x_max = min(x_values), max(x_values)
            matrix = [
                [
                    delta_qx.get((q_index, x_index), 0.0)
                    for x_index in range(x_min, x_max + 1)
                ]
                for q_index in range(q_min, q_max + 1)
            ]
            figure, axis = plt.subplots(figsize=(7, 5.5))
            image = axis.imshow(
                matrix,
                origin="lower",
                aspect="auto",
                extent=(
                    x_min - 0.5,
                    x_max + 0.5,
                    q_min - 0.5,
                    q_max + 0.5,
                ),
            )
            axis.set_xlabel("hard xB bin − observed xB bin")
            axis.set_ylabel("hard Q² bin − observed Q² bin")
            axis.set_title("Cross-section-weighted migration offsets")
            figure.colorbar(
                image, ax=axis, label="sum of trial contributions"
            )
            figure.tight_layout()
            pdf.savefig(figure)
            plt.close(figure)

        figure, axes = plt.subplots(1, 2, figsize=(10, 4.5))
        for axis, values_map, label in (
            (axes[0], delta_t, "hard −t bin − observed −t bin"),
            (axes[1], delta_phi, "wrapped hard φ bin − observed φ bin"),
        ):
            keys = sorted(values_map)
            axis.bar(keys, [values_map[key] for key in keys])
            axis.set_xlabel(label)
            axis.set_ylabel("sum of trial contributions")
            axis.grid(axis="y", alpha=0.25)
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        for stratum_id, _ in records[: args.top_strata]:
            values_map = representative[stratum_id]
            if not values_map:
                continue
            q_values = [key[0] for key in values_map]
            x_values = [key[1] for key in values_map]
            q_min, q_max = min(q_values), max(q_values)
            x_min, x_max = min(x_values), max(x_values)
            matrix = [
                [
                    values_map.get((q_index, x_index), 0.0)
                    for x_index in range(x_min, x_max + 1)
                ]
                for q_index in range(q_min, q_max + 1)
            ]
            figure, axis = plt.subplots(figsize=(7, 5.5))
            image = axis.imshow(
                matrix,
                origin="lower",
                aspect="auto",
                extent=(
                    x_min - 0.5,
                    x_max + 0.5,
                    q_min - 0.5,
                    q_max + 0.5,
                ),
            )
            target = manifest["strata"][stratum_id]["indices"]
            axis.scatter(
                [target["xB"]],
                [target["Q2"]],
                marker="x",
                color="red",
                label="same Born parent",
            )
            axis.set_xlabel("hard xB-bin index")
            axis.set_ylabel("hard Q²-bin index")
            axis.set_title(f"{stratum_id}: held-out hard-parent footprint")
            axis.legend()
            figure.colorbar(
                image, ax=axis, label="sum of trial contributions"
            )
            figure.tight_layout()
            pdf.savefig(figure)
            plt.close(figure)

    guards._write_sha256(pdf_path)
    plot_summary = {
        "schema": PLOT_SCHEMA,
        "created_utc": _now(),
        "manifest": str(manifest_path),
        "manifest_sha256": hashlib.sha256(
            manifest_path.read_bytes()
        ).hexdigest(),
        "validation": str(validation_path),
        "validation_sha256": hashlib.sha256(
            validation_path.read_bytes()
        ).hexdigest(),
        "pdf": pdf_path.name,
        "representative_strata": [
            stratum_id for stratum_id, _ in records[: args.top_strata]
        ],
    }
    summary_path = output / "plot_summary.json"
    summary_path.write_text(
        json.dumps(plot_summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    guards._write_sha256(summary_path)
    return {
        "passed": True,
        "schema": PLOT_SCHEMA,
        "output": str(output),
        "pdf": str(pdf_path),
        "representative_strata": plot_summary["representative_strata"],
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    learn_parser = subparsers.add_parser(
        "learn", help="learn hard-parent footprints from training replicas"
    )
    learn_parser.add_argument("--config", type=Path, required=True)
    learn_parser.add_argument("--survey", type=Path, nargs="+", required=True)
    learn_parser.add_argument("--output", type=Path, required=True)
    learn_parser.add_argument(
        "--target-parent-fraction", type=float, default=0.995
    )
    learn_parser.add_argument("--parent-dilation", type=int, default=0)
    learn_parser.add_argument("--iteration", type=int, default=0)
    learn_parser.add_argument(
        "--minimum-training-rows", type=int, default=10
    )
    learn_parser.add_argument(
        "--minimum-training-ess", type=float, default=5.0
    )
    learn_parser.add_argument(
        "--apply-y-max",
        action="store_true",
        help=(
            "apply phase_space.y_max when assigning final observed strata; "
            "the default applies no y cut"
        ),
    )
    learn_parser.add_argument(
        "--generator-revision",
        help=(
            "Git revision that produced the surveys; defaults to the current "
            "checkout and is recorded as an assumption"
        ),
    )

    validate_parser = subparsers.add_parser(
        "validate", help="validate frozen parent footprints on held-out replicas"
    )
    validate_parser.add_argument("--manifest", type=Path, required=True)
    validate_parser.add_argument(
        "--survey", type=Path, nargs="+", required=True
    )
    validate_parser.add_argument("--output", type=Path, required=True)
    validate_parser.add_argument(
        "--minimum-parent-coverage", type=float, default=0.98
    )

    summarize_parser = subparsers.add_parser(
        "summarize", help="print compact weighted migration coverage"
    )
    summarize_parser.add_argument("--manifest", type=Path, required=True)
    summarize_parser.add_argument("--validation", type=Path)
    summarize_parser.add_argument(
        "--limit",
        type=int,
        default=20,
        help="number of lowest-coverage validation strata to show",
    )

    plot_parser = subparsers.add_parser(
        "plot", help="render migration coverage and hard-parent diagnostics"
    )
    plot_parser.add_argument("--manifest", type=Path, required=True)
    plot_parser.add_argument("--validation", type=Path, required=True)
    plot_parser.add_argument("--output", type=Path, required=True)
    plot_parser.add_argument(
        "--top-strata",
        type=int,
        default=6,
        help="number of highest-cross-section strata to plot individually",
    )
    return parser


def main() -> int:
    args = _build_parser().parse_args()
    try:
        if args.command == "learn":
            result = learn(args)
        elif args.command == "validate":
            result = validate(args)
        elif args.command == "summarize":
            result = summarize(args)
        else:
            result = plot(args)
    except (
        FileNotFoundError,
        FileExistsError,
        KeyError,
        ValueError,
        MigrationError,
    ) as error:
        print(f"ERROR: {error}")
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result.get("passed", True) else 2


if __name__ == "__main__":
    raise SystemExit(main())
