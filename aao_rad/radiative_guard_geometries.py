#!/usr/bin/env python3
"""Compare training-only core-guard geometries for radiative mode 4.

This milestone-2d diagnostic fixes the spatial representation to the
channel-marginalized hard-parent grid selected by milestone 2c.  It compares
independent per-stratum footprints with reproducible neighbor-offset and
boundary-template fallbacks.  It does not activate mode 4, perform
acceptance-rejection, or emit LUND events.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import radiative_guards as guards
import radiative_migrations as migrations
import radiative_survey


COMPARISON_SCHEMA = "aao-rad-guard-geometry-comparison-v1"
RECIPE_SCHEMA = "aao-rad-guard-geometry-recipes-v1"
PLOT_SCHEMA = "aao-rad-guard-geometry-plots-v1"
STATUS_ORDER = (
    "learned",
    "learned_low_support",
    "no_training_contribution",
)
NATIVE_CHANNELS = tuple(range(1, 7))


class GeometryError(RuntimeError):
    """Raised when guard geometries cannot be compared reproducibly."""


@dataclass(frozen=True)
class GeometryCandidate:
    identifier: str
    label: str
    description: str
    construction: str
    dilation_axes: tuple[str, ...] = ()

    def metadata(self) -> dict[str, object]:
        return {
            "identifier": self.identifier,
            "label": self.label,
            "description": self.description,
            "construction": self.construction,
            "dilation_axes": list(self.dilation_axes),
        }


GEOMETRY_CANDIDATES = (
    GeometryCandidate(
        identifier="independent_frozen",
        label="Independent frozen",
        description=(
            "The channel-marginalized seed footprint learned independently "
            "for each final-LUND analysis stratum."
        ),
        construction="own_training_seed_parents",
    ),
    GeometryCandidate(
        identifier="independent_all_axis_dilation",
        label="Independent + all-axis dilation",
        description=(
            "The independent footprint plus one Manhattan step along hard "
            "Q2, xB, minus-t, or periodic phi."
        ),
        construction="own_training_seed_parents_then_axis_dilation",
        dilation_axes=("Q2", "xB", "minus_t", "phi_deg"),
    ),
    GeometryCandidate(
        identifier="hierarchical_neighbor_offsets",
        label="Hierarchical neighbor offsets",
        description=(
            "Supported strata keep their own footprint. Low-support and "
            "empty strata borrow translated hard-minus-observed offsets from "
            "nearby supported strata and use boundary/global templates only "
            "when required."
        ),
        construction="support_adaptive_translated_offset_fallback",
    ),
    GeometryCandidate(
        identifier="hierarchical_neighbor_offsets_qx_dilation",
        label="Hierarchical + Q2/xB dilation",
        description=(
            "The hierarchical fallback followed by one hard-cell Manhattan "
            "step along Q2 or xB only."
        ),
        construction=(
            "support_adaptive_translated_offset_fallback_then_axis_dilation"
        ),
        dilation_axes=("Q2", "xB"),
    ),
)


def _now() -> str:
    return (
        dt.datetime.now(dt.timezone.utc)
        .replace(microsecond=0)
        .isoformat()
    )


def _validate_candidates() -> None:
    identifiers = [candidate.identifier for candidate in GEOMETRY_CANDIDATES]
    if len(identifiers) != len(set(identifiers)):
        raise GeometryError("guard-geometry candidate IDs must be unique")
    allowed_axes = {"Q2", "xB", "minus_t", "phi_deg"}
    for candidate in GEOMETRY_CANDIDATES:
        if not set(candidate.dilation_axes) <= allowed_axes:
            raise GeometryError(
                f"{candidate.identifier}: invalid dilation axis"
            )


def _status(
    total: guards.Moment,
    *,
    minimum_training_rows: int,
    minimum_training_ess: float,
) -> str:
    if total.total <= 0.0:
        return "no_training_contribution"
    supported = (
        total.count >= minimum_training_rows
        and guards._ess(total) >= minimum_training_ess
    )
    return "learned" if supported else "learned_low_support"


def _parent_moments(
    moments: dict[tuple[str, int], guards.Moment],
) -> dict[str, guards.Moment]:
    result: dict[str, guards.Moment] = {}
    for (parent_id, _), moment in moments.items():
        target = result.setdefault(parent_id, guards.Moment())
        migrations._add_moment(target, moment)
    return result


def _select_weighted_keys(
    moments: dict[object, guards.Moment],
    target_fraction: float,
) -> set[object]:
    total = sum(moment.total for moment in moments.values())
    if total <= 0.0:
        return set()
    selected: set[object] = set()
    retained = 0.0
    for key, moment in sorted(
        moments.items(), key=lambda item: (-item[1].total, str(item[0]))
    ):
        if moment.total <= 0.0:
            continue
        selected.add(key)
        retained += moment.total
        if retained >= target_fraction * total:
            break
    return selected


def _target_shape(config: dict) -> tuple[int, int, int, int]:
    binning = config["binning"]
    return (
        len(binning["Q2"]) - 1,
        len(binning["xB"]) - 1,
        len(binning["minus_t"]) - 1,
        len(binning["phi_deg"]) - 1,
    )


def _stratum_key(stratum: guards.Stratum) -> tuple[int, int, int, int]:
    return (stratum.iq2, stratum.ixb, stratum.it, stratum.iphi)


def _observed_neighbor_ids(
    stratum: guards.Stratum,
    lookup: dict[tuple[int, int, int, int], guards.Stratum],
    shape: tuple[int, int, int, int],
    radius: int,
) -> tuple[str, ...]:
    if radius < 0:
        raise GeometryError("observed-neighbor radius must be nonnegative")
    origin = _stratum_key(stratum)
    reached = {origin}
    frontier = {origin}
    for _ in range(radius):
        additions: set[tuple[int, int, int, int]] = set()
        for indices in frontier:
            for dimension in range(4):
                for shift in (-1, 1):
                    neighbor = list(indices)
                    candidate = neighbor[dimension] + shift
                    if dimension == 3:
                        candidate %= shape[dimension]
                    elif not 0 <= candidate < shape[dimension]:
                        continue
                    neighbor[dimension] = candidate
                    additions.add(tuple(neighbor))
        additions -= reached
        reached.update(additions)
        frontier = additions
        if not frontier:
            break
    reached.discard(origin)
    return tuple(
        lookup[key].identifier
        for key in sorted(reached)
        if key in lookup
    )


def _signed_phi_delta(left: int, right: int, bins: int) -> int:
    return (left - right + bins // 2) % bins - bins // 2


def _relative_offset(
    parent_id: str,
    target: guards.Stratum,
    grid: migrations.ParentGrid,
) -> tuple[int, int, int, int]:
    parent = migrations.parent_from_identifier(parent_id, grid)
    return (
        parent.iq2 - target.iq2,
        parent.ixb - target.ixb,
        parent.it - target.it,
        _signed_phi_delta(parent.iphi, target.iphi, grid.shape[3]),
    )


def _translated_parent_id(
    offset: tuple[int, int, int, int],
    target: guards.Stratum,
    grid: migrations.ParentGrid,
) -> str | None:
    nq2 = len(grid.q2_edges) - 1
    nxb = len(grid.xb_edges) - 1
    nt = len(grid.minus_t_edges) - 1
    nphi = len(grid.phi_edges) - 1
    indices = (
        target.iq2 + offset[0],
        target.ixb + offset[1],
        target.it + offset[2],
        (target.iphi + offset[3]) % nphi,
    )
    if not (
        -1 <= indices[0] <= nq2
        and -1 <= indices[1] <= nxb
        and -1 <= indices[2] <= nt
    ):
        return None
    parent = migrations.ParentIndex(*indices)
    return migrations.parent_identifier(parent, grid)


def _translated_parent_ids(
    offsets: Iterable[tuple[int, int, int, int]],
    target: guards.Stratum,
    grid: migrations.ParentGrid,
) -> set[str]:
    result: set[str] = set()
    for offset in offsets:
        parent_id = _translated_parent_id(offset, target, grid)
        if parent_id is not None:
            result.add(parent_id)
    return result


def _dilate_parent_ids(
    seeds: set[str],
    grid: migrations.ParentGrid,
    axes: tuple[str, ...],
    radius: int = 1,
) -> set[str]:
    if radius < 0:
        raise GeometryError("hard-parent dilation must be nonnegative")
    axis_indices = {
        "Q2": 0,
        "xB": 1,
        "minus_t": 2,
        "phi_deg": 3,
    }
    dimensions = tuple(axis_indices[name] for name in axes)
    result = set(seeds)
    frontier = set(seeds)
    regular = (
        len(grid.q2_edges) - 1,
        len(grid.xb_edges) - 1,
        len(grid.minus_t_edges) - 1,
        len(grid.phi_edges) - 1,
    )
    for _ in range(radius):
        additions: set[str] = set()
        for parent_id in frontier:
            parent = migrations.parent_from_identifier(parent_id, grid)
            indices = [parent.iq2, parent.ixb, parent.it, parent.iphi]
            for dimension in dimensions:
                for shift in (-1, 1):
                    neighbor = indices.copy()
                    candidate = neighbor[dimension] + shift
                    if dimension == 3:
                        candidate %= regular[dimension]
                    elif not -1 <= candidate <= regular[dimension]:
                        continue
                    neighbor[dimension] = candidate
                    additions.add(
                        migrations.parent_identifier(
                            migrations.ParentIndex(*neighbor), grid
                        )
                    )
        additions -= result
        result.update(additions)
        frontier = additions
        if not frontier:
            break
    return result


def _boundary_signature(
    stratum: guards.Stratum,
    grid: migrations.ParentGrid,
) -> tuple[str, str, str]:
    return (
        migrations._target_boundary_position(
            stratum.iq2, len(grid.q2_edges) - 1
        ),
        migrations._target_boundary_position(
            stratum.ixb, len(grid.xb_edges) - 1
        ),
        migrations._target_boundary_position(
            stratum.it, len(grid.minus_t_edges) - 1
        ),
    )


def _boundary_identifier(signature: tuple[str, str, str]) -> str:
    return "|".join(
        f"{axis}={position}"
        for axis, position in zip(("Q2", "xB", "minus_t"), signature)
    )


def _serialize_offset(
    offset: tuple[int, int, int, int],
) -> dict[str, int]:
    return {
        "delta_q2_index": offset[0],
        "delta_xb_index": offset[1],
        "delta_minus_t_index": offset[2],
        "delta_phi_index": offset[3],
    }


def _offset_moments(
    catalog: list[guards.Stratum],
    campaign: migrations.MigrationCampaign,
    statuses: dict[str, str],
    grid: migrations.ParentGrid,
    *,
    allowed_statuses: set[str],
) -> tuple[
    dict[tuple[int, int, int, int], guards.Moment],
    dict[
        tuple[str, str, str],
        dict[tuple[int, int, int, int], guards.Moment],
    ],
]:
    global_moments: dict[
        tuple[int, int, int, int], guards.Moment
    ] = {}
    boundary_moments: dict[
        tuple[str, str, str],
        dict[tuple[int, int, int, int], guards.Moment],
    ] = {}
    for stratum in catalog:
        if statuses[stratum.identifier] not in allowed_statuses:
            continue
        signature = _boundary_signature(stratum, grid)
        target_boundary = boundary_moments.setdefault(signature, {})
        for parent_id, moment in _parent_moments(
            campaign.migrations.get(stratum.identifier, {})
        ).items():
            offset = _relative_offset(parent_id, stratum, grid)
            migrations._add_moment(
                global_moments.setdefault(offset, guards.Moment()),
                moment,
            )
            migrations._add_moment(
                target_boundary.setdefault(offset, guards.Moment()),
                moment,
            )
    return global_moments, boundary_moments


def _empty_core_accumulator() -> dict[str, guards.Moment]:
    return {
        "total": guards.Moment(),
        "selected": guards.Moment(),
        "tail": guards.Moment(),
    }


def _accumulate_core(
    accumulator: dict[str, guards.Moment],
    *,
    total: guards.Moment,
    selected: guards.Moment,
    tail: guards.Moment,
) -> None:
    for name, moment in (
        ("total", total),
        ("selected", selected),
        ("tail", tail),
    ):
        migrations._add_moment(accumulator[name], moment)


def _core_metrics(
    accumulator: dict[str, guards.Moment],
    proposals: int,
) -> dict[str, object]:
    total = accumulator["total"]
    selected = accumulator["selected"]
    return {
        "total": guards._compact_metrics(total, proposals),
        "selected_core": guards._compact_metrics(selected, proposals),
        "residual_outside_core": guards._compact_metrics(
            accumulator["tail"], proposals
        ),
        "cross_section_weighted_core_coverage": (
            selected.total / total.total if total.total > 0.0 else None
        ),
        "cross_section_weighted_residual_fraction": (
            accumulator["tail"].total / total.total
            if total.total > 0.0
            else None
        ),
    }


def _material_coverage(
    records: list[dict[str, object]],
    candidate: str,
    *,
    targets: tuple[float, ...] = (0.5, 0.9, 0.95, 0.99),
) -> dict[str, object]:
    contributing = [
        record
        for record in records
        if float(record["validation_total"]) > 0.0
    ]
    contributing.sort(
        key=lambda record: -float(record["validation_total"])
    )
    total = sum(
        float(record["validation_total"]) for record in contributing
    )
    result: dict[str, object] = {}
    for target in targets:
        retained = selected = 0.0
        count = 0
        for record in contributing:
            count += 1
            retained += float(record["validation_total"])
            selected += float(record["candidates"][candidate])
            if retained >= target * total:
                break
        result[f"top_{100.0 * target:g}_percent_cross_section"] = {
            "strata": count,
            "actual_fraction_of_inside_cross_section": (
                retained / total if total > 0.0 else None
            ),
            "cross_section_weighted_core_coverage": (
                selected / retained if retained > 0.0 else None
            ),
        }
    return result


def _native_components(parent_ids: set[str]) -> set[tuple[str, int]]:
    return {
        (parent_id, channel)
        for parent_id in parent_ids
        for channel in NATIVE_CHANNELS
    }


def _candidate_core(
    candidate: GeometryCandidate,
    *,
    stratum: guards.Stratum,
    status: str,
    own_seed_parents: set[str],
    neighbor_offsets: set[tuple[int, int, int, int]],
    boundary_offsets: set[tuple[int, int, int, int]],
    global_offsets: set[tuple[int, int, int, int]],
    grid: migrations.ParentGrid,
) -> tuple[set[str], dict[str, object]]:
    selected = set(own_seed_parents)
    trace: dict[str, object] = {
        "used_neighbor_offsets": False,
        "used_boundary_template": False,
        "used_global_template": False,
        "dilation_axes": list(candidate.dilation_axes),
    }
    if candidate.identifier == "independent_frozen":
        return selected, trace
    if candidate.identifier == "independent_all_axis_dilation":
        return (
            _dilate_parent_ids(
                selected, grid, candidate.dilation_axes, radius=1
            ),
            trace,
        )

    if status != "learned":
        fallback_added = False
        translated_neighbors = _translated_parent_ids(
            neighbor_offsets, stratum, grid
        )
        if translated_neighbors:
            selected.update(translated_neighbors)
            trace["used_neighbor_offsets"] = True
            fallback_added = True
        if status == "no_training_contribution" or not translated_neighbors:
            translated_boundary = _translated_parent_ids(
                boundary_offsets, stratum, grid
            )
            if translated_boundary:
                selected.update(translated_boundary)
                trace["used_boundary_template"] = True
                fallback_added = True
        if not fallback_added:
            translated_global = _translated_parent_ids(
                global_offsets, stratum, grid
            )
            if translated_global:
                selected.update(translated_global)
                trace["used_global_template"] = True

    if candidate.dilation_axes:
        selected = _dilate_parent_ids(
            selected, grid, candidate.dilation_axes, radius=1
        )
    return selected, trace


def _build_geometry_study(
    *,
    config: dict,
    config_path: Path,
    config_sha256: str,
    training: migrations.MigrationCampaign,
    validation: migrations.MigrationCampaign,
    target_parent_fraction: float,
    neighbor_radius: int,
    iteration: int,
    minimum_training_rows: int,
    minimum_training_ess: float,
    minimum_parent_coverage: float,
    apply_y_max: bool,
    generator_revision: str,
    generator_revision_source: str,
    learner_revision: str,
) -> tuple[dict, dict, list[dict[str, object]]]:
    _validate_candidates()
    grid = migrations.ParentGrid.from_config(config)
    catalog = guards.enumerate_strata(config)
    target_lookup = {
        _stratum_key(stratum): stratum for stratum in catalog
    }
    catalog_by_id = {
        stratum.identifier: stratum for stratum in catalog
    }
    shape = _target_shape(config)

    statuses: dict[str, str] = {}
    own_seed_parents: dict[str, set[str]] = {}
    for stratum in catalog:
        total = training.strata.get(
            stratum.identifier, guards.Moment()
        )
        statuses[stratum.identifier] = _status(
            total,
            minimum_training_rows=minimum_training_rows,
            minimum_training_ess=minimum_training_ess,
        )
        parent_moments = _parent_moments(
            training.migrations.get(stratum.identifier, {})
        )
        own_seed_parents[stratum.identifier] = {
            str(value)
            for value in _select_weighted_keys(
                parent_moments, target_parent_fraction
            )
        }

    supported_global, supported_boundary = _offset_moments(
        catalog,
        training,
        statuses,
        grid,
        allowed_statuses={"learned"},
    )
    template_training_source = "well_supported_strata"
    if not supported_global:
        supported_global, supported_boundary = _offset_moments(
            catalog,
            training,
            statuses,
            grid,
            allowed_statuses={"learned_low_support"},
        )
        template_training_source = "low_support_fallback"
    global_offsets = {
        tuple(value)
        for value in _select_weighted_keys(
            supported_global, target_parent_fraction
        )
    }
    boundary_offsets = {
        signature: {
            tuple(value)
            for value in _select_weighted_keys(
                moments, target_parent_fraction
            )
        }
        for signature, moments in supported_boundary.items()
    }

    candidates = {
        candidate.identifier: candidate
        for candidate in GEOMETRY_CANDIDATES
    }
    accumulators: dict[str, dict[str, object]] = {}
    for identifier in candidates:
        accumulators[identifier] = {
            "training": _empty_core_accumulator(),
            "validation": _empty_core_accumulator(),
            "validation_by_status": {
                status: _empty_core_accumulator()
                for status in STATUS_ORDER
            },
            "validation_by_channel": {
                channel: _empty_core_accumulator()
                for channel in NATIVE_CHANNELS
            },
            "training_purity_numerator": 0.0,
            "training_purity_denominator": 0.0,
            "validation_purity_numerator": 0.0,
            "validation_purity_denominator": 0.0,
            "total_hard_cells": 0,
            "total_native_components": 0,
            "nonempty_core_strata": 0,
            "assessed_strata": 0,
            "passed_strata": 0,
            "failed_strata": 0,
            "source_counts": {
                "used_neighbor_offsets": 0,
                "used_boundary_template": 0,
                "used_global_template": 0,
            },
        }

    recipes: dict[str, object] = {}
    rows: list[dict[str, object]] = []
    material_records: list[dict[str, object]] = []
    status_counts = {status: 0 for status in STATUS_ORDER}

    for stratum in catalog:
        stratum_id = stratum.identifier
        status = statuses[stratum_id]
        status_counts[status] += 1
        training_total = training.strata.get(
            stratum_id, guards.Moment()
        )
        validation_total = validation.strata.get(
            stratum_id, guards.Moment()
        )
        training_moments = training.migrations.get(stratum_id, {})
        validation_moments = validation.migrations.get(stratum_id, {})
        neighbor_ids = _observed_neighbor_ids(
            stratum, target_lookup, shape, neighbor_radius
        )
        supported_neighbors = tuple(
            identifier
            for identifier in neighbor_ids
            if statuses[identifier] == "learned"
        )
        fallback_neighbors = (
            supported_neighbors
            if supported_neighbors
            else tuple(
                identifier
                for identifier in neighbor_ids
                if statuses[identifier] == "learned_low_support"
            )
        )
        neighbor_training_source = (
            "well_supported"
            if supported_neighbors
            else (
                "low_support_fallback"
                if fallback_neighbors
                else "none"
            )
        )
        neighbor_offsets: set[tuple[int, int, int, int]] = set()
        for donor_id in fallback_neighbors:
            donor = catalog_by_id[donor_id]
            neighbor_offsets.update(
                _relative_offset(parent_id, donor, grid)
                for parent_id in own_seed_parents[donor_id]
            )
        signature = _boundary_signature(stratum, grid)
        target_boundary_offsets = boundary_offsets.get(signature, set())

        recipe_record: dict[str, object] = {
            **guards._stratum_metadata(stratum),
            "training_status": status,
            "training_total": guards._compact_metrics(
                training_total, training.proposals
            ),
            "own_seed_parent_ids": sorted(
                own_seed_parents[stratum_id]
            ),
            "observed_neighbor_radius": neighbor_radius,
            "supported_neighbor_strata": list(supported_neighbors),
            "fallback_neighbor_strata": list(fallback_neighbors),
            "neighbor_training_source": neighbor_training_source,
            "boundary_template": _boundary_identifier(signature),
            "boundary_template_available": (
                signature in boundary_offsets
            ),
            "candidate_core_counts": {},
            "candidate_construction_trace": {},
        }
        material_record: dict[str, object] = {
            "stratum_id": stratum_id,
            "validation_total": validation_total.total,
            "candidates": {},
        }

        for identifier, candidate in candidates.items():
            selected_parents, trace = _candidate_core(
                candidate,
                stratum=stratum,
                status=status,
                own_seed_parents=own_seed_parents[stratum_id],
                neighbor_offsets=neighbor_offsets,
                boundary_offsets=target_boundary_offsets,
                global_offsets=global_offsets,
                grid=grid,
            )
            selected_native = _native_components(selected_parents)
            training_selected = migrations._subset(
                training_moments, selected_native
            )
            training_tail = migrations._combined(
                moment
                for component, moment in training_moments.items()
                if component not in selected_native
            )
            validation_selected = migrations._subset(
                validation_moments, selected_native
            )
            validation_tail = migrations._combined(
                moment
                for component, moment in validation_moments.items()
                if component not in selected_native
            )

            accumulator = accumulators[identifier]
            _accumulate_core(
                accumulator["training"],
                total=training_total,
                selected=training_selected,
                tail=training_tail,
            )
            _accumulate_core(
                accumulator["validation"],
                total=validation_total,
                selected=validation_selected,
                tail=validation_tail,
            )
            _accumulate_core(
                accumulator["validation_by_status"][status],
                total=validation_total,
                selected=validation_selected,
                tail=validation_tail,
            )
            for channel in NATIVE_CHANNELS:
                channel_total = migrations._combined(
                    moment
                    for (_, native_channel), moment
                    in validation_moments.items()
                    if native_channel == channel
                )
                channel_selected = migrations._combined(
                    moment
                    for (parent_id, native_channel), moment
                    in validation_moments.items()
                    if native_channel == channel
                    and parent_id in selected_parents
                )
                channel_tail = migrations._combined(
                    moment
                    for (parent_id, native_channel), moment
                    in validation_moments.items()
                    if native_channel == channel
                    and parent_id not in selected_parents
                )
                _accumulate_core(
                    accumulator["validation_by_channel"][channel],
                    total=channel_total,
                    selected=channel_selected,
                    tail=channel_tail,
                )

            training_denominator = migrations._purity_denominator(
                training.parent_totals, selected_native
            )
            validation_denominator = migrations._purity_denominator(
                validation.parent_totals, selected_native
            )
            accumulator["training_purity_numerator"] += (
                training_selected.total
            )
            accumulator["training_purity_denominator"] += (
                training_denominator
            )
            accumulator["validation_purity_numerator"] += (
                validation_selected.total
            )
            accumulator["validation_purity_denominator"] += (
                validation_denominator
            )
            accumulator["total_hard_cells"] += len(selected_parents)
            accumulator["total_native_components"] += len(selected_native)
            accumulator["nonempty_core_strata"] += int(
                bool(selected_parents)
            )
            for source_name in accumulator["source_counts"]:
                accumulator["source_counts"][source_name] += int(
                    bool(trace[source_name])
                )

            if validation_total.total > 0.0:
                validation_fraction = (
                    validation_selected.total / validation_total.total
                )
                passed = (
                    validation_fraction >= minimum_parent_coverage
                )
                accumulator["assessed_strata"] += 1
                accumulator["passed_strata"] += int(passed)
                accumulator["failed_strata"] += int(not passed)
            else:
                validation_fraction = None
                passed = None
            training_fraction = (
                training_selected.total / training_total.total
                if training_total.total > 0.0
                else None
            )
            validation_purity = (
                validation_selected.total / validation_denominator
                if validation_denominator > 0.0
                else None
            )

            recipe_record["candidate_core_counts"][identifier] = {
                "hard_cells": len(selected_parents),
                "native_components": len(selected_native),
            }
            recipe_record["candidate_construction_trace"][
                identifier
            ] = trace
            material_record["candidates"][identifier] = (
                validation_selected.total
            )
            rows.append(
                {
                    "stratum_id": stratum_id,
                    "flat_index": stratum.flat_index,
                    "training_status": status,
                    "candidate": identifier,
                    "training_cross_section_microbarn": (
                        training_total.total / training.proposals
                    ),
                    "training_ess": guards._ess(training_total),
                    "validation_cross_section_microbarn": (
                        validation_total.total / validation.proposals
                    ),
                    "validation_ess": guards._ess(validation_total),
                    "own_seed_hard_cells": len(
                        own_seed_parents[stratum_id]
                    ),
                    "neighbor_donor_strata": len(fallback_neighbors),
                    "used_neighbor_offsets": trace[
                        "used_neighbor_offsets"
                    ],
                    "used_boundary_template": trace[
                        "used_boundary_template"
                    ],
                    "used_global_template": trace[
                        "used_global_template"
                    ],
                    "selected_hard_cells": len(selected_parents),
                    "selected_native_components": len(selected_native),
                    "training_core_coverage": training_fraction,
                    "validation_core_coverage": validation_fraction,
                    "validation_core_purity_proxy": validation_purity,
                    "coverage_passed": passed,
                }
            )
        recipes[stratum_id] = recipe_record
        material_records.append(material_record)

    candidate_results: dict[str, object] = {}
    validation_inside = migrations._combined(
        validation.strata.values()
    ).total
    for identifier, candidate in candidates.items():
        accumulator = accumulators[identifier]
        training_metrics = _core_metrics(
            accumulator["training"], training.proposals
        )
        validation_metrics = _core_metrics(
            accumulator["validation"], validation.proposals
        )
        training_purity_denominator = float(
            accumulator["training_purity_denominator"]
        )
        validation_purity_denominator = float(
            accumulator["validation_purity_denominator"]
        )
        validation_by_channel: dict[str, object] = {}
        for channel in NATIVE_CHANNELS:
            channel_metrics = _core_metrics(
                accumulator["validation_by_channel"][channel],
                validation.proposals,
            )
            channel_total = accumulator["validation_by_channel"][
                channel
            ]["total"].total
            channel_metrics[
                "fraction_of_inside_analysis_cross_section"
            ] = (
                channel_total / validation_inside
                if validation_inside > 0.0
                else None
            )
            validation_by_channel[f"intreg_{channel}"] = channel_metrics

        candidate_results[identifier] = {
            **candidate.metadata(),
            "training": {
                **training_metrics,
                "aggregate_core_purity_proxy": (
                    float(accumulator["training_purity_numerator"])
                    / training_purity_denominator
                    if training_purity_denominator > 0.0
                    else None
                ),
            },
            "validation": {
                **validation_metrics,
                "aggregate_core_purity_proxy": (
                    float(accumulator["validation_purity_numerator"])
                    / validation_purity_denominator
                    if validation_purity_denominator > 0.0
                    else None
                ),
            },
            "validation_by_training_status": {
                status: _core_metrics(
                    accumulator["validation_by_status"][status],
                    validation.proposals,
                )
                for status in STATUS_ORDER
            },
            "validation_by_native_intreg": validation_by_channel,
            "validation_material_strata": _material_coverage(
                material_records, identifier
            ),
            "compactness": {
                "analysis_strata": len(catalog),
                "nonempty_core_strata": int(
                    accumulator["nonempty_core_strata"]
                ),
                "total_selected_hard_cells": int(
                    accumulator["total_hard_cells"]
                ),
                "total_selected_native_components": int(
                    accumulator["total_native_components"]
                ),
                "mean_selected_hard_cells_per_analysis_stratum": (
                    int(accumulator["total_hard_cells"]) / len(catalog)
                    if catalog
                    else None
                ),
            },
            "fallback_summary": {
                name: int(value)
                for name, value in accumulator[
                    "source_counts"
                ].items()
            },
            "coverage_summary": {
                "minimum_parent_coverage": minimum_parent_coverage,
                "assessed_strata": int(accumulator["assessed_strata"]),
                "passed_strata": int(accumulator["passed_strata"]),
                "failed_strata": int(accumulator["failed_strata"]),
                "all_assessed_strata_passed": (
                    int(accumulator["assessed_strata"]) > 0
                    and int(accumulator["failed_strata"]) == 0
                ),
                "aggregate_coverage_meets_minimum": (
                    validation_metrics[
                        "cross_section_weighted_core_coverage"
                    ]
                    is not None
                    and float(
                        validation_metrics[
                            "cross_section_weighted_core_coverage"
                        ]
                    )
                    >= minimum_parent_coverage
                ),
            },
        }

    ranking = sorted(
        candidates,
        key=lambda identifier: (
            -float(
                candidate_results[identifier]["validation"][
                    "cross_section_weighted_core_coverage"
                ]
                or 0.0
            ),
            int(
                candidate_results[identifier]["compactness"][
                    "total_selected_native_components"
                ]
            ),
            identifier,
        ),
    )
    training_metadata = migrations._campaign_comparison_metadata(
        training
    )
    validation_metadata = migrations._campaign_comparison_metadata(
        validation
    )
    training_inside = training_metadata["inside_analysis_partition"]
    validation_inside_metrics = validation_metadata[
        "inside_analysis_partition"
    ]
    comparison = {
        "schema": COMPARISON_SCHEMA,
        "created_utc": _now(),
        "passed": True,
        "study_completed": True,
        "production_ready": False,
        "production_readiness_note": (
            "Milestone 2d compares training-only core-guard geometries. "
            "Mode 4 remains disabled, and a nonzero full-support proposal "
            "tail is required before any production unweighting."
        ),
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
        "coordinate_definitions": {
            "hard_parent": migrations.HARD_COORDINATE_DEFINITION,
            "observed_target": migrations.OBSERVED_COORDINATE_DEFINITION,
            "observed_target_interpretation": (
                "final-LUND generator truth, not detector reconstruction"
            ),
        },
        "generator_revision": generator_revision,
        "generator_revision_source": generator_revision_source,
        "geometry_learner_revision": learner_revision,
        "survey_schema": training.norm_reference.get(
            "survey_schema", radiative_survey.SURVEY_SCHEMA
        ),
        "geometry_iteration": iteration,
        "legacy_input_sha256": training.input_signature,
        "legacy_generator_settings": training.legacy_input_settings,
        "legacy_channel_probabilities": training.channel_probabilities,
        "learning": {
            "spatial_channel_representation": "channel_marginalized",
            "native_intreg_sampling_unchanged": True,
            "milestone_2c_baseline_contract": {
                "independent_frozen": (
                    "identical construction to the milestone-2c "
                    "channel_marginalized frozen footprint"
                ),
                "independent_all_axis_dilation": (
                    "identical construction to its one-more-dilation "
                    "diagnostic"
                ),
            },
            "target_parent_fraction": target_parent_fraction,
            "observed_neighbor_manhattan_radius": neighbor_radius,
            "minimum_training_rows_for_supported_label": (
                minimum_training_rows
            ),
            "minimum_training_ess_for_supported_label": (
                minimum_training_ess
            ),
            "minimum_validation_parent_coverage": (
                minimum_parent_coverage
            ),
            "training_status_counts": status_counts,
            "template_training_source": template_training_source,
            "boundary_template_count": len(boundary_offsets),
            "global_template_offset_count": len(global_offsets),
        },
        "proposal_support_requirement": {
            "core_guard_is_not_a_hard_physics_cut": True,
            "future_full_support_tail_required": True,
            "tail_definition": (
                "nonzero proposal density over the complete legacy radiative "
                "proposal domain outside the selected core"
            ),
            "exact_proposal_density_correction_required": True,
        },
        "development_data_note": (
            "The held-out replicas are development data once used to choose "
            "a geometry. A final assessment requires fresh replicas."
        ),
        "training": training_metadata,
        "validation": validation_metadata,
        "inside_training_holdout_relative_difference": (
            (
                float(
                    validation_inside_metrics["cross_section_microbarn"]
                )
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
                float(
                    validation_inside_metrics[
                        "cross_section_microbarn"
                    ]
                ),
                validation_inside_metrics[
                    "cross_section_sem_microbarn"
                ],
            )
        ),
        "candidate_definitions": {
            identifier: candidate.metadata()
            for identifier, candidate in candidates.items()
        },
        "candidates": candidate_results,
        "ranking_by_heldout_coverage_then_compactness": ranking,
        "best_heldout_coverage_candidate": (
            ranking[0] if ranking else None
        ),
    }
    recipe_payload = {
        "schema": RECIPE_SCHEMA,
        "created_utc": comparison["created_utc"],
        "production_ready": False,
        "analysis_config": config,
        "analysis_config_source": str(config_path.resolve()),
        "analysis_config_sha256": config_sha256,
        "generator_revision": generator_revision,
        "geometry_learner_revision": learner_revision,
        "geometry_iteration": iteration,
        "legacy_input_sha256": training.input_signature,
        "training_replica_ids": training_metadata["replica_ids"],
        "target_parent_fraction": target_parent_fraction,
        "observed_neighbor_manhattan_radius": neighbor_radius,
        "minimum_training_rows_for_supported_label": (
            minimum_training_rows
        ),
        "minimum_training_ess_for_supported_label": (
            minimum_training_ess
        ),
        "candidate_definitions": comparison["candidate_definitions"],
        "offset_definition": (
            "hard-parent bin index minus final-LUND target-bin index; phi "
            "uses the signed periodic difference"
        ),
        "template_training_source": template_training_source,
        "global_offset_template": [
            _serialize_offset(offset)
            for offset in sorted(global_offsets)
        ],
        "boundary_offset_templates": {
            _boundary_identifier(signature): [
                _serialize_offset(offset)
                for offset in sorted(offsets)
            ]
            for signature, offsets in sorted(boundary_offsets.items())
        },
        "strata": recipes,
    }
    return comparison, recipe_payload, rows


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    fields = (
        "stratum_id",
        "flat_index",
        "training_status",
        "candidate",
        "training_cross_section_microbarn",
        "training_ess",
        "validation_cross_section_microbarn",
        "validation_ess",
        "own_seed_hard_cells",
        "neighbor_donor_strata",
        "used_neighbor_offsets",
        "used_boundary_template",
        "used_global_template",
        "selected_hard_cells",
        "selected_native_components",
        "training_core_coverage",
        "validation_core_coverage",
        "validation_core_purity_proxy",
        "coverage_passed",
    )
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def compare(args: argparse.Namespace) -> dict[str, object]:
    _validate_candidates()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(
            f"{output} already exists; geometry artifacts are immutable"
        )
    if not 0.0 < args.target_parent_fraction <= 1.0:
        raise GeometryError("--target-parent-fraction must be in (0,1]")
    if args.neighbor_radius < 0:
        raise GeometryError("--neighbor-radius must be nonnegative")
    if args.iteration < 0:
        raise GeometryError("--iteration must be nonnegative")
    if args.minimum_training_rows < 1 or args.minimum_training_ess < 0.0:
        raise GeometryError("training support thresholds are invalid")
    if not 0.0 <= args.minimum_parent_coverage <= 1.0:
        raise GeometryError("--minimum-parent-coverage must be in [0,1]")
    try:
        config, config_sha256 = guards.load_analysis_config(
            args.config.resolve()
        )
    except guards.GuardLearningError as error:
        raise GeometryError(str(error)) from error

    training = migrations.aggregate_surveys(
        args.training_survey,
        config,
        apply_y_max=args.apply_y_max,
    )
    validation = migrations.aggregate_surveys(
        args.validation_survey,
        config,
        apply_y_max=args.apply_y_max,
    )
    if training.input_signature != validation.input_signature:
        raise GeometryError(
            "training and validation surveys have different legacy inputs"
        )
    training_replicas = {
        int(item["replica"]) for item in training.replicas
    }
    validation_replicas = {
        int(item["replica"]) for item in validation.replicas
    }
    overlap = training_replicas & validation_replicas
    if overlap:
        raise GeometryError(
            f"validation replicas overlap training replicas: {sorted(overlap)}"
        )
    for campaign in (training, validation):
        if not math.isclose(
            float(config["beam_energy"]),
            float(campaign.norm_reference["ebeam"]),
            rel_tol=2.0e-7,
        ):
            raise GeometryError("analysis and survey beam energies differ")

    learner_revision = guards._current_revision()
    generator_revision = args.generator_revision or learner_revision
    comparison, recipes, rows = _build_geometry_study(
        config=config,
        config_path=args.config,
        config_sha256=config_sha256,
        training=training,
        validation=validation,
        target_parent_fraction=args.target_parent_fraction,
        neighbor_radius=args.neighbor_radius,
        iteration=args.iteration,
        minimum_training_rows=args.minimum_training_rows,
        minimum_training_ess=args.minimum_training_ess,
        minimum_parent_coverage=args.minimum_parent_coverage,
        apply_y_max=args.apply_y_max,
        generator_revision=generator_revision,
        generator_revision_source=(
            "user_supplied"
            if args.generator_revision
            else "assumed_current_checkout"
        ),
        learner_revision=learner_revision,
    )

    output.mkdir(parents=True)
    recipe_path = output / "guard_geometry_recipes.json"
    recipe_path.write_text(
        json.dumps(recipes, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    guards._write_sha256(recipe_path)
    row_path = output / "guard_geometry_strata.csv"
    _write_rows(row_path, rows)
    guards._write_sha256(row_path)
    comparison["artifacts"] = {
        "frozen_geometry_recipes": recipe_path.name,
        "frozen_geometry_recipes_sha256": hashlib.sha256(
            recipe_path.read_bytes()
        ).hexdigest(),
        "stratum_comparison": row_path.name,
        "stratum_comparison_sha256": hashlib.sha256(
            row_path.read_bytes()
        ).hexdigest(),
        "hash_sidecar_suffix": ".sha256",
    }
    comparison_path = output / "guard_geometry_comparison.json"
    comparison_path.write_text(
        json.dumps(comparison, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    guards._write_sha256(comparison_path)
    return {
        "passed": True,
        "schema": COMPARISON_SCHEMA,
        "comparison": str(comparison_path),
        "recipes": str(recipe_path),
        "stratum_comparison": str(row_path),
        "training_replicas": sorted(training_replicas),
        "validation_replicas": sorted(validation_replicas),
        "ranking_by_heldout_coverage_then_compactness": comparison[
            "ranking_by_heldout_coverage_then_compactness"
        ],
        "candidates": {
            identifier: {
                "heldout_core_coverage": values["validation"][
                    "cross_section_weighted_core_coverage"
                ],
                "heldout_core_purity_proxy": values["validation"][
                    "aggregate_core_purity_proxy"
                ],
                "total_selected_native_components": values[
                    "compactness"
                ]["total_selected_native_components"],
                "passed_strata": values["coverage_summary"][
                    "passed_strata"
                ],
                "failed_strata": values["coverage_summary"][
                    "failed_strata"
                ],
            }
            for identifier, values in comparison["candidates"].items()
        },
    }


def plot(args: argparse.Namespace) -> dict[str, object]:
    comparison_path = args.comparison.resolve()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(
            f"{output} already exists; geometry plots are immutable"
        )
    try:
        comparison = json.loads(
            comparison_path.read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as error:
        raise GeometryError(
            f"cannot read guard-geometry comparison: {error}"
        ) from error
    if comparison.get("schema") != COMPARISON_SCHEMA:
        raise GeometryError("guard-geometry comparison has wrong schema")
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError as error:
        raise GeometryError("plotting requires matplotlib") from error

    identifiers = [
        candidate.identifier for candidate in GEOMETRY_CANDIDATES
    ]
    candidates = comparison["candidates"]
    if set(identifiers) != set(candidates):
        raise GeometryError(
            "guard-geometry comparison lacks expected candidates"
        )
    labels = [str(candidates[name]["label"]) for name in identifiers]
    positions = list(range(len(identifiers)))
    coverages = [
        float(
            candidates[name]["validation"][
                "cross_section_weighted_core_coverage"
            ]
            or 0.0
        )
        for name in identifiers
    ]
    purities = [
        float(
            candidates[name]["validation"][
                "aggregate_core_purity_proxy"
            ]
            or 0.0
        )
        for name in identifiers
    ]
    components = [
        int(
            candidates[name]["compactness"][
                "total_selected_native_components"
            ]
        )
        for name in identifiers
    ]
    minimum = float(
        comparison["learning"]["minimum_validation_parent_coverage"]
    )

    output.mkdir(parents=True)
    pdf_path = output / "guard_geometry_comparison.pdf"
    with PdfPages(pdf_path) as pdf:
        figure, axis = plt.subplots(figsize=(11, 5.5))
        axis.bar(positions, coverages)
        axis.axhline(
            minimum,
            color="black",
            linestyle="--",
            linewidth=1,
            label=f"Requested minimum ({minimum:.1%})",
        )
        axis.set_xticks(positions)
        axis.set_xticklabels(labels, rotation=12, ha="right")
        axis.set_ylim(0.0, 1.03)
        axis.set_ylabel(
            "Cross-section-weighted held-out core coverage"
        )
        axis.set_title("Radiative core-guard geometry comparison")
        axis.grid(axis="y", alpha=0.25)
        axis.legend()
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        figure, axis = plt.subplots(figsize=(10, 5.5))
        axis.scatter(components, purities, s=70)
        for index, label in enumerate(labels):
            axis.annotate(
                label,
                (components[index], purities[index]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize="small",
            )
        axis.set_xscale("log")
        axis.set_xlabel("Total selected native components")
        axis.set_ylabel("Held-out aggregate core purity proxy")
        axis.set_ylim(0.0, max(purities) * 1.2 if purities else 1.0)
        axis.set_title("Core size versus selectivity")
        axis.grid(alpha=0.25)
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        figure, axis = plt.subplots(figsize=(12, 5.5))
        width = 0.2
        status_positions = list(range(len(STATUS_ORDER)))
        for candidate_index, identifier in enumerate(identifiers):
            offset = (
                candidate_index - (len(identifiers) - 1) / 2
            ) * width
            values = [
                float(
                    candidates[identifier][
                        "validation_by_training_status"
                    ][status][
                        "cross_section_weighted_core_coverage"
                    ]
                    or 0.0
                )
                for status in STATUS_ORDER
            ]
            axis.bar(
                [
                    position + offset
                    for position in status_positions
                ],
                values,
                width=width,
                label=labels[candidate_index],
            )
        axis.set_xticks(status_positions)
        axis.set_xticklabels(
            ("Well supported", "Low support", "No training")
        )
        axis.set_ylim(0.0, 1.03)
        axis.set_ylabel("Held-out core coverage")
        axis.set_title("Coverage by training-support class")
        axis.grid(axis="y", alpha=0.25)
        axis.legend(fontsize="small")
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        passed = [
            int(candidates[name]["coverage_summary"]["passed_strata"])
            for name in identifiers
        ]
        failed = [
            int(candidates[name]["coverage_summary"]["failed_strata"])
            for name in identifiers
        ]
        figure, axis = plt.subplots(figsize=(11, 5.5))
        axis.bar(positions, passed, label="At or above threshold")
        axis.bar(
            positions,
            failed,
            bottom=passed,
            label="Below threshold",
        )
        axis.set_xticks(positions)
        axis.set_xticklabels(labels, rotation=12, ha="right")
        axis.set_ylabel("Held-out analysis strata")
        axis.set_title("Per-stratum coverage requirement")
        axis.grid(axis="y", alpha=0.25)
        axis.legend()
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

    guards._write_sha256(pdf_path)
    summary = {
        "schema": PLOT_SCHEMA,
        "created_utc": _now(),
        "comparison": str(comparison_path),
        "comparison_sha256": hashlib.sha256(
            comparison_path.read_bytes()
        ).hexdigest(),
        "pdf": pdf_path.name,
        "candidates": identifiers,
    }
    summary_path = output / "plot_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    guards._write_sha256(summary_path)
    return {
        "passed": True,
        "schema": PLOT_SCHEMA,
        "output": str(output),
        "pdf": str(pdf_path),
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    comparison_parser = subparsers.add_parser(
        "compare",
        help=(
            "compare channel-marginalized core-guard fallback geometries"
        ),
    )
    comparison_parser.add_argument("--config", type=Path, required=True)
    comparison_parser.add_argument(
        "--training-survey", type=Path, nargs="+", required=True
    )
    comparison_parser.add_argument(
        "--validation-survey", type=Path, nargs="+", required=True
    )
    comparison_parser.add_argument("--output", type=Path, required=True)
    comparison_parser.add_argument(
        "--target-parent-fraction", type=float, default=0.995
    )
    comparison_parser.add_argument(
        "--neighbor-radius",
        type=int,
        default=1,
        help=(
            "Manhattan radius in final-LUND analysis-bin indices used to "
            "find offset donors"
        ),
    )
    comparison_parser.add_argument("--iteration", type=int, default=0)
    comparison_parser.add_argument(
        "--minimum-training-rows", type=int, default=10
    )
    comparison_parser.add_argument(
        "--minimum-training-ess", type=float, default=5.0
    )
    comparison_parser.add_argument(
        "--minimum-parent-coverage", type=float, default=0.98
    )
    comparison_parser.add_argument(
        "--apply-y-max",
        action="store_true",
        help=(
            "apply phase_space.y_max in training and validation; the "
            "default applies no y cut"
        ),
    )
    comparison_parser.add_argument(
        "--generator-revision",
        help=(
            "Git revision that produced the surveys; defaults to the current "
            "checkout and is recorded as an assumption"
        ),
    )

    plot_parser = subparsers.add_parser(
        "plot", help="render the milestone-2d geometry comparison"
    )
    plot_parser.add_argument("--comparison", type=Path, required=True)
    plot_parser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> int:
    args = _build_parser().parse_args()
    try:
        if args.command == "compare":
            result = compare(args)
        else:
            result = plot(args)
    except (
        FileNotFoundError,
        FileExistsError,
        KeyError,
        ValueError,
        GeometryError,
    ) as error:
        print(f"ERROR: {error}")
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result.get("passed", True) else 2


if __name__ == "__main__":
    raise SystemExit(main())
