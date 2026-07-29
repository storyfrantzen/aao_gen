#!/usr/bin/env python3
"""Fit continuous native-coordinate guards for radiative mode 4.

This milestone-2e diagnostic replaces whole-cell dilation with weighted,
joint axis-aligned boxes in AAO's normalized proposal coordinates.  Box faces
are optimized at empirical weighted-CDF breakpoints, sparse strata borrow
neighboring final-LUND strata with ESS-controlled regularization, and a
development-only padding scan is evaluated on held-out survey replicas.

The output is not a production mode-4 manifest.  A future proposal must mix
the selected core with a nonzero unrestricted tail and evaluate the exact
mixture density during unweighting.
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
import radiative_survey


COMPARISON_SCHEMA = "aao-rad-continuous-guard-comparison-v1"
RECIPE_SCHEMA = "aao-rad-continuous-guard-recipes-v1"
PLOT_SCHEMA = "aao-rad-continuous-guard-plots-v1"
DEFAULT_PADDINGS = (0.0, 0.0025, 0.005, 0.01, 0.02)
STATUS_ORDER = (
    "learned",
    "learned_low_support",
    "no_training_contribution",
)
NATIVE_CHANNELS = tuple(range(1, 7))


class ContinuousGuardError(RuntimeError):
    """Raised when a continuous guard study cannot be built reproducibly."""


@dataclass(frozen=True)
class Axis:
    name: str
    periodic: bool = False


AXES = (
    Axis("r_u"),
    Axis("r_ep"),
    Axis("u_gamma"),
    Axis("hadron_cosine_base"),
    Axis("hadron_phi_base", periodic=True),
)
PHI_INDEX = len(AXES) - 1


@dataclass(frozen=True)
class Point:
    coordinates: tuple[float, ...]
    weight: float
    channel: int


@dataclass
class ContinuousCampaign:
    base: guards.Campaign
    points: dict[str, list[Point]]

    @property
    def proposals(self) -> int:
        return self.base.proposals


@dataclass(frozen=True)
class Box:
    lower: tuple[float, ...]
    upper: tuple[float, ...]
    phi_origin: float
    fit_iterations: int = 0

    @property
    def volume(self) -> float:
        return math.prod(
            max(0.0, right - left)
            for left, right in zip(self.lower, self.upper)
        )

    def contains(self, coordinates: tuple[float, ...]) -> bool:
        transformed = list(coordinates)
        transformed[PHI_INDEX] = _circular_delta(
            coordinates[PHI_INDEX], self.phi_origin
        )
        tolerance = 2.0e-15
        return all(
            left - tolerance <= value <= right + tolerance
            for value, left, right in zip(
                transformed, self.lower, self.upper
            )
        )

    def padded(self, amount: float) -> "Box":
        if amount < 0.0 or not math.isfinite(amount):
            raise ContinuousGuardError("padding must be finite and nonnegative")
        lower: list[float] = []
        upper: list[float] = []
        for index, (left, right) in enumerate(
            zip(self.lower, self.upper)
        ):
            domain_left, domain_right = (
                (-0.5, 0.5) if index == PHI_INDEX else (0.0, 1.0)
            )
            lower.append(max(domain_left, left - amount))
            upper.append(min(domain_right, right + amount))
        return Box(
            lower=tuple(lower),
            upper=tuple(upper),
            phi_origin=self.phi_origin,
            fit_iterations=self.fit_iterations,
        )

    def metadata(self) -> dict[str, object]:
        axes: list[dict[str, object]] = []
        for index, axis in enumerate(AXES):
            left = self.lower[index]
            right = self.upper[index]
            width = right - left
            if not axis.periodic:
                axes.append(
                    {
                        "name": axis.name,
                        "periodic": False,
                        "lower": left,
                        "upper": right,
                        "width": width,
                    }
                )
                continue
            full = width >= 1.0 - 2.0e-15
            if full:
                start = 0.0
                end = 1.0
                wraps = False
            else:
                start = (self.phi_origin + left) % 1.0
                end = (self.phi_origin + right) % 1.0
                wraps = start > end
            axes.append(
                {
                    "name": axis.name,
                    "periodic": True,
                    "origin": self.phi_origin,
                    "lower_relative_to_origin": left,
                    "upper_relative_to_origin": right,
                    "interval_start": start,
                    "interval_end": end,
                    "wraps": wraps,
                    "full_period": full,
                    "width": width,
                }
            )
        return {
            "representation": "continuous_axis_aligned_native_box",
            "coordinate_domain": [0.0, 1.0],
            "volume": self.volume,
            "fit_iterations": self.fit_iterations,
            "axes": axes,
        }


def _now() -> str:
    return (
        dt.datetime.now(dt.timezone.utc)
        .replace(microsecond=0)
        .isoformat()
    )


def _circular_delta(value: float, origin: float) -> float:
    return (value - origin + 0.5) % 1.0 - 0.5


def _circular_origin(points: list[Point]) -> float:
    cosine = sum(
        point.weight * math.cos(2.0 * math.pi * point.coordinates[PHI_INDEX])
        for point in points
    )
    sine = sum(
        point.weight * math.sin(2.0 * math.pi * point.coordinates[PHI_INDEX])
        for point in points
    )
    scale = sum(point.weight for point in points)
    if scale > 0.0 and math.hypot(cosine, sine) > 1.0e-12 * scale:
        return math.atan2(sine, cosine) / (2.0 * math.pi) % 1.0
    return max(
        points,
        key=lambda point: (
            point.weight,
            -point.coordinates[PHI_INDEX],
        ),
    ).coordinates[PHI_INDEX]


def _expanded_minimum_width(
    left: float,
    right: float,
    minimum_width: float,
    domain_left: float,
    domain_right: float,
) -> tuple[float, float]:
    width = right - left
    if width >= minimum_width:
        return left, right
    center = 0.5 * (left + right)
    new_left = center - 0.5 * minimum_width
    new_right = center + 0.5 * minimum_width
    if new_left < domain_left:
        new_right += domain_left - new_left
        new_left = domain_left
    if new_right > domain_right:
        new_left -= new_right - domain_right
        new_right = domain_right
    return max(domain_left, new_left), min(domain_right, new_right)


def _fit_joint_box(
    points: list[Point],
    target_fraction: float,
    minimum_axis_width: float,
) -> Box:
    """Greedily shrink a unit box at weighted empirical-CDF breakpoints."""
    if not points or sum(point.weight for point in points) <= 0.0:
        raise ContinuousGuardError("cannot fit a box without positive weight")
    if not 0.0 < target_fraction <= 1.0:
        raise ContinuousGuardError("target fraction must be in (0,1]")
    if not 0.0 < minimum_axis_width <= 1.0:
        raise ContinuousGuardError("minimum axis width must be in (0,1]")

    origin = _circular_origin(points)
    transformed: list[tuple[tuple[float, ...], float]] = []
    for point in points:
        coordinates = list(point.coordinates)
        coordinates[PHI_INDEX] = _circular_delta(
            coordinates[PHI_INDEX], origin
        )
        transformed.append((tuple(coordinates), point.weight))

    lower = [0.0] * len(AXES)
    upper = [1.0] * len(AXES)
    lower[PHI_INDEX] = -0.5
    upper[PHI_INDEX] = 0.5
    total = sum(weight for _, weight in transformed)
    required = target_fraction * total
    inside = [True] * len(transformed)
    iterations = 0

    def optimization_volume(
        trial_lower: list[float], trial_upper: list[float]
    ) -> float:
        return math.prod(
            max(right - left, minimum_axis_width)
            for left, right in zip(trial_lower, trial_upper)
        )

    while True:
        current_weight = sum(
            weight
            for keep, (_, weight) in zip(inside, transformed)
            if keep
        )
        slack = max(0.0, current_weight - required)
        current_volume = optimization_volume(lower, upper)
        candidates: list[
            tuple[
                float,
                float,
                int,
                int,
                float,
            ]
        ] = []

        for dimension in range(len(AXES)):
            grouped: dict[float, float] = {}
            for keep, (coordinates, weight) in zip(inside, transformed):
                if keep:
                    grouped[coordinates[dimension]] = (
                        grouped.get(coordinates[dimension], 0.0) + weight
                    )
            ordered = sorted(grouped.items())
            removed = 0.0
            for value, group_weight in ordered:
                if value > lower[dimension] + 2.0e-15:
                    if removed <= slack + 1.0e-13 * total:
                        trial_lower = lower.copy()
                        trial_lower[dimension] = value
                        volume = optimization_volume(trial_lower, upper)
                        saved = current_volume - volume
                        if saved > max(
                            1.0e-300, 1.0e-12 * current_volume
                        ):
                            score = saved / max(removed / total, 1.0e-15)
                            candidates.append(
                                (
                                    score,
                                    saved,
                                    dimension,
                                    0,
                                    value,
                                )
                            )
                removed += group_weight

            removed = 0.0
            for value, group_weight in reversed(ordered):
                if value < upper[dimension] - 2.0e-15:
                    if removed <= slack + 1.0e-13 * total:
                        trial_upper = upper.copy()
                        trial_upper[dimension] = value
                        volume = optimization_volume(lower, trial_upper)
                        saved = current_volume - volume
                        if saved > max(
                            1.0e-300, 1.0e-12 * current_volume
                        ):
                            score = saved / max(removed / total, 1.0e-15)
                            candidates.append(
                                (
                                    score,
                                    saved,
                                    dimension,
                                    1,
                                    value,
                                )
                            )
                removed += group_weight

        if not candidates:
            break
        _, _, dimension, side, value = max(
            candidates,
            key=lambda item: (
                item[0],
                item[1],
                -item[2],
                -item[3],
                -item[4],
            ),
        )
        trial_inside = [
            keep
            and (
                coordinates[dimension] >= value - 2.0e-15
                if side == 0
                else coordinates[dimension] <= value + 2.0e-15
            )
            for keep, (coordinates, _) in zip(inside, transformed)
        ]
        trial_weight = sum(
            weight
            for keep, (_, weight) in zip(trial_inside, transformed)
            if keep
        )
        if trial_weight + 1.0e-13 * total < required:
            raise ContinuousGuardError(
                "internal joint-box optimization violated target coverage"
            )
        if side == 0:
            lower[dimension] = value
        else:
            upper[dimension] = value
        inside = trial_inside
        iterations += 1

    for dimension in range(len(AXES)):
        domain = (
            (-0.5, 0.5) if dimension == PHI_INDEX else (0.0, 1.0)
        )
        lower[dimension], upper[dimension] = _expanded_minimum_width(
            lower[dimension],
            upper[dimension],
            minimum_axis_width,
            *domain,
        )

    box = Box(
        lower=tuple(lower),
        upper=tuple(upper),
        phi_origin=origin,
        fit_iterations=iterations,
    )
    fitted_weight = sum(
        point.weight for point in points if box.contains(point.coordinates)
    )
    if fitted_weight + 1.0e-12 * total < required:
        raise ContinuousGuardError(
            "minimum-width adjustment reduced fitted box coverage"
        )
    return box


def _coordinates(row: dict[str, str]) -> tuple[float, ...]:
    values = tuple(float(row[axis.name]) for axis in AXES)
    if any(
        not math.isfinite(value) or not 0.0 <= value <= 1.0
        for value in values
    ):
        raise ContinuousGuardError(
            "survey base proposal coordinate is outside [0,1]"
        )
    values = list(values)
    values[PHI_INDEX] %= 1.0
    return tuple(values)


def _load_campaign(
    directories: list[Path],
    config: dict,
    *,
    apply_y_max: bool,
) -> ContinuousCampaign:
    base = guards.aggregate_surveys(
        directories,
        config,
        guards.parse_partition(),
        apply_y_max=apply_y_max,
    )
    points: dict[str, list[Point]] = {}
    for requested_directory in directories:
        path = requested_directory.resolve() / radiative_survey.SURVEY_FILENAME
        for row in guards._survey_rows(path):
            stratum_id, _ = guards.assign_stratum(
                row, config, apply_y_max=apply_y_max
            )
            if stratum_id is None:
                continue
            weight = float(row["trial_xsec_observed_microbarn"])
            channel = int(row["intreg"])
            if weight <= 0.0:
                continue
            points.setdefault(stratum_id, []).append(
                Point(
                    coordinates=_coordinates(row),
                    weight=weight,
                    channel=channel,
                )
            )
    for stratum_id, moment in base.strata.items():
        loaded = points.get(stratum_id, [])
        loaded_total = sum(point.weight for point in loaded)
        if len(loaded) != moment.count or not math.isclose(
            loaded_total,
            moment.total,
            rel_tol=2.0e-12,
            abs_tol=1.0e-18,
        ):
            raise ContinuousGuardError(
                f"{stratum_id}: continuous-coordinate reload does not "
                "match validated survey aggregation"
            )
    return ContinuousCampaign(base=base, points=points)


def _status(
    moment: guards.Moment,
    *,
    minimum_training_rows: int,
    minimum_training_ess: float,
) -> str:
    if moment.total <= 0.0:
        return "no_training_contribution"
    if (
        moment.count >= minimum_training_rows
        and guards._ess(moment) >= minimum_training_ess
    ):
        return "learned"
    return "learned_low_support"


def _target_shape(config: dict) -> tuple[int, int, int, int]:
    binning = config["binning"]
    return (
        len(binning["Q2"]) - 1,
        len(binning["xB"]) - 1,
        len(binning["minus_t"]) - 1,
        len(binning["phi_deg"]) - 1,
    )


def _stratum_key(
    stratum: guards.Stratum,
) -> tuple[int, int, int, int]:
    return (stratum.iq2, stratum.ixb, stratum.it, stratum.iphi)


def _neighbor_ids(
    stratum: guards.Stratum,
    lookup: dict[tuple[int, int, int, int], guards.Stratum],
    shape: tuple[int, int, int, int],
    radius: int,
) -> tuple[str, ...]:
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


def _scaled_points(
    points: Iterable[Point], desired_total: float
) -> list[Point]:
    source = list(points)
    total = sum(point.weight for point in source)
    if total <= 0.0 or desired_total <= 0.0:
        return []
    scale = desired_total / total
    return [
        Point(
            coordinates=point.coordinates,
            weight=point.weight * scale,
            channel=point.channel,
        )
        for point in source
    ]


def _bounds_metadata(stratum: guards.Stratum) -> dict[str, object]:
    return {
        "Q2": list(stratum.q2),
        "xB": list(stratum.xb),
        "minus_t": list(stratum.minus_t),
        "phi_deg": list(stratum.phi_deg),
    }


def _candidate_id(padding: float) -> str:
    label = f"{padding:.8g}".replace(".", "p")
    return f"padding_{label}"


def _add_moment(target: guards.Moment, source: guards.Moment) -> None:
    target.count += source.count
    target.total += source.total
    target.square_total += source.square_total
    target.maximum = max(target.maximum, source.maximum)


def _combined(moments: Iterable[guards.Moment]) -> guards.Moment:
    result = guards.Moment()
    for moment in moments:
        _add_moment(result, moment)
    return result


def _fraction_metrics(
    total: guards.Moment,
    selected: guards.Moment,
    tail: guards.Moment,
    proposals: int,
) -> dict[str, object]:
    return {
        "total": guards._compact_metrics(total, proposals),
        "selected_core": guards._compact_metrics(selected, proposals),
        "residual_outside_core": guards._compact_metrics(tail, proposals),
        "cross_section_weighted_core_coverage": (
            selected.total / total.total if total.total > 0.0 else None
        ),
        "cross_section_weighted_residual_fraction": (
            tail.total / total.total if total.total > 0.0 else None
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
            selected += float(record["selected"][candidate])
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


def _campaign_metadata(campaign: ContinuousCampaign) -> dict[str, object]:
    inside = _combined(campaign.base.strata.values())
    return {
        "replica_ids": [
            int(item["replica"]) for item in campaign.base.replicas
        ],
        "replicas": campaign.base.replicas,
        "total_proposals": campaign.proposals,
        "recorded_rows": campaign.base.rows,
        "final_valid_rows": campaign.base.final_valid_rows,
        "global_observed": guards._compact_metrics(
            campaign.base.global_observed, campaign.proposals
        ),
        "inside_analysis_partition": guards._compact_metrics(
            inside, campaign.proposals
        ),
    }


def _build_study(
    *,
    config: dict,
    config_path: Path,
    config_sha256: str,
    training: ContinuousCampaign,
    validation: ContinuousCampaign,
    target_core_fraction: float,
    paddings: tuple[float, ...],
    neighbor_radius: int,
    regularization_ess: float,
    minimum_axis_width: float,
    minimum_training_rows: int,
    minimum_training_ess: float,
    minimum_validation_coverage: float,
    apply_y_max: bool,
    iteration: int,
    generator_revision: str,
    generator_revision_source: str,
    learner_revision: str,
) -> tuple[dict[str, object], dict[str, object], list[dict[str, object]]]:
    catalog = guards.enumerate_strata(config)
    lookup = {_stratum_key(stratum): stratum for stratum in catalog}
    shape = _target_shape(config)
    statuses = {
        stratum.identifier: _status(
            training.base.strata.get(
                stratum.identifier, guards.Moment()
            ),
            minimum_training_rows=minimum_training_rows,
            minimum_training_ess=minimum_training_ess,
        )
        for stratum in catalog
    }
    global_points = [
        point
        for points in training.points.values()
        for point in points
    ]
    if not global_points:
        raise ContinuousGuardError(
            "training surveys contain no positive contribution inside analysis"
        )
    global_box = _fit_joint_box(
        global_points, target_core_fraction, minimum_axis_width
    )

    recipes: dict[str, dict[str, object]] = {}
    boxes: dict[str, Box] = {}
    padding_scales: dict[str, float] = {}
    status_counts = {status: 0 for status in STATUS_ORDER}
    fit_source_counts = {
        "local_only": 0,
        "local_neighbor_blend": 0,
        "neighbor_only": 0,
        "global_fallback": 0,
    }

    for stratum in catalog:
        stratum_id = stratum.identifier
        status = statuses[stratum_id]
        status_counts[status] += 1
        local = training.points.get(stratum_id, [])
        local_moment = training.base.strata.get(
            stratum_id, guards.Moment()
        )
        local_ess = guards._ess(local_moment)
        neighbors = tuple(
            neighbor
            for neighbor in _neighbor_ids(
                stratum, lookup, shape, neighbor_radius
            )
            if training.points.get(neighbor)
        )
        supported_neighbors = tuple(
            neighbor
            for neighbor in neighbors
            if statuses[neighbor] == "learned"
        )
        donors = supported_neighbors or neighbors
        neighbor_points = [
            point
            for neighbor in donors
            for point in training.points.get(neighbor, [])
        ]

        if status == "learned":
            fit_points = local
            local_fraction = 1.0
            source = "local_only"
        elif local and neighbor_points:
            local_fraction = local_ess / (
                local_ess + regularization_ess
            )
            fit_points = (
                _scaled_points(local, local_fraction)
                + _scaled_points(neighbor_points, 1.0 - local_fraction)
            )
            source = "local_neighbor_blend"
        elif local:
            fit_points = local
            local_fraction = 1.0
            source = "local_only"
        elif neighbor_points:
            fit_points = neighbor_points
            local_fraction = 0.0
            source = "neighbor_only"
        else:
            fit_points = []
            local_fraction = 0.0
            source = "global_fallback"

        if source == "global_fallback":
            box = global_box
        else:
            box = _fit_joint_box(
                fit_points, target_core_fraction, minimum_axis_width
            )
        boxes[stratum_id] = box
        fit_source_counts[source] += 1
        padding_scale = (
            1.0
            if status == "learned"
            else 1.0 + (1.0 - local_fraction)
        )
        padding_scales[stratum_id] = padding_scale
        fit_total = sum(point.weight for point in fit_points)
        fit_selected = sum(
            point.weight
            for point in fit_points
            if box.contains(point.coordinates)
        )
        recipes[stratum_id] = {
            "stratum_id": stratum_id,
            "flat_index": stratum.flat_index,
            "indices": {
                "Q2": stratum.iq2,
                "xB": stratum.ixb,
                "minus_t": stratum.it,
                "phi_deg": stratum.iphi,
            },
            "bounds": _bounds_metadata(stratum),
            "training_status": status,
            "training_rows": local_moment.count,
            "training_ess": local_ess,
            "fit_source": source,
            "neighbor_strata": list(donors),
            "local_fit_weight_fraction": local_fraction,
            "neighbor_fit_weight_fraction": (
                1.0 - local_fraction
                if source == "local_neighbor_blend"
                else (1.0 if source == "neighbor_only" else 0.0)
            ),
            "padding_scale": padding_scale,
            "raw_box": box.metadata(),
            "raw_fit_weighted_coverage": (
                fit_selected / fit_total if fit_total > 0.0 else None
            ),
        }

    candidate_ids = tuple(_candidate_id(value) for value in paddings)
    candidate_boxes: dict[str, dict[str, Box]] = {
        candidate: {} for candidate in candidate_ids
    }
    compactness: dict[str, dict[str, object]] = {}
    for candidate, padding in zip(candidate_ids, paddings):
        volumes: list[float] = []
        for stratum in catalog:
            identifier = stratum.identifier
            padded = boxes[identifier].padded(
                padding * padding_scales[identifier]
            )
            candidate_boxes[candidate][identifier] = padded
            volumes.append(padded.volume)
        ordered_volumes = sorted(volumes)
        compactness[candidate] = {
            "analysis_strata": len(catalog),
            "sum_normalized_core_volume": sum(volumes),
            "mean_normalized_core_volume": sum(volumes) / len(volumes),
            "median_normalized_core_volume": ordered_volumes[
                len(ordered_volumes) // 2
            ],
            "maximum_normalized_core_volume": max(volumes),
        }

    def evaluate(
        campaign: ContinuousCampaign,
    ) -> tuple[
        dict[str, dict[str, guards.Moment]],
        dict[str, dict[str, dict[str, guards.Moment]]],
        dict[str, dict[int, dict[str, guards.Moment]]],
        dict[str, dict[str, dict[str, guards.Moment]]],
    ]:
        aggregate = {
            candidate: {
                "total": guards.Moment(),
                "selected": guards.Moment(),
                "tail": guards.Moment(),
            }
            for candidate in candidate_ids
        }
        by_status = {
            candidate: {
                status: {
                    "total": guards.Moment(),
                    "selected": guards.Moment(),
                    "tail": guards.Moment(),
                }
                for status in STATUS_ORDER
            }
            for candidate in candidate_ids
        }
        by_channel = {
            candidate: {
                channel: {
                    "total": guards.Moment(),
                    "selected": guards.Moment(),
                    "tail": guards.Moment(),
                }
                for channel in NATIVE_CHANNELS
            }
            for candidate in candidate_ids
        }
        per_stratum: dict[
            str, dict[str, dict[str, guards.Moment]]
        ] = {}
        for stratum in catalog:
            identifier = stratum.identifier
            points = campaign.points.get(identifier, [])
            per_stratum[identifier] = {}
            for candidate in candidate_ids:
                moments = {
                    "total": guards.Moment(),
                    "selected": guards.Moment(),
                    "tail": guards.Moment(),
                }
                box = candidate_boxes[candidate][identifier]
                for point in points:
                    moments["total"].add(point.weight)
                    target = (
                        moments["selected"]
                        if box.contains(point.coordinates)
                        else moments["tail"]
                    )
                    target.add(point.weight)
                    channel_moments = by_channel[candidate][point.channel]
                    channel_moments["total"].add(point.weight)
                    (
                        channel_moments["selected"]
                        if box.contains(point.coordinates)
                        else channel_moments["tail"]
                    ).add(point.weight)
                per_stratum[identifier][candidate] = moments
                for name in ("total", "selected", "tail"):
                    _add_moment(
                        aggregate[candidate][name], moments[name]
                    )
                    _add_moment(
                        by_status[candidate][statuses[identifier]][name],
                        moments[name],
                    )
        return aggregate, by_status, by_channel, per_stratum

    (
        training_aggregate,
        _,
        _,
        training_per_stratum,
    ) = evaluate(training)
    (
        validation_aggregate,
        validation_by_status,
        validation_by_channel,
        validation_per_stratum,
    ) = evaluate(validation)

    rows: list[dict[str, object]] = []
    material_records: list[dict[str, object]] = []
    coverage_summaries: dict[str, dict[str, object]] = {}
    for candidate, padding in zip(candidate_ids, paddings):
        assessed = passed = failed = 0
        for stratum in catalog:
            identifier = stratum.identifier
            training_moments = training_per_stratum[identifier][candidate]
            validation_moments = validation_per_stratum[identifier][
                candidate
            ]
            training_total = training_moments["total"].total
            validation_total = validation_moments["total"].total
            training_fraction = (
                training_moments["selected"].total / training_total
                if training_total > 0.0
                else None
            )
            validation_fraction = (
                validation_moments["selected"].total / validation_total
                if validation_total > 0.0
                else None
            )
            coverage_passed = (
                validation_fraction >= minimum_validation_coverage
                if validation_fraction is not None
                else None
            )
            if coverage_passed is not None:
                assessed += 1
                passed += int(coverage_passed)
                failed += int(not coverage_passed)
            padded = candidate_boxes[candidate][identifier]
            rows.append(
                {
                    "stratum_id": identifier,
                    "flat_index": stratum.flat_index,
                    "training_status": statuses[identifier],
                    "fit_source": recipes[identifier]["fit_source"],
                    "candidate": candidate,
                    "base_padding": padding,
                    "effective_padding": (
                        padding * padding_scales[identifier]
                    ),
                    "normalized_core_volume": padded.volume,
                    "training_cross_section_microbarn": (
                        training_total / training.proposals
                    ),
                    "training_ess": guards._ess(
                        training_moments["total"]
                    ),
                    "training_core_coverage": training_fraction,
                    "validation_cross_section_microbarn": (
                        validation_total / validation.proposals
                    ),
                    "validation_ess": guards._ess(
                        validation_moments["total"]
                    ),
                    "validation_core_coverage": validation_fraction,
                    "validation_tail_ess": guards._ess(
                        validation_moments["tail"]
                    ),
                    "coverage_passed": coverage_passed,
                }
            )
        coverage_summaries[candidate] = {
            "minimum_validation_coverage": minimum_validation_coverage,
            "assessed_strata": assessed,
            "passed_strata": passed,
            "failed_strata": failed,
            "all_assessed_strata_passed": assessed > 0 and failed == 0,
        }

    for stratum in catalog:
        identifier = stratum.identifier
        material_records.append(
            {
                "validation_total": validation_per_stratum[identifier][
                    candidate_ids[0]
                ]["total"].total,
                "selected": {
                    candidate: validation_per_stratum[identifier][
                        candidate
                    ]["selected"].total
                    for candidate in candidate_ids
                },
            }
        )

    candidate_results: dict[str, object] = {}
    inside_validation = _combined(validation.base.strata.values()).total
    for candidate, padding in zip(candidate_ids, paddings):
        training_metrics = _fraction_metrics(
            **training_aggregate[candidate],
            proposals=training.proposals,
        )
        validation_metrics = _fraction_metrics(
            **validation_aggregate[candidate],
            proposals=validation.proposals,
        )
        candidate_results[candidate] = {
            "identifier": candidate,
            "base_padding": padding,
            "adaptive_padding_rule": (
                "effective_padding = base_padding * "
                "(1 for learned; 1 + (1-local_fit_fraction) otherwise)"
            ),
            "training": training_metrics,
            "validation": validation_metrics,
            "validation_by_training_status": {
                status: _fraction_metrics(
                    **validation_by_status[candidate][status],
                    proposals=validation.proposals,
                )
                for status in STATUS_ORDER
            },
            "validation_by_native_intreg": {
                f"intreg_{channel}": {
                    **_fraction_metrics(
                        **validation_by_channel[candidate][channel],
                        proposals=validation.proposals,
                    ),
                    "fraction_of_inside_analysis_cross_section": (
                        validation_by_channel[candidate][channel][
                            "total"
                        ].total
                        / inside_validation
                        if inside_validation > 0.0
                        else None
                    ),
                }
                for channel in NATIVE_CHANNELS
            },
            "validation_material_strata": _material_coverage(
                material_records, candidate
            ),
            "compactness": compactness[candidate],
            "coverage_summary": {
                **coverage_summaries[candidate],
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
                    >= minimum_validation_coverage
                ),
            },
        }

    passing_candidates = [
        candidate
        for candidate in candidate_ids
        if candidate_results[candidate]["coverage_summary"][
            "aggregate_coverage_meets_minimum"
        ]
    ]
    if passing_candidates:
        recommended = min(
            passing_candidates,
            key=lambda candidate: (
                candidate_results[candidate]["compactness"][
                    "sum_normalized_core_volume"
                ],
                candidate_results[candidate]["base_padding"],
            ),
        )
        recommendation_reason = (
            "smallest total normalized core volume among padding candidates "
            "meeting aggregate held-out coverage"
        )
    else:
        recommended = max(
            candidate_ids,
            key=lambda candidate: (
                candidate_results[candidate]["validation"][
                    "cross_section_weighted_core_coverage"
                ]
                or -1.0,
                -candidate_results[candidate]["compactness"][
                    "sum_normalized_core_volume"
                ],
            ),
        )
        recommendation_reason = (
            "no candidate met aggregate held-out coverage; highest observed "
            "coverage is reported for further development"
        )

    training_inside = _combined(training.base.strata.values())
    validation_inside = _combined(validation.base.strata.values())
    training_sigma = training_inside.total / training.proposals
    validation_sigma = validation_inside.total / validation.proposals
    training_sem = guards._sem(training_inside, training.proposals)
    validation_sem = guards._sem(validation_inside, validation.proposals)
    comparison = {
        "schema": COMPARISON_SCHEMA,
        "created_utc": _now(),
        "passed": True,
        "study_completed": True,
        "production_ready": False,
        "production_readiness_note": (
            "Development-only continuous-core diagnostic. Mode 4 still "
            "requires exact core-plus-global-tail sampling and unweighting, "
            "and the selected padding must be tested on fresh replicas."
        ),
        "continuous_guard_iteration": iteration,
        "generator_revision": generator_revision,
        "generator_revision_source": generator_revision_source,
        "continuous_guard_learner_revision": learner_revision,
        "survey_schema": radiative_survey.SURVEY_SCHEMA,
        "analysis_config_source": str(config_path.resolve()),
        "analysis_config_sha256": config_sha256,
        "analysis_selection": {
            "coordinate_definition": "final_lund_analysis",
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
        "proposal_coordinates": {
            "domain": [0.0, 1.0],
            "axes": [
                {"name": axis.name, "periodic": axis.periodic}
                for axis in AXES
            ],
            "native_intreg_sampling_unchanged": True,
            "unrestricted_in_this_core_study": [
                "photon_side_base",
                "photon_cosine_base",
                "photon_phi_base",
                "incoming_loss_base",
                "outgoing_loss_base",
            ],
        },
        "learning": {
            "algorithm": (
                "greedy_joint_minimum_volume_box_at_weighted_empirical_"
                "cdf_breakpoints"
            ),
            "target_training_core_fraction": target_core_fraction,
            "neighbor_manhattan_radius": neighbor_radius,
            "regularization_ess": regularization_ess,
            "minimum_axis_width": minimum_axis_width,
            "minimum_training_rows_for_supported_label": (
                minimum_training_rows
            ),
            "minimum_training_ess_for_supported_label": (
                minimum_training_ess
            ),
            "minimum_validation_coverage": minimum_validation_coverage,
            "padding_candidates": list(paddings),
            "training_status_counts": status_counts,
            "fit_source_counts": fit_source_counts,
            "global_fallback_box": global_box.metadata(),
        },
        "training": _campaign_metadata(training),
        "validation": _campaign_metadata(validation),
        "inside_training_holdout_relative_difference": (
            (training_sigma - validation_sigma) / validation_sigma
            if validation_sigma != 0.0
            else None
        ),
        "inside_training_holdout_difference_z_score": (
            guards._difference_z_score(
                training_sigma,
                training_sem,
                validation_sigma,
                validation_sem,
            )
            if training_sem is not None and validation_sem is not None
            else None
        ),
        "candidates": candidate_results,
        "ranking_by_heldout_coverage_then_volume": sorted(
            candidate_ids,
            key=lambda candidate: (
                -(
                    candidate_results[candidate]["validation"][
                        "cross_section_weighted_core_coverage"
                    ]
                    or -1.0
                ),
                candidate_results[candidate]["compactness"][
                    "sum_normalized_core_volume"
                ],
            ),
        ),
        "recommended_development_candidate": recommended,
        "recommendation_reason": recommendation_reason,
        "development_data_note": (
            "The held-out replicas used here become development data once "
            "their results influence padding choice. Final assessment must "
            "use fresh survey replicas."
        ),
        "proposal_support_requirement": {
            "continuous_core_is_not_a_hard_physics_cut": True,
            "future_full_support_tail_required": True,
            "exact_mixture_density_correction_required": True,
            "tail_definition": (
                "nonzero legacy proposal density over the complete original "
                "radiative proposal domain"
            ),
        },
    }
    recipe_payload = {
        "schema": RECIPE_SCHEMA,
        "created_utc": comparison["created_utc"],
        "production_ready": False,
        "generator_revision": generator_revision,
        "continuous_guard_learner_revision": learner_revision,
        "analysis_config_source": str(config_path.resolve()),
        "analysis_config_sha256": config_sha256,
        "training_replica_ids": comparison["training"]["replica_ids"],
        "target_training_core_fraction": target_core_fraction,
        "padding_candidates": [
            {
                "identifier": candidate,
                "base_padding": padding,
            }
            for candidate, padding in zip(candidate_ids, paddings)
        ],
        "recommended_development_candidate": recommended,
        "strata": recipes,
    }
    return comparison, recipe_payload, rows


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "stratum_id",
        "flat_index",
        "training_status",
        "fit_source",
        "candidate",
        "base_padding",
        "effective_padding",
        "normalized_core_volume",
        "training_cross_section_microbarn",
        "training_ess",
        "training_core_coverage",
        "validation_cross_section_microbarn",
        "validation_ess",
        "validation_core_coverage",
        "validation_tail_ess",
        "coverage_passed",
    ]
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def compare(args: argparse.Namespace) -> dict[str, object]:
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(
            f"{output} already exists; comparison outputs are immutable"
        )
    config, config_sha256 = guards.load_analysis_config(args.config)
    paddings = tuple(sorted(set(args.padding or DEFAULT_PADDINGS)))
    if any(
        not math.isfinite(value) or not 0.0 <= value <= 0.5
        for value in paddings
    ):
        raise ContinuousGuardError("padding values must lie in [0,0.5]")
    if not paddings:
        raise ContinuousGuardError("at least one padding is required")
    candidate_ids = [_candidate_id(value) for value in paddings]
    if len(candidate_ids) != len(set(candidate_ids)):
        raise ContinuousGuardError(
            "padding values are indistinguishable at manifest precision"
        )
    if args.neighbor_radius < 0:
        raise ContinuousGuardError("neighbor radius must be nonnegative")
    if args.regularization_ess <= 0.0:
        raise ContinuousGuardError("regularization ESS must be positive")
    if args.minimum_training_rows <= 0:
        raise ContinuousGuardError(
            "minimum training rows must be positive"
        )
    if args.minimum_training_ess <= 0.0:
        raise ContinuousGuardError(
            "minimum training ESS must be positive"
        )
    if not 0.0 <= args.minimum_validation_coverage <= 1.0:
        raise ContinuousGuardError(
            "minimum validation coverage must lie in [0,1]"
        )

    training = _load_campaign(
        args.training_survey,
        config,
        apply_y_max=args.apply_y_max,
    )
    validation = _load_campaign(
        args.validation_survey,
        config,
        apply_y_max=args.apply_y_max,
    )
    if training.base.input_signature != validation.base.input_signature:
        raise ContinuousGuardError(
            "training and validation surveys have different legacy inputs"
        )
    guards._compatible_settings(
        training.base.norm_reference, validation.base.norm_reference
    )
    for campaign in (training, validation):
        if not math.isclose(
            float(config["beam_energy"]),
            float(campaign.base.norm_reference["ebeam"]),
            rel_tol=2.0e-7,
        ):
            raise ContinuousGuardError(
                "analysis and survey beam energies differ"
            )
    training_replicas = {
        int(item["replica"]) for item in training.base.replicas
    }
    validation_replicas = {
        int(item["replica"]) for item in validation.base.replicas
    }
    overlap = training_replicas & validation_replicas
    if overlap:
        raise ContinuousGuardError(
            f"validation replicas overlap training replicas: {sorted(overlap)}"
        )

    learner_revision = guards._current_revision()
    generator_revision = args.generator_revision or learner_revision
    comparison, recipes, rows = _build_study(
        config=config,
        config_path=args.config,
        config_sha256=config_sha256,
        training=training,
        validation=validation,
        target_core_fraction=args.target_core_fraction,
        paddings=paddings,
        neighbor_radius=args.neighbor_radius,
        regularization_ess=args.regularization_ess,
        minimum_axis_width=args.minimum_axis_width,
        minimum_training_rows=args.minimum_training_rows,
        minimum_training_ess=args.minimum_training_ess,
        minimum_validation_coverage=args.minimum_validation_coverage,
        apply_y_max=args.apply_y_max,
        iteration=args.iteration,
        generator_revision=generator_revision,
        generator_revision_source=(
            "user_supplied"
            if args.generator_revision
            else "assumed_current_checkout"
        ),
        learner_revision=learner_revision,
    )

    output.mkdir(parents=True)
    recipe_path = output / "continuous_guard_recipes.json"
    recipe_path.write_text(
        json.dumps(recipes, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    guards._write_sha256(recipe_path)
    row_path = output / "continuous_guard_strata.csv"
    _write_rows(row_path, rows)
    guards._write_sha256(row_path)
    comparison["artifacts"] = {
        "recipes": recipe_path.name,
        "recipes_sha256": hashlib.sha256(
            recipe_path.read_bytes()
        ).hexdigest(),
        "stratum_comparison": row_path.name,
        "stratum_comparison_sha256": hashlib.sha256(
            row_path.read_bytes()
        ).hexdigest(),
        "hash_sidecar_suffix": ".sha256",
    }
    comparison_path = output / "continuous_guard_comparison.json"
    comparison_path.write_text(
        json.dumps(comparison, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    guards._write_sha256(comparison_path)
    recommended = comparison["recommended_development_candidate"]
    return {
        "passed": True,
        "schema": COMPARISON_SCHEMA,
        "comparison": str(comparison_path),
        "recipes": str(recipe_path),
        "stratum_comparison": str(row_path),
        "training_replicas": sorted(training_replicas),
        "validation_replicas": sorted(validation_replicas),
        "recommended_development_candidate": recommended,
        "candidates": {
            identifier: {
                "base_padding": record["base_padding"],
                "heldout_core_coverage": record["validation"][
                    "cross_section_weighted_core_coverage"
                ],
                "heldout_tail_ess": record["validation"][
                    "residual_outside_core"
                ]["importance_effective_sample_size"],
                "mean_normalized_core_volume": record["compactness"][
                    "mean_normalized_core_volume"
                ],
                "passed_strata": record["coverage_summary"][
                    "passed_strata"
                ],
                "failed_strata": record["coverage_summary"][
                    "failed_strata"
                ],
            }
            for identifier, record in comparison["candidates"].items()
        },
    }


def plot(args: argparse.Namespace) -> dict[str, object]:
    comparison_path = args.comparison.resolve()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(
            f"{output} already exists; plot outputs are immutable"
        )
    try:
        comparison = json.loads(
            comparison_path.read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as error:
        raise ContinuousGuardError(
            f"cannot read continuous-guard comparison: {error}"
        ) from error
    if comparison.get("schema") != COMPARISON_SCHEMA:
        raise ContinuousGuardError("continuous-guard comparison has wrong schema")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError as error:
        raise ContinuousGuardError(
            "plotting requires matplotlib"
        ) from error

    candidates = comparison["candidates"]
    identifiers = sorted(
        candidates,
        key=lambda item: float(candidates[item]["base_padding"]),
    )
    padding = [
        float(candidates[item]["base_padding"]) for item in identifiers
    ]
    coverage = [
        float(
            candidates[item]["validation"][
                "cross_section_weighted_core_coverage"
            ]
        )
        for item in identifiers
    ]
    volumes = [
        float(
            candidates[item]["compactness"][
                "mean_normalized_core_volume"
            ]
        )
        for item in identifiers
    ]
    tail_ess = [
        float(
            candidates[item]["validation"]["residual_outside_core"][
                "importance_effective_sample_size"
            ]
        )
        for item in identifiers
    ]
    passed = [
        int(candidates[item]["coverage_summary"]["passed_strata"])
        for item in identifiers
    ]
    failed = [
        int(candidates[item]["coverage_summary"]["failed_strata"])
        for item in identifiers
    ]
    minimum = float(
        comparison["learning"]["minimum_validation_coverage"]
    )

    output.mkdir(parents=True)
    pdf_path = output / "continuous_guard_comparison.pdf"
    with PdfPages(pdf_path) as pdf:
        figure, axis = plt.subplots(figsize=(11, 5.5))
        axis.plot(padding, coverage, marker="o")
        axis.axhline(
            minimum,
            color="black",
            linestyle="--",
            label=f"Requested minimum ({100.0 * minimum:.1f}%)",
        )
        axis.set_xlabel("Base outward padding in normalized coordinates")
        axis.set_ylabel("Cross-section-weighted held-out coverage")
        axis.set_title("Continuous native-coordinate guard coverage")
        axis.grid(alpha=0.25)
        axis.legend()
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        figure, axis = plt.subplots(figsize=(11, 5.5))
        axis.plot(volumes, coverage, marker="o")
        axis.text(
            0.02,
            0.98,
            "Padding sequence: "
            + ", ".join(f"{value:g}" for value in padding),
            transform=axis.transAxes,
            va="top",
            bbox={
                "boxstyle": "round",
                "facecolor": "white",
                "edgecolor": "0.8",
                "alpha": 0.9,
            },
        )
        axis.ticklabel_format(
            axis="x", style="sci", scilimits=(-3, 3)
        )
        axis.set_xlabel("Mean normalized core volume per analysis stratum")
        axis.set_ylabel("Cross-section-weighted held-out coverage")
        axis.set_title("Coverage versus continuous core volume")
        axis.grid(alpha=0.25)
        figure.tight_layout(pad=1.4)
        pdf.savefig(figure)
        plt.close(figure)

        figure, axis = plt.subplots(figsize=(11, 5.5))
        axis.plot(padding, tail_ess, marker="o")
        axis.set_xlabel("Base outward padding in normalized coordinates")
        axis.set_ylabel("Held-out residual-tail effective sample size")
        axis.set_title("Statistical support remaining outside the core")
        axis.grid(alpha=0.25)
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        figure, axis = plt.subplots(figsize=(11, 5.5))
        positions = range(len(identifiers))
        labels = [f"{value:g}" for value in padding]
        axis.bar(positions, passed, label="At or above threshold")
        axis.bar(
            positions,
            failed,
            bottom=passed,
            label="Below threshold",
        )
        axis.set_xticks(list(positions), labels)
        axis.set_xlabel("Base outward padding")
        axis.set_ylabel("Held-out analysis strata")
        axis.set_title("Per-stratum continuous-guard coverage")
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

    comparison = subparsers.add_parser(
        "compare",
        help="fit native-coordinate boxes and compare padding candidates",
    )
    comparison.add_argument("--config", type=Path, required=True)
    comparison.add_argument(
        "--training-survey", type=Path, nargs="+", required=True
    )
    comparison.add_argument(
        "--validation-survey", type=Path, nargs="+", required=True
    )
    comparison.add_argument("--output", type=Path, required=True)
    comparison.add_argument(
        "--target-core-fraction", type=float, default=0.995
    )
    comparison.add_argument(
        "--padding",
        type=float,
        action="append",
        help=(
            "repeat to replace the default base-padding scan "
            "0,0.0025,0.005,0.01,0.02"
        ),
    )
    comparison.add_argument("--neighbor-radius", type=int, default=1)
    comparison.add_argument(
        "--regularization-ess", type=float, default=5.0
    )
    comparison.add_argument(
        "--minimum-axis-width", type=float, default=1.0e-4
    )
    comparison.add_argument(
        "--minimum-training-rows", type=int, default=10
    )
    comparison.add_argument(
        "--minimum-training-ess", type=float, default=5.0
    )
    comparison.add_argument(
        "--minimum-validation-coverage", type=float, default=0.98
    )
    comparison.add_argument("--iteration", type=int, default=0)
    comparison.add_argument(
        "--apply-y-max",
        action="store_true",
        help=(
            "apply phase_space.y_max in training and validation; the "
            "default applies no y cut"
        ),
    )
    comparison.add_argument(
        "--generator-revision",
        help=(
            "Git revision that produced the surveys; defaults to the current "
            "checkout and is recorded as an assumption"
        ),
    )

    plotting = subparsers.add_parser(
        "plot", help="render the continuous-guard padding comparison"
    )
    plotting.add_argument("--comparison", type=Path, required=True)
    plotting.add_argument("--output", type=Path, required=True)
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
        guards.GuardLearningError,
        ContinuousGuardError,
    ) as error:
        print(f"ERROR: {error}")
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result.get("passed", True) else 2


if __name__ == "__main__":
    raise SystemExit(main())
