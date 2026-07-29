#!/usr/bin/env python3
"""Learn, validate, and compare hard-parent radiative migration footprints.

The milestone-2b/2c diagnostics relate the Born-like hard-vertex coordinates
of an unrestricted fixed-trial survey to the final-LUND analysis stratum and
compare alternative intreg groupings.  They are diagnostic only: they do not
change AAO sampling, perform unweighting, or emit LUND events.
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
REPRESENTATION_COMPARISON_SCHEMA = (
    "aao-rad-migration-representation-comparison-v2"
)
LEGACY_REPRESENTATION_COMPARISON_SCHEMA = (
    "aao-rad-migration-representation-comparison-v1"
)
REPRESENTATION_FOOTPRINT_SCHEMA = (
    "aao-rad-migration-representation-footprints-v1"
)
REPRESENTATION_PLOT_SCHEMA = (
    "aao-rad-migration-representation-plots-v2"
)
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


@dataclass(frozen=True)
class ChannelGroup:
    identifier: str
    label: str
    channels: tuple[int, ...]


@dataclass(frozen=True)
class ChannelRepresentation:
    identifier: str
    label: str
    description: str
    groups: tuple[ChannelGroup, ...]

    def group_for_channel(self, channel: int) -> str:
        for group in self.groups:
            if channel in group.channels:
                return group.identifier
        raise MigrationError(
            f"{self.identifier}: intreg={channel} is not assigned"
        )

    def channels_for_group(self, identifier: str) -> tuple[int, ...]:
        for group in self.groups:
            if identifier == group.identifier:
                return group.channels
        raise MigrationError(
            f"{self.identifier}: unknown channel group {identifier!r}"
        )

    def metadata(self) -> dict[str, object]:
        return {
            "identifier": self.identifier,
            "label": self.label,
            "description": self.description,
            "channel_groups": [
                {
                    "identifier": group.identifier,
                    "label": group.label,
                    "native_intreg": list(group.channels),
                }
                for group in self.groups
            ],
        }


CHANNEL_REPRESENTATIONS = (
    ChannelRepresentation(
        identifier="six_channel",
        label="Six native intreg channels",
        description=(
            "Milestone-2b baseline: every hard-parent/intreg pair is learned "
            "independently."
        ),
        groups=tuple(
            ChannelGroup(
                identifier=f"intreg_{channel}",
                label=f"intreg {channel}",
                channels=(channel,),
            )
            for channel in range(1, 7)
        ),
    ),
    ChannelRepresentation(
        identifier="four_group",
        label="Four radiative-region groups",
        description=(
            "Combines each narrow photon-angle peak with its surrounding "
            "region, while retaining wide-angle and soft groups."
        ),
        groups=(
            ChannelGroup("peak_a", "peak A: intreg 1 and 3", (1, 3)),
            ChannelGroup("peak_b", "peak B: intreg 2 and 4", (2, 4)),
            ChannelGroup("wide_angle", "wide angle: intreg 5", (5,)),
            ChannelGroup("soft", "soft: intreg 6", (6,)),
        ),
    ),
    ChannelRepresentation(
        identifier="soft_resolved",
        label="Soft versus resolved",
        description=(
            "Learns one spatial footprint for resolved radiation and a "
            "separate footprint for the soft branch."
        ),
        groups=(
            ChannelGroup(
                "resolved",
                "resolved: intreg 1 through 5",
                (1, 2, 3, 4, 5),
            ),
            ChannelGroup("soft", "soft: intreg 6", (6,)),
        ),
    ),
    ChannelRepresentation(
        identifier="channel_marginalized",
        label="Channel marginalized",
        description=(
            "Learns hard-parent geometry after summing over intreg; every "
            "native channel remains eligible in each selected hard cell."
        ),
        groups=(
            ChannelGroup(
                "all_channels",
                "all native intreg channels",
                (1, 2, 3, 4, 5, 6),
            ),
        ),
    ),
)


def _validate_channel_representations() -> None:
    identifiers = [
        representation.identifier
        for representation in CHANNEL_REPRESENTATIONS
    ]
    if len(set(identifiers)) != len(identifiers):
        raise MigrationError("channel representation IDs must be unique")
    expected = list(range(1, 7))
    for representation in CHANNEL_REPRESENTATIONS:
        group_ids = [group.identifier for group in representation.groups]
        if len(set(group_ids)) != len(group_ids):
            raise MigrationError(
                f"{representation.identifier}: channel-group IDs repeat"
            )
        channels = [
            channel
            for group in representation.groups
            for channel in group.channels
        ]
        if sorted(channels) != expected:
            raise MigrationError(
                f"{representation.identifier}: channel groups must partition "
                "native intreg values 1 through 6"
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
    for neighbor in _parent_neighbors(parent_id, config_or_grid):
        yield neighbor, channel


def _parent_neighbors(
    parent_id: str, config_or_grid: dict | ParentGrid
) -> Iterator[str]:
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
            yield parent_identifier(neighbor_parent, config_or_grid)


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


def dilate_parent_groups(
    seeds: set[tuple[str, str]],
    config_or_grid: dict | ParentGrid,
    radius: int,
) -> set[tuple[str, str]]:
    """Dilate spatial parent cells without crossing channel groups."""
    if radius < 0:
        raise MigrationError("parent dilation must be nonnegative")
    result = set(seeds)
    frontier = set(seeds)
    for _ in range(radius):
        additions = {
            (neighbor, group)
            for parent_id, group in frontier
            for neighbor in _parent_neighbors(
                parent_id, config_or_grid
            )
        }
        additions -= result
        result.update(additions)
        frontier = additions
        if not frontier:
            break
    return result


def _grouped_moments(
    moments: dict[tuple[str, int], guards.Moment],
    representation: ChannelRepresentation,
) -> dict[tuple[str, str], guards.Moment]:
    grouped: dict[tuple[str, str], list[guards.Moment]] = {}
    for (parent_id, channel), moment in moments.items():
        key = (parent_id, representation.group_for_channel(channel))
        grouped.setdefault(key, []).append(moment)
    return {
        key: _combined(parts)
        for key, parts in grouped.items()
    }


def _select_seed_parent_groups(
    moments: dict[tuple[str, str], guards.Moment],
    target_fraction: float,
) -> set[tuple[str, str]]:
    total = sum(moment.total for moment in moments.values())
    if total <= 0.0:
        return set()
    selected: set[tuple[str, str]] = set()
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


def _native_components_for_parent_groups(
    parent_groups: set[tuple[str, str]],
    representation: ChannelRepresentation,
) -> set[tuple[str, int]]:
    return {
        (parent_id, channel)
        for parent_id, group in parent_groups
        for channel in representation.channels_for_group(group)
    }


def _serialized_parent_groups(
    parent_groups: set[tuple[str, str]],
) -> list[dict[str, str]]:
    return [
        {"parent_id": parent_id, "channel_group": group}
        for parent_id, group in sorted(parent_groups)
    ]


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


def _add_moment(target: guards.Moment, source: guards.Moment) -> None:
    target.count += source.count
    target.total += source.total
    target.square_total += source.square_total
    target.maximum = max(target.maximum, source.maximum)


def _coverage_from_moments(
    *,
    total: guards.Moment,
    selected: guards.Moment,
    tail: guards.Moment,
    expanded: guards.Moment,
    proposals: int,
) -> dict[str, object]:
    return {
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


def _empty_coverage_accumulator() -> dict[str, guards.Moment]:
    return {
        "total": guards.Moment(),
        "selected": guards.Moment(),
        "tail": guards.Moment(),
        "expanded": guards.Moment(),
    }


def _accumulate_coverage(
    accumulator: dict[str, guards.Moment],
    *,
    total: guards.Moment,
    selected: guards.Moment,
    tail: guards.Moment,
    expanded: guards.Moment,
) -> None:
    for name, moment in (
        ("total", total),
        ("selected", selected),
        ("tail", tail),
        ("expanded", expanded),
    ):
        _add_moment(accumulator[name], moment)


def _serialized_coverage(
    accumulator: dict[str, guards.Moment], proposals: int
) -> dict[str, object]:
    return _coverage_from_moments(
        total=accumulator["total"],
        selected=accumulator["selected"],
        tail=accumulator["tail"],
        expanded=accumulator["expanded"],
        proposals=proposals,
    )


def _purity_denominator(
    parent_totals: dict[tuple[str, int], guards.Moment],
    selected: set[tuple[str, int]],
) -> float:
    return sum(
        parent_totals.get(component, guards.Moment()).total
        for component in selected
    )


def _representation_material_coverage(
    material_records: list[dict[str, object]],
    representation: str,
    *,
    targets: tuple[float, ...] = (0.5, 0.9, 0.95, 0.99),
) -> dict[str, object]:
    contributing = [
        record
        for record in material_records
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
        retained = selected = expanded = 0.0
        count = 0
        for record in contributing:
            count += 1
            retained += float(record["validation_total"])
            values = record["representations"][representation]
            selected += float(values["selected"])
            expanded += float(values["expanded"])
            if retained >= target * total:
                break
        result[f"top_{100.0 * target:g}_percent_cross_section"] = {
            "strata": count,
            "actual_fraction_of_inside_cross_section": (
                retained / total if total > 0.0 else None
            ),
            "cross_section_weighted_selected_parent_fraction": (
                selected / retained if retained > 0.0 else None
            ),
            "cross_section_weighted_one_more_dilation_fraction": (
                expanded / retained if retained > 0.0 else None
            ),
        }
    return result


def _campaign_comparison_metadata(
    campaign: MigrationCampaign,
) -> dict[str, object]:
    inside = _combined(campaign.strata.values())
    return {
        "replica_ids": [item["replica"] for item in campaign.replicas],
        "replicas": campaign.replicas,
        "total_proposals": campaign.proposals,
        "recorded_rows": campaign.rows,
        "final_valid_rows": campaign.final_valid_rows,
        "global_observed": guards._compact_metrics(
            campaign.global_observed, campaign.proposals
        ),
        "inside_analysis_partition": guards._compact_metrics(
            inside, campaign.proposals
        ),
        "hard_parent_relationships_inside_analysis": _fraction_metrics(
            campaign.relationships, campaign.proposals
        ),
        "radiative_channels_inside_analysis": _fraction_metrics(
            campaign.channels, campaign.proposals
        ),
    }


def _target_boundary_position(index: int, bins: int) -> str:
    if bins <= 1:
        return "both_boundaries"
    if index == 0:
        return "lower_boundary"
    if index == bins - 1:
        return "upper_boundary"
    return "interior"


def _migration_offset_key(
    parent_id: str,
    target: guards.Stratum,
    parent_grid: ParentGrid,
) -> tuple[str, str, str, str, int, int, int, int]:
    parent = parent_from_identifier(parent_id, parent_grid)
    metadata = parent_metadata(parent, parent_grid)["coordinates"]
    relationship, _ = parent_relationship(
        parent, target, parent_grid
    )
    phi_bins = len(parent_grid.phi_edges) - 1
    signed_phi_delta = (
        (parent.iphi - target.iphi + phi_bins // 2) % phi_bins
        - phi_bins // 2
    )
    return (
        relationship,
        str(metadata["Q2"]["region"]),
        str(metadata["xB"]["region"]),
        str(metadata["minus_t"]["region"]),
        parent.iq2 - target.iq2,
        parent.ixb - target.ixb,
        parent.it - target.it,
        signed_phi_delta,
    )


def _offset_key_metadata(
    key: tuple[str, str, str, str, int, int, int, int],
) -> dict[str, object]:
    (
        relationship,
        q2_region,
        xb_region,
        t_region,
        delta_q2,
        delta_xb,
        delta_t,
        delta_phi,
    ) = key
    return {
        "relationship": relationship,
        "hard_q2_region": q2_region,
        "hard_xb_region": xb_region,
        "hard_minus_t_region": t_region,
        "delta_q2_index": delta_q2,
        "delta_xb_index": delta_xb,
        "delta_minus_t_index": delta_t,
        "delta_phi_index": delta_phi,
    }


def _empty_offset_accumulator() -> dict[str, guards.Moment]:
    return {
        "total": guards.Moment(),
        "selected": guards.Moment(),
        "frozen_tail": guards.Moment(),
        "expanded": guards.Moment(),
        "recovered": guards.Moment(),
        "residual": guards.Moment(),
    }


def _add_offset_accumulator(
    target: dict[str, guards.Moment],
    source: dict[str, guards.Moment],
) -> None:
    for name in target:
        _add_moment(target[name], source[name])


def _offset_metrics(
    accumulator: dict[str, guards.Moment],
    *,
    proposals: int,
    inside_total: float,
    residual_total: float,
) -> dict[str, object]:
    total = accumulator["total"]
    selected = accumulator["selected"]
    expanded = accumulator["expanded"]
    return {
        "total": guards._compact_metrics(total, proposals),
        "selected_parents": guards._compact_metrics(
            selected, proposals
        ),
        "frozen_tail": guards._compact_metrics(
            accumulator["frozen_tail"], proposals
        ),
        "one_more_dilation": guards._compact_metrics(
            expanded, proposals
        ),
        "recovered_by_one_more_dilation": guards._compact_metrics(
            accumulator["recovered"], proposals
        ),
        "residual_after_one_more_dilation": guards._compact_metrics(
            accumulator["residual"], proposals
        ),
        "frozen_parent_fraction": (
            selected.total / total.total if total.total > 0.0 else None
        ),
        "one_more_dilation_fraction": (
            expanded.total / total.total if total.total > 0.0 else None
        ),
        "fraction_of_inside_analysis_cross_section": (
            total.total / inside_total if inside_total > 0.0 else None
        ),
        "fraction_of_residual_after_one_more_dilation": (
            accumulator["residual"].total / residual_total
            if residual_total > 0.0
            else None
        ),
    }


def _residual_offset_summary(
    detailed: dict[
        tuple[
            str,
            int,
            str,
            str,
            str,
            str,
            str,
            str,
            str,
            int,
            int,
            int,
            int,
        ],
        dict[str, guards.Moment],
    ],
    *,
    proposals: int,
    limit: int = 50,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    pooled: dict[
        tuple[str, str, str, str, int, int, int, int],
        dict[str, guards.Moment],
    ] = {}
    by_status = {
        status: _empty_offset_accumulator()
        for status in (
            "learned",
            "learned_low_support",
            "no_training_contribution",
        )
    }
    by_channel = {
        channel: _empty_offset_accumulator()
        for channel in range(1, 7)
    }
    by_boundary: dict[
        str, dict[str, dict[str, guards.Moment]]
    ] = {
        "Q2": {},
        "xB": {},
        "minus_t": {},
    }
    rows: list[dict[str, object]] = []
    total_accumulator = _empty_offset_accumulator()

    for key, accumulator in detailed.items():
        (
            status,
            channel,
            q2_position,
            xb_position,
            t_position,
            *offset_values,
        ) = key
        offset_key = tuple(offset_values)
        pooled_entry = pooled.setdefault(
            offset_key, _empty_offset_accumulator()
        )
        _add_offset_accumulator(pooled_entry, accumulator)
        status_entry = by_status[status]
        _add_offset_accumulator(status_entry, accumulator)
        channel_entry = by_channel[channel]
        _add_offset_accumulator(channel_entry, accumulator)
        for axis, position in (
            ("Q2", q2_position),
            ("xB", xb_position),
            ("minus_t", t_position),
        ):
            boundary_entry = by_boundary[axis].setdefault(
                position, _empty_offset_accumulator()
            )
            _add_offset_accumulator(boundary_entry, accumulator)
        _add_offset_accumulator(total_accumulator, accumulator)

    inside_total = total_accumulator["total"].total
    residual_total = total_accumulator["residual"].total
    for key, accumulator in sorted(detailed.items()):
        (
            status,
            channel,
            q2_position,
            xb_position,
            t_position,
            *offset_values,
        ) = key
        metrics = _offset_metrics(
            accumulator,
            proposals=proposals,
            inside_total=inside_total,
            residual_total=residual_total,
        )
        row = {
            "training_status": status,
            "native_intreg": channel,
            "target_q2_position": q2_position,
            "target_xb_position": xb_position,
            "target_minus_t_position": t_position,
            **_offset_key_metadata(tuple(offset_values)),
            "contributing_rows": metrics["total"]["contributing_rows"],
            "total_cross_section_microbarn": metrics["total"][
                "cross_section_microbarn"
            ],
            "selected_cross_section_microbarn": metrics[
                "selected_parents"
            ]["cross_section_microbarn"],
            "frozen_tail_cross_section_microbarn": metrics[
                "frozen_tail"
            ]["cross_section_microbarn"],
            "one_more_dilation_cross_section_microbarn": metrics[
                "one_more_dilation"
            ]["cross_section_microbarn"],
            "recovered_cross_section_microbarn": metrics[
                "recovered_by_one_more_dilation"
            ]["cross_section_microbarn"],
            "residual_cross_section_microbarn": metrics[
                "residual_after_one_more_dilation"
            ]["cross_section_microbarn"],
            "frozen_parent_fraction": metrics["frozen_parent_fraction"],
            "one_more_dilation_fraction": metrics[
                "one_more_dilation_fraction"
            ],
            "fraction_of_inside_analysis_cross_section": metrics[
                "fraction_of_inside_analysis_cross_section"
            ],
            "fraction_of_residual_after_one_more_dilation": metrics[
                "fraction_of_residual_after_one_more_dilation"
            ],
        }
        rows.append(row)

    pooled_metrics = [
        {
            **_offset_key_metadata(key),
            **_offset_metrics(
                accumulator,
                proposals=proposals,
                inside_total=inside_total,
                residual_total=residual_total,
            ),
        }
        for key, accumulator in pooled.items()
    ]
    pooled_metrics.sort(
        key=lambda item: (
            -float(
                item["residual_after_one_more_dilation"][
                    "cross_section_microbarn"
                ]
            ),
            item["relationship"],
            item["delta_q2_index"],
            item["delta_xb_index"],
            item["delta_minus_t_index"],
            item["delta_phi_index"],
        )
    )

    marginal_by_axis: dict[str, list[dict[str, object]]] = {}
    axis_positions = {
        "Q2": (1, 4),
        "xB": (2, 5),
        "minus_t": (3, 6),
        "phi_deg": (None, 7),
    }
    for axis, (region_position, delta_position) in axis_positions.items():
        marginal: dict[
            tuple[str, int], dict[str, guards.Moment]
        ] = {}
        for key, accumulator in pooled.items():
            region = (
                "analysis_bin"
                if region_position is None
                else str(key[region_position])
            )
            delta = int(key[delta_position])
            entry = marginal.setdefault(
                (region, delta), _empty_offset_accumulator()
            )
            _add_offset_accumulator(entry, accumulator)
        items = [
            {
                "hard_region": region,
                "delta_index": delta,
                **_offset_metrics(
                    accumulator,
                    proposals=proposals,
                    inside_total=inside_total,
                    residual_total=residual_total,
                ),
            }
            for (region, delta), accumulator in marginal.items()
        ]
        items.sort(
            key=lambda item: (
                -float(
                    item["residual_after_one_more_dilation"][
                        "cross_section_microbarn"
                    ]
                ),
                item["hard_region"],
                item["delta_index"],
            )
        )
        marginal_by_axis[axis] = items

    residual_values = [
        float(
            item["residual_after_one_more_dilation"][
                "cross_section_microbarn"
            ]
        )
        for item in pooled_metrics
    ]
    residual_cross_section = residual_total / proposals
    concentration = {
        f"top_{count}_offsets": (
            sum(residual_values[:count]) / residual_cross_section
            if residual_cross_section > 0.0
            else None
        )
        for count in (1, 5, 10, 20, 50)
    }
    summary = {
        "all_offsets": _offset_metrics(
            total_accumulator,
            proposals=proposals,
            inside_total=inside_total,
            residual_total=residual_total,
        ),
        "top_residual_offsets": [
            item
            for item in pooled_metrics
            if float(
                item["residual_after_one_more_dilation"][
                    "cross_section_microbarn"
                ]
            )
            > 0.0
        ][:limit],
        "residual_concentration": concentration,
        "by_training_status": {
            status: _offset_metrics(
                accumulator,
                proposals=proposals,
                inside_total=inside_total,
                residual_total=residual_total,
            )
            for status, accumulator in sorted(by_status.items())
        },
        "by_native_intreg": {
            f"intreg_{channel}": _offset_metrics(
                accumulator,
                proposals=proposals,
                inside_total=inside_total,
                residual_total=residual_total,
            )
            for channel, accumulator in sorted(by_channel.items())
        },
        "by_target_boundary_position": {
            axis: {
                position: _offset_metrics(
                    accumulator,
                    proposals=proposals,
                    inside_total=inside_total,
                    residual_total=residual_total,
                )
                for position, accumulator in sorted(values.items())
            }
            for axis, values in by_boundary.items()
        },
        "marginal_by_axis_offset": marginal_by_axis,
    }
    return summary, rows


def _build_representation_study(
    *,
    config: dict,
    config_path: Path,
    config_sha256: str,
    training: MigrationCampaign,
    validation: MigrationCampaign,
    target_parent_fraction: float,
    parent_dilation: int,
    iteration: int,
    minimum_training_rows: int,
    minimum_training_ess: float,
    minimum_parent_coverage: float,
    apply_y_max: bool,
    generator_revision: str,
    generator_revision_source: str,
    learner_revision: str,
) -> tuple[
    dict,
    dict,
    list[dict[str, object]],
    list[dict[str, object]],
]:
    parent_grid = ParentGrid.from_config(config)
    catalog = guards.enumerate_strata(config)
    accumulators: dict[str, dict[str, object]] = {}
    for representation in CHANNEL_REPRESENTATIONS:
        accumulators[representation.identifier] = {
            "training": _empty_coverage_accumulator(),
            "validation": _empty_coverage_accumulator(),
            "validation_by_status": {
                status: _empty_coverage_accumulator()
                for status in (
                    "learned",
                    "learned_low_support",
                    "no_training_contribution",
                )
            },
            "validation_by_channel": {
                channel: _empty_coverage_accumulator()
                for channel in range(1, 7)
            },
            "training_purity_numerator": 0.0,
            "training_purity_denominator": 0.0,
            "training_expanded_purity_numerator": 0.0,
            "training_expanded_purity_denominator": 0.0,
            "validation_purity_numerator": 0.0,
            "validation_purity_denominator": 0.0,
            "validation_expanded_purity_numerator": 0.0,
            "validation_expanded_purity_denominator": 0.0,
            "validation_offsets": {},
            "seed_parent_groups": 0,
            "selected_parent_groups": 0,
            "selected_hard_cells": 0,
            "selected_native_components": 0,
            "expanded_parent_groups": 0,
            "expanded_hard_cells": 0,
            "expanded_native_components": 0,
            "footprints": 0,
            "assessed_strata": 0,
            "passed_strata": 0,
            "failed_strata": 0,
            "expanded_passed_strata": 0,
            "expanded_failed_strata": 0,
            "training_strata_without_holdout_contribution": 0,
        }

    footprints: dict[str, object] = {}
    rows: list[dict[str, object]] = []
    material_records: list[dict[str, object]] = []
    learned = low_support = empty = 0

    for stratum in catalog:
        stratum_id = stratum.identifier
        training_total = training.strata.get(
            stratum_id, guards.Moment()
        )
        validation_total = validation.strata.get(
            stratum_id, guards.Moment()
        )
        training_moments = training.migrations.get(stratum_id, {})
        validation_moments = validation.migrations.get(stratum_id, {})
        if training_total.total <= 0.0:
            status = "no_training_contribution"
            empty += 1
        else:
            supported = (
                training_total.count >= minimum_training_rows
                and guards._ess(training_total) >= minimum_training_ess
            )
            status = "learned" if supported else "learned_low_support"
            learned += int(supported)
            low_support += int(not supported)

        footprint_record: dict[str, object] = {
            **guards._stratum_metadata(stratum),
            "training_status": status,
            "training_total": guards._compact_metrics(
                training_total, training.proposals
            ),
            "representations": {},
        }
        material_record: dict[str, object] = {
            "stratum_id": stratum_id,
            "validation_total": validation_total.total,
            "representations": {},
        }
        q2_position = _target_boundary_position(
            stratum.iq2, len(parent_grid.q2_edges) - 1
        )
        xb_position = _target_boundary_position(
            stratum.ixb, len(parent_grid.xb_edges) - 1
        )
        t_position = _target_boundary_position(
            stratum.it, len(parent_grid.minus_t_edges) - 1
        )

        for representation in CHANNEL_REPRESENTATIONS:
            identifier = representation.identifier
            accumulator = accumulators[identifier]
            grouped = _grouped_moments(
                training_moments, representation
            )
            seeds = _select_seed_parent_groups(
                grouped, target_parent_fraction
            )
            selected_groups = dilate_parent_groups(
                seeds, parent_grid, parent_dilation
            )
            expanded_groups = dilate_parent_groups(
                seeds, parent_grid, parent_dilation + 1
            )
            selected_native = _native_components_for_parent_groups(
                selected_groups, representation
            )
            expanded_native = _native_components_for_parent_groups(
                expanded_groups, representation
            )

            training_selected = _subset(
                training_moments, selected_native
            )
            training_tail = _combined(
                moment
                for component, moment in training_moments.items()
                if component not in selected_native
            )
            training_expanded = _subset(
                training_moments, expanded_native
            )
            validation_selected = _subset(
                validation_moments, selected_native
            )
            validation_tail = _combined(
                moment
                for component, moment in validation_moments.items()
                if component not in selected_native
            )
            validation_expanded = _subset(
                validation_moments, expanded_native
            )

            _accumulate_coverage(
                accumulator["training"],
                total=training_total,
                selected=training_selected,
                tail=training_tail,
                expanded=training_expanded,
            )
            _accumulate_coverage(
                accumulator["validation"],
                total=validation_total,
                selected=validation_selected,
                tail=validation_tail,
                expanded=validation_expanded,
            )
            _accumulate_coverage(
                accumulator["validation_by_status"][status],
                total=validation_total,
                selected=validation_selected,
                tail=validation_tail,
                expanded=validation_expanded,
            )

            training_purity_denominator = _purity_denominator(
                training.parent_totals, selected_native
            )
            validation_purity_denominator = _purity_denominator(
                validation.parent_totals, selected_native
            )
            training_expanded_purity_denominator = _purity_denominator(
                training.parent_totals, expanded_native
            )
            validation_expanded_purity_denominator = _purity_denominator(
                validation.parent_totals, expanded_native
            )
            accumulator["training_purity_numerator"] += (
                training_selected.total
            )
            accumulator["training_purity_denominator"] += (
                training_purity_denominator
            )
            accumulator["training_expanded_purity_numerator"] += (
                training_expanded.total
            )
            accumulator["training_expanded_purity_denominator"] += (
                training_expanded_purity_denominator
            )
            accumulator["validation_purity_numerator"] += (
                validation_selected.total
            )
            accumulator["validation_purity_denominator"] += (
                validation_purity_denominator
            )
            accumulator["validation_expanded_purity_numerator"] += (
                validation_expanded.total
            )
            accumulator["validation_expanded_purity_denominator"] += (
                validation_expanded_purity_denominator
            )

            if training_total.total > 0.0:
                accumulator["footprints"] += 1
                accumulator["seed_parent_groups"] += len(seeds)
                accumulator["selected_parent_groups"] += len(
                    selected_groups
                )
                accumulator["selected_hard_cells"] += len(
                    {parent_id for parent_id, _ in selected_groups}
                )
                accumulator["selected_native_components"] += len(
                    selected_native
                )
                accumulator["expanded_parent_groups"] += len(
                    expanded_groups
                )
                accumulator["expanded_hard_cells"] += len(
                    {parent_id for parent_id, _ in expanded_groups}
                )
                accumulator["expanded_native_components"] += len(
                    expanded_native
                )

            if validation_total.total > 0.0:
                selected_fraction = (
                    validation_selected.total / validation_total.total
                )
                coverage_passed = (
                    selected_fraction >= minimum_parent_coverage
                )
                expanded_fraction = (
                    validation_expanded.total / validation_total.total
                )
                expanded_coverage_passed = (
                    expanded_fraction >= minimum_parent_coverage
                )
                accumulator["assessed_strata"] += 1
                accumulator["passed_strata"] += int(coverage_passed)
                accumulator["failed_strata"] += int(
                    not coverage_passed
                )
                accumulator["expanded_passed_strata"] += int(
                    expanded_coverage_passed
                )
                accumulator["expanded_failed_strata"] += int(
                    not expanded_coverage_passed
                )
            else:
                selected_fraction = None
                expanded_fraction = None
                coverage_passed = None
                expanded_coverage_passed = None
                if training_total.total > 0.0:
                    accumulator[
                        "training_strata_without_holdout_contribution"
                    ] += 1

            validation_offsets = accumulator["validation_offsets"]
            for component, moment in validation_moments.items():
                parent_id, channel = component
                offset_key = _migration_offset_key(
                    parent_id, stratum, parent_grid
                )
                detailed_key = (
                    status,
                    channel,
                    q2_position,
                    xb_position,
                    t_position,
                    *offset_key,
                )
                offset_accumulator = validation_offsets.setdefault(
                    detailed_key, _empty_offset_accumulator()
                )
                _add_moment(offset_accumulator["total"], moment)
                if component in selected_native:
                    _add_moment(offset_accumulator["selected"], moment)
                else:
                    _add_moment(
                        offset_accumulator["frozen_tail"], moment
                    )
                if component in expanded_native:
                    _add_moment(offset_accumulator["expanded"], moment)
                    if component not in selected_native:
                        _add_moment(
                            offset_accumulator["recovered"], moment
                        )
                else:
                    _add_moment(offset_accumulator["residual"], moment)

            for channel in range(1, 7):
                channel_total = _combined(
                    moment
                    for (_, native_channel), moment
                    in validation_moments.items()
                    if native_channel == channel
                )
                channel_selected = _combined(
                    moment
                    for component, moment in validation_moments.items()
                    if component[1] == channel
                    and component in selected_native
                )
                channel_tail = _combined(
                    moment
                    for component, moment in validation_moments.items()
                    if component[1] == channel
                    and component not in selected_native
                )
                channel_expanded = _combined(
                    moment
                    for component, moment in validation_moments.items()
                    if component[1] == channel
                    and component in expanded_native
                )
                _accumulate_coverage(
                    accumulator["validation_by_channel"][channel],
                    total=channel_total,
                    selected=channel_selected,
                    tail=channel_tail,
                    expanded=channel_expanded,
                )

            training_selected_fraction = (
                training_selected.total / training_total.total
                if training_total.total > 0.0
                else None
            )
            training_expanded_fraction = (
                training_expanded.total / training_total.total
                if training_total.total > 0.0
                else None
            )
            validation_expanded_fraction = (
                expanded_fraction
            )
            training_purity = (
                training_selected.total / training_purity_denominator
                if training_purity_denominator > 0.0
                else None
            )
            validation_purity = (
                validation_selected.total / validation_purity_denominator
                if validation_purity_denominator > 0.0
                else None
            )
            training_expanded_purity = (
                training_expanded.total
                / training_expanded_purity_denominator
                if training_expanded_purity_denominator > 0.0
                else None
            )
            validation_expanded_purity = (
                validation_expanded.total
                / validation_expanded_purity_denominator
                if validation_expanded_purity_denominator > 0.0
                else None
            )

            footprint_record["representations"][identifier] = {
                "seed_parent_groups": _serialized_parent_groups(seeds),
                "dilation_axis_steps": parent_dilation,
                "selected_parent_group_count": len(selected_groups),
                "selected_hard_cell_count": len(
                    {parent_id for parent_id, _ in selected_groups}
                ),
                "selected_native_component_count": len(selected_native),
                "one_more_dilation_parent_group_count": len(
                    expanded_groups
                ),
                "one_more_dilation_hard_cell_count": len(
                    {parent_id for parent_id, _ in expanded_groups}
                ),
                "one_more_dilation_native_component_count": len(
                    expanded_native
                ),
            }
            material_record["representations"][identifier] = {
                "selected": validation_selected.total,
                "expanded": validation_expanded.total,
            }
            rows.append(
                {
                    "stratum_id": stratum_id,
                    "flat_index": stratum.flat_index,
                    "training_status": status,
                    "representation": identifier,
                    "training_cross_section_microbarn": (
                        training_total.total / training.proposals
                    ),
                    "training_ess": guards._ess(training_total),
                    "validation_cross_section_microbarn": (
                        validation_total.total / validation.proposals
                    ),
                    "validation_ess": guards._ess(validation_total),
                    "seed_parent_groups": len(seeds),
                    "selected_parent_groups": len(selected_groups),
                    "selected_hard_cells": len(
                        {parent_id for parent_id, _ in selected_groups}
                    ),
                    "selected_native_components": len(selected_native),
                    "one_more_dilation_parent_groups": len(
                        expanded_groups
                    ),
                    "one_more_dilation_hard_cells": len(
                        {parent_id for parent_id, _ in expanded_groups}
                    ),
                    "one_more_dilation_native_components": len(
                        expanded_native
                    ),
                    "training_selected_parent_fraction": (
                        training_selected_fraction
                    ),
                    "training_one_more_dilation_fraction": (
                        training_expanded_fraction
                    ),
                    "training_selected_parent_purity": training_purity,
                    "training_one_more_dilation_purity": (
                        training_expanded_purity
                    ),
                    "validation_selected_parent_fraction": (
                        selected_fraction
                    ),
                    "validation_one_more_dilation_fraction": (
                        validation_expanded_fraction
                    ),
                    "validation_selected_parent_purity": (
                        validation_purity
                    ),
                    "validation_one_more_dilation_purity": (
                        validation_expanded_purity
                    ),
                    "coverage_passed": coverage_passed,
                    "one_more_dilation_coverage_passed": (
                        expanded_coverage_passed
                    ),
                }
            )

        footprints[stratum_id] = footprint_record
        material_records.append(material_record)

    representation_results: dict[str, object] = {}
    residual_offset_rows: list[dict[str, object]] = []
    for representation in CHANNEL_REPRESENTATIONS:
        identifier = representation.identifier
        accumulator = accumulators[identifier]
        training_coverage = _serialized_coverage(
            accumulator["training"], training.proposals
        )
        validation_coverage = _serialized_coverage(
            accumulator["validation"], validation.proposals
        )
        footprints_count = int(accumulator["footprints"])
        compactness = {
            "contributing_training_strata": footprints_count,
        }
        for name in (
            "seed_parent_groups",
            "selected_parent_groups",
            "selected_hard_cells",
            "selected_native_components",
            "expanded_parent_groups",
            "expanded_hard_cells",
            "expanded_native_components",
        ):
            value = int(accumulator[name])
            compactness[f"total_{name}"] = value
            compactness[f"mean_{name}_per_contributing_stratum"] = (
                value / footprints_count
                if footprints_count > 0
                else None
            )

        training_purity_denominator = float(
            accumulator["training_purity_denominator"]
        )
        validation_purity_denominator = float(
            accumulator["validation_purity_denominator"]
        )
        training_expanded_purity_denominator = float(
            accumulator["training_expanded_purity_denominator"]
        )
        validation_expanded_purity_denominator = float(
            accumulator["validation_expanded_purity_denominator"]
        )
        per_channel: dict[str, object] = {}
        validation_inside_total = accumulator["validation"]["total"].total
        for channel in range(1, 7):
            channel_coverage = _serialized_coverage(
                accumulator["validation_by_channel"][channel],
                validation.proposals,
            )
            channel_total = accumulator["validation_by_channel"][
                channel
            ]["total"].total
            channel_coverage[
                "fraction_of_inside_analysis_cross_section"
            ] = (
                channel_total / validation_inside_total
                if validation_inside_total > 0.0
                else None
            )
            per_channel[f"intreg_{channel}"] = channel_coverage

        residual_summary, representation_offset_rows = (
            _residual_offset_summary(
                accumulator["validation_offsets"],
                proposals=validation.proposals,
            )
        )
        residual_offset_rows.extend(
            {
                "representation": identifier,
                **row,
            }
            for row in representation_offset_rows
        )

        representation_results[identifier] = {
            **representation.metadata(),
            "training": {
                **training_coverage,
                "aggregate_selected_parent_purity_proxy": (
                    float(accumulator["training_purity_numerator"])
                    / training_purity_denominator
                    if training_purity_denominator > 0.0
                    else None
                ),
                "aggregate_one_more_dilation_purity_proxy": (
                    float(
                        accumulator[
                            "training_expanded_purity_numerator"
                        ]
                    )
                    / training_expanded_purity_denominator
                    if training_expanded_purity_denominator > 0.0
                    else None
                ),
            },
            "validation": {
                **validation_coverage,
                "aggregate_selected_parent_purity_proxy": (
                    float(accumulator["validation_purity_numerator"])
                    / validation_purity_denominator
                    if validation_purity_denominator > 0.0
                    else None
                ),
                "aggregate_one_more_dilation_purity_proxy": (
                    float(
                        accumulator[
                            "validation_expanded_purity_numerator"
                        ]
                    )
                    / validation_expanded_purity_denominator
                    if validation_expanded_purity_denominator > 0.0
                    else None
                ),
            },
            "compactness": compactness,
            "validation_by_training_status": {
                status: _serialized_coverage(
                    accumulator["validation_by_status"][status],
                    validation.proposals,
                )
                for status in (
                    "learned",
                    "learned_low_support",
                    "no_training_contribution",
                )
            },
            "validation_by_native_intreg": per_channel,
            "validation_material_strata": (
                _representation_material_coverage(
                    material_records, identifier
                )
            ),
            "validation_residual_offsets": residual_summary,
            "coverage_summary": {
                "minimum_parent_coverage": minimum_parent_coverage,
                "assessed_strata": int(accumulator["assessed_strata"]),
                "passed_strata": int(accumulator["passed_strata"]),
                "failed_strata": int(accumulator["failed_strata"]),
                "one_more_dilation_passed_strata": int(
                    accumulator["expanded_passed_strata"]
                ),
                "one_more_dilation_failed_strata": int(
                    accumulator["expanded_failed_strata"]
                ),
                "training_strata_without_holdout_contribution": int(
                    accumulator[
                        "training_strata_without_holdout_contribution"
                    ]
                ),
                "all_assessed_strata_passed": (
                    int(accumulator["assessed_strata"]) > 0
                    and int(accumulator["failed_strata"]) == 0
                ),
                "all_assessed_strata_passed_after_one_more_dilation": (
                    int(accumulator["assessed_strata"]) > 0
                    and int(accumulator["expanded_failed_strata"]) == 0
                ),
                "aggregate_coverage_meets_minimum": (
                    validation_coverage[
                        "cross_section_weighted_selected_parent_fraction"
                    ]
                    is not None
                    and float(
                        validation_coverage[
                            "cross_section_weighted_selected_parent_fraction"
                        ]
                    )
                    >= minimum_parent_coverage
                ),
                "aggregate_one_more_dilation_coverage_meets_minimum": (
                    validation_coverage[
                        "cross_section_weighted_one_more_dilation_fraction"
                    ]
                    is not None
                    and float(
                        validation_coverage[
                            "cross_section_weighted_one_more_dilation_fraction"
                        ]
                    )
                    >= minimum_parent_coverage
                ),
            },
        }

    ranked = sorted(
        representation_results,
        key=lambda identifier: (
            -float(
                representation_results[identifier]["validation"][
                    "cross_section_weighted_selected_parent_fraction"
                ]
                or 0.0
            ),
            int(
                representation_results[identifier]["compactness"][
                    "total_selected_native_components"
                ]
            ),
            identifier,
        ),
    )
    training_metadata = _campaign_comparison_metadata(training)
    validation_metadata = _campaign_comparison_metadata(validation)
    training_inside = training_metadata["inside_analysis_partition"]
    validation_inside = validation_metadata["inside_analysis_partition"]
    comparison = {
        "schema": REPRESENTATION_COMPARISON_SCHEMA,
        "created_utc": _now(),
        "passed": True,
        "study_completed": True,
        "production_ready": False,
        "production_readiness_note": (
            "Milestone 2c compares frozen guard representations only. It "
            "does not activate mode 4 or change AAO channel sampling."
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
            "hard_parent": HARD_COORDINATE_DEFINITION,
            "hard_columns": HARD_COLUMNS,
            "observed_target": OBSERVED_COORDINATE_DEFINITION,
        },
        "generator_revision": generator_revision,
        "generator_revision_source": generator_revision_source,
        "migration_learner_revision": learner_revision,
        "survey_schema": radiative_survey.SURVEY_SCHEMA,
        "representation_iteration": iteration,
        "legacy_input_sha256": training.input_signature,
        "legacy_generator_settings": training.legacy_input_settings,
        "legacy_channel_probabilities": training.channel_probabilities,
        "learning": {
            "target_parent_fraction_before_dilation": (
                target_parent_fraction
            ),
            "parent_dilation_axis_steps": parent_dilation,
            "minimum_training_rows_for_supported_label": (
                minimum_training_rows
            ),
            "minimum_training_ess_for_supported_label": (
                minimum_training_ess
            ),
            "minimum_validation_parent_coverage": (
                minimum_parent_coverage
            ),
            "learned_strata": learned,
            "learned_low_support_strata": low_support,
            "no_training_contribution_strata": empty,
        },
        "representation_interpretation": {
            "channel_grouping_scope": (
                "Grouping changes only how spatial footprints are learned. "
                "Native intreg probabilities and Jacobians remain unchanged."
            ),
            "purity_proxy": (
                "Sum of selected target contributions divided by the sum of "
                "global observed contributions from each target's selected "
                "parents. Parent contributions are intentionally counted once "
                "per target footprint."
            ),
            "selection_ranking": (
                "Representations are ranked first by held-out selected-parent "
                "coverage, then by fewer selected native components."
            ),
        },
        "training": training_metadata,
        "validation": validation_metadata,
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
        "representations": representation_results,
        "ranking_by_heldout_coverage_then_compactness": ranked,
        "best_heldout_coverage_representation": (
            ranked[0] if ranked else None
        ),
    }
    footprint_payload = {
        "schema": REPRESENTATION_FOOTPRINT_SCHEMA,
        "created_utc": comparison["created_utc"],
        "production_ready": False,
        "analysis_config": config,
        "analysis_config_source": str(config_path.resolve()),
        "analysis_config_sha256": config_sha256,
        "generator_revision": generator_revision,
        "migration_learner_revision": learner_revision,
        "representation_iteration": iteration,
        "legacy_input_sha256": training.input_signature,
        "legacy_generator_settings": training.legacy_input_settings,
        "legacy_channel_probabilities": training.channel_probabilities,
        "training_replica_ids": training_metadata["replica_ids"],
        "target_parent_fraction_before_dilation": target_parent_fraction,
        "parent_dilation_axis_steps": parent_dilation,
        "representations": {
            representation.identifier: representation.metadata()
            for representation in CHANNEL_REPRESENTATIONS
        },
        "strata": footprints,
    }
    return comparison, footprint_payload, rows, residual_offset_rows


def _write_representation_rows(
    path: Path, rows: list[dict[str, object]]
) -> None:
    fields = (
        "stratum_id",
        "flat_index",
        "training_status",
        "representation",
        "training_cross_section_microbarn",
        "training_ess",
        "validation_cross_section_microbarn",
        "validation_ess",
        "seed_parent_groups",
        "selected_parent_groups",
        "selected_hard_cells",
        "selected_native_components",
        "one_more_dilation_parent_groups",
        "one_more_dilation_hard_cells",
        "one_more_dilation_native_components",
        "training_selected_parent_fraction",
        "training_one_more_dilation_fraction",
        "training_selected_parent_purity",
        "training_one_more_dilation_purity",
        "validation_selected_parent_fraction",
        "validation_one_more_dilation_fraction",
        "validation_selected_parent_purity",
        "validation_one_more_dilation_purity",
        "coverage_passed",
        "one_more_dilation_coverage_passed",
    )
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _write_residual_offset_rows(
    path: Path, rows: list[dict[str, object]]
) -> None:
    fields = (
        "representation",
        "training_status",
        "native_intreg",
        "target_q2_position",
        "target_xb_position",
        "target_minus_t_position",
        "relationship",
        "hard_q2_region",
        "hard_xb_region",
        "hard_minus_t_region",
        "delta_q2_index",
        "delta_xb_index",
        "delta_minus_t_index",
        "delta_phi_index",
        "contributing_rows",
        "total_cross_section_microbarn",
        "selected_cross_section_microbarn",
        "frozen_tail_cross_section_microbarn",
        "one_more_dilation_cross_section_microbarn",
        "recovered_cross_section_microbarn",
        "residual_cross_section_microbarn",
        "frozen_parent_fraction",
        "one_more_dilation_fraction",
        "fraction_of_inside_analysis_cross_section",
        "fraction_of_residual_after_one_more_dilation",
    )
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def compare_representations(args: argparse.Namespace) -> dict:
    _validate_channel_representations()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(
            f"{output} already exists; representation artifacts are immutable"
        )
    if not 0.0 < args.target_parent_fraction <= 1.0:
        raise MigrationError("--target-parent-fraction must be in (0,1]")
    if args.parent_dilation < 0:
        raise MigrationError("--parent-dilation must be nonnegative")
    if args.iteration < 0:
        raise MigrationError("--iteration must be nonnegative")
    if args.minimum_training_rows < 1 or args.minimum_training_ess < 0.0:
        raise MigrationError("training support thresholds are invalid")
    if not 0.0 <= args.minimum_parent_coverage <= 1.0:
        raise MigrationError("--minimum-parent-coverage must be in [0,1]")
    try:
        config, config_sha256 = guards.load_analysis_config(
            args.config.resolve()
        )
    except guards.GuardLearningError as error:
        raise MigrationError(str(error)) from error

    training = aggregate_surveys(
        args.training_survey,
        config,
        apply_y_max=args.apply_y_max,
    )
    validation = aggregate_surveys(
        args.validation_survey,
        config,
        apply_y_max=args.apply_y_max,
    )
    if training.input_signature != validation.input_signature:
        raise MigrationError(
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
        raise MigrationError(
            f"validation replicas overlap training replicas: {sorted(overlap)}"
        )
    for campaign in (training, validation):
        if not math.isclose(
            float(config["beam_energy"]),
            float(campaign.norm_reference["ebeam"]),
            rel_tol=2.0e-7,
        ):
            raise MigrationError("analysis and survey beam energies differ")

    learner_revision = guards._current_revision()
    generator_revision = args.generator_revision or learner_revision
    (
        comparison,
        footprints,
        rows,
        residual_offset_rows,
    ) = _build_representation_study(
        config=config,
        config_path=args.config,
        config_sha256=config_sha256,
        training=training,
        validation=validation,
        target_parent_fraction=args.target_parent_fraction,
        parent_dilation=args.parent_dilation,
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
    footprint_path = output / "representation_footprints.json"
    footprint_path.write_text(
        json.dumps(footprints, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    guards._write_sha256(footprint_path)
    comparison["artifacts"] = {
        "frozen_footprints": footprint_path.name,
        "frozen_footprints_sha256": hashlib.sha256(
            footprint_path.read_bytes()
        ).hexdigest(),
        "stratum_comparison": "representation_strata.csv",
        "residual_offset_comparison": (
            "representation_residual_offsets.csv"
        ),
        "hash_sidecar_suffix": ".sha256",
    }
    comparison_path = output / "representation_comparison.json"
    comparison_path.write_text(
        json.dumps(comparison, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    row_path = output / "representation_strata.csv"
    _write_representation_rows(row_path, rows)
    residual_offset_path = (
        output / "representation_residual_offsets.csv"
    )
    _write_residual_offset_rows(
        residual_offset_path, residual_offset_rows
    )
    guards._write_sha256(comparison_path)
    guards._write_sha256(row_path)
    guards._write_sha256(residual_offset_path)
    return {
        "passed": True,
        "schema": REPRESENTATION_COMPARISON_SCHEMA,
        "comparison": str(comparison_path),
        "footprints": str(footprint_path),
        "stratum_comparison": str(row_path),
        "residual_offset_comparison": str(residual_offset_path),
        "training_replicas": sorted(training_replicas),
        "validation_replicas": sorted(validation_replicas),
        "ranking_by_heldout_coverage_then_compactness": comparison[
            "ranking_by_heldout_coverage_then_compactness"
        ],
        "representations": {
            identifier: {
                "heldout_selected_parent_fraction": values["validation"][
                    "cross_section_weighted_selected_parent_fraction"
                ],
                "heldout_one_more_dilation_fraction": values["validation"][
                    "cross_section_weighted_one_more_dilation_fraction"
                ],
                "aggregate_selected_parent_purity_proxy": values[
                    "validation"
                ]["aggregate_selected_parent_purity_proxy"],
                "aggregate_one_more_dilation_purity_proxy": values[
                    "validation"
                ]["aggregate_one_more_dilation_purity_proxy"],
                "one_more_dilation_passed_strata": values[
                    "coverage_summary"
                ]["one_more_dilation_passed_strata"],
                "one_more_dilation_failed_strata": values[
                    "coverage_summary"
                ]["one_more_dilation_failed_strata"],
                "total_selected_hard_cells": values["compactness"][
                    "total_selected_hard_cells"
                ],
                "total_selected_native_components": values["compactness"][
                    "total_selected_native_components"
                ],
            }
            for identifier, values in comparison[
                "representations"
            ].items()
        },
    }


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


def plot_representations(args: argparse.Namespace) -> dict:
    _validate_channel_representations()
    comparison_path = args.comparison.resolve()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(
            f"{output} already exists; representation plots are immutable"
        )
    try:
        comparison = json.loads(
            comparison_path.read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as error:
        raise MigrationError(
            f"cannot read representation comparison: {error}"
        ) from error
    if comparison.get("schema") not in {
        REPRESENTATION_COMPARISON_SCHEMA,
        LEGACY_REPRESENTATION_COMPARISON_SCHEMA,
    }:
        raise MigrationError("representation comparison has the wrong schema")
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError as error:
        raise MigrationError(
            "plotting requires matplotlib"
        ) from error

    identifiers = [
        representation.identifier
        for representation in CHANNEL_REPRESENTATIONS
    ]
    representations = comparison["representations"]
    if set(identifiers) != set(representations):
        raise MigrationError(
            "representation comparison does not contain the expected models"
        )
    labels = [
        str(representations[identifier]["label"])
        for identifier in identifiers
    ]
    positions = list(range(len(identifiers)))
    selected = [
        float(
            representations[identifier]["validation"][
                "cross_section_weighted_selected_parent_fraction"
            ]
            or 0.0
        )
        for identifier in identifiers
    ]
    expanded = [
        float(
            representations[identifier]["validation"][
                "cross_section_weighted_one_more_dilation_fraction"
            ]
            or 0.0
        )
        for identifier in identifiers
    ]
    minimum = float(
        comparison["learning"]["minimum_validation_parent_coverage"]
    )

    output.mkdir(parents=True)
    pdf_path = output / "representation_comparison.pdf"
    with PdfPages(pdf_path) as pdf:
        figure, axis = plt.subplots(figsize=(11, 5.5))
        width = 0.36
        axis.bar(
            [position - width / 2 for position in positions],
            selected,
            width=width,
            label="Frozen footprint",
        )
        axis.bar(
            [position + width / 2 for position in positions],
            expanded,
            width=width,
            label="One more spatial dilation",
        )
        axis.axhline(
            minimum,
            color="black",
            linestyle="--",
            linewidth=1,
            label=f"Requested minimum ({minimum:.1%})",
        )
        axis.set_xticks(positions, labels, rotation=12, ha="right")
        axis.set_ylim(0.0, 1.03)
        axis.set_ylabel(
            "Cross-section-weighted held-out parent coverage"
        )
        axis.set_title("Native-intreg representation comparison")
        axis.grid(axis="y", alpha=0.25)
        axis.legend()
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        hard_cells = [
            int(
                representations[identifier]["compactness"][
                    "total_selected_hard_cells"
                ]
            )
            for identifier in identifiers
        ]
        native_components = [
            int(
                representations[identifier]["compactness"][
                    "total_selected_native_components"
                ]
            )
            for identifier in identifiers
        ]
        purities = [
            float(
                representations[identifier]["validation"][
                    "aggregate_selected_parent_purity_proxy"
                ]
                or 0.0
            )
            for identifier in identifiers
        ]
        expanded_purities = [
            float(
                representations[identifier]["validation"].get(
                    "aggregate_one_more_dilation_purity_proxy",
                    representations[identifier]["validation"][
                        "aggregate_selected_parent_purity_proxy"
                    ],
                )
                or 0.0
            )
            for identifier in identifiers
        ]
        figure, axes = plt.subplots(1, 2, figsize=(12, 5))
        axes[0].bar(
            [position - width / 2 for position in positions],
            hard_cells,
            width=width,
            label="Hard cells",
        )
        axes[0].bar(
            [position + width / 2 for position in positions],
            native_components,
            width=width,
            label="Native components",
        )
        axes[0].set_yscale("log")
        axes[0].set_xticks(positions, labels, rotation=15, ha="right")
        axes[0].set_ylabel("Total selected across training strata")
        axes[0].set_title("Footprint compactness")
        axes[0].grid(axis="y", alpha=0.25)
        axes[0].legend()
        axes[1].bar(
            [position - width / 2 for position in positions],
            purities,
            width=width,
            label="Frozen footprint",
        )
        axes[1].bar(
            [position + width / 2 for position in positions],
            expanded_purities,
            width=width,
            label="One more spatial dilation",
        )
        axes[1].set_xticks(positions, labels, rotation=15, ha="right")
        axes[1].set_ylim(0.0, 1.0)
        axes[1].set_ylabel("Held-out aggregate purity proxy")
        axes[1].set_title("Coverage-versus-purity cost of dilation")
        axes[1].grid(axis="y", alpha=0.25)
        axes[1].legend()
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        best_identifier = comparison.get(
            "best_heldout_coverage_representation"
        )
        if best_identifier in representations:
            residual_offsets = representations[best_identifier].get(
                "validation_residual_offsets", {}
            ).get("top_residual_offsets", [])
            residual_offsets = [
                item
                for item in residual_offsets[:10]
                if float(
                    item["residual_after_one_more_dilation"][
                        "cross_section_microbarn"
                    ]
                )
                > 0.0
            ]
            if residual_offsets:
                residual_offsets.reverse()
                offset_labels = []
                residual_fractions = []
                for item in residual_offsets:
                    label = (
                        "Δ=("
                        f"{item['delta_q2_index']},"
                        f"{item['delta_xb_index']},"
                        f"{item['delta_minus_t_index']},"
                        f"{item['delta_phi_index']})"
                    )
                    if item["relationship"] == (
                        "hard_parent_underflow_or_overflow"
                    ):
                        regions = "/".join(
                            str(item[name])
                            for name in (
                                "hard_q2_region",
                                "hard_xb_region",
                                "hard_minus_t_region",
                            )
                            if item[name] != "analysis_bin"
                        )
                        label = f"{label} [{regions}]"
                    offset_labels.append(label)
                    residual_fractions.append(
                        float(
                            item[
                                "fraction_of_residual_after_one_more_dilation"
                            ]
                            or 0.0
                        )
                    )
                figure, axis = plt.subplots(figsize=(11, 6))
                axis.barh(offset_labels, residual_fractions)
                axis.set_xlim(
                    0.0,
                    max(residual_fractions) * 1.08
                    if residual_fractions
                    else 1.0,
                )
                axis.set_xlabel(
                    "Fraction of residual cross section after one more dilation"
                )
                axis.set_title(
                    "Largest remaining hard-minus-observed offsets\n"
                    f"{representations[best_identifier]['label']}"
                )
                axis.grid(axis="x", alpha=0.25)
                figure.tight_layout()
                pdf.savefig(figure)
                plt.close(figure)

        figure, axis = plt.subplots(figsize=(12, 5.5))
        channel_positions = list(range(1, 7))
        channel_width = 0.18
        for representation_index, identifier in enumerate(identifiers):
            offset = (
                representation_index - (len(identifiers) - 1) / 2
            ) * channel_width
            channel_coverages = [
                float(
                    representations[identifier][
                        "validation_by_native_intreg"
                    ][f"intreg_{channel}"][
                        "cross_section_weighted_selected_parent_fraction"
                    ]
                    or 0.0
                )
                for channel in channel_positions
            ]
            axis.bar(
                [channel + offset for channel in channel_positions],
                channel_coverages,
                width=channel_width,
                label=labels[representation_index],
            )
        axis.set_xticks(channel_positions)
        axis.set_xlabel("Native intreg")
        axis.set_ylabel("Held-out selected-parent coverage")
        axis.set_ylim(0.0, 1.03)
        axis.set_title("Coverage retained separately by native intreg")
        axis.grid(axis="y", alpha=0.25)
        axis.legend(fontsize="small")
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        material_keys = (
            "top_50_percent_cross_section",
            "top_90_percent_cross_section",
            "top_95_percent_cross_section",
            "top_99_percent_cross_section",
        )
        material_labels = ("Top 50%", "Top 90%", "Top 95%", "Top 99%")
        figure, axis = plt.subplots(figsize=(12, 5.5))
        material_width = 0.18
        material_positions = list(range(len(material_keys)))
        for representation_index, identifier in enumerate(identifiers):
            offset = (
                representation_index - (len(identifiers) - 1) / 2
            ) * material_width
            values = [
                float(
                    representations[identifier][
                        "validation_material_strata"
                    ][key][
                        "cross_section_weighted_selected_parent_fraction"
                    ]
                    or 0.0
                )
                for key in material_keys
            ]
            axis.bar(
                [
                    position + offset
                    for position in material_positions
                ],
                values,
                width=material_width,
                label=labels[representation_index],
            )
        axis.set_xticks(material_positions, material_labels)
        axis.set_ylim(0.0, 1.03)
        axis.set_ylabel("Held-out selected-parent coverage")
        axis.set_title("Coverage in cross-section-dominant strata")
        axis.grid(axis="y", alpha=0.25)
        axis.legend(fontsize="small")
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

    guards._write_sha256(pdf_path)
    summary = {
        "schema": REPRESENTATION_PLOT_SCHEMA,
        "created_utc": _now(),
        "comparison": str(comparison_path),
        "comparison_sha256": hashlib.sha256(
            comparison_path.read_bytes()
        ).hexdigest(),
        "pdf": pdf_path.name,
        "representations": identifiers,
    }
    summary_path = output / "plot_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    guards._write_sha256(summary_path)
    return {
        "passed": True,
        "schema": REPRESENTATION_PLOT_SCHEMA,
        "output": str(output),
        "pdf": str(pdf_path),
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

    comparison_parser = subparsers.add_parser(
        "compare-representations",
        help=(
            "compare native-intreg groupings using frozen training and "
            "held-out replicas"
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
        "--parent-dilation", type=int, default=0
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
            "apply phase_space.y_max in both training and validation; the "
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
    representation_plot_parser = subparsers.add_parser(
        "plot-representations",
        help="render the milestone-2c representation comparison",
    )
    representation_plot_parser.add_argument(
        "--comparison", type=Path, required=True
    )
    representation_plot_parser.add_argument(
        "--output", type=Path, required=True
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
        elif args.command == "compare-representations":
            result = compare_representations(args)
        elif args.command == "plot-representations":
            result = plot_representations(args)
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
