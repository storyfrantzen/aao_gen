#!/usr/bin/env python3
"""Learn and validate radiative proposal guards from fixed-trial AAO surveys.

This is milestone 2 of the radiative bin-conditional workflow.  It does not
generate LUND events and it does not change AAO's proposal.  It converts
unrestricted milestone-1 surveys into a frozen, versioned description of a
high-contribution core and its nonzero global complement for every observed
analysis stratum.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import datetime as dt
import hashlib
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import radiative_survey


MANIFEST_SCHEMA = "aao-rad-guard-v1"
VALIDATION_SCHEMA = "aao-rad-guard-validation-v1"
COORDINATE_DEFINITION = "final_lund_analysis"
TAIL_DEFINITION = "global_unit_hypercube_minus_core_cells"

DEFAULT_PARTITION = (
    ("r_u", 8),
    ("r_ep", 8),
    ("u_gamma", 6),
    ("hadron_cosine_base", 6),
    ("hadron_phi_base", 8),
)
ALLOWED_PARTITION_COLUMNS = {
    "r_u",
    "r_ep",
    "u_gamma",
    "photon_cosine_base",
    "photon_phi_base",
    "hadron_cosine_base",
    "hadron_phi_base",
}
PERIODIC_PARTITION_COLUMNS = {"photon_phi_base", "hadron_phi_base"}
COMPATIBLE_NORM_KEYS = (
    "ebeam",
    "q2_min",
    "q2_max",
    "ep_min",
    "ep_max_effective",
    "delta",
    "epirea",
    "th_opt",
    "res_opt",
    "survey_phase_volume",
)


class GuardLearningError(RuntimeError):
    """Raised when inputs cannot define one reproducible guard campaign."""


@dataclass(frozen=True)
class PartitionAxis:
    name: str
    bins: int
    periodic: bool


@dataclass
class Moment:
    count: int = 0
    total: float = 0.0
    square_total: float = 0.0
    maximum: float = 0.0

    def add(self, value: float) -> None:
        self.count += 1
        self.total += value
        self.square_total += value * value
        self.maximum = max(self.maximum, value)


@dataclass(frozen=True)
class Stratum:
    identifier: str
    flat_index: int
    iq2: int
    ixb: int
    it: int
    iphi: int
    q2: tuple[float, float]
    xb: tuple[float, float]
    minus_t: tuple[float, float]
    phi_deg: tuple[float, float]


@dataclass
class Campaign:
    proposals: int
    rows: int
    final_valid_rows: int
    replicas: list[dict[str, object]]
    global_observed: Moment
    outside: dict[str, Moment]
    strata: dict[str, Moment]
    cells: dict[str, dict[tuple[int, int], Moment]]
    channels: dict[str, dict[int, Moment]]
    norm_reference: dict[str, str]
    input_signature: str
    legacy_input_settings: dict[str, object]
    channel_probabilities: dict[str, float]


def _strict_edges(values: Iterable[float], name: str) -> tuple[float, ...]:
    edges = tuple(float(value) for value in values)
    if len(edges) < 2 or any(not math.isfinite(value) for value in edges):
        raise GuardLearningError(f"{name} must contain at least two finite edges")
    if any(right <= left for left, right in zip(edges, edges[1:])):
        raise GuardLearningError(f"{name} edges must be strictly increasing")
    return edges


def load_analysis_config(path: Path) -> tuple[dict, str]:
    raw = path.read_bytes()
    try:
        config = json.loads(raw)
    except json.JSONDecodeError as error:
        raise GuardLearningError(f"{path}: invalid JSON: {error}") from error
    if "beam_energy" not in config:
        raise GuardLearningError("analysis config is missing beam_energy")
    phase_space = config.get("phase_space", {})
    for name in ("W_min",):
        if name not in phase_space:
            raise GuardLearningError(f"analysis config is missing phase_space.{name}")
    binning = config.get("binning", {})
    for name in ("Q2", "xB", "minus_t", "phi_deg"):
        if name not in binning:
            raise GuardLearningError(f"analysis config is missing binning.{name}")
        _strict_edges(binning[name], name)
    phi_edges = _strict_edges(binning["phi_deg"], "phi_deg")
    if not math.isclose(phi_edges[-1] - phi_edges[0], 360.0, abs_tol=1.0e-9):
        raise GuardLearningError("phi_deg binning must span one 360-degree period")
    return config, hashlib.sha256(raw).hexdigest()


def legacy_flat_index(
    iq2: int,
    ixb: int,
    it: int,
    iphi: int,
    *,
    nq2: int,
    nt: int,
    nphi: int,
) -> int:
    """Match the analysis response order: xB, Q2, phi, then -t fastest."""
    return ixb * nq2 * nphi * nt + iq2 * nphi * nt + iphi * nt + it


def enumerate_strata(config: dict) -> list[Stratum]:
    binning = config["binning"]
    q2_edges = _strict_edges(binning["Q2"], "Q2")
    xb_edges = _strict_edges(binning["xB"], "xB")
    t_edges = _strict_edges(binning["minus_t"], "minus_t")
    phi_edges = _strict_edges(binning["phi_deg"], "phi_deg")
    nq2 = len(q2_edges) - 1
    nt = len(t_edges) - 1
    nphi = len(phi_edges) - 1
    strata: list[Stratum] = []
    for iq2 in range(nq2):
        for ixb in range(len(xb_edges) - 1):
            for it in range(nt):
                for iphi in range(nphi):
                    flat = legacy_flat_index(
                        iq2, ixb, it, iphi, nq2=nq2, nt=nt, nphi=nphi
                    )
                    strata.append(
                        Stratum(
                            identifier=f"s{flat:05d}",
                            flat_index=flat,
                            iq2=iq2,
                            ixb=ixb,
                            it=it,
                            iphi=iphi,
                            q2=(q2_edges[iq2], q2_edges[iq2 + 1]),
                            xb=(xb_edges[ixb], xb_edges[ixb + 1]),
                            minus_t=(t_edges[it], t_edges[it + 1]),
                            phi_deg=(phi_edges[iphi], phi_edges[iphi + 1]),
                        )
                    )
    return sorted(strata, key=lambda item: item.flat_index)


def _bin_index(value: float, edges: tuple[float, ...]) -> int | None:
    if not math.isfinite(value) or value < edges[0] or value >= edges[-1]:
        return None
    index = bisect.bisect_right(edges, value) - 1
    return index if 0 <= index < len(edges) - 1 else None


def assign_stratum(
    row: dict[str, str], config: dict, *, apply_y_max: bool = False
) -> tuple[str | None, str | None]:
    """Return the canonical observed stratum or an explicit outside category."""
    if int(row["final_valid"]) != 1:
        return None, "invalid_final_candidate"
    values = {
        "Q2": float(row["q2_observed"]),
        "xB": float(row["xb_observed"]),
        "minus_t": float(row["minus_t_observed"]),
        "phi_deg": float(row["phi_observed_deg"]),
        "W": float(row["w_observed"]),
        "y": float(row["y_observed"]),
    }
    if any(not math.isfinite(value) for value in values.values()):
        return None, "nonfinite_observed_coordinates"

    phase_space = config["phase_space"]
    q2_min = float(phase_space.get("Q2_min", config["binning"]["Q2"][0]))
    if values["Q2"] < q2_min or values["W"] < float(phase_space["W_min"]):
        return None, "failed_analysis_phase_space"
    if apply_y_max:
        if "y_max" not in phase_space:
            raise GuardLearningError(
                "a y_max cut was requested, but phase_space.y_max is absent"
            )
        if values["y"] > float(phase_space["y_max"]):
            return None, "failed_analysis_phase_space"

    binning = config["binning"]
    q2_edges = _strict_edges(binning["Q2"], "Q2")
    xb_edges = _strict_edges(binning["xB"], "xB")
    t_edges = _strict_edges(binning["minus_t"], "minus_t")
    phi_edges = _strict_edges(binning["phi_deg"], "phi_deg")
    phi = (values["phi_deg"] - phi_edges[0]) % 360.0 + phi_edges[0]
    iq2 = _bin_index(values["Q2"], q2_edges)
    ixb = _bin_index(values["xB"], xb_edges)
    it = _bin_index(values["minus_t"], t_edges)
    iphi = _bin_index(phi, phi_edges)
    if None in (iq2, ixb, it, iphi):
        return None, "outside_analysis_binning"
    flat = legacy_flat_index(
        int(iq2),
        int(ixb),
        int(it),
        int(iphi),
        nq2=len(q2_edges) - 1,
        nt=len(t_edges) - 1,
        nphi=len(phi_edges) - 1,
    )
    return f"s{flat:05d}", None


def parse_partition(specifications: list[str] | None) -> tuple[PartitionAxis, ...]:
    raw = DEFAULT_PARTITION if not specifications else tuple(
        _parse_axis_specification(item) for item in specifications
    )
    axes = tuple(
        PartitionAxis(
            name=name,
            bins=bins,
            periodic=name in PERIODIC_PARTITION_COLUMNS,
        )
        for name, bins in raw
    )
    names = [axis.name for axis in axes]
    if len(names) != len(set(names)):
        raise GuardLearningError("partition dimensions must be unique")
    return axes


def _parse_axis_specification(value: str) -> tuple[str, int]:
    try:
        name, count = value.split("=", 1)
        bins = int(count)
    except (ValueError, TypeError) as error:
        raise GuardLearningError(
            f"invalid partition {value!r}; expected survey_column=number_of_bins"
        ) from error
    if name not in ALLOWED_PARTITION_COLUMNS:
        raise GuardLearningError(
            f"cannot partition {name!r}; allowed: "
            f"{', '.join(sorted(ALLOWED_PARTITION_COLUMNS))}"
        )
    if bins < 2 or bins > 256:
        raise GuardLearningError("each partition dimension must have 2 to 256 bins")
    return name, bins


def cell_indices(
    row: dict[str, str], partition: tuple[PartitionAxis, ...]
) -> tuple[int, ...]:
    indices: list[int] = []
    for axis in partition:
        value = float(row[axis.name])
        if not math.isfinite(value) or value < 0.0 or value > 1.0:
            raise GuardLearningError(
                f"{axis.name}={value!r} is outside its survey base domain [0,1]"
            )
        if axis.periodic:
            value %= 1.0
        index = min(int(value * axis.bins), axis.bins - 1)
        indices.append(index)
    return tuple(indices)


def flatten_cell(
    indices: tuple[int, ...], partition: tuple[PartitionAxis, ...]
) -> int:
    cell_id = 0
    for index, axis in zip(indices, partition):
        if not 0 <= index < axis.bins:
            raise GuardLearningError(f"cell index {index} is outside {axis.name}")
        cell_id = cell_id * axis.bins + index
    return cell_id


def unflatten_cell(
    cell_id: int, partition: tuple[PartitionAxis, ...]
) -> tuple[int, ...]:
    if cell_id < 0:
        raise GuardLearningError("cell ID must be nonnegative")
    indices = [0] * len(partition)
    remainder = cell_id
    for position in range(len(partition) - 1, -1, -1):
        indices[position] = remainder % partition[position].bins
        remainder //= partition[position].bins
    if remainder:
        raise GuardLearningError(f"cell ID {cell_id} exceeds the partition")
    return tuple(indices)


def dilate_cells(
    seeds: set[tuple[int, int]],
    partition: tuple[PartitionAxis, ...],
    radius: int,
) -> set[tuple[int, int]]:
    """Dilate by repeated axis-neighbor steps, wrapping periodic axes."""
    if radius < 0:
        raise GuardLearningError("dilation radius must be nonnegative")
    result = set(seeds)
    frontier = set(seeds)
    for _ in range(radius):
        additions: set[tuple[int, int]] = set()
        for channel, cell_id in frontier:
            indices = list(unflatten_cell(cell_id, partition))
            for dimension, axis in enumerate(partition):
                for shift in (-1, 1):
                    neighbor = indices.copy()
                    candidate = neighbor[dimension] + shift
                    if axis.periodic:
                        candidate %= axis.bins
                    elif not 0 <= candidate < axis.bins:
                        continue
                    neighbor[dimension] = candidate
                    additions.add(
                        (channel, flatten_cell(tuple(neighbor), partition))
                    )
        additions -= result
        result.update(additions)
        frontier = additions
        if not frontier:
            break
    return result


def _survey_rows(path: Path) -> Iterator[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as source:
        first = source.readline().strip()
        if first != f"# schema={radiative_survey.SURVEY_SCHEMA}":
            raise GuardLearningError(
                f"{path}: expected '# schema={radiative_survey.SURVEY_SCHEMA}'"
            )
        reader = csv.DictReader(source)
        if reader.fieldnames != radiative_survey.SURVEY_COLUMNS:
            raise GuardLearningError(f"{path}: unexpected survey columns")
        yield from reader


def _validate_guard_row(
    row: dict[str, str],
    *,
    directory: Path,
    phase_volume: float,
    beam_energy: float,
) -> None:
    """Recheck the survey invariants on which guard classification depends."""
    try:
        valid = int(row["final_valid"])
        status = int(row["candidate_status"])
        internal_integrand = float(row["integrand_internal"])
        observed_integrand = float(row["integrand_observed"])
        internal_contribution = float(row["trial_xsec_internal_microbarn"])
        observed_contribution = float(row["trial_xsec_observed_microbarn"])
    except (ValueError, TypeError) as error:
        raise GuardLearningError(f"{directory}: malformed survey row") from error
    if valid not in (0, 1) or status not in (0, 1, 2, 3):
        raise GuardLearningError(f"{directory}: invalid candidate state")
    if valid == 1 and status != 0:
        raise GuardLearningError(f"{directory}: valid candidate has nonzero status")
    if valid == 0 and status == 0:
        raise GuardLearningError(f"{directory}: invalid candidate has zero status")
    for name in (
        "r_u",
        "r_ep",
        "u_gamma",
        "channel_selector",
        "vertex_fraction",
        "incoming_loss_base",
        "incoming_loss_accept_test",
        "photon_side_base",
        "photon_cosine_base",
        "photon_phi_base",
        "hadron_cosine_base",
        "hadron_phi_base",
        "outgoing_loss_base",
        "outgoing_loss_accept_test",
    ):
        value = float(row[name])
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            raise GuardLearningError(f"{directory}: {name} is outside [0,1]")
    if internal_integrand <= 0.0 or not math.isclose(
        internal_contribution,
        phase_volume * internal_integrand,
        rel_tol=3.0e-6,
        abs_tol=1.0e-18,
    ):
        raise GuardLearningError(f"{directory}: invalid internal contribution")
    if valid == 0:
        if observed_integrand != 0.0 or observed_contribution != 0.0:
            raise GuardLearningError(
                f"{directory}: invalid candidate has observed contribution"
            )
        return
    if not math.isclose(
        observed_integrand,
        internal_integrand,
        rel_tol=3.0e-6,
        abs_tol=1.0e-18,
    ) or not math.isclose(
        observed_contribution,
        phase_volume * observed_integrand,
        rel_tol=3.0e-6,
        abs_tol=1.0e-18,
    ):
        raise GuardLearningError(f"{directory}: invalid observed contribution")
    try:
        calculated = radiative_survey.observed_coordinates(
            beam_energy,
            tuple(
                float(row[name])
                for name in (
                    "final_e_px",
                    "final_e_py",
                    "final_e_pz",
                    "final_e_energy",
                )
            ),
            tuple(
                float(row[name])
                for name in (
                    "final_p_px",
                    "final_p_py",
                    "final_p_pz",
                    "final_p_energy",
                )
            ),
        )
    except (ValueError, radiative_survey.SurveyValidationError) as error:
        raise GuardLearningError(
            f"{directory}: cannot reconstruct observed coordinates"
        ) from error
    for name, tolerance in radiative_survey.COORDINATE_TOLERANCES.items():
        recorded = float(row[name])
        if name == "phi_observed_deg":
            difference = abs((calculated[name] - recorded + 180.0) % 360.0 - 180.0)
        else:
            difference = abs(calculated[name] - recorded)
        if difference > tolerance:
            raise GuardLearningError(
                f"{directory}: observed-coordinate mismatch in {name}"
            )


def _legacy_input_metadata(
    path: Path,
) -> tuple[str, dict[str, float], dict[str, object]]:
    if not path.is_file():
        raise GuardLearningError(f"missing survey input snapshot {path}")
    records = [
        line.split("!", 1)[0].strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.split("!", 1)[0].strip()
    ]
    if len(records) < 21 or records[-4] != "1":
        raise GuardLearningError(f"{path}: not a milestone-1 survey input")
    legacy = records[:-4]
    try:
        probabilities = [float(value) for value in legacy[2].split()]
    except (ValueError, IndexError) as error:
        raise GuardLearningError(f"{path}: cannot parse channel probabilities") from error
    if len(probabilities) != 4 or any(value < 0.0 for value in probabilities):
        raise GuardLearningError(f"{path}: expected four nonnegative channel sizes")
    fifth = 1.0 - sum(probabilities)
    if fifth < -1.0e-7:
        raise GuardLearningError(f"{path}: radiative channel sizes sum above one")
    channel_probabilities = {
        str(index + 1): value for index, value in enumerate(probabilities)
    }
    channel_probabilities["5"] = max(0.0, fifth)
    try:
        settings = {
            "physics_model": int(legacy[0].split()[0]),
            "electron_helicity_flag": int(legacy[1].split()[0]),
            "npart": int(legacy[3].split()[0]),
            "epirea": int(legacy[4].split()[0]),
            "missing_mass_squared_tolerance": float(legacy[5].split()[0]),
            "target": {
                "length_cm": float(legacy[6].split()[0]),
                "radius_cm": float(legacy[7].split()[0]),
                "beam_position_cm": [
                    float(legacy[index].split()[0]) for index in (8, 9, 10)
                ],
            },
            "beam_energy_GeV": float(legacy[11].split()[0]),
            "input_q2_range_GeV2": [
                float(value) for value in legacy[12].split()
            ],
            "input_scattered_electron_energy_range_GeV": [
                float(value) for value in legacy[13].split()
            ],
            "minimum_internal_photon_energy_GeV": float(
                legacy[14].split()[0]
            ),
        }
    except (ValueError, IndexError) as error:
        raise GuardLearningError(f"{path}: cannot parse legacy settings") from error
    signature = hashlib.sha256(("\n".join(legacy) + "\n").encode()).hexdigest()
    return signature, channel_probabilities, settings


def _validated_norm(directory: Path) -> dict[str, str]:
    path = directory / radiative_survey.NORM_FILENAME
    if not path.is_file():
        raise GuardLearningError(f"missing {path}")
    norm = radiative_survey.parse_norm(path)
    expected = {
        "generator": "aao_rad",
        "sampling_mode": "1",
        "fixed_trial_survey": "1",
        "survey_schema": radiative_survey.SURVEY_SCHEMA,
        "survey_emits_lund": "0",
    }
    for key, value in expected.items():
        if norm.get(key) != value:
            raise GuardLearningError(f"{path}: {key}={norm.get(key)!r}, expected {value!r}")
    return norm


def _compatible_settings(reference: dict[str, str], candidate: dict[str, str]) -> None:
    for key in COMPATIBLE_NORM_KEYS:
        try:
            left = float(reference[key])
            right = float(candidate[key])
        except (KeyError, ValueError) as error:
            raise GuardLearningError(f"normalization metadata lacks {key}") from error
        if not math.isclose(left, right, rel_tol=2.0e-7, abs_tol=1.0e-12):
            raise GuardLearningError(
                f"survey setting {key} differs across replicas: {left} versus {right}"
            )


def aggregate_surveys(
    directories: list[Path],
    config: dict,
    partition: tuple[PartitionAxis, ...],
    *,
    apply_y_max: bool = False,
) -> Campaign:
    if not directories:
        raise GuardLearningError("at least one survey directory is required")
    proposals = rows = final_valid_rows = 0
    replicas: list[dict[str, object]] = []
    replica_ids: set[int] = set()
    global_observed = Moment()
    outside: dict[str, Moment] = {}
    strata: dict[str, Moment] = {}
    cells: dict[str, dict[tuple[int, int], Moment]] = {}
    channels: dict[str, dict[int, Moment]] = {}
    norm_reference: dict[str, str] | None = None
    input_signature: str | None = None
    legacy_input_settings: dict[str, object] | None = None
    channel_probabilities: dict[str, float] | None = None

    for requested_directory in directories:
        directory = requested_directory.resolve()
        norm = _validated_norm(directory)
        if norm_reference is None:
            norm_reference = norm
        else:
            _compatible_settings(norm_reference, norm)
        signature, probabilities, settings = _legacy_input_metadata(
            directory / "survey_input.inp"
        )
        if input_signature is None:
            input_signature = signature
            channel_probabilities = probabilities
            legacy_input_settings = settings
        elif signature != input_signature:
            raise GuardLearningError(
                "survey legacy inputs differ; training/validation replicas must "
                "have identical physical and proposal settings"
            )

        try:
            ntrials = int(norm["ntries"])
            requested = int(norm["survey_ntrials_requested"])
            expected_rows = int(norm["survey_rows"])
            expected_valid = int(norm["survey_final_valid"])
            replica = int(norm["survey_replica"])
        except (KeyError, ValueError) as error:
            raise GuardLearningError(
                f"{directory}: incomplete survey normalization metadata"
            ) from error
        if ntrials != requested or ntrials <= 0:
            raise GuardLearningError(f"{directory}: incomplete fixed-trial survey")
        if replica in replica_ids:
            raise GuardLearningError(f"duplicate survey replica ID {replica}")
        replica_ids.add(replica)

        local = Moment()
        local_rows = local_valid = 0
        previous_trial = 0
        phase_volume = float(norm["survey_phase_volume"])
        beam_energy = float(norm["ebeam"])
        for row in _survey_rows(directory / radiative_survey.SURVEY_FILENAME):
            local_rows += 1
            _validate_guard_row(
                row,
                directory=directory,
                phase_volume=phase_volume,
                beam_energy=beam_energy,
            )
            try:
                trial = int(row["trial"])
                row_replica = int(row["replica"])
                valid = int(row["final_valid"])
                channel = int(row["intreg"])
                weight = float(row["trial_xsec_observed_microbarn"])
            except (ValueError, TypeError) as error:
                raise GuardLearningError(
                    f"{directory}: invalid numeric survey row {local_rows}"
                ) from error
            if trial <= previous_trial or trial > ntrials:
                raise GuardLearningError(
                    f"{directory}: trial IDs must increase within [1,ntries]"
                )
            previous_trial = trial
            if row_replica != replica:
                raise GuardLearningError(f"{directory}: CSV replica ID mismatch")
            if valid not in (0, 1) or channel not in range(1, 7):
                raise GuardLearningError(f"{directory}: invalid status or channel")
            if not math.isfinite(weight) or weight < 0.0:
                raise GuardLearningError(f"{directory}: invalid observed contribution")
            if valid == 0 and weight != 0.0:
                raise GuardLearningError(
                    f"{directory}: invalid final candidate has nonzero contribution"
                )

            rows += 1
            if valid == 1:
                final_valid_rows += 1
                local_valid += 1
                local.add(weight)
                global_observed.add(weight)
            stratum_id, outside_reason = assign_stratum(
                row, config, apply_y_max=apply_y_max
            )
            if stratum_id is None:
                outside.setdefault(outside_reason or "unknown", Moment()).add(weight)
                continue
            strata.setdefault(stratum_id, Moment()).add(weight)
            channel_moment = channels.setdefault(stratum_id, {}).setdefault(
                channel, Moment()
            )
            channel_moment.add(weight)
            key = (
                channel,
                flatten_cell(cell_indices(row, partition), partition),
            )
            cells.setdefault(stratum_id, {}).setdefault(key, Moment()).add(weight)

        if local_rows != expected_rows or local_valid != expected_valid:
            raise GuardLearningError(
                f"{directory}: CSV row counts disagree with normalization metadata"
            )
        observed_estimate = local.total / ntrials
        expected_estimate = float(norm["survey_observed_sig_sum"])
        if not math.isclose(
            observed_estimate, expected_estimate, rel_tol=5.0e-6, abs_tol=1.0e-15
        ):
            raise GuardLearningError(
                f"{directory}: CSV observed integral {observed_estimate:.12g} "
                f"does not match norm {expected_estimate:.12g}"
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
                "observed_cross_section_sem_microbarn": _sem(local, ntrials),
            }
        )

    return Campaign(
        proposals=proposals,
        rows=rows,
        final_valid_rows=final_valid_rows,
        replicas=sorted(replicas, key=lambda item: int(item["replica"])),
        global_observed=global_observed,
        outside=outside,
        strata=strata,
        cells=cells,
        channels=channels,
        norm_reference=norm_reference or {},
        input_signature=input_signature or "",
        legacy_input_settings=legacy_input_settings or {},
        channel_probabilities=channel_probabilities or {},
    )


def _sem(moment: Moment, proposals: int) -> float | None:
    if proposals <= 1:
        return None
    variance = max(
        0.0, moment.square_total - moment.total * moment.total / proposals
    ) / (proposals - 1)
    return math.sqrt(variance / proposals)


def _ess(moment: Moment) -> float:
    if moment.square_total <= 0.0:
        return 0.0
    return moment.total * moment.total / moment.square_total


def _metrics(moment: Moment, proposals: int) -> dict[str, float | int | None]:
    cross_section = moment.total / proposals
    sem = _sem(moment, proposals)
    return {
        "contributing_rows": moment.count,
        "contributing_row_fraction_of_all_proposals": moment.count / proposals,
        "sum_trial_contributions_microbarn": moment.total,
        "sum_squared_trial_contributions_microbarn2": moment.square_total,
        "largest_trial_contribution_microbarn": moment.maximum,
        "cross_section_microbarn": cross_section,
        "cross_section_sem_microbarn": sem,
        "cross_section_relative_standard_error": (
            sem / cross_section if sem is not None and cross_section > 0.0 else None
        ),
        "importance_effective_sample_size": _ess(moment),
    }


def _difference_z_score(
    first_value: float,
    first_sem: float | None,
    second_value: float,
    second_sem: float | None,
) -> float | None:
    if first_sem is None or second_sem is None:
        return None
    uncertainty = math.hypot(first_sem, second_sem)
    if uncertainty <= 0.0:
        return None
    return (second_value - first_value) / uncertainty


def _subset_moment(
    cell_moments: dict[tuple[int, int], Moment],
    selected: set[tuple[int, int]],
) -> Moment:
    result = Moment()
    for key, moment in cell_moments.items():
        if key in selected:
            result.count += moment.count
            result.total += moment.total
            result.square_total += moment.square_total
            result.maximum = max(result.maximum, moment.maximum)
    return result


def _combined_moment(moments: Iterable[Moment]) -> Moment:
    result = Moment()
    for moment in moments:
        result.count += moment.count
        result.total += moment.total
        result.square_total += moment.square_total
        result.maximum = max(result.maximum, moment.maximum)
    return result


def _select_seed_cells(
    cell_moments: dict[tuple[int, int], Moment], target_fraction: float
) -> set[tuple[int, int]]:
    total = sum(moment.total for moment in cell_moments.values())
    if total <= 0.0:
        return set()
    ordered = sorted(
        cell_moments.items(),
        key=lambda item: (-item[1].total, item[0][0], item[0][1]),
    )
    selected: set[tuple[int, int]] = set()
    retained = 0.0
    for key, moment in ordered:
        if moment.total <= 0.0:
            continue
        selected.add(key)
        retained += moment.total
        if retained / total >= target_fraction:
            break
    return selected


def _group_cell_ids(cells: set[tuple[int, int]]) -> dict[str, list[int]]:
    grouped: dict[str, list[int]] = {}
    for channel, cell_id in sorted(cells):
        grouped.setdefault(str(channel), []).append(cell_id)
    return grouped


def _core_sampling_model(
    core: set[tuple[int, int]],
    cell_moments: dict[tuple[int, int], Moment],
    uniform_floor_fraction: float,
) -> tuple[dict, dict[str, float]]:
    weighted = {
        key: moment.total
        for key, moment in cell_moments.items()
        if key in core and moment.total > 0.0
    }
    weighted_total = sum(weighted.values())
    cell_count = len(core)
    probabilities: dict[tuple[int, int], float] = {}
    for key in core:
        probability = uniform_floor_fraction / cell_count
        if weighted_total > 0.0:
            probability += (
                (1.0 - uniform_floor_fraction) * weighted.get(key, 0.0) / weighted_total
            )
        probabilities[key] = probability
    channel_probabilities: dict[str, float] = {}
    for (channel, _), probability in probabilities.items():
        key = str(channel)
        channel_probabilities[key] = channel_probabilities.get(key, 0.0) + probability
    weighted_by_channel: dict[str, list[dict[str, float | int]]] = {}
    for (channel, cell_id), contribution in sorted(weighted.items()):
        weighted_by_channel.setdefault(str(channel), []).append(
            {
                "cell_id": cell_id,
                "training_sum_trial_contributions_microbarn": contribution,
            }
        )
    model = {
        "description": (
            "(1-uniform_floor_fraction) times normalized training contribution "
            "plus uniform_floor_fraction uniformly over all derived core cells"
        ),
        "uniform_floor_fraction": uniform_floor_fraction,
        "weighted_cells_by_channel": weighted_by_channel,
    }
    return model, channel_probabilities


def _stratum_metadata(stratum: Stratum) -> dict[str, object]:
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


def _manifest_global_domain(
    norm: dict[str, str], channel_probabilities: dict[str, float]
) -> dict[str, object]:
    return {
        "base_coordinate_domain": {
            name: [0.0, 1.0]
            for name in (
                "r_u",
                "r_ep",
                "u_gamma",
                "channel_selector",
                "vertex_fraction",
                "incoming_loss_base",
                "incoming_loss_accept_test",
                "photon_side_base",
                "photon_cosine_base",
                "photon_phi_base",
                "hadron_cosine_base",
                "hadron_phi_base",
                "outgoing_loss_base",
                "outgoing_loss_accept_test",
            )
        },
        "transforms": {
            "r_u": {
                "quantity": "inverse_q2_leptonic",
                "minimum": 1.0 / float(norm["q2_max"]),
                "maximum": 1.0 / float(norm["q2_min"]),
            },
            "r_ep": {
                "quantity": "energy_e_pre_external_GeV",
                "at_zero": float(norm["ep_max_effective"]),
                "at_one": float(norm["ep_min"]),
            },
            "u_gamma": {
                "quantity": "energy_gamma_GeV",
                "formula": "-ln(u_gamma)/5",
            },
            "hadron_cosine_base": {
                "quantity": "cos_theta_cm",
                "formula": "-1+2*hadron_cosine_base",
            },
            "hadron_phi_base": {
                "quantity": "phi_cm_deg",
                "formula": "360*hadron_phi_base",
            },
        },
        "radiative_angular_channel_base_probabilities": channel_probabilities,
        "soft_channel_note": (
            "Survey intreg=6 is the soft-photon branch selected after u_gamma; "
            "the listed probabilities describe the five angular-channel draws."
        ),
    }


def build_manifest(
    *,
    config: dict,
    config_path: Path,
    config_sha256: str,
    campaign: Campaign,
    partition: tuple[PartitionAxis, ...],
    target_core_fraction: float,
    core_mixture_probability: float,
    uniform_floor_fraction: float,
    dilation_radius: int,
    guard_iteration: int,
    generator_revision: str,
    generator_revision_source: str,
    learner_revision: str,
    minimum_training_rows: int,
    minimum_training_ess: float,
    apply_y_max: bool,
) -> dict:
    catalog = enumerate_strata(config)
    strata_output: dict[str, dict] = {}
    learned = low_support = empty = 0
    for stratum in catalog:
        total = campaign.strata.get(stratum.identifier, Moment())
        cell_moments = campaign.cells.get(stratum.identifier, {})
        record = _stratum_metadata(stratum)
        record["training_total"] = _metrics(total, campaign.proposals)
        if total.total <= 0.0:
            record.update(
                {
                    "status": "no_training_contribution",
                    "core_cells": {
                        "representation": "axis_dilation_of_seed_cells",
                        "seed_cell_ids_by_channel": {},
                        "dilation_axis_steps": dilation_radius,
                        "derived_cell_count": 0,
                    },
                    "estimated_core_fraction": 0.0,
                    "estimated_tail_fraction": 0.0,
                }
            )
            empty += 1
            strata_output[stratum.identifier] = record
            continue

        seeds = _select_seed_cells(cell_moments, target_core_fraction)
        core = dilate_cells(seeds, partition, dilation_radius)
        expanded = dilate_cells(seeds, partition, dilation_radius + 1)
        core_moment = _subset_moment(cell_moments, core)
        tail_moment = _subset_moment(
            cell_moments, set(cell_moments).difference(core)
        )
        expanded_moment = _subset_moment(cell_moments, expanded)
        sampling_model, channel_probabilities = _core_sampling_model(
            core, cell_moments, uniform_floor_fraction
        )
        sufficient = total.count >= minimum_training_rows and _ess(
            total
        ) >= minimum_training_ess
        status = "learned" if sufficient else "learned_low_support"
        learned += int(sufficient)
        low_support += int(not sufficient)
        record.update(
            {
                "status": status,
                "core_cells": {
                    "representation": "axis_dilation_of_seed_cells",
                    "neighborhood": "repeated_axis_neighbors",
                    "seed_cell_ids_by_channel": _group_cell_ids(seeds),
                    "dilation_axis_steps": dilation_radius,
                    "derived_cell_count": len(core),
                },
                "core_sampling_model": sampling_model,
                "core_channel_probabilities": channel_probabilities,
                "training_channel_cross_section_fractions": {
                    str(channel): moment.total / total.total
                    for channel, moment in sorted(
                        campaign.channels.get(stratum.identifier, {}).items()
                    )
                },
                "training_core": _metrics(core_moment, campaign.proposals),
                "training_tail": _metrics(tail_moment, campaign.proposals),
                "training_one_more_dilation": _metrics(
                    expanded_moment, campaign.proposals
                ),
                "estimated_core_fraction": core_moment.total / total.total,
                "estimated_tail_fraction": tail_moment.total / total.total,
                "one_more_dilation_fraction": expanded_moment.total / total.total,
            }
        )
        strata_output[stratum.identifier] = record

    created = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    inside = _combined_moment(campaign.strata.values())
    outside_total = _combined_moment(campaign.outside.values())
    return {
        "schema": MANIFEST_SCHEMA,
        "created_utc": created,
        "production_ready": False,
        "production_readiness_note": (
            "Milestone 2 learning artifact; mode-4 proposal correction and "
            "unweighting are implemented in milestone 3."
        ),
        "coordinate_definition": COORDINATE_DEFINITION,
        "artifacts": {
            "training_coverage": "training_coverage.csv",
            "training_cell_statistics": "training_cells.csv",
            "hash_sidecar_suffix": ".sha256",
        },
        "guard_iteration": guard_iteration,
        "generator_revision": generator_revision,
        "generator_revision_source": generator_revision_source,
        "guard_learner_revision": learner_revision,
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
        "global_proposal_domain": _manifest_global_domain(
            campaign.norm_reference, campaign.channel_probabilities
        ),
        "legacy_generator_settings": campaign.legacy_input_settings,
        "partition_definition": {
            "coordinate_domain": [0.0, 1.0],
            "cell_id_order": "row_major_in_listed_axis_order",
            "cell_base_coordinate_volume": math.prod(
                1.0 / axis.bins for axis in partition
            ),
            "axes": [
                {
                    "name": axis.name,
                    "bins": axis.bins,
                    "periodic": axis.periodic,
                }
                for axis in partition
            ],
        },
        "tail_definition": TAIL_DEFINITION,
        "proposal_mixture": {
            "core_probability": core_mixture_probability,
            "tail_probability": 1.0 - core_mixture_probability,
            "tail_has_nonzero_probability": core_mixture_probability < 1.0,
        },
        "learning": {
            "target_core_fraction_before_dilation": target_core_fraction,
            "dilation_axis_steps": dilation_radius,
            "uniform_core_cell_probability_floor_fraction": uniform_floor_fraction,
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
            "global_observed": _metrics(
                campaign.global_observed, campaign.proposals
            ),
            "inside_analysis_partition": _metrics(inside, campaign.proposals),
            "outside_analysis_partition": {
                name: _metrics(moment, campaign.proposals)
                for name, moment in sorted(campaign.outside.items())
            },
            "cross_section_partition_closure_microbarn": (
                campaign.global_observed.total
                - inside.total
                - outside_total.total
            )
            / campaign.proposals,
        },
        "strata": strata_output,
    }


def _partition_from_manifest(manifest: dict) -> tuple[PartitionAxis, ...]:
    try:
        return tuple(
            PartitionAxis(
                name=str(axis["name"]),
                bins=int(axis["bins"]),
                periodic=bool(axis["periodic"]),
            )
            for axis in manifest["partition_definition"]["axes"]
        )
    except (KeyError, TypeError, ValueError) as error:
        raise GuardLearningError("manifest has an invalid partition definition") from error


def _core_from_record(
    record: dict, partition: tuple[PartitionAxis, ...], extra_steps: int = 0
) -> set[tuple[int, int]]:
    definition = record.get("core_cells", {})
    seeds = {
        (int(channel), int(cell_id))
        for channel, cell_ids in definition.get(
            "seed_cell_ids_by_channel", {}
        ).items()
        for cell_id in cell_ids
    }
    radius = int(definition.get("dilation_axis_steps", 0)) + extra_steps
    return dilate_cells(seeds, partition, radius)


def validate_manifest(
    manifest: dict,
    manifest_path: Path,
    campaign: Campaign,
    partition: tuple[PartitionAxis, ...],
    minimum_core_fraction: float,
) -> dict:
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise GuardLearningError(
            f"manifest schema is {manifest.get('schema')!r}, expected {MANIFEST_SCHEMA}"
        )
    expected_signature = manifest["training"]["legacy_input_sha256"]
    if campaign.input_signature != expected_signature:
        raise GuardLearningError(
            "validation surveys do not match the manifest's legacy input settings"
        )
    training_ids = set(int(value) for value in manifest["training"]["replica_ids"])
    validation_ids = {int(item["replica"]) for item in campaign.replicas}
    overlap = training_ids & validation_ids
    if overlap:
        raise GuardLearningError(
            f"validation replicas overlap training replicas: {sorted(overlap)}"
        )

    strata_output: dict[str, dict] = {}
    assessed = passed = failed = no_holdout = 0
    for stratum_id, record in manifest["strata"].items():
        training_total = float(
            record["training_total"]["sum_trial_contributions_microbarn"]
        )
        holdout_total = campaign.strata.get(stratum_id, Moment())
        if training_total <= 0.0 and holdout_total.total <= 0.0:
            continue
        cells = campaign.cells.get(stratum_id, {})
        core = _core_from_record(record, partition)
        expanded = _core_from_record(record, partition, extra_steps=1)
        core_moment = _subset_moment(cells, core)
        tail_moment = _subset_moment(cells, set(cells).difference(core))
        expanded_moment = _subset_moment(cells, expanded)
        tail_cells = sorted(
            (
                (key, moment)
                for key, moment in cells.items()
                if key not in core and moment.total > 0.0
            ),
            key=lambda item: -item[1].total,
        )
        if holdout_total.total > 0.0:
            core_fraction = core_moment.total / holdout_total.total
            expanded_fraction = expanded_moment.total / holdout_total.total
            coverage_passed = core_fraction >= minimum_core_fraction
            assessed += 1
            passed += int(coverage_passed)
            failed += int(not coverage_passed)
        else:
            core_fraction = expanded_fraction = None
            coverage_passed = None
            no_holdout += 1
        channel_fractions = {
            str(channel): moment.total / holdout_total.total
            for channel, moment in sorted(
                campaign.channels.get(stratum_id, {}).items()
            )
            if holdout_total.total > 0.0
        }
        largest_tail_cells = []
        for (channel, cell_id), moment in tail_cells[:5]:
            largest_tail_cells.append(
                {
                    "channel": channel,
                    "cell_id": cell_id,
                    **_metrics(moment, campaign.proposals),
                    "fraction_of_holdout_cross_section": (
                        moment.total / holdout_total.total
                        if holdout_total.total > 0.0
                        else None
                    ),
                }
            )
        training_metrics = record["training_total"]
        holdout_metrics = _metrics(holdout_total, campaign.proposals)
        strata_output[stratum_id] = {
            "training_status": record["status"],
            "training_core_fraction": record.get("estimated_core_fraction", 0.0),
            "holdout_total": holdout_metrics,
            "holdout_core": _metrics(core_moment, campaign.proposals),
            "holdout_tail": _metrics(tail_moment, campaign.proposals),
            "holdout_one_more_dilation": _metrics(
                expanded_moment, campaign.proposals
            ),
            "holdout_core_fraction": core_fraction,
            "holdout_tail_fraction": (
                1.0 - core_fraction if core_fraction is not None else None
            ),
            "holdout_one_more_dilation_fraction": expanded_fraction,
            "channel_cross_section_fractions": channel_fractions,
            "largest_tail_cells": largest_tail_cells,
            "training_holdout_difference_z_score": _difference_z_score(
                float(training_metrics["cross_section_microbarn"]),
                training_metrics["cross_section_sem_microbarn"],
                float(holdout_metrics["cross_section_microbarn"]),
                holdout_metrics["cross_section_sem_microbarn"],
            ),
            "coverage_passed": coverage_passed,
        }

    inside = _combined_moment(campaign.strata.values())
    outside_total = _combined_moment(campaign.outside.values())
    training_global = manifest["training"]["global_observed"]
    validation_global = _metrics(campaign.global_observed, campaign.proposals)
    return {
        "schema": VALIDATION_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc)
        .replace(microsecond=0)
        .isoformat(),
        "manifest": str(manifest_path.resolve()),
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "coordinate_definition": COORDINATE_DEFINITION,
        "artifacts": {
            "validation_coverage": "validation_coverage.csv",
            "validation_cell_statistics": "validation_cells.csv",
            "hash_sidecar_suffix": ".sha256",
        },
        "minimum_core_fraction": minimum_core_fraction,
        "passed": failed == 0 and assessed > 0,
        "global_training_holdout_difference_z_score": _difference_z_score(
            float(training_global["cross_section_microbarn"]),
            training_global["cross_section_sem_microbarn"],
            float(validation_global["cross_section_microbarn"]),
            validation_global["cross_section_sem_microbarn"],
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
            "inside_analysis_partition": _metrics(inside, campaign.proposals),
            "outside_analysis_partition": {
                name: _metrics(moment, campaign.proposals)
                for name, moment in sorted(campaign.outside.items())
            },
            "cross_section_partition_closure_microbarn": (
                campaign.global_observed.total
                - inside.total
                - outside_total.total
            )
            / campaign.proposals,
        },
        "strata": strata_output,
    }


def _write_json_artifact(output: Path, filename: str, payload: dict) -> Path:
    if output.exists():
        raise FileExistsError(
            f"{output} already exists; guard artifacts are immutable"
        )
    output.mkdir(parents=True)
    path = output / filename
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    path.write_text(serialized, encoding="utf-8")
    digest = hashlib.sha256(serialized.encode()).hexdigest()
    (output / f"{filename}.sha256").write_text(
        f"{digest}  {filename}\n", encoding="utf-8"
    )
    return path


def _write_learning_csv(path: Path, manifest: dict) -> None:
    fields = (
        "stratum_id",
        "flat_index",
        "status",
        "training_cross_section_microbarn",
        "training_sem_microbarn",
        "training_ess",
        "training_rows",
        "core_fraction",
        "tail_fraction",
        "one_more_dilation_fraction",
        "seed_cells",
        "derived_core_cells",
    )
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        for stratum_id, record in manifest["strata"].items():
            total = record["training_total"]
            seeds = sum(
                len(ids)
                for ids in record["core_cells"][
                    "seed_cell_ids_by_channel"
                ].values()
            )
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
                    "training_ess": total["importance_effective_sample_size"],
                    "training_rows": total["contributing_rows"],
                    "core_fraction": record["estimated_core_fraction"],
                    "tail_fraction": record["estimated_tail_fraction"],
                    "one_more_dilation_fraction": record.get(
                        "one_more_dilation_fraction", ""
                    ),
                    "seed_cells": seeds,
                    "derived_core_cells": record["core_cells"][
                        "derived_cell_count"
                    ],
                }
            )


def _write_validation_csv(path: Path, validation: dict) -> None:
    fields = (
        "stratum_id",
        "training_status",
        "training_core_fraction",
        "holdout_cross_section_microbarn",
        "holdout_sem_microbarn",
        "holdout_ess",
        "holdout_core_fraction",
        "holdout_tail_fraction",
        "holdout_one_more_dilation_fraction",
        "coverage_passed",
        "largest_tail_fraction",
    )
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        for stratum_id, record in validation["strata"].items():
            total = record["holdout_total"]
            largest_cells = record["largest_tail_cells"]
            largest = largest_cells[0] if largest_cells else None
            writer.writerow(
                {
                    "stratum_id": stratum_id,
                    "training_status": record["training_status"],
                    "training_core_fraction": record["training_core_fraction"],
                    "holdout_cross_section_microbarn": total[
                        "cross_section_microbarn"
                    ],
                    "holdout_sem_microbarn": total[
                        "cross_section_sem_microbarn"
                    ],
                    "holdout_ess": total["importance_effective_sample_size"],
                    "holdout_core_fraction": record["holdout_core_fraction"],
                    "holdout_tail_fraction": record["holdout_tail_fraction"],
                    "holdout_one_more_dilation_fraction": record[
                        "holdout_one_more_dilation_fraction"
                    ],
                    "coverage_passed": record["coverage_passed"],
                    "largest_tail_fraction": (
                        largest["fraction_of_holdout_cross_section"]
                        if largest
                        else ""
                    ),
                }
            )


def _write_cell_statistics(
    path: Path,
    *,
    manifest: dict,
    campaign: Campaign,
    partition: tuple[PartitionAxis, ...],
    validation: bool,
) -> None:
    prefix = "holdout" if validation else "training"
    fields = [
        "stratum_id",
        "channel",
        "cell_id",
        *(f"{axis.name}_index" for axis in partition),
        "contributing_rows",
        f"{prefix}_sum_trial_contributions_microbarn",
        f"{prefix}_sum_squared_trial_contributions_microbarn2",
        f"{prefix}_largest_trial_contribution_microbarn",
        f"{prefix}_cross_section_microbarn",
        "selected_seed",
        "in_core",
        "in_one_more_dilation",
    ]
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        for stratum_id, cell_moments in sorted(campaign.cells.items()):
            record = manifest["strata"][stratum_id]
            definition = record["core_cells"]
            seeds = {
                (int(channel), int(cell_id))
                for channel, cell_ids in definition[
                    "seed_cell_ids_by_channel"
                ].items()
                for cell_id in cell_ids
            }
            core = _core_from_record(record, partition)
            expanded = _core_from_record(record, partition, extra_steps=1)
            for (channel, cell_id), moment in sorted(cell_moments.items()):
                indices = unflatten_cell(cell_id, partition)
                row: dict[str, object] = {
                    "stratum_id": stratum_id,
                    "channel": channel,
                    "cell_id": cell_id,
                    "contributing_rows": moment.count,
                    f"{prefix}_sum_trial_contributions_microbarn": moment.total,
                    f"{prefix}_sum_squared_trial_contributions_microbarn2": (
                        moment.square_total
                    ),
                    f"{prefix}_largest_trial_contribution_microbarn": (
                        moment.maximum
                    ),
                    f"{prefix}_cross_section_microbarn": (
                        moment.total / campaign.proposals
                    ),
                    "selected_seed": int((channel, cell_id) in seeds),
                    "in_core": int((channel, cell_id) in core),
                    "in_one_more_dilation": int(
                        (channel, cell_id) in expanded
                    ),
                }
                row.update(
                    {
                        f"{axis.name}_index": index
                        for axis, index in zip(partition, indices)
                    }
                )
                writer.writerow(row)


def _write_sha256(path: Path) -> None:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    (path.parent / f"{path.name}.sha256").write_text(
        f"{digest.hexdigest()}  {path.name}\n", encoding="utf-8"
    )


def _current_revision() -> str:
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).resolve().parent,
            text=True,
            capture_output=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return completed.stdout.strip()


def learn_guards(args: argparse.Namespace) -> dict:
    if args.output.resolve().exists():
        raise FileExistsError(
            f"{args.output.resolve()} already exists; guard artifacts are immutable"
        )
    if not 0.0 < args.target_core_fraction <= 1.0:
        raise GuardLearningError("--target-core-fraction must be in (0,1]")
    if not 0.0 < args.core_probability < 1.0:
        raise GuardLearningError("--core-probability must be in (0,1)")
    if not 0.0 < args.uniform_cell_floor < 1.0:
        raise GuardLearningError("--uniform-cell-floor must be in (0,1)")
    if args.minimum_training_rows < 1 or args.minimum_training_ess < 0.0:
        raise GuardLearningError("training support thresholds are invalid")
    config, config_sha = load_analysis_config(args.config.resolve())
    partition = parse_partition(args.partition)
    campaign = aggregate_surveys(
        args.survey, config, partition, apply_y_max=args.apply_y_max
    )
    if not math.isclose(
        float(config["beam_energy"]),
        float(campaign.norm_reference["ebeam"]),
        rel_tol=2.0e-7,
    ):
        raise GuardLearningError("analysis and survey beam energies differ")
    learner_revision = _current_revision()
    generator_revision = args.generator_revision or learner_revision
    manifest = build_manifest(
        config=config,
        config_path=args.config,
        config_sha256=config_sha,
        campaign=campaign,
        partition=partition,
        target_core_fraction=args.target_core_fraction,
        core_mixture_probability=args.core_probability,
        uniform_floor_fraction=args.uniform_cell_floor,
        dilation_radius=args.dilation,
        guard_iteration=args.iteration,
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
    path = _write_json_artifact(args.output.resolve(), "guard_manifest.json", manifest)
    coverage_path = path.parent / "training_coverage.csv"
    cells_path = path.parent / "training_cells.csv"
    _write_learning_csv(coverage_path, manifest)
    _write_cell_statistics(
        cells_path,
        manifest=manifest,
        campaign=campaign,
        partition=partition,
        validation=False,
    )
    _write_sha256(coverage_path)
    _write_sha256(cells_path)
    return {
        "passed": True,
        "schema": MANIFEST_SCHEMA,
        "manifest": str(path),
        "manifest_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "training_replicas": manifest["training"]["replica_ids"],
        "total_proposals": campaign.proposals,
        **manifest["learning"],
    }


def validate_guards(args: argparse.Namespace) -> dict:
    if args.output.resolve().exists():
        raise FileExistsError(
            f"{args.output.resolve()} already exists; guard artifacts are immutable"
        )
    if not 0.0 <= args.minimum_core_fraction <= 1.0:
        raise GuardLearningError("--minimum-core-fraction must be in [0,1]")
    manifest_path = args.manifest.resolve()
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GuardLearningError(f"cannot read manifest {manifest_path}: {error}") from error
    partition = _partition_from_manifest(manifest)
    config = manifest["analysis_config"]
    campaign = aggregate_surveys(
        args.survey,
        config,
        partition,
        apply_y_max=bool(
            manifest.get("analysis_selection", {}).get("apply_y_max", False)
        ),
    )
    validation = validate_manifest(
        manifest,
        manifest_path,
        campaign,
        partition,
        args.minimum_core_fraction,
    )
    path = _write_json_artifact(
        args.output.resolve(), "guard_validation.json", validation
    )
    coverage_path = path.parent / "validation_coverage.csv"
    cells_path = path.parent / "validation_cells.csv"
    _write_validation_csv(coverage_path, validation)
    _write_cell_statistics(
        cells_path,
        manifest=manifest,
        campaign=campaign,
        partition=partition,
        validation=True,
    )
    _write_sha256(coverage_path)
    _write_sha256(cells_path)
    return {
        "passed": validation["passed"],
        "schema": VALIDATION_SCHEMA,
        "validation": str(path),
        "validation_replicas": validation["validation"]["replica_ids"],
        **validation["coverage_summary"],
    }


def _aggregate_record_moments(records: Iterable[dict], key: str) -> Moment:
    result = Moment()
    for record in records:
        metrics = record.get(key)
        if not metrics:
            continue
        result.count += int(metrics["contributing_rows"])
        result.total += float(metrics["sum_trial_contributions_microbarn"])
        result.square_total += float(
            metrics["sum_squared_trial_contributions_microbarn2"]
        )
        result.maximum = max(
            result.maximum,
            float(metrics["largest_trial_contribution_microbarn"]),
        )
    return result


def _compact_metrics(moment: Moment, proposals: int) -> dict[str, object]:
    metrics = _metrics(moment, proposals)
    return {
        "contributing_rows": metrics["contributing_rows"],
        "cross_section_microbarn": metrics["cross_section_microbarn"],
        "cross_section_sem_microbarn": metrics[
            "cross_section_sem_microbarn"
        ],
        "cross_section_relative_standard_error": metrics[
            "cross_section_relative_standard_error"
        ],
        "importance_effective_sample_size": metrics[
            "importance_effective_sample_size"
        ],
        "largest_trial_contribution_microbarn": metrics[
            "largest_trial_contribution_microbarn"
        ],
    }


def _weighted_coverage(
    records: list[dict],
    *,
    proposals: int,
    total_key: str,
    core_key: str,
    tail_key: str,
    expanded_key: str,
) -> dict[str, object]:
    total = _aggregate_record_moments(records, total_key)
    core = _aggregate_record_moments(records, core_key)
    tail = _aggregate_record_moments(records, tail_key)
    expanded = _aggregate_record_moments(records, expanded_key)
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
        "total": _compact_metrics(total, proposals),
        "core": _compact_metrics(core, proposals),
        "tail": _compact_metrics(tail, proposals),
        "one_more_dilation": _compact_metrics(expanded, proposals),
        "cross_section_weighted_core_fraction": (
            core.total / total.total if total.total > 0.0 else None
        ),
        "cross_section_weighted_tail_fraction": (
            tail.total / total.total if total.total > 0.0 else None
        ),
        "cross_section_weighted_one_more_dilation_fraction": (
            expanded.total / total.total if total.total > 0.0 else None
        ),
        "extra_cross_section_fraction_recovered_by_one_more_dilation": (
            (expanded.total - core.total) / total.total
            if total.total > 0.0
            else None
        ),
    }


def _material_coverage(
    records: list[dict],
    *,
    proposals: int,
    targets: tuple[float, ...] = (0.5, 0.9, 0.95, 0.99),
) -> dict[str, object]:
    contributing = [
        record
        for record in records
        if float(
            record["holdout_total"]["sum_trial_contributions_microbarn"]
        )
        > 0.0
    ]
    contributing.sort(
        key=lambda record: -float(
            record["holdout_total"]["sum_trial_contributions_microbarn"]
        )
    )
    total = sum(
        float(record["holdout_total"]["sum_trial_contributions_microbarn"])
        for record in contributing
    )
    result: dict[str, object] = {}
    for target in targets:
        retained = 0.0
        stop = 0
        for stop, record in enumerate(contributing, 1):
            retained += float(
                record["holdout_total"][
                    "sum_trial_contributions_microbarn"
                ]
            )
            if retained >= target * total:
                break
        selected = contributing[:stop]
        coverage = _weighted_coverage(
            selected,
            proposals=proposals,
            total_key="holdout_total",
            core_key="holdout_core",
            tail_key="holdout_tail",
            expanded_key="holdout_one_more_dilation",
        )
        result[f"top_{100.0 * target:g}_percent_cross_section"] = {
            "strata": stop,
            "actual_fraction_of_inside_cross_section": (
                retained / total if total > 0.0 else None
            ),
            "minimum_stratum_cross_section_microbarn": (
                float(
                    selected[-1]["holdout_total"][
                        "cross_section_microbarn"
                    ]
                )
                if selected
                else None
            ),
            "cross_section_weighted_core_fraction": coverage[
                "cross_section_weighted_core_fraction"
            ],
            "cross_section_weighted_one_more_dilation_fraction": coverage[
                "cross_section_weighted_one_more_dilation_fraction"
            ],
        }
    return result


def _channel_cross_section_fractions(records: list[dict]) -> dict[str, float]:
    totals: dict[str, float] = {}
    total_cross_section = 0.0
    for record in records:
        cross_section = float(
            record["holdout_total"]["cross_section_microbarn"]
        )
        total_cross_section += cross_section
        for channel, channel_fraction in record[
            "channel_cross_section_fractions"
        ].items():
            totals[channel] = totals.get(channel, 0.0) + (
                cross_section * float(channel_fraction)
            )
    return {
        channel: value / total_cross_section
        for channel, value in sorted(totals.items(), key=lambda item: int(item[0]))
        if total_cross_section > 0.0
    }


def summarize_coverage(args: argparse.Namespace) -> dict:
    if args.limit < 0:
        raise GuardLearningError("--limit must be nonnegative")
    manifest = json.loads(args.manifest.resolve().read_text(encoding="utf-8"))
    manifest_sha256 = hashlib.sha256(
        args.manifest.resolve().read_bytes()
    ).hexdigest()
    training_records = list(manifest["strata"].values())
    training_proposals = int(manifest["training"]["total_proposals"])
    summary = {
        "manifest_schema": manifest["schema"],
        "manifest_sha256": manifest_sha256,
        "analysis_selection": manifest["analysis_selection"],
        "partition_axes": manifest["partition_definition"]["axes"],
        "training_replicas": manifest["training"]["replica_ids"],
        "training_global_cross_section_microbarn": manifest["training"][
            "global_observed"
        ]["cross_section_microbarn"],
        "weighted_coverage": {
            "training_inside_analysis_partition": _weighted_coverage(
                training_records,
                proposals=training_proposals,
                total_key="training_total",
                core_key="training_core",
                tail_key="training_tail",
                expanded_key="training_one_more_dilation",
            )
        },
        **manifest["learning"],
    }
    if args.validation:
        validation = json.loads(
            args.validation.resolve().read_text(encoding="utf-8")
        )
        if validation.get("manifest_sha256") != manifest_sha256:
            raise GuardLearningError(
                "validation artifact does not reference the supplied manifest"
            )
        validation_records = list(validation["strata"].values())
        validation_proposals = int(
            validation["validation"]["total_proposals"]
        )
        assessed = [
            (stratum_id, record)
            for stratum_id, record in validation["strata"].items()
            if record["holdout_core_fraction"] is not None
        ]
        worst = sorted(
            assessed,
            key=lambda item: (
                item[1]["holdout_core_fraction"],
                item[0],
            ),
        )[: args.limit]
        validation_weighted = _weighted_coverage(
            validation_records,
            proposals=validation_proposals,
            total_key="holdout_total",
            core_key="holdout_core",
            tail_key="holdout_tail",
            expanded_key="holdout_one_more_dilation",
        )
        by_training_status = {
            status: _weighted_coverage(
                [
                    record
                    for record in validation_records
                    if record["training_status"] == status
                ],
                proposals=validation_proposals,
                total_key="holdout_total",
                core_key="holdout_core",
                tail_key="holdout_tail",
                expanded_key="holdout_one_more_dilation",
            )
            for status in (
                "learned",
                "learned_low_support",
                "no_training_contribution",
            )
        }
        by_coverage_result = {
            label: _weighted_coverage(
                [
                    record
                    for record in validation_records
                    if record["coverage_passed"] is expected
                ],
                proposals=validation_proposals,
                total_key="holdout_total",
                core_key="holdout_core",
                tail_key="holdout_tail",
                expanded_key="holdout_one_more_dilation",
            )
            for label, expected in (("passed", True), ("failed", False))
        }
        training_inside = manifest["training"]["inside_analysis_partition"]
        validation_inside = validation["validation"][
            "inside_analysis_partition"
        ]
        summary["weighted_coverage"].update(
            {
                "validation_inside_analysis_partition": validation_weighted,
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
                    _difference_z_score(
                        float(training_inside["cross_section_microbarn"]),
                        training_inside["cross_section_sem_microbarn"],
                        float(validation_inside["cross_section_microbarn"]),
                        validation_inside["cross_section_sem_microbarn"],
                    )
                ),
                "validation_by_training_status": by_training_status,
                "validation_by_coverage_result": by_coverage_result,
                "validation_material_strata": _material_coverage(
                    validation_records, proposals=validation_proposals
                ),
                "validation_radiative_channel_cross_section_fractions": (
                    _channel_cross_section_fractions(validation_records)
                ),
            }
        )
        summary.update(
            {
                "validation_schema": validation["schema"],
                "validation_passed": validation["passed"],
                "validation_replicas": validation["validation"]["replica_ids"],
                "validation_global_cross_section_microbarn": validation[
                    "validation"
                ]["global_observed"]["cross_section_microbarn"],
                "global_training_holdout_difference_z_score": validation[
                    "global_training_holdout_difference_z_score"
                ],
                "worst_holdout_core_fractions": [
                    {
                        "stratum_id": stratum_id,
                        "holdout_core_fraction": record[
                            "holdout_core_fraction"
                        ],
                        "holdout_cross_section_microbarn": record[
                            "holdout_total"
                        ]["cross_section_microbarn"],
                        "holdout_ess": record["holdout_total"][
                            "importance_effective_sample_size"
                        ],
                        "largest_tail_fraction": (
                            record["largest_tail_cells"][0][
                                "fraction_of_holdout_cross_section"
                            ]
                            if record["largest_tail_cells"]
                            else None
                        ),
                    }
                    for stratum_id, record in worst
                ],
                **validation["coverage_summary"],
            }
        )
    return summary


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    learn = subparsers.add_parser(
        "learn-guards", help="learn and freeze guards from training replicas"
    )
    learn.add_argument("--config", type=Path, required=True)
    learn.add_argument("--survey", type=Path, nargs="+", required=True)
    learn.add_argument("--output", type=Path, required=True)
    learn.add_argument("--target-core-fraction", type=float, default=0.995)
    learn.add_argument("--core-probability", type=float, default=0.98)
    learn.add_argument("--uniform-cell-floor", type=float, default=0.05)
    learn.add_argument("--dilation", type=int, default=1)
    learn.add_argument("--iteration", type=int, default=0)
    learn.add_argument("--minimum-training-rows", type=int, default=10)
    learn.add_argument("--minimum-training-ess", type=float, default=5.0)
    learn.add_argument(
        "--apply-y-max",
        action="store_true",
        help=(
            "apply phase_space.y_max when assigning observed strata; the "
            "default applies no y cut"
        ),
    )
    learn.add_argument(
        "--partition",
        action="append",
        help="replace defaults with repeated survey_column=number_of_bins",
    )
    learn.add_argument(
        "--generator-revision",
        help=(
            "Git revision that produced the surveys; defaults to the current "
            "checkout and is recorded as an assumption"
        ),
    )

    validate = subparsers.add_parser(
        "validate-guards", help="measure a frozen guard on held-out replicas"
    )
    validate.add_argument("--manifest", type=Path, required=True)
    validate.add_argument("--survey", type=Path, nargs="+", required=True)
    validate.add_argument("--output", type=Path, required=True)
    validate.add_argument("--minimum-core-fraction", type=float, default=0.98)

    summarize = subparsers.add_parser(
        "summarize-coverage", help="print compact training/validation coverage"
    )
    summarize.add_argument("--manifest", type=Path, required=True)
    summarize.add_argument("--validation", type=Path)
    summarize.add_argument(
        "--limit",
        type=int,
        default=20,
        help="number of lowest-coverage validation strata to show",
    )
    return parser


def main() -> int:
    args = _build_parser().parse_args()
    try:
        if args.command == "learn-guards":
            result = learn_guards(args)
        elif args.command == "validate-guards":
            result = validate_guards(args)
        else:
            result = summarize_coverage(args)
    except (
        FileNotFoundError,
        FileExistsError,
        GuardLearningError,
        KeyError,
        ValueError,
    ) as error:
        print(f"ERROR: {error}")
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result.get("passed", True) else 2


if __name__ == "__main__":
    raise SystemExit(main())
