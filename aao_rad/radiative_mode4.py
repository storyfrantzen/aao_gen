#!/usr/bin/env python3
"""Calibrate, generate, validate, and finalize radiative AAO mode-4 strata.

Mode 4 generates one final-LUND analysis stratum per invocation from an exact
mixture of a learned continuous native-coordinate core and the unrestricted
legacy proposal.  The legacy component has strictly positive probability, so
an imperfect learned guard can reduce efficiency but cannot remove physical
support.  Its fixed-trial calibration operation measures a practical
acceptance envelope without emitting LUND events.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import shutil
import statistics
import struct
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import radiative_guards
import radiative_survey


MANIFEST_SCHEMA = "aao-rad-mode4-manifest-v4"
LEGACY_MANIFEST_SCHEMAS = ("aao-rad-mode4-manifest-v3",)
RUN_SCHEMA = "aao-rad-mode4-run-v3"
WEIGHTS_SCHEMA = "aao-rad-mode4-weights-v2"
CALIBRATION_SCHEMA = "aao-rad-mode4-envelope-calibration-v4"
PILOT_VALIDATION_SCHEMA = "aao-rad-mode4-pilot-validation-v1"
RECIPE_SCHEMA = "aao-rad-continuous-guard-recipes-v1"
REFINEMENT_SCHEMA = "aao-rad-mode4-guard-refinements-v1"
REFINEMENT_COORDINATE_SPACE = "post_padding_guard_bounds"
SAMPLING_MODE = 4
MODE4_KINEMATICS_FILENAME = "aao_rad.mode4.csv"
MODE4_KINEMATICS_SCHEMA = "aao-rad-mode4-events-v1"
MODE4_KINEMATICS_COLUMNS = (
    "event",
    "proposal_component",
    "proposal_density_ratio",
    "r_u",
    "r_ep",
    "u_gamma",
    "hadron_cosine_base",
    "hadron_phi_base",
    "q2_observed",
    "xb_observed",
    "minus_t_observed",
    "phi_observed_deg",
    "w_observed",
    "y_observed",
    "integrand_corrected",
)
MODE4_HEARTBEAT_FILENAME = "aao_rad.mode4.heartbeat.csv"
MODE4_HEARTBEAT_SCHEMA = "aao-rad-mode4-heartbeat-v3"
MODE4_HEARTBEAT_COLUMNS = (
    "proposals",
    "events",
    "core_trials",
    "noncore_trials",
    "internal_rows",
    "final_candidates",
    "target_candidates",
    "core_targets",
    "noncore_targets",
    "core_events",
    "noncore_events",
    "emitting_candidates",
    "duplicate_events",
    "mcall_max",
)
MODE4_CALIBRATION_FILENAME = "aao_rad.mode4.calibration.csv"
MODE4_CALIBRATION_EVENT_SCHEMA = "aao-rad-mode4-calibration-events-v2"
MODE4_CALIBRATION_COLUMNS = (
    "trial",
    "proposal_component",
    "inside_core",
    "proposal_density_ratio",
    "integrand_corrected",
    "component_importance_weight",
    "r_u",
    "r_ep",
    "u_gamma",
    "hadron_cosine_base",
    "hadron_phi_base",
    "q2_observed",
    "xb_observed",
    "minus_t_observed",
    "phi_observed_deg",
    "w_observed",
    "y_observed",
)
AXES = (
    "r_u",
    "r_ep",
    "u_gamma",
    "hadron_cosine_base",
    "hadron_phi_base",
)


class Mode4Error(RuntimeError):
    """Raised when a mode-4 artifact violates a required invariant."""


@dataclass(frozen=True)
class GuardBox:
    nonperiodic: dict[str, tuple[float, float]]
    phi_origin: float
    phi_relative: tuple[float, float]

    @property
    def volume(self) -> float:
        widths = [
            self.nonperiodic[name][1] - self.nonperiodic[name][0]
            for name in AXES[:-1]
        ]
        widths.append(self.phi_relative[1] - self.phi_relative[0])
        return math.prod(widths)

    def contains(self, coordinates: dict[str, float]) -> bool:
        for name in AXES[:-1]:
            lower, upper = self.nonperiodic[name]
            if not lower <= coordinates[name] <= upper:
                return False
        relative = (
            coordinates["hadron_phi_base"] - self.phi_origin + 0.5
        ) % 1.0 - 0.5
        return self.phi_relative[0] <= relative <= self.phi_relative[1]

    def manifest_record(self) -> dict[str, object]:
        return {
            "axes": {
                name: list(self.nonperiodic[name]) for name in AXES[:-1]
            },
            "hadron_phi_base": {
                "origin": self.phi_origin,
                "lower_relative_to_origin": self.phi_relative[0],
                "upper_relative_to_origin": self.phi_relative[1],
                "width": self.phi_relative[1] - self.phi_relative[0],
            },
            "normalized_volume": self.volume,
            "u_gamma_soft_endpoint_anchored": True,
        }


def _fortran_real32(value: float) -> float:
    """Round one value exactly as the generator's default REAL does."""
    return struct.unpack("=f", struct.pack("=f", float(value)))[0]


def _fortran_guard_contains(
    box: GuardBox, coordinates: dict[str, float]
) -> bool:
    """Reproduce mode4_proposal_ratio's single-precision guard test.

    Guard bounds are serialized from Python doubles but read by AAO into
    default Fortran REAL variables.  Diagnostic coordinates are also default
    REAL values.  Replaying the comparison in Python double precision can
    disagree on the one-float shell at a rounded face, so artifact validation
    must use the generator's actual numerical guard.
    """
    for name in AXES[:-1]:
        value = _fortran_real32(coordinates[name])
        lower = _fortran_real32(box.nonperiodic[name][0])
        upper = _fortran_real32(box.nonperiodic[name][1])
        if value < lower or value > upper:
            return False
    phi = _fortran_real32(coordinates["hadron_phi_base"])
    origin = _fortran_real32(box.phi_origin)
    relative = _fortran_real32(phi - origin)
    relative = _fortran_real32(relative + _fortran_real32(0.5))
    relative = _fortran_real32(relative % _fortran_real32(1.0))
    relative = _fortran_real32(relative - _fortran_real32(0.5))
    lower = _fortran_real32(box.phi_relative[0])
    upper = _fortran_real32(box.phi_relative[1])
    return lower <= relative <= upper


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_revision() -> str:
    """Best-effort revision of the repository containing this wrapper."""
    completed = subprocess.run(
        [
            "git",
            "-C",
            str(Path(__file__).resolve().parent),
            "rev-parse",
            "HEAD",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    revision = completed.stdout.strip()
    return revision if completed.returncode == 0 and revision else "UNKNOWN"


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _supported_manifest_schema(payload: dict) -> bool:
    return payload.get("schema") in (MANIFEST_SCHEMA, *LEGACY_MANIFEST_SCHEMAS)


def _run_sigr_max(manifest: dict, record: dict) -> float:
    """Resolve a run envelope, including legacy shared-envelope manifests."""
    raw = record.get("sigr_max", manifest.get("sigr_max"))
    try:
        value = float(raw)
    except (TypeError, ValueError) as error:
        raise Mode4Error("manifest run lacks a valid sigr_max") from error
    if not math.isfinite(value) or value <= 0.0:
        raise Mode4Error("manifest run sigr_max must be finite and positive")
    return value


def _candidate_padding(recipes: dict, candidate: str) -> float:
    matches = [
        item
        for item in recipes.get("padding_candidates", [])
        if item.get("identifier") == candidate
    ]
    if len(matches) != 1:
        raise ValueError(
            f"guard recipes contain {len(matches)} entries for {candidate!r}"
        )
    try:
        padding = float(matches[0]["base_padding"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{candidate}: malformed base padding") from error
    if not math.isfinite(padding) or padding < 0.0:
        raise ValueError(f"{candidate}: padding must be finite and nonnegative")
    return padding


def reconstruct_guard_box(recipe: dict, base_padding: float) -> GuardBox:
    """Reconstruct the exact support-adaptively padded recipe box."""
    try:
        scale = float(recipe["padding_scale"])
        raw_axes = {
            item["name"]: item for item in recipe["raw_box"]["axes"]
        }
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("malformed continuous guard recipe") from error
    if not math.isfinite(scale) or not 1.0 <= scale <= 2.0:
        raise ValueError("guard padding scale must lie in [1,2]")
    amount = base_padding * scale
    nonperiodic: dict[str, tuple[float, float]] = {}
    for name in AXES[:-1]:
        try:
            axis = raw_axes[name]
            if bool(axis["periodic"]):
                raise ValueError(f"{name} is unexpectedly periodic")
            lower = max(0.0, float(axis["lower"]) - amount)
            upper = (
                1.0
                if name == "u_gamma"
                else min(1.0, float(axis["upper"]) + amount)
            )
        except (KeyError, TypeError, ValueError) as error:
            if isinstance(error, ValueError) and "unexpectedly" in str(error):
                raise
            raise ValueError(f"malformed {name} guard axis") from error
        if not 0.0 <= lower < upper <= 1.0:
            raise ValueError(f"{name} guard bounds are invalid")
        nonperiodic[name] = (lower, upper)
    try:
        phi = raw_axes["hadron_phi_base"]
        if not bool(phi["periodic"]):
            raise ValueError("hadron_phi_base must be periodic")
        origin = float(phi["origin"]) % 1.0
        lower_relative = max(
            -0.5, float(phi["lower_relative_to_origin"]) - amount
        )
        upper_relative = min(
            0.5, float(phi["upper_relative_to_origin"]) + amount
        )
    except (KeyError, TypeError, ValueError) as error:
        if isinstance(error, ValueError) and "must be periodic" in str(error):
            raise
        raise ValueError("malformed periodic guard axis") from error
    if not -0.5 <= lower_relative < upper_relative <= 0.5:
        raise ValueError("periodic guard bounds are invalid")
    box = GuardBox(
        nonperiodic=nonperiodic,
        phi_origin=origin,
        phi_relative=(lower_relative, upper_relative),
    )
    if not 0.0 < box.volume <= 1.0:
        raise ValueError("guard volume must lie in (0,1]")
    return box


def apply_guard_refinement(
    box: GuardBox,
    specification: dict,
    *,
    stratum_id: str,
) -> tuple[GuardBox, dict[str, object]]:
    """Expand selected final guard faces and describe every applied change."""
    try:
        rationale = str(specification["rationale"]).strip()
        evidence = specification["evidence"]
        faces = specification["faces"]
    except (KeyError, TypeError) as error:
        raise ValueError(
            f"{stratum_id}: malformed guard refinement"
        ) from error
    if not rationale:
        raise ValueError(f"{stratum_id}: refinement rationale is empty")
    if not isinstance(evidence, list) or not evidence:
        raise ValueError(
            f"{stratum_id}: refinement must cite at least one evidence artifact"
        )
    for item in evidence:
        if not isinstance(item, dict):
            raise ValueError(f"{stratum_id}: malformed refinement evidence")
        try:
            evidence_path = str(item["path"]).strip()
            evidence_hash = str(item["sha256"]).strip().lower()
        except (KeyError, TypeError) as error:
            raise ValueError(
                f"{stratum_id}: malformed refinement evidence"
            ) from error
        if not evidence_path or len(evidence_hash) != 64 or any(
            character not in "0123456789abcdef" for character in evidence_hash
        ):
            raise ValueError(
                f"{stratum_id}: malformed refinement evidence provenance"
            )
    if not isinstance(faces, dict) or not faces:
        raise ValueError(f"{stratum_id}: refinement faces are empty")
    unknown_axes = set(faces).difference(AXES)
    if unknown_axes:
        raise ValueError(
            f"{stratum_id}: unknown refinement axes "
            f"{sorted(unknown_axes)}"
        )

    nonperiodic = dict(box.nonperiodic)
    phi_relative = list(box.phi_relative)
    changes: list[dict[str, object]] = []
    for axis, requested in faces.items():
        if not isinstance(requested, dict) or not requested:
            raise ValueError(
                f"{stratum_id}: {axis} refinement faces are empty"
            )
        unknown_faces = set(requested).difference(("lower", "upper"))
        if unknown_faces:
            raise ValueError(
                f"{stratum_id}: {axis} has unknown faces "
                f"{sorted(unknown_faces)}"
            )
        if axis == "hadron_phi_base":
            old_bounds = box.phi_relative
            domain = (-0.5, 0.5)
        else:
            old_bounds = box.nonperiodic[axis]
            domain = (0.0, 1.0)
        new_bounds = list(old_bounds)
        for face, raw_value in requested.items():
            try:
                value = float(raw_value)
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"{stratum_id}: {axis} {face} is not numeric"
                ) from error
            if not math.isfinite(value) or not domain[0] <= value <= domain[1]:
                raise ValueError(
                    f"{stratum_id}: {axis} {face} lies outside "
                    f"[{domain[0]},{domain[1]}]"
                )
            index = 0 if face == "lower" else 1
            original = old_bounds[index]
            if face == "lower" and value > original:
                raise ValueError(
                    f"{stratum_id}: {axis} lower refinement would contract "
                    "the guard"
                )
            if face == "upper" and value < original:
                raise ValueError(
                    f"{stratum_id}: {axis} upper refinement would contract "
                    "the guard"
                )
            new_bounds[index] = value
            if not _close(value, original):
                changes.append(
                    {
                        "axis": axis,
                        "face": face,
                        "original": original,
                        "refined": value,
                        "signed_change": value - original,
                    }
                )
        if not new_bounds[0] < new_bounds[1]:
            raise ValueError(
                f"{stratum_id}: refined {axis} bounds are invalid"
            )
        if axis == "hadron_phi_base":
            phi_relative = new_bounds
        else:
            nonperiodic[axis] = (new_bounds[0], new_bounds[1])
    if not changes:
        raise ValueError(f"{stratum_id}: refinement does not move any face")
    if not _close(nonperiodic["u_gamma"][1], 1.0):
        raise ValueError(
            f"{stratum_id}: refinement breaks the u_gamma endpoint anchor"
        )
    refined = GuardBox(
        nonperiodic=nonperiodic,
        phi_origin=box.phi_origin,
        phi_relative=(phi_relative[0], phi_relative[1]),
    )
    if not box.volume < refined.volume <= 1.0:
        raise ValueError(
            f"{stratum_id}: refinement must increase guard volume within "
            "the native domain"
        )
    record = {
        "coordinate_space": REFINEMENT_COORDINATE_SPACE,
        "rationale": rationale,
        "evidence": evidence,
        "applied_face_changes": changes,
        "original_normalized_volume": box.volume,
        "refined_normalized_volume": refined.volume,
        "volume_ratio": refined.volume / box.volume,
    }
    return refined, record


def proposal_density_ratio(
    coordinates: dict[str, float],
    box: GuardBox,
    core_fraction: float,
) -> float:
    """Return q_legacy/q_mix for the core-plus-legacy proposal."""
    return _proposal_density_ratio_for_membership(
        box, core_fraction, box.contains(coordinates)
    )


def _proposal_density_ratio_for_membership(
    box: GuardBox,
    core_fraction: float,
    inside: bool,
) -> float:
    """Return the proposal correction for an established guard class."""
    if not 0.0 < core_fraction < 1.0:
        raise ValueError("core fraction must lie strictly between zero and one")
    tail_fraction = 1.0 - core_fraction
    denominator = tail_fraction
    if inside:
        denominator += core_fraction / box.volume
    return 1.0 / denominator


def _guard_box_from_manifest(record: dict) -> GuardBox:
    try:
        guard = record["guard"]
        axes = guard["axes"]
        phi = guard["hadron_phi_base"]
        box = GuardBox(
            nonperiodic={
                name: (float(axes[name][0]), float(axes[name][1]))
                for name in AXES[:-1]
            },
            phi_origin=float(phi["origin"]),
            phi_relative=(
                float(phi["lower_relative_to_origin"]),
                float(phi["upper_relative_to_origin"]),
            ),
        )
    except (KeyError, IndexError, TypeError, ValueError) as error:
        raise Mode4Error("manifest contains a malformed guard box") from error
    if not _close(box.volume, float(guard["normalized_volume"])):
        raise Mode4Error("manifest guard volume differs from its bounds")
    if not bool(guard.get("u_gamma_soft_endpoint_anchored")):
        raise Mode4Error("manifest guard lacks the soft-photon endpoint anchor")
    if not _close(box.nonperiodic["u_gamma"][1], 1.0):
        raise Mode4Error("manifest guard does not end at u_gamma=1")
    return box


def _legacy_mode4_input(
    source: str,
    source_path: Path,
    *,
    events: int,
    sigr_max: float,
) -> tuple[str, list[str]]:
    records = radiative_survey._records(source)
    count = radiative_survey._legacy_record_count(records, source_path)
    if len(records) != count:
        raise ValueError(
            f"{source_path}: pass a legacy input without a sampling trailer"
        )
    try:
        theory = int(records[0].split()[0])
        original_fmcall = float(records[16].split()[0])
    except (ValueError, IndexError) as error:
        raise ValueError(f"{source_path}: malformed legacy input") from error
    position = 17
    if original_fmcall == 0.0:
        position += 1
    w_minimum_record = records[position] if theory > 10 else None
    rewritten = list(records[:17])
    rewritten[15] = str(events)
    rewritten[16] = "0"
    rewritten.append(f"{sigr_max:.17g}")
    if w_minimum_record is not None:
        rewritten.append(w_minimum_record)
    text = "\n".join(rewritten) + "\n"
    radiative_survey._validate_legacy_input_shape(text, source_path)
    return text, rewritten


def _mode4_trailer(
    *,
    seed: int,
    replica: int,
    stratum: radiative_guards.Stratum,
    w_minimum: float,
    apply_y_max: bool,
    y_maximum: float,
    core_fraction: float,
    base_padding: float,
    box: GuardBox,
    operation: str,
    trial_limit: int,
    heartbeat_interval: int,
    component_core_fraction: float,
    calibration_proposal: int,
) -> str:
    axes = box.nonperiodic
    return "\n".join(
        [
            str(SAMPLING_MODE),
            str(seed),
            str(replica),
            (
                f"{stratum.flat_index} {stratum.iq2} {stratum.ixb} "
                f"{stratum.it} {stratum.iphi}"
            ),
            f"{stratum.q2[0]:.17g} {stratum.q2[1]:.17g}",
            f"{stratum.xb[0]:.17g} {stratum.xb[1]:.17g}",
            f"{stratum.minus_t[0]:.17g} {stratum.minus_t[1]:.17g}",
            f"{stratum.phi_deg[0]:.17g} {stratum.phi_deg[1]:.17g}",
            f"{w_minimum:.17g}",
            f"{int(apply_y_max)} {y_maximum:.17g}",
            f"{core_fraction:.17g}",
            f"{base_padding:.17g}",
            f"{box.volume:.17g}",
            f"{axes['r_u'][0]:.17g} {axes['r_u'][1]:.17g}",
            f"{axes['r_ep'][0]:.17g} {axes['r_ep'][1]:.17g}",
            f"{axes['u_gamma'][0]:.17g} {axes['u_gamma'][1]:.17g}",
            (
                f"{axes['hadron_cosine_base'][0]:.17g} "
                f"{axes['hadron_cosine_base'][1]:.17g}"
            ),
            (
                f"{box.phi_origin:.17g} {box.phi_relative[0]:.17g} "
                f"{box.phi_relative[1]:.17g}"
            ),
            str(1 if operation == "calibration" else 0),
            str(trial_limit),
            str(heartbeat_interval),
            f"{component_core_fraction:.17g}",
            str(calibration_proposal),
        ]
    ) + "\n"


def _load_config_and_recipes(
    config_path: Path,
    recipes_path: Path,
) -> tuple[dict, dict, str, str]:
    config_raw = config_path.read_bytes()
    recipes_raw = recipes_path.read_bytes()
    try:
        config = json.loads(config_raw)
        recipes = json.loads(recipes_raw)
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid JSON input: {error}") from error
    if recipes.get("schema") != RECIPE_SCHEMA:
        raise ValueError(
            f"unsupported recipe schema: {recipes.get('schema')!r}"
        )
    config_sha256 = hashlib.sha256(config_raw).hexdigest()
    if recipes.get("analysis_config_sha256") != config_sha256:
        raise ValueError("guard recipes and analysis config hashes differ")
    return (
        config,
        recipes,
        config_sha256,
        hashlib.sha256(recipes_raw).hexdigest(),
    )


def _load_configuration(
    config_path: Path,
    recipes_path: Path,
    input_path: Path,
) -> tuple[dict, dict, str, str, str]:
    (
        config,
        recipes,
        config_sha256,
        recipes_sha256,
    ) = _load_config_and_recipes(config_path, recipes_path)
    legacy = input_path.read_text(encoding="utf-8")
    radiative_survey._load_balanced_config(
        config_path,
        legacy_input=legacy,
        legacy_path=input_path,
    )
    try:
        phase_space = config["phase_space"]
        w_minimum = float(phase_space["W_min"])
        y_maximum = float(phase_space.get("y_max", 1.0))
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("analysis config lacks valid phase-space metadata") from error
    if w_minimum <= radiative_survey.PROTON_MASS_GEV:
        raise ValueError("W_min is not physical")
    if not 0.0 < y_maximum <= 1.0:
        raise ValueError("y_max must lie in (0,1]")
    return (
        config,
        recipes,
        config_sha256,
        recipes_sha256,
        legacy,
    )


def _load_guard_refinements(
    path: Path,
    *,
    config_sha256: str,
    recipes_sha256: str,
    candidate: str,
    recipes: dict,
    base_padding: float,
) -> tuple[dict, str]:
    raw = path.read_bytes()
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as error:
        raise ValueError(f"{path}: invalid refinement JSON: {error}") from error
    if payload.get("schema") != REFINEMENT_SCHEMA:
        raise ValueError(
            f"{path}: unsupported refinement schema "
            f"{payload.get('schema')!r}"
        )
    if payload.get("coordinate_space") != REFINEMENT_COORDINATE_SPACE:
        raise ValueError(
            f"{path}: refinements must use final post-padding guard bounds"
        )
    expected = {
        "analysis_config_sha256": config_sha256,
        "guard_recipes_sha256": recipes_sha256,
        "guard_candidate": candidate,
    }
    for name, value in expected.items():
        if payload.get(name) != value:
            raise ValueError(
                f"{path}: {name} does not match the selected guard inputs"
            )
    specifications = payload.get("strata")
    if not isinstance(specifications, dict) or not specifications:
        raise ValueError(f"{path}: refinement strata are empty")
    for stratum_id, specification in specifications.items():
        if stratum_id not in recipes.get("strata", {}):
            raise ValueError(
                f"{path}: refinement names unknown stratum {stratum_id}"
            )
        original = reconstruct_guard_box(
            recipes["strata"][stratum_id], base_padding
        )
        apply_guard_refinement(
            original, specification, stratum_id=stratum_id
        )
    return payload, hashlib.sha256(raw).hexdigest()


def _analysis_selection(config: dict, *, apply_y_max: bool) -> dict[str, object]:
    return {
        "coordinate_definition": "final_lund_analysis",
        "w_minimum": float(config["phase_space"]["W_min"]),
        "apply_y_max": bool(apply_y_max),
        "y_maximum": (
            float(config["phase_space"].get("y_max", 1.0))
            if apply_y_max
            else None
        ),
        "no_implicit_y_minimum": True,
    }


def _selected_strata(
    strata: list[radiative_guards.Stratum], args: argparse.Namespace
) -> tuple[list[radiative_guards.Stratum], dict[str, object]]:
    raw_sparse = getattr(args, "flat_indices", None) or []
    raw_index_file = getattr(args, "flat_index_file", None)
    start = int(getattr(args, "bin_start", 0) or 0)
    stop = getattr(args, "bin_stop", None)
    source_path: Optional[Path] = None
    if raw_index_file is not None:
        if raw_sparse:
            raise ValueError(
                "--flat-index-file cannot be combined with --flat-index"
            )
        source_path = Path(raw_index_file).expanduser().resolve()
        if not source_path.is_file():
            raise FileNotFoundError(source_path)
        parsed: list[int] = []
        for line_number, raw_line in enumerate(
            source_path.read_text(encoding="utf-8").splitlines(), start=1
        ):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                parsed.append(int(line))
            except ValueError as error:
                raise ValueError(
                    f"{source_path}:{line_number}: expected one flat index"
                ) from error
        if not parsed:
            raise ValueError(f"{source_path}: no flat indices")
        raw_sparse = parsed
    if raw_sparse:
        if start != 0 or stop is not None:
            raise ValueError(
                "sparse flat-index selection cannot be combined with "
                "--bin-start or --bin-stop"
            )
        indices = [int(value) for value in raw_sparse]
        if len(set(indices)) != len(indices):
            raise ValueError("--flat-index values must be unique")
        if any(index < 0 or index >= len(strata) for index in indices):
            raise ValueError(
                f"--flat-index values must lie in [0,{len(strata)})"
            )
        ordered = sorted(indices)
        metadata: dict[str, object] = {
            "mode": (
                "flat_index_file"
                if source_path is not None
                else "sparse_flat_indices"
            ),
            "flat_indices": ordered,
        }
        if source_path is not None:
            metadata.update(
                {
                    "source": str(source_path),
                    "source_sha256": _sha256(source_path),
                    "snapshot": "flat_index_selection.txt",
                }
            )
        return [strata[index] for index in ordered], metadata
    resolved_stop = len(strata) if stop is None else int(stop)
    if start < 0 or resolved_stop < start or resolved_stop > len(strata):
        raise ValueError(f"invalid stratum range [{start},{resolved_stop})")
    return (
        strata[start:resolved_stop],
        {
            "mode": "contiguous_range",
            "bin_start": start,
            "bin_stop": resolved_stop,
        },
    )


def _load_generation_envelopes(
    path: Path,
    *,
    config_sha256: str,
    recipes_sha256: str,
    refinements_sha256: Optional[str],
    candidate: str,
    core_fraction: float,
    selection: dict[str, object],
    generator_revision: str,
    allow_revision_mismatch: bool,
    revision_compatibility_rationale: str,
) -> tuple[dict, str, dict[str, dict], dict[str, object]]:
    raw = path.read_bytes()
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as error:
        raise ValueError(f"{path}: invalid envelope JSON: {error}") from error
    if payload.get("schema") != CALIBRATION_SCHEMA:
        raise ValueError(
            f"{path}: expected calibration schema {CALIBRATION_SCHEMA}"
        )
    expected = {
        "analysis_config_sha256": config_sha256,
        "guard_recipes_sha256": recipes_sha256,
        "guard_refinements_sha256": refinements_sha256,
        "guard_candidate": candidate,
    }
    for name, value in expected.items():
        if payload.get(name) != value:
            raise ValueError(
                f"{path}: {name} does not match the selected generation inputs"
            )
    try:
        report_core_fraction = float(payload["core_fraction"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}: invalid calibration core fraction") from error
    if not _close(report_core_fraction, core_fraction):
        raise ValueError(
            f"{path}: core_fraction does not match the generation mixture"
        )
    report_selection, revisions, report_guards = (
        _calibration_validation_metadata(path, payload)
    )
    if report_selection != selection:
        raise ValueError(
            f"{path}: analysis selection does not match generation"
        )
    revision_matches = generator_revision in revisions
    if not revision_matches and not allow_revision_mismatch:
        raise ValueError(
            f"{path}: generator revision {generator_revision!r} is absent "
            "from the calibration report; use the explicit audited override"
        )
    if allow_revision_mismatch and not revision_compatibility_rationale:
        raise ValueError(
            "--envelope-revision-compatibility-rationale is required with "
            "--allow-envelope-revision-mismatch"
        )
    records: dict[str, dict] = {}
    for item in payload.get("strata", []):
        stratum_id = str(item.get("stratum_id", ""))
        if not stratum_id or stratum_id in records:
            raise ValueError(f"{path}: duplicate or empty calibration stratum")
        readiness = item.get("pilot_readiness")
        recommendation = item.get("recommended_envelope")
        if readiness not in ("ready", "ready_provisional_zero_complement"):
            continue
        if not isinstance(recommendation, dict):
            raise ValueError(
                f"{path}: ready stratum {stratum_id} lacks an envelope"
            )
        try:
            envelope = float(recommendation["sigr_max"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(
                f"{path}: {stratum_id} has an invalid envelope"
            ) from error
        if not math.isfinite(envelope) or envelope <= 0.0:
            raise ValueError(
                f"{path}: {stratum_id} envelope must be finite and positive"
            )
        if stratum_id not in report_guards:
            raise ValueError(f"{path}: {stratum_id} lacks a calibrated guard")
        records[stratum_id] = {
            **item,
            "guard": report_guards[stratum_id],
            "resolved_sigr_max": envelope,
        }
    override = {
        "enabled": bool(not revision_matches and allow_revision_mismatch),
        "calibration_generator_revisions": sorted(revisions),
        "generation_generator_revision": generator_revision,
        "rationale": (
            revision_compatibility_rationale
            if not revision_matches and allow_revision_mismatch
            else None
        ),
    }
    return payload, hashlib.sha256(raw).hexdigest(), records, override


def create_refinement(args: argparse.Namespace) -> Path:
    """Create one validated, evidence-hashed stratum refinement artifact."""
    config_path = args.config.resolve()
    recipes_path = args.recipes.resolve()
    output = args.output.resolve()
    for path in (config_path, recipes_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    (
        _config,
        recipes,
        config_sha256,
        recipes_sha256,
    ) = _load_config_and_recipes(config_path, recipes_path)
    base_padding = _candidate_padding(recipes, args.candidate)
    if args.stratum not in recipes.get("strata", {}):
        raise ValueError(f"guard recipes lack {args.stratum}")
    faces: dict[str, dict[str, float]] = {}
    for specification in args.face:
        pieces = specification.split(":")
        if len(pieces) != 3:
            raise ValueError(
                f"{specification!r}: expected AXIS:FACE:VALUE"
            )
        axis, face, raw_value = pieces
        if axis not in AXES or face not in ("lower", "upper"):
            raise ValueError(
                f"{specification!r}: invalid axis or face"
            )
        if face in faces.setdefault(axis, {}):
            raise ValueError(
                f"{args.stratum}: duplicate {axis} {face} refinement"
            )
        try:
            faces[axis][face] = float(raw_value)
        except ValueError as error:
            raise ValueError(
                f"{specification!r}: refinement value is not numeric"
            ) from error
    evidence: list[dict[str, object]] = []
    for raw_path in args.evidence:
        path = raw_path.resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        evidence.append(
            {
                "path": str(path),
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
            }
        )
    stratum_refinement = {
        "rationale": args.rationale,
        "evidence": evidence,
        "faces": faces,
    }
    original = reconstruct_guard_box(
        recipes["strata"][args.stratum], base_padding
    )
    refined, applied = apply_guard_refinement(
        original, stratum_refinement, stratum_id=args.stratum
    )
    payload: dict[str, object] = {
        "schema": REFINEMENT_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "coordinate_space": REFINEMENT_COORDINATE_SPACE,
        "analysis_config_source": str(config_path),
        "analysis_config_sha256": config_sha256,
        "guard_recipes_source": str(recipes_path),
        "guard_recipes_sha256": recipes_sha256,
        "guard_candidate": args.candidate,
        "strata": {args.stratum: stratum_refinement},
        "preview": {
            args.stratum: {
                "original_guard": original.manifest_record(),
                "refined_guard": refined.manifest_record(),
                **applied,
            }
        },
    }
    if output.exists() and not args.overwrite:
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    _write_json(output, payload)
    return output


def prepare(args: argparse.Namespace) -> Path:
    operation = (
        "calibration"
        if getattr(args, "command", "prepare") == "prepare-calibration"
        else "generation"
    )
    events_per_stratum = int(getattr(args, "events_per_stratum", 0))
    calibration_trials = int(getattr(args, "trials", 0))
    heartbeat_interval = int(getattr(args, "heartbeat_interval", 100000))
    component_core_fraction = float(
        getattr(
            args,
            "calibration_inside_guard_fraction",
            getattr(args, "calibration_core_fraction", args.core_fraction),
        )
    )
    calibration_proposal = 1 if operation == "calibration" else 0
    if operation == "generation" and events_per_stratum <= 0:
        raise ValueError("--events-per-stratum must be positive")
    if operation == "calibration" and calibration_trials <= 0:
        raise ValueError("--trials must be positive")
    if heartbeat_interval <= 0:
        raise ValueError("--heartbeat-interval must be positive")
    if not 0.0 < component_core_fraction < 1.0:
        raise ValueError(
            "--inside-guard-trial-fraction must lie strictly between zero "
            "and one"
        )
    if args.replicas <= 0:
        raise ValueError("--replicas must be positive")
    if args.seed_base <= 0:
        raise ValueError("--seed-base must be positive")
    if not 0.0 < args.core_fraction < 1.0:
        raise ValueError("--core-fraction must lie strictly between zero and one")
    raw_shared_sigr_max = getattr(args, "sigr_max", None)
    raw_envelope_path = getattr(args, "envelope_report", None)
    if operation == "generation":
        if (raw_shared_sigr_max is None) == (raw_envelope_path is None):
            raise ValueError(
                "generation requires exactly one of --sigr-max or "
                "--envelope-report"
            )
        shared_sigr_max = (
            float(raw_shared_sigr_max)
            if raw_shared_sigr_max is not None
            else None
        )
        if shared_sigr_max is not None and (
            not math.isfinite(shared_sigr_max) or shared_sigr_max <= 0.0
        ):
            raise ValueError("--sigr-max must be finite and positive")
    else:
        shared_sigr_max = 1.0
        if raw_envelope_path is not None:
            raise ValueError(
                "--envelope-report is valid only for generation"
            )

    config_path = args.config.resolve()
    recipes_path = args.recipes.resolve()
    input_path = args.input.resolve()
    for path in (config_path, recipes_path, input_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    (
        config,
        recipes,
        config_sha256,
        recipes_sha256,
        legacy_source,
    ) = _load_configuration(config_path, recipes_path, input_path)
    base_padding = _candidate_padding(recipes, args.candidate)
    refinements_path = getattr(args, "refinements", None)
    refinements: Optional[dict] = None
    refinements_sha256: Optional[str] = None
    if refinements_path is not None:
        refinements_path = refinements_path.resolve()
        if not refinements_path.is_file():
            raise FileNotFoundError(refinements_path)
        refinements, refinements_sha256 = _load_guard_refinements(
            refinements_path,
            config_sha256=config_sha256,
            recipes_sha256=recipes_sha256,
            candidate=args.candidate,
            recipes=recipes,
            base_padding=base_padding,
        )
    selection = _analysis_selection(
        config, apply_y_max=bool(args.apply_y_max)
    )
    envelope_path: Optional[Path] = None
    envelope_sha256: Optional[str] = None
    envelope_records: dict[str, dict] = {}
    envelope_revision_override: Optional[dict[str, object]] = None
    if raw_envelope_path is not None:
        envelope_path = Path(raw_envelope_path).expanduser().resolve()
        if not envelope_path.is_file():
            raise FileNotFoundError(envelope_path)
        revision_rationale = str(
            getattr(
                args,
                "envelope_revision_compatibility_rationale",
                "",
            )
            or ""
        ).strip()
        (
            _envelope_payload,
            envelope_sha256,
            envelope_records,
            envelope_revision_override,
        ) = _load_generation_envelopes(
            envelope_path,
            config_sha256=config_sha256,
            recipes_sha256=recipes_sha256,
            refinements_sha256=refinements_sha256,
            candidate=args.candidate,
            core_fraction=float(args.core_fraction),
            selection=selection,
            generator_revision=str(args.generator_revision),
            allow_revision_mismatch=bool(
                getattr(args, "allow_envelope_revision_mismatch", False)
            ),
            revision_compatibility_rationale=revision_rationale,
        )
    elif bool(getattr(args, "allow_envelope_revision_mismatch", False)):
        raise ValueError(
            "--allow-envelope-revision-mismatch requires --envelope-report"
        )
    strata = radiative_guards.enumerate_strata(config)
    selected_strata, stratum_selection = _selected_strata(strata, args)
    output = args.output.resolve()
    manifest_path = output / "manifest.json"
    if output.exists() and any(output.iterdir()) and not args.overwrite:
        raise FileExistsError(
            f"{output} is not empty; pass --overwrite to replace prepared files"
        )

    records: list[dict[str, object]] = []
    prepared_inputs: dict[Path, str] = {}
    shared_legacy_records: Optional[list[str]] = None
    used_envelope_strata: set[str] = set()
    for stratum in selected_strata:
        try:
            recipe = recipes["strata"][stratum.identifier]
        except KeyError as error:
            raise ValueError(
                f"guard recipes lack {stratum.identifier}"
            ) from error
        expected_bounds = {
            "Q2": list(stratum.q2),
            "xB": list(stratum.xb),
            "minus_t": list(stratum.minus_t),
            "phi_deg": list(stratum.phi_deg),
        }
        if recipe.get("bounds") != expected_bounds:
            raise ValueError(
                f"{stratum.identifier}: recipe bounds differ from config"
            )
        original_box = reconstruct_guard_box(recipe, base_padding)
        refinement_record: Optional[dict[str, object]] = None
        specification = (
            refinements["strata"].get(stratum.identifier)
            if refinements is not None
            else None
        )
        if specification is not None:
            box, refinement_record = apply_guard_refinement(
                original_box,
                specification,
                stratum_id=stratum.identifier,
            )
        else:
            box = original_box
        envelope_record = envelope_records.get(stratum.identifier)
        if operation == "generation" and envelope_path is not None:
            if envelope_record is None:
                raise ValueError(
                    f"{envelope_path}: selected stratum {stratum.identifier} "
                    "is absent or not pilot-ready"
                )
            if int(envelope_record["flat_index"]) != stratum.flat_index:
                raise ValueError(
                    f"{envelope_path}: {stratum.identifier} flat index differs"
                )
            if envelope_record.get("bounds") != expected_bounds:
                raise ValueError(
                    f"{envelope_path}: {stratum.identifier} bounds differ"
                )
            if envelope_record.get("guard") != box.manifest_record():
                raise ValueError(
                    f"{envelope_path}: {stratum.identifier} guard differs"
                )
            stratum_sigr_max = float(
                envelope_record["resolved_sigr_max"]
            )
            used_envelope_strata.add(stratum.identifier)
        else:
            stratum_sigr_max = float(shared_sigr_max)
        if operation == "calibration" and box.volume >= 1.0 - 1.0e-7:
            raise ValueError(
                f"{stratum.identifier}: the anchored guard fills the native "
                "hypercube, so its complement cannot be calibrated"
            )
        legacy_input, legacy_records = _legacy_mode4_input(
            legacy_source,
            input_path,
            events=(events_per_stratum if operation == "generation" else 1),
            sigr_max=stratum_sigr_max,
        )
        if envelope_path is None:
            shared_legacy_records = legacy_records
        for replica_index in range(args.replicas):
            seed = (
                args.seed_base + 1000 * stratum.flat_index + replica_index
            )
            generation_id = f"g{replica_index:04d}"
            stem = f"{stratum.identifier}__{generation_id}"
            input_relative = Path("inputs") / f"{stem}.inp"
            output_stem = Path("runs") / stratum.identifier / stem
            text = legacy_input + _mode4_trailer(
                seed=seed,
                replica=replica_index,
                stratum=stratum,
                w_minimum=float(config["phase_space"]["W_min"]),
                apply_y_max=args.apply_y_max,
                y_maximum=float(config["phase_space"].get("y_max", 1.0)),
                core_fraction=args.core_fraction,
                base_padding=base_padding,
                box=box,
                operation=operation,
                trial_limit=(
                    calibration_trials if operation == "calibration" else 0
                ),
                heartbeat_interval=heartbeat_interval,
                component_core_fraction=component_core_fraction,
                calibration_proposal=calibration_proposal,
            )
            prepared_inputs[input_relative] = text
            records.append(
                {
                    "stratum_id": stratum.identifier,
                    "flat_index": stratum.flat_index,
                    "indices": {
                        "iq2": stratum.iq2,
                        "ixb": stratum.ixb,
                        "it": stratum.it,
                        "iphi": stratum.iphi,
                    },
                    "bounds": expected_bounds,
                    "training_status": recipe["training_status"],
                    "fit_source": recipe["fit_source"],
                    "generation_id": generation_id,
                    "replica_index": replica_index,
                    "seed": seed,
                    "operation": operation,
                    "sigr_max": (
                        stratum_sigr_max
                        if operation == "generation"
                        else None
                    ),
                    "envelope_recommendation_status": (
                        envelope_record.get("recommendation_status")
                        if envelope_record is not None
                        else None
                    ),
                    "envelope_pilot_readiness": (
                        envelope_record.get("pilot_readiness")
                        if envelope_record is not None
                        else None
                    ),
                    "envelope_recommendation_basis": (
                        envelope_record.get("recommendation_basis")
                        if envelope_record is not None
                        else None
                    ),
                    "events_requested": (
                        events_per_stratum if operation == "generation" else 0
                    ),
                    "trials_requested": (
                        calibration_trials if operation == "calibration" else 0
                    ),
                    "input_file": str(input_relative),
                    "input_sha256": hashlib.sha256(
                        text.encode("utf-8")
                    ).hexdigest(),
                    "legacy_input_records": legacy_records,
                    "output_stem": str(output_stem),
                    "guard_original": original_box.manifest_record(),
                    "guard": box.manifest_record(),
                    "guard_refinement": refinement_record,
                }
            )
    if not records:
        raise ValueError("selected strata produced no runs")
    manifest: dict[str, object] = {
        "schema": MANIFEST_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "development_only": True,
        "production_ready": False,
        "sampling_mode": SAMPLING_MODE,
        "operation": operation,
        "generator_revision": args.generator_revision,
        "stratum_selection": stratum_selection,
        "analysis_config_source": str(config_path),
        "analysis_config_sha256": config_sha256,
        "guard_recipes_source": str(recipes_path),
        "guard_recipes_sha256": recipes_sha256,
        "guard_refinements_source": (
            str(refinements_path) if refinements_path is not None else None
        ),
        "guard_refinements_sha256": refinements_sha256,
        "guard_refinement_schema": (
            REFINEMENT_SCHEMA if refinements is not None else None
        ),
        "guard_refined_strata": (
            sorted(refinements["strata"]) if refinements is not None else []
        ),
        "guard_learner_revision": recipes.get(
            "continuous_guard_learner_revision"
        ),
        "guard_training_generator_revision": recipes.get(
            "generator_revision"
        ),
        "guard_candidate": args.candidate,
        "base_padding": base_padding,
        "core_fraction": args.core_fraction,
        "legacy_tail_fraction": 1.0 - args.core_fraction,
        "full_support_guaranteed": True,
        "events_per_stratum": (
            events_per_stratum if operation == "generation" else 0
        ),
        "calibration_trials_per_replica": (
            calibration_trials if operation == "calibration" else 0
        ),
        "calibration_proposal": (
            "guard_partition" if operation == "calibration" else None
        ),
        "calibration_component_mapping": (
            {
                "1": "inside_guard",
                "0": "guard_complement",
            }
            if operation == "calibration"
            else None
        ),
        "component_core_fraction": component_core_fraction,
        "calibration_inside_guard_fraction": (
            component_core_fraction if operation == "calibration" else None
        ),
        "heartbeat_interval": heartbeat_interval,
        "replicas_per_stratum": args.replicas,
        "envelope_mode": (
            "per_stratum_calibration"
            if envelope_path is not None
            else "shared_scalar"
            if operation == "generation"
            else None
        ),
        "sigr_max": (
            shared_sigr_max
            if operation == "generation" and envelope_path is None
            else None
        ),
        "envelope_calibration_source": (
            str(envelope_path) if envelope_path is not None else None
        ),
        "envelope_calibration_sha256": envelope_sha256,
        "envelope_calibration_schema": (
            CALIBRATION_SCHEMA if envelope_path is not None else None
        ),
        "envelope_calibration_snapshot": (
            "envelope_calibration.json" if envelope_path is not None else None
        ),
        "envelope_revision_compatibility_override": (
            envelope_revision_override
        ),
        "envelope_strata_used": sorted(used_envelope_strata),
        "analysis_selection": selection,
        "legacy_input_records": shared_legacy_records,
        "runs": records,
    }
    input_directory = output / "inputs"
    input_directory.mkdir(parents=True, exist_ok=True)
    shutil.copy2(config_path, output / "analysis_config.json")
    shutil.copy2(recipes_path, output / "continuous_guard_recipes.json")
    if stratum_selection.get("mode") == "flat_index_file":
        shutil.copy2(
            Path(str(stratum_selection["source"])),
            output / str(stratum_selection["snapshot"]),
        )
    if refinements_path is not None:
        shutil.copy2(
            refinements_path, output / "guard_refinements.json"
        )
    if envelope_path is not None:
        shutil.copy2(
            envelope_path, output / "envelope_calibration.json"
        )
    (output / "legacy_input.inp").write_text(
        legacy_source, encoding="utf-8"
    )
    for relative, text in prepared_inputs.items():
        (output / relative).write_text(text, encoding="utf-8")
    _write_json(manifest_path, manifest)
    return manifest_path


def _norm_int(norm: dict[str, str], name: str) -> int:
    try:
        return int(norm[name])
    except (KeyError, ValueError) as error:
        raise Mode4Error(f"normalization lacks integer {name}") from error


def _norm_float(norm: dict[str, str], name: str) -> float:
    try:
        value = float(norm[name])
    except (KeyError, ValueError) as error:
        raise Mode4Error(f"normalization lacks float {name}") from error
    if not math.isfinite(value):
        raise Mode4Error(f"normalization {name} is not finite")
    return value


def _close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=3.0e-6, abs_tol=3.0e-8)


def _validate_norm(
    norm: dict[str, str],
    manifest: dict,
    record: dict,
) -> None:
    expected_ints = {
        "sampling_mode": SAMPLING_MODE,
        "fixed_trial_survey": 0,
        "mode4_replica": int(record["replica_index"]),
        "mode4_seed": int(record["seed"]),
        "stratum_flat_index": int(record["flat_index"]),
        "stratum_iq2": int(record["indices"]["iq2"]),
        "stratum_ixb": int(record["indices"]["ixb"]),
        "stratum_it": int(record["indices"]["it"]),
        "stratum_iphi": int(record["indices"]["iphi"]),
        "mode4_apply_y_max": int(
            manifest["analysis_selection"]["apply_y_max"]
        ),
        "mode4_operation": int(manifest["operation"] == "calibration"),
        "mode4_trial_limit": int(record["trials_requested"]),
        "mode4_heartbeat_interval": int(manifest["heartbeat_interval"]),
        "mode4_calibration_proposal": int(
            manifest["operation"] == "calibration"
        ),
    }
    for name, expected in expected_ints.items():
        actual = _norm_int(norm, name)
        if actual != expected:
            raise Mode4Error(f"{name}={actual}, expected {expected}")
    expected_floats = {
        "mode4_core_fraction": float(manifest["core_fraction"]),
        "mode4_legacy_fraction": float(manifest["legacy_tail_fraction"]),
        "mode4_base_padding": float(manifest["base_padding"]),
        "mode4_core_volume": float(record["guard"]["normalized_volume"]),
        "mode4_w_minimum": float(
            manifest["analysis_selection"]["w_minimum"]
        ),
        "mode4_component_core_fraction": float(
            manifest["component_core_fraction"]
        ),
    }
    for name, expected in expected_floats.items():
        actual = _norm_float(norm, name)
        if not _close(actual, expected):
            raise Mode4Error(f"{name}={actual}, expected {expected}")
    if _norm_int(norm, "ntries") <= 0:
        raise Mode4Error("mode 4 produced no proposals")
    if manifest["operation"] == "generation":
        if _norm_int(norm, "nevent") <= 0:
            raise Mode4Error("mode-4 generation produced no events")
        if _norm_float(norm, "sig_sum") <= 0.0:
            raise Mode4Error("mode 4 produced a nonpositive stratum integral")
    else:
        if _norm_int(norm, "nevent") != 0:
            raise Mode4Error("mode-4 calibration unexpectedly emitted events")
        if _norm_int(norm, "ntries") != int(record["trials_requested"]):
            raise Mode4Error("mode-4 calibration stopped at the wrong trial count")
    if _norm_int(norm, "mode4_core_trials") + _norm_int(
        norm, "mode4_legacy_trials"
    ) != _norm_int(norm, "ntries"):
        raise Mode4Error("mode-4 component trial counts do not sum to ntries")
    if _norm_int(norm, "mode4_core_targets") + _norm_int(
        norm, "mode4_legacy_targets"
    ) != _norm_int(norm, "mode4_target_candidates"):
        raise Mode4Error("mode-4 component target counts are inconsistent")
    if _norm_int(norm, "mode4_core_events") + _norm_int(
        norm, "mode4_legacy_events"
    ) != _norm_int(norm, "nevent"):
        raise Mode4Error("mode-4 component event counts are inconsistent")
    if _norm_int(norm, "mode4_emitting_candidates") + _norm_int(
        norm, "mode4_duplicate_events"
    ) != _norm_int(norm, "nevent"):
        raise Mode4Error(
            "mode-4 emitting and duplicate counts do not sum to nevent"
        )
    if manifest["operation"] == "calibration" and (
        _norm_int(norm, "mode4_emitting_candidates") != 0
        or _norm_int(norm, "mode4_duplicate_events") != 0
    ):
        raise Mode4Error(
            "mode-4 calibration unexpectedly recorded emitted events"
        )


def _inside(value: float, bounds: list[float], tolerance: float) -> bool:
    return bounds[0] - tolerance <= value <= bounds[1] + tolerance


def _validate_event_diagnostics(
    path: Path,
    manifest: dict,
    record: dict,
    expected_events: int,
) -> dict[str, object]:
    with path.open(encoding="utf-8", newline="") as source:
        first = source.readline().strip()
        if first != f"# schema={MODE4_KINEMATICS_SCHEMA}":
            raise Mode4Error(f"{path}: unexpected event schema")
        reader = csv.DictReader(source)
        if tuple(reader.fieldnames or ()) != MODE4_KINEMATICS_COLUMNS:
            raise Mode4Error(f"{path}: unexpected event columns")
        rows = list(reader)
    if len(rows) != expected_events:
        raise Mode4Error(
            f"{path}: {len(rows)} rows for {expected_events} events"
        )
    bounds = record["bounds"]
    core_events = legacy_events = 0
    maximum_ratio = 0.0
    guard_box = _guard_box_from_manifest(record)
    for index, row in enumerate(rows, start=1):
        if int(row["event"]) != index:
            raise Mode4Error(f"{path}: event numbers are not sequential")
        component = int(row["proposal_component"])
        if component == 1:
            core_events += 1
        elif component == 0:
            legacy_events += 1
        else:
            raise Mode4Error(f"{path}: invalid proposal component")
        values = {
            name: float(row[name])
            for name in MODE4_KINEMATICS_COLUMNS[2:]
        }
        if any(not math.isfinite(value) for value in values.values()):
            raise Mode4Error(f"{path}: nonfinite event diagnostic")
        if not _inside(values["q2_observed"], bounds["Q2"], 4.0e-6):
            raise Mode4Error(f"{path}: Q2 outside requested stratum")
        if not _inside(values["xb_observed"], bounds["xB"], 4.0e-7):
            raise Mode4Error(f"{path}: xB outside requested stratum")
        if not _inside(
            values["minus_t_observed"], bounds["minus_t"], 4.0e-6
        ):
            raise Mode4Error(f"{path}: -t outside requested stratum")
        phi = values["phi_observed_deg"] % 360.0
        if not _inside(phi, bounds["phi_deg"], 4.0e-4):
            raise Mode4Error(f"{path}: phi outside requested stratum")
        if values["w_observed"] + 2.0e-6 < float(
            manifest["analysis_selection"]["w_minimum"]
        ):
            raise Mode4Error(f"{path}: W below requested selection")
        if manifest["analysis_selection"]["apply_y_max"] and values[
            "y_observed"
        ] > float(manifest["analysis_selection"]["y_maximum"]) + 1.0e-7:
            raise Mode4Error(f"{path}: y above requested selection")
        if values["proposal_density_ratio"] <= 0.0:
            raise Mode4Error(f"{path}: nonpositive density correction")
        proposal_coordinates = {
            name: values[name] for name in AXES
        }
        generator_inside = _fortran_guard_contains(
            guard_box,
            proposal_coordinates,
        )
        expected_ratio = _proposal_density_ratio_for_membership(
            guard_box,
            float(manifest["core_fraction"]),
            generator_inside,
        )
        if not math.isclose(
            values["proposal_density_ratio"],
            expected_ratio,
            rel_tol=3.0e-5,
            abs_tol=3.0e-7,
        ):
            raise Mode4Error(
                f"{path}: proposal density correction differs from manifest"
            )
        maximum_ratio = max(
            maximum_ratio, values["proposal_density_ratio"]
        )
    return {
        "core_events": core_events,
        "legacy_events": legacy_events,
        "maximum_event_density_ratio": maximum_ratio,
    }


def _validate_heartbeat(
    path: Path,
    norm: dict[str, str],
) -> dict[str, int]:
    with path.open(encoding="utf-8", newline="") as source:
        first = source.readline().strip()
        if first != f"# schema={MODE4_HEARTBEAT_SCHEMA}":
            raise Mode4Error(f"{path}: unexpected heartbeat schema")
        reader = csv.DictReader(source)
        if tuple(reader.fieldnames or ()) != MODE4_HEARTBEAT_COLUMNS:
            raise Mode4Error(f"{path}: unexpected heartbeat columns")
        rows = [{name: int(row[name]) for name in MODE4_HEARTBEAT_COLUMNS}
                for row in reader]
    if not rows:
        raise Mode4Error(f"{path}: heartbeat contains no snapshots")
    previous = -1
    for row in rows:
        if row["proposals"] < previous:
            raise Mode4Error(f"{path}: heartbeat proposals are not monotonic")
        previous = row["proposals"]
    final = rows[-1]
    expected = {
        "proposals": _norm_int(norm, "ntries"),
        "events": _norm_int(norm, "nevent"),
        "core_trials": _norm_int(norm, "mode4_core_trials"),
        "noncore_trials": _norm_int(norm, "mode4_legacy_trials"),
        "internal_rows": _norm_int(norm, "mode4_internal_rows"),
        "final_candidates": _norm_int(norm, "mode4_final_candidates"),
        "target_candidates": _norm_int(norm, "mode4_target_candidates"),
        "core_targets": _norm_int(norm, "mode4_core_targets"),
        "noncore_targets": _norm_int(norm, "mode4_legacy_targets"),
        "core_events": _norm_int(norm, "mode4_core_events"),
        "noncore_events": _norm_int(norm, "mode4_legacy_events"),
        "emitting_candidates": _norm_int(
            norm, "mode4_emitting_candidates"
        ),
        "duplicate_events": _norm_int(norm, "mode4_duplicate_events"),
        "mcall_max": _norm_int(norm, "mcall_max"),
    }
    if final != expected:
        raise Mode4Error(f"{path}: final heartbeat differs from normalization")
    return final


def _validate_calibration_diagnostics(
    path: Path,
    manifest: dict,
    record: dict,
    norm: dict[str, str],
) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", newline="") as source:
        first = source.readline().strip()
        if first != f"# schema={MODE4_CALIBRATION_EVENT_SCHEMA}":
            raise Mode4Error(f"{path}: unexpected calibration schema")
        reader = csv.DictReader(source)
        if tuple(reader.fieldnames or ()) != MODE4_CALIBRATION_COLUMNS:
            raise Mode4Error(f"{path}: unexpected calibration columns")
        raw_rows = list(reader)
    expected_rows = _norm_int(norm, "mode4_target_candidates")
    if len(raw_rows) != expected_rows:
        raise Mode4Error(
            f"{path}: {len(raw_rows)} target rows for {expected_rows} candidates"
        )
    alpha = float(manifest["core_fraction"])
    beta = float(manifest["component_core_fraction"])
    box = _guard_box_from_manifest(record)
    inside_mass = alpha + (1.0 - alpha) * box.volume
    complement_mass = (1.0 - alpha) * (1.0 - box.volume)
    bounds = record["bounds"]
    rows: list[dict[str, object]] = []
    previous_trial = 0
    for raw in raw_rows:
        trial = int(raw["trial"])
        component = int(raw["proposal_component"])
        inside = int(raw["inside_core"])
        if trial <= previous_trial or trial > _norm_int(norm, "ntries"):
            raise Mode4Error(f"{path}: invalid calibration trial ordering")
        previous_trial = trial
        if component not in (0, 1) or inside not in (0, 1):
            raise Mode4Error(f"{path}: invalid component or core indicator")
        values = {
            name: float(raw[name])
            for name in MODE4_CALIBRATION_COLUMNS[3:]
        }
        if any(not math.isfinite(value) for value in values.values()):
            raise Mode4Error(f"{path}: nonfinite calibration diagnostic")
        if values["integrand_corrected"] <= 0.0:
            raise Mode4Error(f"{path}: nonpositive target integrand")
        coordinates = {name: values[name] for name in AXES}
        expected_inside = int(_fortran_guard_contains(box, coordinates))
        if (
            inside != expected_inside
            or (component == 1 and inside != 1)
            or (component == 0 and inside != 0)
        ):
            raise Mode4Error(f"{path}: inconsistent learned-core membership")
        expected_ratio = _proposal_density_ratio_for_membership(
            box, alpha, bool(expected_inside)
        )
        if not math.isclose(
            values["proposal_density_ratio"],
            expected_ratio,
            rel_tol=3.0e-5,
            abs_tol=3.0e-7,
        ):
            raise Mode4Error(f"{path}: inconsistent proposal correction")
        expected_weight = inside_mass / beta if component == 1 else (
            complement_mass / (1.0 - beta)
        )
        if not _close(
            values["component_importance_weight"], expected_weight
        ):
            raise Mode4Error(f"{path}: inconsistent component weight")
        if not _inside(values["q2_observed"], bounds["Q2"], 4.0e-6):
            raise Mode4Error(f"{path}: Q2 outside requested stratum")
        if not _inside(values["xb_observed"], bounds["xB"], 4.0e-7):
            raise Mode4Error(f"{path}: xB outside requested stratum")
        if not _inside(
            values["minus_t_observed"], bounds["minus_t"], 4.0e-6
        ):
            raise Mode4Error(f"{path}: -t outside requested stratum")
        if not _inside(
            values["phi_observed_deg"] % 360.0,
            bounds["phi_deg"],
            4.0e-4,
        ):
            raise Mode4Error(f"{path}: phi outside requested stratum")
        if values["w_observed"] + 2.0e-6 < float(
            manifest["analysis_selection"]["w_minimum"]
        ):
            raise Mode4Error(f"{path}: W below requested selection")
        if manifest["analysis_selection"]["apply_y_max"] and values[
            "y_observed"
        ] > float(manifest["analysis_selection"]["y_maximum"]) + 1.0e-7:
            raise Mode4Error(f"{path}: y above requested selection")
        rows.append(
            {
                "trial": trial,
                "proposal_component": component,
                "inside_core": inside,
                **values,
            }
        )
    core_rows = sum(int(row["proposal_component"] == 1) for row in rows)
    if core_rows != _norm_int(norm, "mode4_core_targets"):
        raise Mode4Error(f"{path}: core target count differs from normalization")
    if len(rows) - core_rows != _norm_int(norm, "mode4_legacy_targets"):
        raise Mode4Error(f"{path}: tail target count differs from normalization")
    reconstructed_sigma = (
        _norm_float(norm, "mode4_phase_volume")
        / _norm_int(norm, "ntries")
        * sum(
            float(row["integrand_corrected"])
            * float(row["component_importance_weight"])
            for row in rows
        )
    )
    if not math.isclose(
        reconstructed_sigma,
        _norm_float(norm, "sig_sum"),
        rel_tol=5.0e-5,
        abs_tol=1.0e-12,
    ):
        raise Mode4Error(
            f"{path}: calibration rows do not reproduce sig_sum"
        )
    return rows


def _validate_lund(path: Path, events: int) -> None:
    lines = [
        line.split()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(lines) != 5 * events:
        raise Mode4Error(
            f"{path}: {len(lines)} nonempty lines for {events} events"
        )
    for offset in range(0, len(lines), 5):
        if lines[offset][0] != "4":
            raise Mode4Error(f"{path}: LUND header does not declare 4 particles")
        for particle, row in enumerate(lines[offset + 1 : offset + 5], start=1):
            if int(row[0]) != particle:
                raise Mode4Error(f"{path}: malformed LUND particle order")


def run(args: argparse.Namespace) -> Path:
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not _supported_manifest_schema(manifest):
        raise ValueError(f"unsupported manifest schema: {manifest.get('schema')}")
    matches = [
        record
        for record in manifest["runs"]
        if int(record["flat_index"]) == args.flat_index
        and int(record["replica_index"]) == args.replica_index
    ]
    if len(matches) != 1:
        raise ValueError(
            f"manifest has {len(matches)} matching mode-4 run records"
        )
    record = matches[0]
    root = manifest_path.parent
    stratum_selection = manifest.get("stratum_selection") or {}
    selection_sha256 = stratum_selection.get("source_sha256")
    if selection_sha256 is not None:
        selection_snapshot = root / str(
            stratum_selection.get("snapshot", "flat_index_selection.txt")
        )
        if (
            not selection_snapshot.is_file()
            or _sha256(selection_snapshot) != selection_sha256
        ):
            raise Mode4Error(
                "flat-index selection snapshot is missing or differs from "
                "the manifest"
            )
    run_sigr_max = (
        _run_sigr_max(manifest, record)
        if manifest["operation"] == "generation"
        else None
    )
    refinements_sha256 = manifest.get("guard_refinements_sha256")
    if refinements_sha256 is not None:
        refinement_snapshot = root / "guard_refinements.json"
        if (
            not refinement_snapshot.is_file()
            or _sha256(refinement_snapshot) != refinements_sha256
        ):
            raise Mode4Error(
                "guard-refinement snapshot is missing or differs from "
                "the manifest"
            )
    envelope_sha256 = manifest.get("envelope_calibration_sha256")
    if envelope_sha256 is not None:
        envelope_snapshot = root / str(
            manifest.get(
                "envelope_calibration_snapshot",
                "envelope_calibration.json",
            )
        )
        if (
            not envelope_snapshot.is_file()
            or _sha256(envelope_snapshot) != envelope_sha256
        ):
            raise Mode4Error(
                "envelope-calibration snapshot is missing or differs from "
                "the manifest"
            )
    input_path = root / record["input_file"]
    if record.get("input_sha256") is not None and _sha256(input_path) != record[
        "input_sha256"
    ]:
        raise Mode4Error("prepared mode-4 input differs from the manifest")
    output_stem = root / record["output_stem"]
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    run_record_path = output_stem.with_suffix(".json")
    if run_record_path.exists() and not args.overwrite:
        raise FileExistsError(run_record_path)
    executable = args.executable.resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise FileNotFoundError(executable)
    with tempfile.TemporaryDirectory(
        prefix=f"{record['stratum_id']}_mode4_"
    ) as temporary:
        work = Path(temporary)
        print(f"mode4_work_directory={work}", flush=True)
        print(
            f"mode4_live_heartbeat={work / MODE4_HEARTBEAT_FILENAME}",
            flush=True,
        )
        completed = subprocess.run(
            [str(executable)],
            input=input_path.read_text(encoding="utf-8"),
            text=True,
            cwd=work,
            capture_output=True,
            check=False,
        )
        (work / "generator.stdout.txt").write_text(
            completed.stdout, encoding="utf-8"
        )
        (work / "generator.stderr.txt").write_text(
            completed.stderr, encoding="utf-8"
        )
        if completed.returncode != 0:
            raise Mode4Error(
                f"AAO failed with exit code {completed.returncode}: "
                f"{completed.stderr[-2000:]}"
            )
        norm_path = work / radiative_survey.NORM_FILENAME
        lund_path = work / radiative_survey.LUND_FILENAME
        heartbeat_path = work / MODE4_HEARTBEAT_FILENAME
        required = [norm_path, lund_path, heartbeat_path]
        if manifest["operation"] == "generation":
            diagnostic_path = work / MODE4_KINEMATICS_FILENAME
        else:
            diagnostic_path = work / MODE4_CALIBRATION_FILENAME
        required.append(diagnostic_path)
        for path in required:
            if not path.is_file():
                raise Mode4Error(f"generator did not write {path.name}")
        norm = radiative_survey.parse_norm(norm_path)
        _validate_norm(norm, manifest, record)
        heartbeat_summary = _validate_heartbeat(heartbeat_path, norm)
        events = _norm_int(norm, "nevent")
        calibration_rows: list[dict[str, object]] = []
        if manifest["operation"] == "generation":
            event_summary = _validate_event_diagnostics(
                diagnostic_path, manifest, record, events
            )
            _validate_lund(lund_path, events)
            if event_summary["core_events"] != _norm_int(
                norm, "mode4_core_events"
            ) or event_summary["legacy_events"] != _norm_int(
                norm, "mode4_legacy_events"
            ):
                raise Mode4Error(
                    "event diagnostics and norm component counts differ"
                )
        else:
            event_summary = {
                "core_events": 0,
                "legacy_events": 0,
                "maximum_event_density_ratio": None,
            }
            if lund_path.stat().st_size != 0:
                raise Mode4Error("mode-4 calibration emitted LUND content")
            calibration_rows = _validate_calibration_diagnostics(
                diagnostic_path, manifest, record, norm
            )
        with norm_path.open("a", encoding="utf-8") as norm_output:
            norm_output.write(
                f"mode4_manifest_sha256={_sha256(manifest_path)}\n"
                f"mode4_analysis_config_sha256="
                f"{manifest['analysis_config_sha256']}\n"
                f"mode4_guard_recipes_sha256="
                f"{manifest['guard_recipes_sha256']}\n"
                f"mode4_guard_candidate={manifest['guard_candidate']}\n"
                f"mode4_stratum_id={record['stratum_id']}\n"
                f"mode4_generator_revision={manifest['generator_revision']}\n"
            )
            if manifest.get("guard_refinements_sha256") is not None:
                norm_output.write(
                    "mode4_guard_refinements_sha256="
                    f"{manifest['guard_refinements_sha256']}\n"
                )
            if run_sigr_max is not None:
                norm_output.write(f"mode4_sigr_max={run_sigr_max:.17g}\n")
            if envelope_sha256 is not None:
                norm_output.write(
                    "mode4_envelope_calibration_sha256="
                    f"{envelope_sha256}\n"
                )
        products = [
            (".norm", radiative_survey.NORM_FILENAME),
            (".heartbeat.csv", MODE4_HEARTBEAT_FILENAME),
            (".sum", "aao_rad.sum"),
            (".out", "aao_rad.out"),
            (".stdout", "generator.stdout.txt"),
            (".stderr", "generator.stderr.txt"),
        ]
        if manifest["operation"] == "generation":
            products.extend(
                [
                    (".lund", radiative_survey.LUND_FILENAME),
                    (".mode4.csv", MODE4_KINEMATICS_FILENAME),
                ]
            )
        else:
            products.append(
                (".calibration.csv", MODE4_CALIBRATION_FILENAME)
            )
        for suffix, source_name in products:
            source = work / source_name
            if source.exists():
                destination = Path(str(output_stem) + suffix)
                if destination.exists():
                    destination.unlink()
                shutil.move(str(source), destination)
        shutil.copy2(input_path, Path(str(output_stem) + ".inp"))
    sig_sum = _norm_float(norm, "sig_sum")
    run_record: dict[str, object] = {
        **record,
        "schema": RUN_SCHEMA,
        "source_manifest": str(manifest_path),
        "source_manifest_sha256": _sha256(manifest_path),
        "sigr_max": run_sigr_max,
        "sig_sum_microbarn": sig_sum,
        "sig_int_microbarn": _norm_float(norm, "sig_int"),
        "operation": manifest["operation"],
        "events": events,
        "event_overshoot": events - int(record["events_requested"]),
        "ntries": _norm_int(norm, "ntries"),
        "mcall_max": _norm_int(norm, "mcall_max"),
        "multiplicity_correction_used": _norm_int(norm, "mcall_max") > 1,
        "emitting_candidates": _norm_int(
            norm, "mode4_emitting_candidates"
        ),
        "duplicate_events": _norm_int(norm, "mode4_duplicate_events"),
        "duplicate_event_fraction": (
            _norm_int(norm, "mode4_duplicate_events") / events
            if events
            else 0.0
        ),
        "event_weight_microbarn": (
            sig_sum / events if events else None
        ),
        "event_yield_per_proposal": events / _norm_int(norm, "ntries"),
        "core_trials": _norm_int(norm, "mode4_core_trials"),
        "noncore_trials": _norm_int(norm, "mode4_legacy_trials"),
        "legacy_trials": _norm_int(norm, "mode4_legacy_trials"),
        "target_candidates": _norm_int(norm, "mode4_target_candidates"),
        "core_targets": _norm_int(norm, "mode4_core_targets"),
        "noncore_targets": _norm_int(norm, "mode4_legacy_targets"),
        "legacy_targets": _norm_int(norm, "mode4_legacy_targets"),
        "core_events": event_summary["core_events"],
        "noncore_events": event_summary["legacy_events"],
        "legacy_events": event_summary["legacy_events"],
        "noncore_component": (
            "guard_complement"
            if manifest["operation"] == "calibration"
            else "legacy"
        ),
        "maximum_event_density_ratio": event_summary[
            "maximum_event_density_ratio"
        ],
        "final_heartbeat": heartbeat_summary,
    }
    if manifest["operation"] == "calibration":
        run_record["calibration_target_rows"] = len(calibration_rows)
        run_record["calibration_component_core_fraction"] = manifest[
            "component_core_fraction"
        ]
        run_record["calibration_proposal"] = manifest[
            "calibration_proposal"
        ]
    _write_json(run_record_path, run_record)
    return run_record_path


def _quantile(values: list[float], probability: float) -> Optional[float]:
    if not values:
        return None
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def _component_summary(
    values: list[float],
    trials: int,
    production_probability_mass: float,
    phase_volume: float,
) -> dict[str, object]:
    if trials <= 0:
        raise Mode4Error("calibration component has no proposals")
    total = sum(values)
    total_squared = sum(value * value for value in values)
    mean = total / trials
    variance_of_mean = 0.0
    if trials > 1:
        variance = max(
            0.0,
            (total_squared - total * total / trials) / (trials - 1),
        )
        variance_of_mean = variance / trials
    integral = phase_volume * production_probability_mass * mean
    sem = (
        phase_volume
        * production_probability_mass
        * math.sqrt(variance_of_mean)
    )
    ess = total * total / total_squared if total_squared > 0.0 else 0.0
    return {
        "trials": trials,
        "production_probability_mass": production_probability_mass,
        "target_candidates": len(values),
        "target_candidate_rate": len(values) / trials,
        "sum_corrected_integrand": total,
        "integrated_cross_section_microbarn": integral,
        "integrated_cross_section_sem_microbarn": sem,
        "effective_sample_size": ess,
        "positive_integrand_quantiles": {
            "p50": _quantile(values, 0.50),
            "p90": _quantile(values, 0.90),
            "p99": _quantile(values, 0.99),
            "p999": _quantile(values, 0.999),
            "maximum": max(values) if values else None,
        },
    }


def _envelope_evaluation(
    envelope: float,
    inside_values: list[float],
    complement_values: list[float],
    inside_trials: int,
    complement_trials: int,
    inside_probability_mass: float,
    complement_probability_mass: float,
) -> dict[str, float]:
    def component(values: list[float], trials: int) -> tuple[float, float, float]:
        ratios = [value / envelope for value in values]
        emitted = sum(ratios) / trials
        emitting = sum(min(1.0, ratio) for ratio in ratios) / trials
        duplicates = sum(max(0.0, ratio - 1.0) for ratio in ratios) / trials
        return emitted, emitting, duplicates

    inside_emitted, inside_emitting, inside_duplicates = component(
        inside_values, inside_trials
    )
    complement_emitted, complement_emitting, complement_duplicates = component(
        complement_values, complement_trials
    )
    emitted = (
        inside_probability_mass * inside_emitted
        + complement_probability_mass * complement_emitted
    )
    emitting = (
        inside_probability_mass * inside_emitting
        + complement_probability_mass * complement_emitting
    )
    duplicates = (
        inside_probability_mass * inside_duplicates
        + complement_probability_mass * complement_duplicates
    )
    maximum = max(inside_values + complement_values, default=0.0)
    return {
        "sigr_max": envelope,
        "expected_events_per_proposal": emitted,
        "expected_emitting_proposals_per_proposal": emitting,
        "expected_duplicate_events_per_proposal": duplicates,
        "expected_duplicate_event_fraction": (
            duplicates / emitted if emitted > 0.0 else 0.0
        ),
        "expected_events_per_emitting_proposal": (
            emitted / emitting if emitting > 0.0 else 0.0
        ),
        "observed_maximum_mcall_ratio": (
            maximum / envelope if envelope > 0.0 else math.inf
        ),
        "expected_inside_guard_event_fraction": (
            inside_probability_mass * inside_emitted / emitted
            if emitted > 0.0
            else 0.0
        ),
        "expected_guard_complement_event_fraction": (
            complement_probability_mass * complement_emitted / emitted
            if emitted > 0.0
            else 0.0
        ),
    }


def _calibration_rows(path: Path) -> list[tuple[int, int, float]]:
    with path.open(encoding="utf-8", newline="") as source:
        schema = source.readline().strip()
        if schema != f"# schema={MODE4_CALIBRATION_EVENT_SCHEMA}":
            raise Mode4Error(f"{path}: unexpected calibration schema")
        reader = csv.DictReader(source)
        if tuple(reader.fieldnames or ()) != MODE4_CALIBRATION_COLUMNS:
            raise Mode4Error(f"{path}: unexpected calibration columns")
        return [
            (
                int(row["proposal_component"]),
                int(row["inside_core"]),
                float(row["integrand_corrected"]),
            )
            for row in reader
        ]


def _additional_inside_pilot_observation(path: Path) -> dict[str, object]:
    run_path = path.resolve()
    run = json.loads(run_path.read_text(encoding="utf-8"))
    if run.get("schema") != RUN_SCHEMA:
        raise Mode4Error(f"{run_path}: unsupported mode-4 run schema")
    if run.get("operation") != "generation":
        raise Mode4Error(f"{run_path}: expected a generation pilot")
    if int(run.get("noncore_events", -1)) != 0:
        raise Mode4Error(
            f"{run_path}: additional envelope evidence must contain only "
            "inside-guard events"
        )
    event_path = run_path.with_suffix(".mode4.csv")
    with event_path.open(encoding="utf-8", newline="") as source:
        schema = source.readline().strip()
        if schema != f"# schema={MODE4_KINEMATICS_SCHEMA}":
            raise Mode4Error(f"{event_path}: unexpected event schema")
        reader = csv.DictReader(source)
        if tuple(reader.fieldnames or ()) != MODE4_KINEMATICS_COLUMNS:
            raise Mode4Error(f"{event_path}: unexpected event columns")
        rows = list(reader)
    if len(rows) != int(run["events"]):
        raise Mode4Error(
            f"{event_path}: event count differs from {run_path}"
        )
    if not rows:
        raise Mode4Error(f"{event_path}: no pilot events")
    if any(int(row["proposal_component"]) != 1 for row in rows):
        raise Mode4Error(
            f"{event_path}: additional envelope evidence contains a "
            "noncore proposal component"
        )
    values = [float(row["integrand_corrected"]) for row in rows]
    if any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise Mode4Error(f"{event_path}: invalid corrected integrand")
    return {
        "stratum_id": str(run["stratum_id"]),
        "flat_index": int(run["flat_index"]),
        "indices": run["indices"],
        "bounds": run["bounds"],
        "guard": run["guard"],
        "run_path": str(run_path),
        "run_sha256": _sha256(run_path),
        "event_path": str(event_path),
        "event_sha256": _sha256(event_path),
        "events": len(rows),
        "observed_maximum_integrand_corrected": max(values),
        "used_as_envelope_floor_only": True,
        "excluded_from_fixed_trial_calibration_statistics": True,
    }


def _zero_success_upper_rate(trials: int, confidence: float) -> float:
    """Exact one-sided binomial upper rate after zero observed successes."""
    if trials <= 0:
        raise ValueError("zero-success bound requires positive trials")
    if not 0.0 < confidence < 1.0:
        raise ValueError("zero-success confidence must lie in (0,1)")
    return -math.expm1(math.log1p(-confidence) / trials)


def _finalize_calibration(
    args: argparse.Namespace,
    sources: list[tuple[Path, dict]],
    grouped: dict[str, list[tuple[Path, dict, dict, Path]]],
) -> Path:
    first_manifest_path, manifest = sources[0]
    root = first_manifest_path.parent
    safety_factor = float(getattr(args, "envelope_safety_factor", 1.2))
    maximum_duplicate_fraction = float(
        getattr(args, "maximum_duplicate_fraction", 0.05)
    )
    minimum_targets = int(
        getattr(args, "minimum_component_targets", 20)
    )
    minimum_provisional_inside_targets = int(
        getattr(args, "minimum_provisional_inside_targets", 1000)
    )
    allow_zero_complement = bool(
        getattr(args, "allow_zero_complement", False)
    )
    zero_complement_confidence = float(
        getattr(args, "zero_complement_confidence", 0.95)
    )
    maximum_zero_complement_target_rate = float(
        getattr(args, "maximum_zero_complement_target_rate", 1.0e-6)
    )
    allow_revision_mismatch = bool(
        getattr(args, "allow_calibration_revision_mismatch", False)
    )
    revision_compatibility_rationale = str(
        getattr(args, "revision_compatibility_rationale", "") or ""
    ).strip()
    generator_revisions = sorted(
        {str(source["generator_revision"]) for _, source in sources}
    )
    raw_pilot_paths = (
        getattr(args, "additional_inside_pilot_run", None) or []
    )
    if len({Path(path).resolve() for path in raw_pilot_paths}) != len(
        raw_pilot_paths
    ):
        raise ValueError(
            "the same --additional-inside-pilot-run was supplied more "
            "than once"
        )
    pilot_observations: dict[str, list[dict[str, object]]] = {}
    for path in raw_pilot_paths:
        observation = _additional_inside_pilot_observation(Path(path))
        pilot_observations.setdefault(
            str(observation["stratum_id"]), []
        ).append(observation)
    if safety_factor < 1.0:
        raise ValueError("--envelope-safety-factor must be at least one")
    if not 0.0 <= maximum_duplicate_fraction < 1.0:
        raise ValueError("--maximum-duplicate-fraction must lie in [0,1)")
    if minimum_targets < 1:
        raise ValueError("--minimum-component-targets must be positive")
    if minimum_provisional_inside_targets < minimum_targets:
        raise ValueError(
            "--minimum-provisional-inside-targets must be at least "
            "--minimum-component-targets"
        )
    if not 0.0 < zero_complement_confidence < 1.0:
        raise ValueError(
            "--zero-complement-confidence must lie strictly between zero "
            "and one"
        )
    if not 0.0 < maximum_zero_complement_target_rate < 1.0:
        raise ValueError(
            "--maximum-zero-complement-target-rate must lie in (0,1)"
        )
    alpha = float(manifest["core_fraction"])
    strata: list[dict[str, object]] = []
    for stratum_id, items in sorted(
        grouped.items(), key=lambda item: int(item[1][0][1]["flat_index"])
    ):
        inside_values: list[float] = []
        complement_values: list[float] = []
        for item_root, record, _run, _source_path in items:
            csv_path = item_root / (
                str(record["output_stem"]) + ".calibration.csv"
            )
            for component, inside, value in _calibration_rows(csv_path):
                if component == 1:
                    if inside != 1:
                        raise Mode4Error(
                            f"{csv_path}: inside-guard component escaped guard"
                        )
                    inside_values.append(value)
                elif component == 0:
                    if inside != 0:
                        raise Mode4Error(
                            f"{csv_path}: complement component entered guard"
                        )
                    complement_values.append(value)
                else:
                    raise Mode4Error(f"{csv_path}: invalid component")
        inside_trials = sum(int(run["core_trials"]) for _, _, run, _ in items)
        complement_trials = sum(
            int(run["noncore_trials"]) for _, _, run, _ in items
        )
        phase_volumes = [
            float(
                radiative_survey.parse_norm(
                    item_root / (str(record["output_stem"]) + ".norm")
                )["mode4_phase_volume"]
            )
            for item_root, record, _run, _source_path in items
        ]
        phase_volume = phase_volumes[0]
        if any(not _close(value, phase_volume) for value in phase_volumes[1:]):
            raise Mode4Error(
                f"{stratum_id}: pooled runs disagree on phase-space volume"
            )
        first = items[0][1]
        additional_observations = pilot_observations.pop(stratum_id, [])
        for observation in additional_observations:
            frozen = {
                name: observation[name]
                for name in ("flat_index", "indices", "bounds", "guard")
            }
            expected = {
                name: first[name]
                for name in ("flat_index", "indices", "bounds", "guard")
            }
            if frozen != expected:
                raise Mode4Error(
                    f"{observation['run_path']}: pilot stratum geometry "
                    f"differs from calibration {stratum_id}"
                )
        guard_volume = float(first["guard"]["normalized_volume"])
        inside_mass = alpha + (1.0 - alpha) * guard_volume
        complement_mass = (1.0 - alpha) * (1.0 - guard_volume)
        inside_summary = _component_summary(
            inside_values, inside_trials, inside_mass, phase_volume
        )
        complement_summary = _component_summary(
            complement_values,
            complement_trials,
            complement_mass,
            phase_volume,
        )
        sigma = (
            float(inside_summary["integrated_cross_section_microbarn"])
            + float(complement_summary["integrated_cross_section_microbarn"])
        )
        sem = math.hypot(
            float(inside_summary["integrated_cross_section_sem_microbarn"]),
            float(
                complement_summary[
                    "integrated_cross_section_sem_microbarn"
                ]
            ),
        )
        envelope_sources: list[tuple[str, float]] = []
        for prefix, values in (
            ("inside_guard", inside_values),
            ("guard_complement", complement_values),
            ("all", inside_values + complement_values),
        ):
            for label, probability in (
                ("p90", 0.90),
                ("p99", 0.99),
                ("p999", 0.999),
                ("maximum", 1.0),
            ):
                value = _quantile(values, probability)
                if value is not None and value > 0.0:
                    envelope_sources.append(
                        (f"{prefix}_{label}", safety_factor * value)
                    )
        evaluations: list[dict[str, object]] = []
        seen: list[float] = []
        for source, envelope in sorted(
            envelope_sources, key=lambda item: item[1]
        ):
            if any(math.isclose(envelope, old, rel_tol=1.0e-12)
                   for old in seen):
                continue
            seen.append(envelope)
            evaluations.append(
                {
                    "source": source,
                    **_envelope_evaluation(
                        envelope,
                        inside_values,
                        complement_values,
                        inside_trials,
                        complement_trials,
                        inside_mass,
                        complement_mass,
                    ),
                }
            )
        inside_support = len(inside_values) >= minimum_targets
        complement_support = len(complement_values) >= minimum_targets
        enough_support = inside_support and complement_support
        complement_zero = not complement_values
        zero_upper_rate = (
            _zero_success_upper_rate(
                complement_trials, zero_complement_confidence
            )
            if complement_zero
            else None
        )
        zero_rate_passes = (
            zero_upper_rate is not None
            and zero_upper_rate <= maximum_zero_complement_target_rate
        )
        acceptable = [
            item
            for item in evaluations
            if float(item["expected_duplicate_event_fraction"])
            <= maximum_duplicate_fraction
        ]
        provisional_inside_support = (
            len(inside_values) >= minimum_provisional_inside_targets
        )
        inside_observed_maximum = (
            max(inside_values) if inside_values else None
        )
        additional_observed_maximum = (
            max(
                float(
                    observation[
                        "observed_maximum_integrand_corrected"
                    ]
                )
                for observation in additional_observations
            )
            if additional_observations
            else None
        )
        combined_observed_maximum = max(
            value
            for value in (
                inside_observed_maximum,
                additional_observed_maximum,
            )
            if value is not None
        ) if (
            inside_observed_maximum is not None
            or additional_observed_maximum is not None
        ) else None
        provisional_envelope = (
            safety_factor * combined_observed_maximum
            if combined_observed_maximum is not None
            else None
        )
        provisional_source = (
            "additional_inside_pilot_observed_maximum"
            if (
                additional_observed_maximum is not None
                and (
                    inside_observed_maximum is None
                    or additional_observed_maximum
                    > inside_observed_maximum
                )
            )
            else "inside_guard_maximum"
        )
        provisional_candidate = (
            {
                "source": provisional_source,
                **_envelope_evaluation(
                    provisional_envelope,
                    inside_values,
                    complement_values,
                    inside_trials,
                    complement_trials,
                    inside_mass,
                    complement_mass,
                ),
            }
            if provisional_envelope is not None
            else None
        )
        provisional_candidate_passes = (
            provisional_candidate is not None
            and float(
                provisional_candidate[
                    "expected_duplicate_event_fraction"
                ]
            )
            <= maximum_duplicate_fraction
        )
        strict_recommendation = enough_support and bool(acceptable)
        provisional_recommendation = (
            allow_zero_complement
            and inside_support
            and provisional_inside_support
            and complement_zero
            and zero_rate_passes
            and provisional_candidate_passes
        )
        recommendation: Optional[dict[str, object]] = None
        recommendation_status: str
        recommendation_basis: Optional[str] = None
        pilot_readiness = "not_ready"
        if strict_recommendation:
            recommendation = {
                **acceptable[0],
                "provisional": False,
                "basis": "both_calibration_components_observed",
            }
            recommendation_status = "recommended"
            recommendation_basis = str(recommendation["basis"])
            pilot_readiness = "ready"
        elif provisional_recommendation:
            recommendation = {
                **provisional_candidate,
                "provisional": True,
                "basis": (
                    "safety_scaled_inside_guard_observed_maximum_with_"
                    "zero_complement_rate_bound"
                ),
                "unobserved_complement_warning": (
                    "Expected yield and duplicate metrics do not include an "
                    "unobserved complement contribution; stochastic "
                    "multiplicity preserves correctness if one appears."
                ),
            }
            recommendation_status = "provisional_zero_complement"
            recommendation_basis = str(recommendation["basis"])
            pilot_readiness = "ready_provisional_zero_complement"
        elif not inside_support and not complement_support:
            recommendation_status = (
                "insufficient_both_component_target_support"
            )
        elif not inside_support:
            recommendation_status = "insufficient_inside_guard_target_support"
        elif (
            complement_zero
            and allow_zero_complement
            and not zero_rate_passes
        ):
            recommendation_status = "insufficient_zero_complement_exposure"
        elif (
            complement_zero
            and allow_zero_complement
            and zero_rate_passes
            and not provisional_inside_support
        ):
            recommendation_status = (
                "insufficient_provisional_inside_envelope_support"
            )
        elif (
            complement_zero
            and allow_zero_complement
            and zero_rate_passes
            and not provisional_candidate_passes
        ):
            recommendation_status = "no_candidate_meets_duplicate_limit"
        elif not complement_support:
            recommendation_status = (
                "insufficient_guard_complement_target_support"
            )
        else:
            recommendation_status = "no_candidate_meets_duplicate_limit"
        strata.append(
            {
                "stratum_id": stratum_id,
                "flat_index": first["flat_index"],
                "indices": first["indices"],
                "bounds": first["bounds"],
                "guard": first["guard"],
                "calibration_runs": len(items),
                "source_campaigns": len(
                    {str(source_path) for *_, source_path in items}
                ),
                "guard_volume": guard_volume,
                "inside_guard": inside_summary,
                "guard_complement": complement_summary,
                "integrated_cross_section_microbarn": sigma,
                "integrated_cross_section_sem_microbarn": sem,
                "estimated_guard_complement_cross_section_fraction": (
                    float(
                        complement_summary[
                            "integrated_cross_section_microbarn"
                        ]
                    ) / sigma
                    if sigma > 0.0
                    else None
                ),
                "envelope_candidates": evaluations,
                "additional_inside_pilot_observations": (
                    additional_observations
                ),
                "provisional_inside_envelope_support": {
                    "observed_targets": len(inside_values),
                    "minimum_required_targets": (
                        minimum_provisional_inside_targets
                    ),
                    "passes_target_threshold": (
                        provisional_inside_support
                    ),
                    "empirical_next_target_rank_resolution": (
                        1.0 / (len(inside_values) + 1.0)
                    ),
                    "observed_maximum_integrand_corrected": (
                        inside_observed_maximum
                    ),
                    "additional_pilot_observed_maximum_"
                    "integrand_corrected": (
                        additional_observed_maximum
                    ),
                    "combined_observed_maximum_integrand_corrected": (
                        combined_observed_maximum
                    ),
                    "safety_scaled_observed_maximum_sigr_max": (
                        provisional_envelope
                    ),
                    "rank_resolution_is_not_a_confidence_bound": True,
                },
                "zero_complement_stopping_test": {
                    "policy_enabled": allow_zero_complement,
                    "observed_complement_targets": len(complement_values),
                    "complement_trials": complement_trials,
                    "confidence_level": zero_complement_confidence,
                    "one_sided_upper_target_rate": zero_upper_rate,
                    "maximum_allowed_target_rate": (
                        maximum_zero_complement_target_rate
                    ),
                    "passes_rate_threshold": zero_rate_passes,
                    "bounds_occurrence_rate_only": True,
                    "does_not_bound_cross_section_or_integrand_magnitude": True,
                },
                "recommendation_status": recommendation_status,
                "recommendation_basis": recommendation_basis,
                "pilot_readiness": pilot_readiness,
                "recommended_envelope": recommendation,
            }
        )
    if pilot_observations:
        unknown = ", ".join(sorted(pilot_observations))
        raise Mode4Error(
            "additional inside-pilot evidence does not match a finalized "
            f"stratum: {unknown}"
        )
    payload: dict[str, object] = {
        "schema": CALIBRATION_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "finalizer_revision": _source_revision(),
        "finalizer_source_sha256": _sha256(Path(__file__).resolve()),
        "source_manifests": [
            {
                "path": str(path),
                "sha256": _sha256(path),
                "inside_guard_trial_fraction": source[
                    "calibration_inside_guard_fraction"
                ],
                "trials_per_replica": source[
                    "calibration_trials_per_replica"
                ],
                "replicas_per_stratum": source["replicas_per_stratum"],
                "generator_revision": source["generator_revision"],
            }
            for path, source in sources
        ],
        "generator_revision": (
            generator_revisions[0]
            if len(generator_revisions) == 1
            else None
        ),
        "generator_revisions": generator_revisions,
        "revision_compatibility_override": {
            "enabled": allow_revision_mismatch,
            "rationale": (
                revision_compatibility_rationale
                if allow_revision_mismatch
                else None
            ),
        },
        "analysis_config_sha256": manifest["analysis_config_sha256"],
        "guard_recipes_sha256": manifest["guard_recipes_sha256"],
        "guard_refinements_sha256": manifest.get(
            "guard_refinements_sha256"
        ),
        "guard_candidate": manifest["guard_candidate"],
        "core_fraction": alpha,
        "analysis_selection": manifest["analysis_selection"],
        "calibration_proposal": "guard_partition",
        "envelope_safety_factor": safety_factor,
        "maximum_duplicate_fraction": maximum_duplicate_fraction,
        "minimum_component_targets": minimum_targets,
        "minimum_provisional_inside_targets": (
            minimum_provisional_inside_targets
        ),
        "zero_complement_policy": {
            "enabled": allow_zero_complement,
            "confidence_level": zero_complement_confidence,
            "maximum_target_rate": maximum_zero_complement_target_rate,
            "bounds_occurrence_rate_only": True,
            "does_not_bound_cross_section_or_integrand_magnitude": True,
            "strict_observed_complement_rule_retained": True,
        },
        "stratum_count": len(strata),
        "strata": strata,
    }
    output = (
        args.output.resolve()
        if args.output is not None
        else root / "envelope_calibration.json"
    )
    _write_calibration_report(output, payload)
    return output


def _write_calibration_report(output: Path, payload: dict) -> None:
    """Write one calibration JSON and its stable tabular readiness view."""
    _write_json(output, payload)
    with output.with_suffix(".tsv").open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        writer = csv.writer(destination, delimiter="\t")
        writer.writerow(
            [
                "flat_index",
                "stratum_id",
                "recommendation_status",
                "pilot_readiness",
                "recommendation_basis",
                "recommended_sigr_max",
                "expected_events_per_proposal",
                "expected_duplicate_event_fraction",
                "inside_guard_targets",
                "inside_guard_targets_required_for_provisional",
                "inside_empirical_next_target_rank_resolution",
                "guard_complement_targets",
                "zero_complement_upper_target_rate",
                "zero_complement_rate_threshold",
                "estimated_guard_complement_cross_section_fraction",
            ]
        )
        for stratum in payload["strata"]:
            recommendation = stratum["recommended_envelope"]
            writer.writerow(
                [
                    stratum["flat_index"],
                    stratum["stratum_id"],
                    stratum["recommendation_status"],
                    stratum["pilot_readiness"],
                    stratum["recommendation_basis"] or "",
                    (
                        recommendation["sigr_max"]
                        if recommendation is not None
                        else ""
                    ),
                    (
                        recommendation["expected_events_per_proposal"]
                        if recommendation is not None
                        else ""
                    ),
                    (
                        recommendation["expected_duplicate_event_fraction"]
                        if recommendation is not None
                        else ""
                    ),
                    stratum["inside_guard"]["target_candidates"],
                    stratum["provisional_inside_envelope_support"][
                        "minimum_required_targets"
                    ],
                    stratum["provisional_inside_envelope_support"][
                        "empirical_next_target_rank_resolution"
                    ],
                    stratum["guard_complement"]["target_candidates"],
                    stratum["zero_complement_stopping_test"][
                        "one_sided_upper_target_rate"
                    ],
                    stratum["zero_complement_stopping_test"][
                        "maximum_allowed_target_rate"
                    ],
                    stratum[
                        "estimated_guard_complement_cross_section_fraction"
                    ],
                ]
            )


def finalize(args: argparse.Namespace) -> Path:
    raw_paths = getattr(args, "manifests", None)
    if raw_paths is None:
        raw_paths = [args.manifest]
    manifest_paths = [path.resolve() for path in raw_paths]
    if len(set(manifest_paths)) != len(manifest_paths):
        raise ValueError("the same manifest was supplied more than once")
    sources: list[tuple[Path, dict]] = []
    for manifest_path in manifest_paths:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not _supported_manifest_schema(manifest):
            raise ValueError(
                f"{manifest_path}: unsupported manifest schema "
                f"{manifest.get('schema')}"
            )
        sources.append((manifest_path, manifest))
    manifest_path, manifest = sources[0]
    operation = manifest["operation"]
    raw_stratum_ids = getattr(args, "stratum_ids", None)
    selected_stratum_ids = (
        {str(identifier) for identifier in raw_stratum_ids}
        if raw_stratum_ids is not None
        else None
    )
    if selected_stratum_ids is not None:
        if operation != "calibration":
            raise ValueError(
                "stratum-subset finalization is valid only for calibration"
            )
        if not selected_stratum_ids:
            raise ValueError("stratum-subset finalization cannot be empty")
    allow_revision_mismatch = bool(
        getattr(args, "allow_calibration_revision_mismatch", False)
    )
    revision_compatibility_rationale = str(
        getattr(args, "revision_compatibility_rationale", "") or ""
    ).strip()
    if allow_revision_mismatch:
        if operation != "calibration":
            raise ValueError(
                "--allow-calibration-revision-mismatch is valid only for "
                "calibration manifests"
            )
        if not revision_compatibility_rationale:
            raise ValueError(
                "--revision-compatibility-rationale is required with "
                "--allow-calibration-revision-mismatch"
            )
    compatibility_keys = (
        "operation",
        "analysis_config_sha256",
        "guard_recipes_sha256",
        "guard_refinements_sha256",
        "guard_candidate",
        "base_padding",
        "core_fraction",
        "analysis_selection",
        "calibration_proposal",
    )
    if not allow_revision_mismatch:
        compatibility_keys = (
            "generator_revision",
            *compatibility_keys,
        )
    for path, candidate in sources[1:]:
        for name in compatibility_keys:
            if candidate.get(name) != manifest.get(name):
                raise Mode4Error(
                    f"{path}: {name} is incompatible with "
                    f"{manifest_path}"
                )
    if len(sources) > 1 and operation != "calibration":
        raise ValueError(
            "multiple-manifest pooling is supported only for calibration"
        )

    if operation == "calibration":
        grouped_calibration: dict[
            str, list[tuple[Path, dict, dict, Path]]
        ] = {}
        frozen_records: dict[str, dict[str, object]] = {}
        seen_stratum_seeds: set[tuple[str, int]] = set()
        for source_path, source in sources:
            source_root = source_path.parent
            for record in source["runs"]:
                stratum_id = str(record["stratum_id"])
                if (
                    selected_stratum_ids is not None
                    and stratum_id not in selected_stratum_ids
                ):
                    continue
                stratum_seed = (
                    stratum_id,
                    int(record["seed"]),
                )
                if stratum_seed in seen_stratum_seeds:
                    raise Mode4Error(
                        f"{source_path}: repeated seed {record['seed']} for "
                        f"{record['stratum_id']} would double count or "
                        "correlate calibration trials"
                    )
                seen_stratum_seeds.add(stratum_seed)
                frozen = {
                    name: record[name]
                    for name in ("flat_index", "indices", "bounds", "guard")
                }
                old = frozen_records.setdefault(record["stratum_id"], frozen)
                if old != frozen:
                    raise Mode4Error(
                        f"{source_path}: {record['stratum_id']} guard or "
                        "analysis bounds differ across campaigns"
                    )
                run_path = source_root / (
                    str(record["output_stem"]) + ".json"
                )
                if not run_path.is_file():
                    raise FileNotFoundError(run_path)
                completed = json.loads(run_path.read_text(encoding="utf-8"))
                if completed.get("schema") != RUN_SCHEMA:
                    raise Mode4Error(f"{run_path}: unsupported run schema")
                if completed.get("source_manifest_sha256") != _sha256(
                    source_path
                ):
                    raise Mode4Error(
                        f"{run_path}: source manifest hash differs from "
                        "the supplied manifest"
                    )
                for name in ("stratum_id", "flat_index", "replica_index"):
                    if completed[name] != record[name]:
                        raise Mode4Error(
                            f"{run_path}: {name} differs from manifest"
                        )
                if completed.get("operation") != operation:
                    raise Mode4Error(
                        f"{run_path}: operation differs from manifest"
                    )
                grouped_calibration.setdefault(
                    stratum_id, []
                ).append((source_root, record, completed, source_path))
        if selected_stratum_ids is not None:
            missing = sorted(
                selected_stratum_ids - set(grouped_calibration)
            )
            if missing:
                raise Mode4Error(
                    "requested calibration strata have no pooled runs: "
                    + ", ".join(missing)
                )
        return _finalize_calibration(args, sources, grouped_calibration)

    root = manifest_path.parent
    grouped: dict[str, list[tuple[dict, dict]]] = {}
    for record in manifest["runs"]:
        run_path = root / (str(record["output_stem"]) + ".json")
        if not run_path.is_file():
            raise FileNotFoundError(run_path)
        completed = json.loads(run_path.read_text(encoding="utf-8"))
        if completed.get("schema") != RUN_SCHEMA:
            raise Mode4Error(f"{run_path}: unsupported run schema")
        if completed.get("source_manifest_sha256") != _sha256(manifest_path):
            raise Mode4Error(
                f"{run_path}: source manifest hash differs from the "
                "supplied manifest"
            )
        for name in ("stratum_id", "flat_index", "replica_index"):
            if completed[name] != record[name]:
                raise Mode4Error(f"{run_path}: {name} differs from manifest")
        if completed.get("operation") != operation:
            raise Mode4Error(f"{run_path}: operation differs from manifest")
        expected_envelope = _run_sigr_max(manifest, record)
        if completed.get("sigr_max") is not None and not _close(
            float(completed["sigr_max"]), expected_envelope
        ):
            raise Mode4Error(f"{run_path}: sigr_max differs from manifest")
        grouped.setdefault(record["stratum_id"], []).append(
            (record, completed)
        )
    strata: list[dict[str, object]] = []
    for stratum_id, items in sorted(
        grouped.items(), key=lambda item: int(item[1][0][0]["flat_index"])
    ):
        envelopes = [
            _run_sigr_max(manifest, record) for record, _run in items
        ]
        if any(not _close(value, envelopes[0]) for value in envelopes[1:]):
            raise Mode4Error(
                f"{stratum_id}: generation replicas use different envelopes"
            )
        total_proposals = sum(int(run["ntries"]) for _, run in items)
        total_events = sum(int(run["events"]) for _, run in items)
        combined_sigma = sum(
            int(run["ntries"]) * float(run["sig_sum_microbarn"])
            for _, run in items
        ) / total_proposals
        first = items[0][0]
        strata.append(
            {
                "stratum_id": stratum_id,
                "flat_index": first["flat_index"],
                "indices": first["indices"],
                "bounds": first["bounds"],
                "sigr_max": envelopes[0],
                "replicas": len(items),
                "total_proposals": total_proposals,
                "total_events": total_events,
                "combined_sig_sum_microbarn": combined_sigma,
                "pooled_event_weight_microbarn": (
                    combined_sigma / total_events
                ),
                "event_yield_per_proposal": total_events / total_proposals,
                "mcall_max": max(int(run["mcall_max"]) for _, run in items),
                "emitting_candidates": sum(
                    int(run["emitting_candidates"]) for _, run in items
                ),
                "duplicate_events": sum(
                    int(run["duplicate_events"]) for _, run in items
                ),
                "duplicate_event_fraction": (
                    sum(int(run["duplicate_events"]) for _, run in items)
                    / total_events
                ),
                "core_events": sum(
                    int(run["core_events"]) for _, run in items
                ),
                "legacy_events": sum(
                    int(run["legacy_events"]) for _, run in items
                ),
            }
        )
    payload: dict[str, object] = {
        "schema": WEIGHTS_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "finalizer_revision": _source_revision(),
        "finalizer_source_sha256": _sha256(Path(__file__).resolve()),
        "source_manifest": str(manifest_path),
        "source_manifest_sha256": _sha256(manifest_path),
        "generator_revision": manifest["generator_revision"],
        "analysis_config_sha256": manifest["analysis_config_sha256"],
        "guard_recipes_sha256": manifest["guard_recipes_sha256"],
        "guard_refinements_sha256": manifest.get(
            "guard_refinements_sha256"
        ),
        "envelope_mode": manifest.get("envelope_mode", "shared_scalar"),
        "envelope_calibration_sha256": manifest.get(
            "envelope_calibration_sha256"
        ),
        "guard_candidate": manifest["guard_candidate"],
        "core_fraction": manifest["core_fraction"],
        "legacy_tail_fraction": manifest["legacy_tail_fraction"],
        "stratum_count": len(strata),
        "strata": strata,
    }
    output = (
        args.output.resolve()
        if args.output is not None
        else root / "campaign_weights.json"
    )
    _write_json(output, payload)
    tsv = output.with_suffix(".tsv")
    with tsv.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.writer(destination, delimiter="\t")
        writer.writerow(
            [
                "flat_index",
                "stratum_id",
                "replicas",
                "total_events",
                "total_proposals",
                "sigr_max",
                "combined_sig_sum_microbarn",
                "pooled_event_weight_microbarn",
                "event_yield_per_proposal",
                "emitting_candidates",
                "duplicate_events",
                "duplicate_event_fraction",
                "core_events",
                "legacy_events",
            ]
        )
        for stratum in strata:
            writer.writerow(
                [
                    stratum[name]
                    for name in (
                        "flat_index",
                        "stratum_id",
                        "replicas",
                        "total_events",
                        "total_proposals",
                        "sigr_max",
                        "combined_sig_sum_microbarn",
                        "pooled_event_weight_microbarn",
                        "event_yield_per_proposal",
                        "emitting_candidates",
                        "duplicate_events",
                        "duplicate_event_fraction",
                        "core_events",
                        "legacy_events",
                    )
                ]
            )
    return output


def _wilson_upper_fraction(
    successes: int, trials: int, confidence: float
) -> Optional[float]:
    """One-sided Wilson upper diagnostic for an observed event fraction."""
    if trials <= 0:
        return None
    if successes < 0 or successes > trials:
        raise ValueError("Wilson successes must lie in [0,trials]")
    if not 0.5 < confidence < 1.0:
        raise ValueError("Wilson confidence must lie in (0.5,1)")
    z = statistics.NormalDist().inv_cdf(confidence)
    observed = successes / trials
    z_squared = z * z
    denominator = 1.0 + z_squared / trials
    numerator = (
        observed
        + z_squared / (2.0 * trials)
        + z
        * math.sqrt(
            observed * (1.0 - observed) / trials
            + z_squared / (4.0 * trials * trials)
        )
    )
    return min(1.0, numerator / denominator)


def _guard_face_violations(
    box: GuardBox, coordinates: dict[str, float]
) -> list[dict[str, object]]:
    violations: list[dict[str, object]] = []
    for axis in AXES[:-1]:
        value = coordinates[axis]
        lower, upper = box.nonperiodic[axis]
        if value < lower:
            violations.append(
                {
                    "axis": axis,
                    "face": "lower",
                    "value": value,
                    "boundary": lower,
                    "excursion": lower - value,
                }
            )
        elif value > upper:
            violations.append(
                {
                    "axis": axis,
                    "face": "upper",
                    "value": value,
                    "boundary": upper,
                    "excursion": value - upper,
                }
            )
    phi = coordinates["hadron_phi_base"]
    relative = (phi - box.phi_origin + 0.5) % 1.0 - 0.5
    lower, upper = box.phi_relative
    if relative < lower:
        violations.append(
            {
                "axis": "hadron_phi_base",
                "face": "lower",
                "value": phi,
                "relative_value": relative,
                "boundary": lower,
                "excursion": lower - relative,
            }
        )
    elif relative > upper:
        violations.append(
            {
                "axis": "hadron_phi_base",
                "face": "upper",
                "value": phi,
                "relative_value": relative,
                "boundary": upper,
                "excursion": relative - upper,
            }
        )
    return violations


def _pilot_source_manifest(
    run_path: Path, run: dict
) -> tuple[Path, dict]:
    candidates = [Path(str(run.get("source_manifest", "")))]
    if len(run_path.parents) >= 3:
        candidates.append(run_path.parents[2] / "manifest.json")
    expected_hash = str(run.get("source_manifest_sha256", ""))
    visited: set[Path] = set()
    for candidate in candidates:
        if not str(candidate):
            continue
        resolved = candidate.expanduser().resolve()
        if resolved in visited or not resolved.is_file():
            continue
        visited.add(resolved)
        if _sha256(resolved) != expected_hash:
            continue
        manifest = json.loads(resolved.read_text(encoding="utf-8"))
        if not _supported_manifest_schema(manifest):
            raise Mode4Error(f"{resolved}: unsupported manifest schema")
        if manifest.get("operation") != "generation":
            raise Mode4Error(f"{resolved}: expected generation manifest")
        envelope_sha256 = manifest.get("envelope_calibration_sha256")
        if envelope_sha256 is not None:
            snapshot = resolved.parent / str(
                manifest.get(
                    "envelope_calibration_snapshot",
                    "envelope_calibration.json",
                )
            )
            if not snapshot.is_file() or _sha256(snapshot) != envelope_sha256:
                raise Mode4Error(
                    f"{resolved}: envelope-calibration snapshot changed"
                )
        return resolved, manifest
    raise Mode4Error(
        f"{run_path}: cannot locate the hashed source generation manifest"
    )


def _pilot_maximum_record(
    candidate: dict[str, object], run_path: Path
) -> dict[str, object]:
    return {
        "run_path": str(run_path),
        "proposal_component": candidate["proposal_component"],
        "proposal_component_name": candidate["proposal_component_name"],
        "inside_geometric_guard": candidate["inside_geometric_guard"],
        "integrand_corrected": candidate["integrand_corrected"],
        "ratio_to_run_sigr_max": candidate["ratio_to_run_sigr_max"],
        "multiplicity": candidate["multiplicity"],
        "coordinates": candidate["coordinates"],
        "guard_face_violations": candidate["guard_face_violations"],
    }


def _load_pilot_run(path: Path) -> dict[str, object]:
    run_path = path.expanduser().resolve()
    run = json.loads(run_path.read_text(encoding="utf-8"))
    if run.get("schema") != RUN_SCHEMA:
        raise Mode4Error(f"{run_path}: unsupported mode-4 run schema")
    if run.get("operation") != "generation":
        raise Mode4Error(f"{run_path}: expected a generation pilot")
    manifest_path, manifest = _pilot_source_manifest(run_path, run)
    matching = [
        record
        for record in manifest["runs"]
        if (
            record["stratum_id"] == run["stratum_id"]
            and int(record["flat_index"]) == int(run["flat_index"])
            and int(record["replica_index"]) == int(run["replica_index"])
        )
    ]
    if len(matching) != 1:
        raise Mode4Error(
            f"{run_path}: source manifest does not identify one run"
        )
    record = matching[0]
    for name in (
        "stratum_id",
        "flat_index",
        "replica_index",
        "indices",
        "bounds",
        "guard",
        "seed",
    ):
        if run[name] != record[name]:
            raise Mode4Error(f"{run_path}: {name} differs from manifest")
    event_path = run_path.with_suffix(".mode4.csv")
    event_summary = _validate_event_diagnostics(
        event_path, manifest, record, int(run["events"])
    )
    if (
        int(event_summary["core_events"]) != int(run["core_events"])
        or int(event_summary["legacy_events"]) != int(run["legacy_events"])
    ):
        raise Mode4Error(
            f"{event_path}: proposal-component counts differ from {run_path}"
        )
    with event_path.open(encoding="utf-8", newline="") as source:
        schema = source.readline().strip()
        if schema != f"# schema={MODE4_KINEMATICS_SCHEMA}":
            raise Mode4Error(f"{event_path}: unexpected event schema")
        reader = csv.DictReader(source)
        if tuple(reader.fieldnames or ()) != MODE4_KINEMATICS_COLUMNS:
            raise Mode4Error(f"{event_path}: unexpected event columns")
        rows = list(reader)
    if len(rows) != int(run["events"]):
        raise Mode4Error(
            f"{event_path}: event count differs from {run_path}"
        )
    if not rows:
        raise Mode4Error(f"{event_path}: pilot has no events")
    envelope = _run_sigr_max(manifest, record)
    if run.get("sigr_max") is not None and not _close(
        float(run["sigr_max"]), envelope
    ):
        raise Mode4Error(f"{run_path}: sigr_max differs from manifest")
    groups: dict[tuple[str, ...], list[dict[str, str]]] = {}
    for expected_event, row in enumerate(rows, start=1):
        if int(row["event"]) != expected_event:
            raise Mode4Error(f"{event_path}: nonsequential events")
        key = tuple(row[name].strip() for name in MODE4_KINEMATICS_COLUMNS[1:])
        groups.setdefault(key, []).append(row)
    if len(groups) != int(run["emitting_candidates"]):
        raise Mode4Error(
            f"{event_path}: distinct candidates differ from {run_path}"
        )
    duplicate_events = sum(len(items) - 1 for items in groups.values())
    if duplicate_events != int(run["duplicate_events"]):
        raise Mode4Error(
            f"{event_path}: duplicate count differs from {run_path}"
        )
    box = _guard_box_from_manifest(record)
    classifications = {
        name: {
            "emitting_candidates": 0,
            "events": 0,
            "duplicate_events": 0,
        }
        for name in (
            "guard_focused_inside",
            "guard_focused_outside",
            "legacy_inside",
            "legacy_outside",
        )
    }
    maxima: dict[str, dict[str, object]] = {}
    outside_faces: dict[str, dict[str, object]] = {}
    duplicate_candidates = 0
    envelope_exceeding_candidates = 0
    maximum_multiplicity = 0
    for items in groups.values():
        raw = items[0]
        multiplicity = len(items)
        component = int(raw["proposal_component"])
        if component not in (0, 1):
            raise Mode4Error(f"{event_path}: invalid proposal component")
        coordinates = {name: float(raw[name]) for name in AXES}
        if any(not math.isfinite(value) for value in coordinates.values()):
            raise Mode4Error(f"{event_path}: nonfinite proposal coordinate")
        integrand = float(raw["integrand_corrected"])
        if not math.isfinite(integrand) or integrand <= 0.0:
            raise Mode4Error(f"{event_path}: invalid corrected integrand")
        inside_guard = _fortran_guard_contains(box, coordinates)
        if component == 1 and not inside_guard:
            raise Mode4Error(
                f"{event_path}: guard-focused component escaped guard"
            )
        component_name = "guard_focused" if component == 1 else "legacy"
        region_name = "inside" if inside_guard else "outside"
        classification_name = f"{component_name}_{region_name}"
        classification = classifications[classification_name]
        classification["emitting_candidates"] += 1
        classification["events"] += multiplicity
        classification["duplicate_events"] += multiplicity - 1
        if multiplicity > 1:
            duplicate_candidates += 1
        ratio = integrand / envelope
        if ratio > 1.0:
            envelope_exceeding_candidates += 1
        maximum_multiplicity = max(maximum_multiplicity, multiplicity)
        violations = (
            [] if inside_guard else _guard_face_violations(box, coordinates)
        )
        candidate: dict[str, object] = {
            "proposal_component": component,
            "proposal_component_name": component_name,
            "inside_geometric_guard": inside_guard,
            "integrand_corrected": integrand,
            "ratio_to_run_sigr_max": ratio,
            "multiplicity": multiplicity,
            "coordinates": coordinates,
            "guard_face_violations": violations,
        }
        for category in (
            "all",
            component_name,
            "inside_guard" if inside_guard else "guard_complement",
            classification_name,
        ):
            old = maxima.get(category)
            if (
                old is None
                or integrand > float(old["integrand_corrected"])
            ):
                maxima[category] = _pilot_maximum_record(
                    candidate, run_path
                )
        for violation in violations:
            key = f"{violation['axis']}:{violation['face']}"
            summary = outside_faces.setdefault(
                key,
                {
                    "axis": violation["axis"],
                    "face": violation["face"],
                    "boundary": violation["boundary"],
                    "emitting_candidates": 0,
                    "events": 0,
                    "duplicate_events": 0,
                    "maximum_excursion": 0.0,
                    "most_extreme_value": None,
                },
            )
            summary["emitting_candidates"] += 1
            summary["events"] += multiplicity
            summary["duplicate_events"] += multiplicity - 1
            if float(violation["excursion"]) > float(
                summary["maximum_excursion"]
            ):
                summary["maximum_excursion"] = violation["excursion"]
                summary["most_extreme_value"] = violation["value"]
    if maximum_multiplicity > int(run["mcall_max"]):
        raise Mode4Error(
            f"{event_path}: observed multiplicity exceeds recorded maximum"
        )
    if sum(
        int(item["events"]) for item in classifications.values()
    ) != int(run["events"]):
        raise Mode4Error(f"{event_path}: classification count is incomplete")
    return {
        "stratum_id": str(run["stratum_id"]),
        "flat_index": int(run["flat_index"]),
        "indices": run["indices"],
        "bounds": run["bounds"],
        "guard": run["guard"],
        "seed": int(run["seed"]),
        "replica_index": int(run["replica_index"]),
        "generator_revision": manifest["generator_revision"],
        "analysis_config_sha256": manifest["analysis_config_sha256"],
        "guard_recipes_sha256": manifest["guard_recipes_sha256"],
        "guard_refinements_sha256": manifest.get(
            "guard_refinements_sha256"
        ),
        "guard_candidate": manifest["guard_candidate"],
        "core_fraction": float(manifest["core_fraction"]),
        "analysis_selection": manifest["analysis_selection"],
        "sigr_max": envelope,
        "events": int(run["events"]),
        "ntries": int(run["ntries"]),
        "sig_sum_microbarn": float(run["sig_sum_microbarn"]),
        "event_overshoot": int(run["event_overshoot"]),
        "emitting_candidates": int(run["emitting_candidates"]),
        "duplicate_events": int(run["duplicate_events"]),
        "duplicate_candidates": duplicate_candidates,
        "mcall_max": int(run["mcall_max"]),
        "envelope_exceeding_candidates": envelope_exceeding_candidates,
        "classifications": classifications,
        "maxima": maxima,
        "guard_complement_faces": sorted(
            outside_faces.values(),
            key=lambda item: (
                str(item["axis"]),
                str(item["face"]),
            ),
        ),
        "run_path": str(run_path),
        "run_sha256": _sha256(run_path),
        "event_path": str(event_path),
        "event_sha256": _sha256(event_path),
        "manifest_path": str(manifest_path),
        "manifest_sha256": _sha256(manifest_path),
    }


def _merge_pilot_maxima(
    runs: list[dict[str, object]]
) -> dict[str, dict[str, object]]:
    merged: dict[str, dict[str, object]] = {}
    for run in runs:
        for category, candidate in run["maxima"].items():
            old = merged.get(category)
            if (
                old is None
                or float(candidate["integrand_corrected"])
                > float(old["integrand_corrected"])
            ):
                merged[category] = candidate
    return merged


def _merge_guard_complement_faces(
    runs: list[dict[str, object]]
) -> list[dict[str, object]]:
    merged: dict[str, dict[str, object]] = {}
    for run in runs:
        for face in run["guard_complement_faces"]:
            key = f"{face['axis']}:{face['face']}"
            summary = merged.setdefault(
                key,
                {
                    "axis": face["axis"],
                    "face": face["face"],
                    "boundary": face["boundary"],
                    "emitting_candidates": 0,
                    "events": 0,
                    "duplicate_events": 0,
                    "maximum_excursion": 0.0,
                    "most_extreme_value": None,
                },
            )
            summary["emitting_candidates"] += int(
                face["emitting_candidates"]
            )
            summary["events"] += int(face["events"])
            summary["duplicate_events"] += int(face["duplicate_events"])
            if float(face["maximum_excursion"]) > float(
                summary["maximum_excursion"]
            ):
                summary["maximum_excursion"] = face["maximum_excursion"]
                summary["most_extreme_value"] = face[
                    "most_extreme_value"
                ]
    return sorted(
        merged.values(),
        key=lambda item: (str(item["axis"]), str(item["face"])),
    )


def _calibration_validation_metadata(
    calibration_path: Path, calibration: dict
) -> tuple[dict, set[str], dict[str, dict]]:
    selection = calibration.get("analysis_selection")
    revisions = {
        str(value)
        for value in (
            calibration.get("generator_revisions")
            or (
                [calibration["generator_revision"]]
                if calibration.get("generator_revision") is not None
                else []
            )
        )
    }
    guards = {
        str(item["stratum_id"]): item["guard"]
        for item in calibration["strata"]
        if item.get("guard") is not None
    }
    required_strata = {
        str(item["stratum_id"]) for item in calibration["strata"]
    }
    if (
        selection is not None
        and revisions
        and required_strata.issubset(guards)
    ):
        return selection, revisions, guards
    loaded = 0
    for source_record in calibration.get("source_manifests", []):
        if isinstance(source_record, str):
            raw_path = source_record
            expected_hash = None
        else:
            raw_path = str(source_record["path"])
            expected_hash = source_record.get("sha256")
        candidate = Path(raw_path).expanduser()
        candidates = [candidate]
        if not candidate.is_absolute():
            candidates.append(calibration_path.parent / candidate)
        source_path = next(
            (
                path.resolve()
                for path in candidates
                if path.resolve().is_file()
            ),
            None,
        )
        if source_path is None:
            continue
        if expected_hash is not None and _sha256(source_path) != expected_hash:
            raise Mode4Error(
                f"{source_path}: calibration source-manifest hash changed"
            )
        manifest = json.loads(source_path.read_text(encoding="utf-8"))
        if not _supported_manifest_schema(manifest):
            raise Mode4Error(
                f"{source_path}: unsupported calibration source manifest"
            )
        if manifest.get("operation") != "calibration":
            raise Mode4Error(
                f"{source_path}: expected calibration source manifest"
            )
        loaded += 1
        candidate_selection = manifest["analysis_selection"]
        if selection is None:
            selection = candidate_selection
        elif selection != candidate_selection:
            raise Mode4Error(
                f"{source_path}: calibration selections disagree"
            )
        revisions.add(str(manifest["generator_revision"]))
        for record in manifest["runs"]:
            stratum_id = str(record["stratum_id"])
            old = guards.setdefault(stratum_id, record["guard"])
            if old != record["guard"]:
                raise Mode4Error(
                    f"{source_path}: calibration guards disagree for "
                    f"{stratum_id}"
                )
    if loaded == 0 and (
        selection is None
        or not revisions
        or not required_strata.issubset(guards)
    ):
        raise Mode4Error(
            f"{calibration_path}: validation metadata are incomplete and "
            "hashed source manifests are unavailable"
        )
    if selection is None or not revisions:
        raise Mode4Error(
            f"{calibration_path}: incomplete calibration provenance"
        )
    missing_guards = sorted(required_strata - set(guards))
    if missing_guards:
        raise Mode4Error(
            f"{calibration_path}: missing guards for "
            + ", ".join(missing_guards)
        )
    return selection, revisions, guards


def validate_pilots(args: argparse.Namespace) -> Path:
    calibration_path = args.calibration.expanduser().resolve()
    calibration = json.loads(
        calibration_path.read_text(encoding="utf-8")
    )
    if calibration.get("schema") != CALIBRATION_SCHEMA:
        raise Mode4Error(
            f"{calibration_path}: expected {CALIBRATION_SCHEMA}"
        )
    (
        calibration_selection,
        calibration_revisions,
        calibration_guards,
    ) = _calibration_validation_metadata(calibration_path, calibration)
    run_paths = [path.expanduser().resolve() for path in args.runs]
    if len(set(run_paths)) != len(run_paths):
        raise ValueError("the same pilot run was supplied more than once")
    minimum_runs = int(args.minimum_runs)
    minimum_events = int(args.minimum_events)
    maximum_duplicate_fraction = float(args.maximum_duplicate_fraction)
    maximum_guard_complement_fraction = float(
        args.maximum_guard_complement_fraction
    )
    maximum_relative_cross_section_difference = float(
        args.maximum_relative_cross_section_difference
    )
    maximum_cross_section_z_score = float(
        args.maximum_cross_section_z_score
    )
    confidence = float(args.confidence)
    allow_revision_mismatch = bool(
        getattr(args, "allow_pilot_revision_mismatch", False)
    )
    revision_compatibility_rationale = str(
        getattr(args, "revision_compatibility_rationale", "") or ""
    ).strip()
    if minimum_runs < 1 or minimum_events < 1:
        raise ValueError("minimum pilot support must be positive")
    for name, value in (
        ("maximum duplicate fraction", maximum_duplicate_fraction),
        (
            "maximum guard-complement fraction",
            maximum_guard_complement_fraction,
        ),
        (
            "maximum relative cross-section difference",
            maximum_relative_cross_section_difference,
        ),
    ):
        if not 0.0 < value < 1.0:
            raise ValueError(f"{name} must lie in (0,1)")
    if maximum_cross_section_z_score <= 0.0:
        raise ValueError("maximum cross-section z score must be positive")
    if not 0.5 < confidence < 1.0:
        raise ValueError("confidence must lie in (0.5,1)")
    if allow_revision_mismatch and not revision_compatibility_rationale:
        raise ValueError(
            "--revision-compatibility-rationale is required with "
            "--allow-pilot-revision-mismatch"
        )
    loaded = [_load_pilot_run(path) for path in run_paths]
    pilot_revisions = {
        str(run["generator_revision"]) for run in loaded
    }
    revisions_match = (
        len(pilot_revisions) == 1
        and pilot_revisions.issubset(calibration_revisions)
    )
    if not revisions_match and not allow_revision_mismatch:
        raise Mode4Error(
            "pilot generator revisions differ from one another or from "
            "the audited calibration revisions"
        )
    grouped: dict[str, list[dict[str, object]]] = {}
    for run in loaded:
        grouped.setdefault(str(run["stratum_id"]), []).append(run)
    calibration_strata = {
        str(item["stratum_id"]): item for item in calibration["strata"]
    }
    unknown = sorted(set(grouped) - set(calibration_strata))
    if unknown:
        raise Mode4Error(
            "pilot strata are absent from the calibration report: "
            + ", ".join(unknown)
        )
    strata: list[dict[str, object]] = []
    for stratum_id, runs in sorted(
        grouped.items(), key=lambda item: int(item[1][0]["flat_index"])
    ):
        calibration_stratum = calibration_strata[stratum_id]
        first = runs[0]
        if first["bounds"] != calibration_stratum["bounds"]:
            raise Mode4Error(
                f"{stratum_id}: pilot and calibration bounds differ"
            )
        if first["guard"] != calibration_guards[stratum_id]:
            raise Mode4Error(
                f"{stratum_id}: pilot and calibration guards differ"
            )
        compatibility = {
            "analysis_config_sha256": calibration[
                "analysis_config_sha256"
            ],
            "guard_recipes_sha256": calibration["guard_recipes_sha256"],
            "guard_refinements_sha256": calibration.get(
                "guard_refinements_sha256"
            ),
            "guard_candidate": calibration["guard_candidate"],
            "core_fraction": float(calibration["core_fraction"]),
            "analysis_selection": calibration_selection,
        }
        for run in runs:
            for name, expected in compatibility.items():
                if run[name] != expected:
                    raise Mode4Error(
                        f"{run['run_path']}: {name} differs from calibration"
                    )
            for name in ("flat_index", "indices", "bounds", "guard"):
                if run[name] != first[name]:
                    raise Mode4Error(
                        f"{run['run_path']}: pilot {name} differs across runs"
                    )
        seeds = [int(run["seed"]) for run in runs]
        if len(set(seeds)) != len(seeds):
            raise Mode4Error(f"{stratum_id}: pilot seeds are not independent")
        recommended = calibration_stratum.get("recommended_envelope")
        target_envelope = (
            float(recommended["sigr_max"])
            if recommended is not None
            else None
        )
        qualifying = [
            run
            for run in runs
            if (
                target_envelope is not None
                and float(run["sigr_max"])
                <= target_envelope * (1.0 + 1.0e-12)
            )
        ]
        total_proposals = sum(int(run["ntries"]) for run in runs)
        total_events = sum(int(run["events"]) for run in runs)
        pooled_sigma = sum(
            int(run["ntries"]) * float(run["sig_sum_microbarn"])
            for run in runs
        ) / total_proposals
        run_sigmas = [float(run["sig_sum_microbarn"]) for run in runs]
        run_to_run_sem = (
            statistics.stdev(run_sigmas) / math.sqrt(len(run_sigmas))
            if len(run_sigmas) >= 2
            else None
        )
        calibration_sigma = float(
            calibration_stratum["integrated_cross_section_microbarn"]
        )
        calibration_sem = float(
            calibration_stratum[
                "integrated_cross_section_sem_microbarn"
            ]
        )
        difference = pooled_sigma - calibration_sigma
        relative_difference = (
            difference / calibration_sigma
            if calibration_sigma != 0.0
            else None
        )
        combined_sem = (
            math.hypot(calibration_sem, run_to_run_sem)
            if run_to_run_sem is not None
            else calibration_sem
        )
        difference_z_score = (
            difference / combined_sem if combined_sem > 0.0 else None
        )
        qualifying_runs = len(qualifying)
        qualifying_events = sum(int(run["events"]) for run in qualifying)
        qualifying_emitting = sum(
            int(run["emitting_candidates"]) for run in qualifying
        )
        qualifying_duplicate_events = sum(
            int(run["duplicate_events"]) for run in qualifying
        )
        qualifying_duplicate_candidates = sum(
            int(run["duplicate_candidates"]) for run in qualifying
        )
        duplicate_event_fraction = (
            qualifying_duplicate_events / qualifying_events
            if qualifying_events
            else None
        )
        duplicate_candidate_fraction = (
            qualifying_duplicate_candidates / qualifying_emitting
            if qualifying_emitting
            else None
        )
        duplicate_candidate_upper = _wilson_upper_fraction(
            qualifying_duplicate_candidates,
            qualifying_emitting,
            confidence,
        )
        classifications = {
            name: {
                "emitting_candidates": sum(
                    int(run["classifications"][name]["emitting_candidates"])
                    for run in runs
                ),
                "events": sum(
                    int(run["classifications"][name]["events"])
                    for run in runs
                ),
                "duplicate_events": sum(
                    int(run["classifications"][name]["duplicate_events"])
                    for run in runs
                ),
            }
            for name in first["classifications"]
        }
        outside_candidates = sum(
            int(item["emitting_candidates"])
            for name, item in classifications.items()
            if name.endswith("_outside")
        )
        outside_events = sum(
            int(item["events"])
            for name, item in classifications.items()
            if name.endswith("_outside")
        )
        total_emitting = sum(
            int(run["emitting_candidates"]) for run in runs
        )
        guard_complement_event_fraction = outside_events / total_events
        guard_complement_candidate_fraction = (
            outside_candidates / total_emitting
        )
        guard_complement_candidate_upper = _wilson_upper_fraction(
            outside_candidates, total_emitting, confidence
        )
        support_passes = (
            qualifying_runs >= minimum_runs
            and qualifying_events >= minimum_events
        )
        duplicate_passes = (
            duplicate_event_fraction is not None
            and duplicate_event_fraction <= maximum_duplicate_fraction
            and duplicate_candidate_upper is not None
            and duplicate_candidate_upper <= maximum_duplicate_fraction
        )
        guard_passes = (
            guard_complement_event_fraction
            <= maximum_guard_complement_fraction
            and guard_complement_candidate_upper is not None
            and guard_complement_candidate_upper
            <= maximum_guard_complement_fraction
        )
        closure_passes = (
            relative_difference is not None
            and abs(relative_difference)
            <= maximum_relative_cross_section_difference
            and (
                difference_z_score is None
                or abs(difference_z_score)
                <= maximum_cross_section_z_score
            )
        )
        overshoot_passes = all(
            int(run["event_overshoot"]) == 0 for run in runs
        )
        calibration_ready = (
            calibration_stratum.get("pilot_readiness")
            in ("ready", "ready_provisional_zero_complement")
            and recommended is not None
        )
        if not calibration_ready:
            recommendation = "calibration_not_ready"
        elif not support_passes:
            recommendation = "collect_more_pilot_support"
        elif not duplicate_passes:
            recommendation = "increase_envelope_and_revalidate"
        elif not guard_passes:
            recommendation = "review_guard_complement_geometry"
        elif not closure_passes:
            recommendation = "review_cross_section_closure"
        elif not overshoot_passes:
            recommendation = "review_event_overshoot"
        else:
            recommendation = "ready_for_multi_stratum_pilot"
        passed = recommendation == "ready_for_multi_stratum_pilot"
        strata.append(
            {
                "stratum_id": stratum_id,
                "flat_index": first["flat_index"],
                "indices": first["indices"],
                "bounds": first["bounds"],
                "calibration_recommendation_status": calibration_stratum[
                    "recommendation_status"
                ],
                "target_sigr_max": target_envelope,
                "pilot_runs": len(runs),
                "qualifying_conservative_runs": qualifying_runs,
                "pilot_sigr_max_range": [
                    min(float(run["sigr_max"]) for run in runs),
                    max(float(run["sigr_max"]) for run in runs),
                ],
                "total_proposals": total_proposals,
                "total_events": total_events,
                "qualifying_events": qualifying_events,
                "event_yield_per_proposal": (
                    total_events / total_proposals
                ),
                "pooled_sig_sum_microbarn": pooled_sigma,
                "pilot_run_to_run_sem_microbarn": run_to_run_sem,
                "calibration_cross_section_microbarn": calibration_sigma,
                "calibration_sem_microbarn": calibration_sem,
                "cross_section_difference_microbarn": difference,
                "relative_cross_section_difference": relative_difference,
                "combined_calibration_and_run_sem_microbarn": combined_sem,
                "cross_section_difference_z_score": difference_z_score,
                "emitting_candidates": total_emitting,
                "qualifying_emitting_candidates": qualifying_emitting,
                "duplicate_events": sum(
                    int(run["duplicate_events"]) for run in runs
                ),
                "qualifying_duplicate_events": (
                    qualifying_duplicate_events
                ),
                "qualifying_duplicate_event_fraction": (
                    duplicate_event_fraction
                ),
                "qualifying_duplicate_candidates": (
                    qualifying_duplicate_candidates
                ),
                "qualifying_duplicate_candidate_fraction": (
                    duplicate_candidate_fraction
                ),
                "duplicate_candidate_wilson_upper_fraction": (
                    duplicate_candidate_upper
                ),
                "maximum_mcall": max(
                    int(run["mcall_max"]) for run in runs
                ),
                "envelope_exceeding_candidates": sum(
                    int(run["envelope_exceeding_candidates"])
                    for run in runs
                ),
                "event_classification": classifications,
                "guard_complement_events": outside_events,
                "guard_complement_event_fraction": (
                    guard_complement_event_fraction
                ),
                "guard_complement_emitting_candidates": (
                    outside_candidates
                ),
                "guard_complement_candidate_fraction": (
                    guard_complement_candidate_fraction
                ),
                "guard_complement_candidate_wilson_upper_fraction": (
                    guard_complement_candidate_upper
                ),
                "guard_complement_faces": (
                    _merge_guard_complement_faces(runs)
                ),
                "maximum_integrands": _merge_pilot_maxima(runs),
                "support_passes": support_passes,
                "duplicate_overhead_passes": duplicate_passes,
                "guard_complement_passes": guard_passes,
                "cross_section_closure_passes": closure_passes,
                "event_overshoot_passes": overshoot_passes,
                "envelope_decision": (
                    "hold"
                    if duplicate_passes
                    else "increase_and_revalidate"
                ),
                "guard_decision": (
                    "hold"
                    if guard_passes
                    else "review_complement_geometry"
                ),
                "recommendation": recommendation,
                "passed": passed,
                "runs": [
                    {
                        name: run[name]
                        for name in (
                            "run_path",
                            "run_sha256",
                            "event_path",
                            "event_sha256",
                            "manifest_path",
                            "manifest_sha256",
                            "seed",
                            "replica_index",
                            "sigr_max",
                            "events",
                            "ntries",
                            "sig_sum_microbarn",
                            "event_overshoot",
                            "emitting_candidates",
                            "duplicate_events",
                            "duplicate_candidates",
                            "mcall_max",
                            "envelope_exceeding_candidates",
                            "classifications",
                            "maxima",
                            "guard_complement_faces",
                        )
                    }
                    for run in runs
                ],
            }
        )
    passed = bool(strata) and all(bool(item["passed"]) for item in strata)
    payload: dict[str, object] = {
        "schema": PILOT_VALIDATION_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "validator_revision": _source_revision(),
        "validator_source_sha256": _sha256(Path(__file__).resolve()),
        "calibration_report": str(calibration_path),
        "calibration_report_sha256": _sha256(calibration_path),
        "calibration_generator_revisions": sorted(
            calibration_revisions
        ),
        "pilot_generator_revisions": sorted(pilot_revisions),
        "pilot_revision_compatibility_override": {
            "enabled": allow_revision_mismatch,
            "rationale": (
                revision_compatibility_rationale
                if allow_revision_mismatch
                else None
            ),
        },
        "confidence_level": confidence,
        "thresholds": {
            "minimum_runs": minimum_runs,
            "minimum_events": minimum_events,
            "maximum_duplicate_fraction": maximum_duplicate_fraction,
            "maximum_guard_complement_fraction": (
                maximum_guard_complement_fraction
            ),
            "maximum_relative_cross_section_difference": (
                maximum_relative_cross_section_difference
            ),
            "maximum_cross_section_z_score": (
                maximum_cross_section_z_score
            ),
        },
        "wilson_bounds_are_diagnostics_not_formal_iid_guarantees": True,
        "lower_envelope_runs_are_conservative_for_target_envelope": True,
        "passed": passed,
        "recommendation": (
            "ready_for_multi_stratum_pilot"
            if passed
            else "review_stratum_recommendations"
        ),
        "stratum_count": len(strata),
        "passed_strata": sum(bool(item["passed"]) for item in strata),
        "strata": strata,
    }
    output = (
        args.output.expanduser().resolve()
        if args.output is not None
        else calibration_path.with_name("pilot_validation.json")
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    _write_json(output, payload)
    with output.with_suffix(".tsv").open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        writer = csv.writer(destination, delimiter="\t")
        writer.writerow(
            [
                "flat_index",
                "stratum_id",
                "passed",
                "recommendation",
                "target_sigr_max",
                "pilot_runs",
                "qualifying_runs",
                "total_events",
                "qualifying_events",
                "pooled_sig_sum_microbarn",
                "calibration_cross_section_microbarn",
                "relative_cross_section_difference",
                "cross_section_difference_z_score",
                "qualifying_duplicate_event_fraction",
                "duplicate_candidate_wilson_upper_fraction",
                "guard_complement_event_fraction",
                "guard_complement_candidate_wilson_upper_fraction",
                "maximum_mcall",
                "envelope_decision",
                "guard_decision",
            ]
        )
        for stratum in strata:
            writer.writerow(
                [
                    stratum["flat_index"],
                    stratum["stratum_id"],
                    stratum["passed"],
                    stratum["recommendation"],
                    stratum["target_sigr_max"],
                    stratum["pilot_runs"],
                    stratum["qualifying_conservative_runs"],
                    stratum["total_events"],
                    stratum["qualifying_events"],
                    stratum["pooled_sig_sum_microbarn"],
                    stratum["calibration_cross_section_microbarn"],
                    stratum["relative_cross_section_difference"],
                    stratum["cross_section_difference_z_score"],
                    stratum["qualifying_duplicate_event_fraction"],
                    stratum[
                        "duplicate_candidate_wilson_upper_fraction"
                    ],
                    stratum["guard_complement_event_fraction"],
                    stratum[
                        "guard_complement_candidate_wilson_upper_fraction"
                    ],
                    stratum["maximum_mcall"],
                    stratum["envelope_decision"],
                    stratum["guard_decision"],
                ]
            )
    return output


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    refinement_parser = subparsers.add_parser(
        "create-refinement",
        help=(
            "write one evidence-hashed, post-padding guard-face refinement"
        ),
    )
    refinement_parser.add_argument("--config", type=Path, required=True)
    refinement_parser.add_argument("--recipes", type=Path, required=True)
    refinement_parser.add_argument(
        "--candidate", default="padding_0p035"
    )
    refinement_parser.add_argument("--output", type=Path, required=True)
    refinement_parser.add_argument("--stratum", required=True)
    refinement_parser.add_argument(
        "--face",
        action="append",
        required=True,
        metavar="AXIS:FACE:VALUE",
        help=(
            "expand one final guard face; repeat for additional faces"
        ),
    )
    refinement_parser.add_argument("--rationale", required=True)
    refinement_parser.add_argument(
        "--evidence",
        type=Path,
        action="append",
        required=True,
        help="evidence artifact to hash; repeat for additional artifacts",
    )
    refinement_parser.add_argument("--overwrite", action="store_true")

    def add_prepare_common(target: argparse.ArgumentParser) -> None:
        target.add_argument("--config", type=Path, required=True)
        target.add_argument("--recipes", type=Path, required=True)
        target.add_argument(
            "--refinements",
            type=Path,
            help=(
                "optional evidence-hashed post-padding guard refinements"
            ),
        )
        target.add_argument("--input", type=Path, required=True)
        target.add_argument("--output", type=Path, required=True)
        target.add_argument("--candidate", default="padding_0p035")
        target.add_argument("--core-fraction", type=float, default=0.9)
        target.add_argument("--replicas", type=int, default=1)
        target.add_argument("--seed-base", type=int, default=481001)
        target.add_argument("--bin-start", type=int, default=0)
        target.add_argument("--bin-stop", type=int)
        target.add_argument(
            "--flat-index",
            dest="flat_indices",
            type=int,
            action="append",
            help=(
                "select one sparse analysis stratum; repeat as needed and "
                "do not combine with --bin-start/--bin-stop"
            ),
        )
        target.add_argument(
            "--flat-index-file",
            type=Path,
            help=(
                "file containing one flat index per line; comments and "
                "blank lines are ignored, and it cannot be combined with "
                "other selection options"
            ),
        )
        target.add_argument("--apply-y-max", action="store_true")
        target.add_argument("--heartbeat-interval", type=int, default=100000)
        target.add_argument("--generator-revision", default="UNKNOWN")
        target.add_argument("--overwrite", action="store_true")

    prepare_parser = subparsers.add_parser(
        "prepare", help="write one radiative mode-4 input per stratum"
    )
    add_prepare_common(prepare_parser)
    envelope_group = prepare_parser.add_mutually_exclusive_group(
        required=True
    )
    envelope_group.add_argument(
        "--sigr-max",
        type=float,
        help="one shared generation envelope for every selected stratum",
    )
    envelope_group.add_argument(
        "--envelope-report",
        type=Path,
        help=(
            "finalized calibration report supplying one audited envelope "
            "per selected stratum"
        ),
    )
    prepare_parser.add_argument(
        "--allow-envelope-revision-mismatch",
        action="store_true",
        help=(
            "use an envelope report from a different generator revision "
            "after an explicit compatibility audit"
        ),
    )
    prepare_parser.add_argument(
        "--envelope-revision-compatibility-rationale",
        help=(
            "required audit note when the envelope calibration and "
            "generation revisions differ"
        ),
    )
    prepare_parser.add_argument("--events-per-stratum", type=int, default=5000)
    calibration_parser = subparsers.add_parser(
        "prepare-calibration",
        help="write fixed-trial mode-4 envelope-calibration inputs",
    )
    add_prepare_common(calibration_parser)
    calibration_parser.add_argument("--trials", type=int, required=True)
    calibration_parser.add_argument(
        "--inside-guard-trial-fraction",
        "--calibration-core-fraction",
        dest="calibration_inside_guard_fraction",
        type=float,
        default=0.5,
        help=(
            "fraction of calibration trials drawn uniformly inside the "
            "guard; the remainder are drawn from its exact complement"
        ),
    )

    run_parser = subparsers.add_parser(
        "run", help="execute and validate one prepared mode-4 stratum"
    )
    run_parser.add_argument("manifest", type=Path)
    run_parser.add_argument("--flat-index", type=int, required=True)
    run_parser.add_argument("--replica-index", type=int, default=0)
    run_parser.add_argument("--executable", type=Path, required=True)
    run_parser.add_argument("--overwrite", action="store_true")

    finalize_parser = subparsers.add_parser(
        "finalize",
        help=(
            "pool generation replicas, or compatible calibration campaigns, "
            "by stratum"
        ),
    )
    finalize_parser.add_argument("manifests", type=Path, nargs="+")
    finalize_parser.add_argument("--output", type=Path)
    finalize_parser.add_argument(
        "--envelope-safety-factor", type=float, default=1.2
    )
    finalize_parser.add_argument(
        "--maximum-duplicate-fraction", type=float, default=0.05
    )
    finalize_parser.add_argument(
        "--minimum-component-targets", type=int, default=20
    )
    finalize_parser.add_argument(
        "--minimum-provisional-inside-targets",
        type=int,
        default=1000,
        help=(
            "minimum observed inside-guard targets required before a "
            "zero-complement envelope may be used for a pilot"
        ),
    )
    finalize_parser.add_argument(
        "--additional-inside-pilot-run",
        type=Path,
        action="append",
        help=(
            "validated inside-only mode-4 generation run JSON whose "
            "observed maximum raises the provisional envelope floor; "
            "repeat for more pilots"
        ),
    )
    finalize_parser.add_argument(
        "--allow-zero-complement",
        action="store_true",
        help=(
            "allow a provisional inside-derived envelope after zero observed "
            "complement targets pass the configured occurrence-rate bound"
        ),
    )
    finalize_parser.add_argument(
        "--zero-complement-confidence",
        type=float,
        default=0.95,
        help=(
            "one-sided binomial confidence level for zero complement targets"
        ),
    )
    finalize_parser.add_argument(
        "--maximum-zero-complement-target-rate",
        type=float,
        default=1.0e-6,
        help=(
            "largest allowed upper occurrence-rate bound for a provisional "
            "zero-complement envelope"
        ),
    )
    finalize_parser.add_argument(
        "--allow-calibration-revision-mismatch",
        action="store_true",
        help=(
            "pool calibration manifests with different recorded generator "
            "revisions after an explicit compatibility audit"
        ),
    )
    finalize_parser.add_argument(
        "--revision-compatibility-rationale",
        help=(
            "required audit note explaining why differing generator "
            "revisions have identical calibration physics and proposals"
        ),
    )
    validation_parser = subparsers.add_parser(
        "validate-pilots",
        help=(
            "pool independent generation pilots and decide whether each "
            "calibrated stratum is ready for a broader pilot"
        ),
    )
    validation_parser.add_argument(
        "--calibration", type=Path, required=True
    )
    validation_parser.add_argument(
        "--run",
        dest="runs",
        type=Path,
        action="append",
        required=True,
        help="completed mode-4 generation run JSON; repeat for more runs",
    )
    validation_parser.add_argument("--output", type=Path)
    validation_parser.add_argument("--minimum-runs", type=int, default=2)
    validation_parser.add_argument("--minimum-events", type=int, default=400)
    validation_parser.add_argument(
        "--maximum-duplicate-fraction", type=float, default=0.05
    )
    validation_parser.add_argument(
        "--maximum-guard-complement-fraction",
        type=float,
        default=0.02,
    )
    validation_parser.add_argument(
        "--maximum-relative-cross-section-difference",
        type=float,
        default=0.10,
    )
    validation_parser.add_argument(
        "--maximum-cross-section-z-score", type=float, default=3.0
    )
    validation_parser.add_argument(
        "--confidence", type=float, default=0.95
    )
    validation_parser.add_argument(
        "--allow-pilot-revision-mismatch",
        action="store_true",
        help=(
            "accept differing pilot or calibration generator revisions "
            "after an explicit compatibility audit"
        ),
    )
    validation_parser.add_argument(
        "--revision-compatibility-rationale",
        help=(
            "required audit note when pilot generator revisions differ "
            "from one another or the calibration"
        ),
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.command == "create-refinement":
        result = create_refinement(args)
    elif args.command in ("prepare", "prepare-calibration"):
        result = prepare(args)
    elif args.command == "run":
        result = run(args)
    elif args.command == "validate-pilots":
        result = validate_pilots(args)
    else:
        result = finalize(args)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
