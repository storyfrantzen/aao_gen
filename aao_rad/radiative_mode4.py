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
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import radiative_guards
import radiative_survey


MANIFEST_SCHEMA = "aao-rad-mode4-manifest-v3"
RUN_SCHEMA = "aao-rad-mode4-run-v3"
WEIGHTS_SCHEMA = "aao-rad-mode4-weights-v1"
CALIBRATION_SCHEMA = "aao-rad-mode4-envelope-calibration-v2"
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
MODE4_HEARTBEAT_SCHEMA = "aao-rad-mode4-heartbeat-v2"
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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


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
    if not 0.0 < core_fraction < 1.0:
        raise ValueError("core fraction must lie strictly between zero and one")
    tail_fraction = 1.0 - core_fraction
    denominator = tail_fraction
    if box.contains(coordinates):
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
    sigr_max = (
        float(args.sigr_max) if operation == "generation" else 1.0
    )
    if not math.isfinite(sigr_max) or sigr_max <= 0.0:
        raise ValueError("--sigr-max must be finite and positive")

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
    legacy_input, legacy_records = _legacy_mode4_input(
        legacy_source,
        input_path,
        events=events_per_stratum if operation == "generation" else 1,
        sigr_max=sigr_max,
    )
    output = args.output.resolve()
    manifest_path = output / "manifest.json"
    if output.exists() and any(output.iterdir()) and not args.overwrite:
        raise FileExistsError(
            f"{output} is not empty; pass --overwrite to replace prepared files"
        )
    input_directory = output / "inputs"
    input_directory.mkdir(parents=True, exist_ok=True)
    shutil.copy2(config_path, output / "analysis_config.json")
    shutil.copy2(recipes_path, output / "continuous_guard_recipes.json")
    if refinements_path is not None:
        shutil.copy2(
            refinements_path, output / "guard_refinements.json"
        )
    (output / "legacy_input.inp").write_text(
        legacy_source, encoding="utf-8"
    )

    strata = radiative_guards.enumerate_strata(config)
    stop = len(strata) if args.bin_stop is None else args.bin_stop
    if args.bin_start < 0 or stop < args.bin_start or stop > len(strata):
        raise ValueError(f"invalid stratum range [{args.bin_start},{stop})")
    records: list[dict[str, object]] = []
    for stratum in strata:
        if not args.bin_start <= stratum.flat_index < stop:
            continue
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
        if operation == "calibration" and box.volume >= 1.0 - 1.0e-7:
            raise ValueError(
                f"{stratum.identifier}: the anchored guard fills the native "
                "hypercube, so its complement cannot be calibrated"
            )
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
            (output / input_relative).write_text(text, encoding="utf-8")
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
                    "events_requested": (
                        events_per_stratum if operation == "generation" else 0
                    ),
                    "trials_requested": (
                        calibration_trials if operation == "calibration" else 0
                    ),
                    "input_file": str(input_relative),
                    "output_stem": str(output_stem),
                    "guard_original": original_box.manifest_record(),
                    "guard": box.manifest_record(),
                    "guard_refinement": refinement_record,
                }
            )
    if not records:
        raise ValueError("selected stratum range produced no runs")
    manifest: dict[str, object] = {
        "schema": MANIFEST_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "development_only": True,
        "production_ready": False,
        "sampling_mode": SAMPLING_MODE,
        "operation": operation,
        "generator_revision": args.generator_revision,
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
        "sigr_max": sigr_max if operation == "generation" else None,
        "analysis_selection": {
            "coordinate_definition": "final_lund_analysis",
            "w_minimum": float(config["phase_space"]["W_min"]),
            "apply_y_max": bool(args.apply_y_max),
            "y_maximum": (
                float(config["phase_space"].get("y_max", 1.0))
                if args.apply_y_max
                else None
            ),
            "no_implicit_y_minimum": True,
        },
        "legacy_input_records": legacy_records,
        "runs": records,
    }
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
        expected_ratio = proposal_density_ratio(
            proposal_coordinates,
            guard_box,
            float(manifest["core_fraction"]),
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
        expected_inside = int(box.contains(coordinates))
        if (
            inside != expected_inside
            or (component == 1 and inside != 1)
            or (component == 0 and inside != 0)
        ):
            raise Mode4Error(f"{path}: inconsistent learned-core membership")
        expected_ratio = proposal_density_ratio(coordinates, box, alpha)
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
    if manifest.get("schema") != MANIFEST_SCHEMA:
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
    input_path = root / record["input_file"]
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
        "sig_sum_microbarn": sig_sum,
        "sig_int_microbarn": _norm_float(norm, "sig_int"),
        "operation": manifest["operation"],
        "events": events,
        "event_overshoot": events - int(record["events_requested"]),
        "ntries": _norm_int(norm, "ntries"),
        "mcall_max": _norm_int(norm, "mcall_max"),
        "multiplicity_correction_used": _norm_int(norm, "mcall_max") > 1,
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
    if safety_factor < 1.0:
        raise ValueError("--envelope-safety-factor must be at least one")
    if not 0.0 <= maximum_duplicate_fraction < 1.0:
        raise ValueError("--maximum-duplicate-fraction must lie in [0,1)")
    if minimum_targets < 1:
        raise ValueError("--minimum-component-targets must be positive")
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
        acceptable = [
            item
            for item in evaluations
            if float(item["expected_duplicate_event_fraction"])
            <= maximum_duplicate_fraction
        ]
        recommendation = (
            acceptable[0] if enough_support and acceptable else None
        )
        strata.append(
            {
                "stratum_id": stratum_id,
                "flat_index": first["flat_index"],
                "indices": first["indices"],
                "bounds": first["bounds"],
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
                "recommendation_status": (
                    "recommended"
                    if recommendation is not None
                    else (
                        "insufficient_both_component_target_support"
                        if not inside_support and not complement_support
                        else "insufficient_inside_guard_target_support"
                        if not inside_support
                        else "insufficient_guard_complement_target_support"
                        if not complement_support
                        else "no_candidate_meets_duplicate_limit"
                    )
                ),
                "recommended_envelope": recommendation,
            }
        )
    payload: dict[str, object] = {
        "schema": CALIBRATION_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
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
            }
            for path, source in sources
        ],
        "generator_revision": manifest["generator_revision"],
        "analysis_config_sha256": manifest["analysis_config_sha256"],
        "guard_recipes_sha256": manifest["guard_recipes_sha256"],
        "guard_refinements_sha256": manifest.get(
            "guard_refinements_sha256"
        ),
        "guard_candidate": manifest["guard_candidate"],
        "core_fraction": alpha,
        "calibration_proposal": "guard_partition",
        "envelope_safety_factor": safety_factor,
        "maximum_duplicate_fraction": maximum_duplicate_fraction,
        "minimum_component_targets": minimum_targets,
        "stratum_count": len(strata),
        "strata": strata,
    }
    output = (
        args.output.resolve()
        if args.output is not None
        else root / "envelope_calibration.json"
    )
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
                "recommended_sigr_max",
                "expected_events_per_proposal",
                "expected_duplicate_event_fraction",
                "inside_guard_targets",
                "guard_complement_targets",
                "estimated_guard_complement_cross_section_fraction",
            ]
        )
        for stratum in strata:
            recommendation = stratum["recommended_envelope"]
            writer.writerow(
                [
                    stratum["flat_index"],
                    stratum["stratum_id"],
                    stratum["recommendation_status"],
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
                    stratum["guard_complement"]["target_candidates"],
                    stratum[
                        "estimated_guard_complement_cross_section_fraction"
                    ],
                ]
            )
    return output


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
        if manifest.get("schema") != MANIFEST_SCHEMA:
            raise ValueError(
                f"{manifest_path}: unsupported manifest schema "
                f"{manifest.get('schema')}"
            )
        sources.append((manifest_path, manifest))
    manifest_path, manifest = sources[0]
    operation = manifest["operation"]
    compatibility_keys = (
        "operation",
        "generator_revision",
        "analysis_config_sha256",
        "guard_recipes_sha256",
        "guard_refinements_sha256",
        "guard_candidate",
        "base_padding",
        "core_fraction",
        "analysis_selection",
        "calibration_proposal",
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
                stratum_seed = (
                    str(record["stratum_id"]),
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
                    record["stratum_id"], []
                ).append((source_root, record, completed, source_path))
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
        grouped.setdefault(record["stratum_id"], []).append(
            (record, completed)
        )
    strata: list[dict[str, object]] = []
    for stratum_id, items in sorted(
        grouped.items(), key=lambda item: int(item[1][0][0]["flat_index"])
    ):
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
                "replicas": len(items),
                "total_proposals": total_proposals,
                "total_events": total_events,
                "combined_sig_sum_microbarn": combined_sigma,
                "pooled_event_weight_microbarn": (
                    combined_sigma / total_events
                ),
                "event_yield_per_proposal": total_events / total_proposals,
                "mcall_max": max(int(run["mcall_max"]) for _, run in items),
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
        "source_manifest": str(manifest_path),
        "source_manifest_sha256": _sha256(manifest_path),
        "generator_revision": manifest["generator_revision"],
        "analysis_config_sha256": manifest["analysis_config_sha256"],
        "guard_recipes_sha256": manifest["guard_recipes_sha256"],
        "guard_refinements_sha256": manifest.get(
            "guard_refinements_sha256"
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
                "combined_sig_sum_microbarn",
                "pooled_event_weight_microbarn",
                "event_yield_per_proposal",
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
                        "combined_sig_sum_microbarn",
                        "pooled_event_weight_microbarn",
                        "event_yield_per_proposal",
                        "core_events",
                        "legacy_events",
                    )
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
        target.add_argument("--apply-y-max", action="store_true")
        target.add_argument("--heartbeat-interval", type=int, default=100000)
        target.add_argument("--generator-revision", default="UNKNOWN")
        target.add_argument("--overwrite", action="store_true")

    prepare_parser = subparsers.add_parser(
        "prepare", help="write one radiative mode-4 input per stratum"
    )
    add_prepare_common(prepare_parser)
    prepare_parser.add_argument("--sigr-max", type=float, required=True)
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
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.command == "create-refinement":
        result = create_refinement(args)
    elif args.command in ("prepare", "prepare-calibration"):
        result = prepare(args)
    elif args.command == "run":
        result = run(args)
    else:
        result = finalize(args)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
