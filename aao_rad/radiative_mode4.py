#!/usr/bin/env python3
"""Prepare, run, validate, and finalize radiative AAO mode-4 strata.

Mode 4 generates one final-LUND analysis stratum per invocation from an exact
mixture of a learned continuous native-coordinate core and the unrestricted
legacy proposal.  The legacy component has strictly positive probability, so
an imperfect learned guard can reduce efficiency but cannot remove physical
support.
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

import radiative_guards
import radiative_survey


MANIFEST_SCHEMA = "aao-rad-mode4-manifest-v1"
RUN_SCHEMA = "aao-rad-mode4-run-v1"
WEIGHTS_SCHEMA = "aao-rad-mode4-weights-v1"
RECIPE_SCHEMA = "aao-rad-continuous-guard-recipes-v1"
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
            upper = min(1.0, float(axis["upper"]) + amount)
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
        ]
    ) + "\n"


def _load_configuration(
    config_path: Path,
    recipes_path: Path,
    input_path: Path,
) -> tuple[dict, dict, str, str, str]:
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
        hashlib.sha256(recipes_raw).hexdigest(),
        legacy,
    )


def prepare(args: argparse.Namespace) -> Path:
    if args.events_per_stratum <= 0:
        raise ValueError("--events-per-stratum must be positive")
    if args.replicas <= 0:
        raise ValueError("--replicas must be positive")
    if args.seed_base <= 0:
        raise ValueError("--seed-base must be positive")
    if not 0.0 < args.core_fraction < 1.0:
        raise ValueError("--core-fraction must lie strictly between zero and one")
    if not math.isfinite(args.sigr_max) or args.sigr_max <= 0.0:
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
    legacy_input, legacy_records = _legacy_mode4_input(
        legacy_source,
        input_path,
        events=args.events_per_stratum,
        sigr_max=args.sigr_max,
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
        box = reconstruct_guard_box(recipe, base_padding)
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
                    "events_requested": args.events_per_stratum,
                    "input_file": str(input_relative),
                    "output_stem": str(output_stem),
                    "guard": box.manifest_record(),
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
        "generator_revision": args.generator_revision,
        "analysis_config_source": str(config_path),
        "analysis_config_sha256": config_sha256,
        "guard_recipes_source": str(recipes_path),
        "guard_recipes_sha256": recipes_sha256,
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
        "events_per_stratum": args.events_per_stratum,
        "replicas_per_stratum": args.replicas,
        "sigr_max": args.sigr_max,
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
    }
    for name, expected in expected_floats.items():
        actual = _norm_float(norm, name)
        if not _close(actual, expected):
            raise Mode4Error(f"{name}={actual}, expected {expected}")
    if _norm_int(norm, "nevent") <= 0 or _norm_int(norm, "ntries") <= 0:
        raise Mode4Error("mode 4 produced no events or proposals")
    if _norm_float(norm, "sig_sum") <= 0.0:
        raise Mode4Error("mode 4 produced a nonpositive stratum integral")
    if _norm_int(norm, "mode4_core_trials") + _norm_int(
        norm, "mode4_legacy_trials"
    ) != _norm_int(norm, "ntries"):
        raise Mode4Error("mode-4 component trial counts do not sum to ntries")


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
        diagnostic_path = work / MODE4_KINEMATICS_FILENAME
        for path in (norm_path, lund_path, diagnostic_path):
            if not path.is_file():
                raise Mode4Error(f"generator did not write {path.name}")
        norm = radiative_survey.parse_norm(norm_path)
        _validate_norm(norm, manifest, record)
        events = _norm_int(norm, "nevent")
        event_summary = _validate_event_diagnostics(
            diagnostic_path, manifest, record, events
        )
        _validate_lund(lund_path, events)
        if event_summary["core_events"] != _norm_int(
            norm, "mode4_core_events"
        ) or event_summary["legacy_events"] != _norm_int(
            norm, "mode4_legacy_events"
        ):
            raise Mode4Error("event diagnostics and norm component counts differ")
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
        products = (
            (".lund", radiative_survey.LUND_FILENAME),
            (".norm", radiative_survey.NORM_FILENAME),
            (".mode4.csv", MODE4_KINEMATICS_FILENAME),
            (".sum", "aao_rad.sum"),
            (".out", "aao_rad.out"),
            (".stdout", "generator.stdout.txt"),
            (".stderr", "generator.stderr.txt"),
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
        "events": events,
        "event_overshoot": events - int(record["events_requested"]),
        "ntries": _norm_int(norm, "ntries"),
        "mcall_max": _norm_int(norm, "mcall_max"),
        "multiplicity_correction_used": _norm_int(norm, "mcall_max") > 1,
        "event_weight_microbarn": sig_sum / events,
        "event_yield_per_proposal": events / _norm_int(norm, "ntries"),
        "core_trials": _norm_int(norm, "mode4_core_trials"),
        "legacy_trials": _norm_int(norm, "mode4_legacy_trials"),
        "target_candidates": _norm_int(norm, "mode4_target_candidates"),
        "core_events": event_summary["core_events"],
        "legacy_events": event_summary["legacy_events"],
        "maximum_event_density_ratio": event_summary[
            "maximum_event_density_ratio"
        ],
    }
    _write_json(run_record_path, run_record)
    return run_record_path


def finalize(args: argparse.Namespace) -> Path:
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError(f"unsupported manifest schema: {manifest.get('schema')}")
    root = manifest_path.parent
    grouped: dict[str, list[tuple[dict, dict]]] = {}
    for record in manifest["runs"]:
        run_path = root / (str(record["output_stem"]) + ".json")
        if not run_path.is_file():
            raise FileNotFoundError(run_path)
        completed = json.loads(run_path.read_text(encoding="utf-8"))
        if completed.get("schema") != RUN_SCHEMA:
            raise Mode4Error(f"{run_path}: unsupported run schema")
        for name in ("stratum_id", "flat_index", "replica_index"):
            if completed[name] != record[name]:
                raise Mode4Error(f"{run_path}: {name} differs from manifest")
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
    prepare_parser = subparsers.add_parser(
        "prepare", help="write one radiative mode-4 input per stratum"
    )
    prepare_parser.add_argument("--config", type=Path, required=True)
    prepare_parser.add_argument("--recipes", type=Path, required=True)
    prepare_parser.add_argument("--input", type=Path, required=True)
    prepare_parser.add_argument("--output", type=Path, required=True)
    prepare_parser.add_argument("--candidate", default="padding_0p035")
    prepare_parser.add_argument("--core-fraction", type=float, default=0.9)
    prepare_parser.add_argument("--sigr-max", type=float, required=True)
    prepare_parser.add_argument("--events-per-stratum", type=int, default=5000)
    prepare_parser.add_argument("--replicas", type=int, default=1)
    prepare_parser.add_argument("--seed-base", type=int, default=481001)
    prepare_parser.add_argument("--bin-start", type=int, default=0)
    prepare_parser.add_argument("--bin-stop", type=int)
    prepare_parser.add_argument("--apply-y-max", action="store_true")
    prepare_parser.add_argument(
        "--generator-revision", default="UNKNOWN"
    )
    prepare_parser.add_argument("--overwrite", action="store_true")

    run_parser = subparsers.add_parser(
        "run", help="execute and validate one prepared mode-4 stratum"
    )
    run_parser.add_argument("manifest", type=Path)
    run_parser.add_argument("--flat-index", type=int, required=True)
    run_parser.add_argument("--replica-index", type=int, default=0)
    run_parser.add_argument("--executable", type=Path, required=True)
    run_parser.add_argument("--overwrite", action="store_true")

    finalize_parser = subparsers.add_parser(
        "finalize", help="pool replicas into one weight per stratum"
    )
    finalize_parser.add_argument("manifest", type=Path)
    finalize_parser.add_argument("--output", type=Path)
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.command == "prepare":
        result = prepare(args)
    elif args.command == "run":
        result = run(args)
    else:
        result = finalize(args)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
