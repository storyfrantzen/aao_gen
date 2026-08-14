#!/usr/bin/env python3
"""Calibrate, generate, and validate global radiative AAO mode-3 samples.

Mode 3 mixes a direct proposal in ``(1/Q_l^2, x_l, -t_h, phi_h)`` with
AAO's unrestricted legacy radiative proposal.  It emits equal-weight LUND
events in one padded final-LUND analysis volume.  This driver freezes inputs,
splits large productions into immutable jobs, validates every artifact, and
pools run normalization without loading LUND files into memory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import socket
import stat
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path


MANIFEST_SCHEMA = "aao-rad-mode3-campaign-v1"
TASK_SCHEMA = "aao-rad-mode3-task-v1"
RUN_SCHEMA = "aao-rad-mode3-run-v1"
CALIBRATION_SCHEMA = "aao-rad-mode3-envelope-v1"
WEIGHTS_SCHEMA = "aao-rad-mode3-weights-v1"
MODE3_SCHEMA = "aao-rad-mode3-v1"
CALIBRATION_COLUMNS = [
    "trial",
    "proposal_component",
    "inside_direct",
    "proposal_density_ratio",
    "integrand_corrected",
    "xb_leptonic",
    "minus_t_hard",
    "phi_cm_deg",
    "q2_observed",
    "xb_observed",
    "minus_t_observed",
    "phi_observed_deg",
    "w_observed",
    "y_observed",
]


class Mode3Error(RuntimeError):
    """Raised when a mode-3 campaign or generated artifact is inconsistent."""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _load_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise Mode3Error(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise Mode3Error(f"{path}: expected a JSON object")
    return value


def _records(text: str) -> list[str]:
    return [
        line.split("!", 1)[0].strip()
        for line in text.splitlines()
        if line.split("!", 1)[0].strip()
    ]


def _legacy_records(path: Path) -> list[str]:
    records = _records(path.read_text(encoding="utf-8"))
    if len(records) < 17:
        raise Mode3Error(f"{path}: expected at least 17 legacy AAO records")
    try:
        theory = int(records[0].split()[0])
        fmcall = float(records[16].split()[0])
    except (ValueError, IndexError) as error:
        raise Mode3Error(f"{path}: cannot parse theory/fmcall") from error
    count = 17 + int(fmcall == 0.0) + int(theory > 10)
    if len(records) != count:
        raise Mode3Error(
            f"{path}: pass a legacy input without a sampling-mode trailer"
        )
    return records


def _strict_edges(config: dict[str, object], name: str) -> list[float]:
    try:
        values = config["binning"][name]  # type: ignore[index]
        edges = [float(value) for value in values]  # type: ignore[arg-type]
    except (KeyError, TypeError, ValueError) as error:
        raise Mode3Error(f"analysis config lacks numeric binning.{name}") from error
    if (
        len(edges) < 2
        or any(not math.isfinite(value) for value in edges)
        or any(right <= left for left, right in zip(edges, edges[1:]))
    ):
        raise Mode3Error(f"binning.{name} must contain increasing finite edges")
    return edges


def _padded_range(
    edges: list[float], fraction: float, *, nonnegative: bool = False
) -> tuple[float, float]:
    width = edges[-1] - edges[0]
    lower = edges[0] - fraction * width
    upper = edges[-1] + fraction * width
    if nonnegative:
        lower = max(0.0, lower)
    return lower, upper


def _configuration(path: Path, padding_fraction: float) -> dict[str, object]:
    config = _load_json(path)
    if not 0.0 <= padding_fraction <= 0.5:
        raise ValueError("--padding-fraction must lie in [0, 0.5]")
    try:
        beam = float(config["beam_energy"])
        mass = float(config["target_mass"])
        phase = config["phase_space"]
        w_min = float(phase["W_min"])  # type: ignore[index]
        electron_p_min = float(phase["electron_p_min"])  # type: ignore[index]
    except (KeyError, TypeError, ValueError) as error:
        raise Mode3Error(
            "mode 3 requires beam_energy, target_mass, phase_space.W_min, "
            "and phase_space.electron_p_min"
        ) from error
    if beam <= 0.0 or mass <= 0.0 or w_min <= 0.0 or electron_p_min <= 0.0:
        raise Mode3Error("beam, mass, W minimum, and electron momentum must be positive")
    q2 = _padded_range(_strict_edges(config, "Q2"), padding_fraction)
    xb = _padded_range(_strict_edges(config, "xB"), padding_fraction)
    minus_t = _padded_range(
        _strict_edges(config, "minus_t"), padding_fraction, nonnegative=True
    )
    phi_edges = _strict_edges(config, "phi_deg")
    if not (
        math.isclose(phi_edges[0], 0.0, abs_tol=1.0e-9)
        and math.isclose(phi_edges[-1], 360.0, abs_tol=1.0e-9)
    ):
        raise Mode3Error("mode 3 currently requires phi binning spanning [0,360]")
    y_value = phase.get("y_max")  # type: ignore[union-attr]
    y_max = None if y_value is None else float(y_value)
    if y_max is not None and not 0.0 < y_max <= 1.0:
        raise Mode3Error("phase_space.y_max must be null or lie in (0,1]")
    if xb[0] <= 0.0:
        raise Mode3Error(
            "padding makes the xB lower bound nonpositive; reduce the padding"
        )
    return {
        "beam_energy": beam,
        "target_mass": mass,
        "padding_fraction": padding_fraction,
        "bounds": {
            "Q2": list(q2),
            "xB": list(xb),
            "minus_t": list(minus_t),
            "phi_deg": [0.0, 360.0],
        },
        "W_min": w_min,
        "electron_p_min": electron_p_min,
        "apply_y_max": y_max is not None,
        "y_max": 1.0 if y_max is None else y_max,
    }


def proposal_density_ratio(
    *,
    q2: float,
    xb: float,
    minus_t: float,
    phi_deg: float,
    t_jacobian: float,
    ep_range: float,
    direct_fraction: float,
    xb_bounds: tuple[float, float],
    t_bounds: tuple[float, float],
    phi_bounds: tuple[float, float],
    target_mass: float,
) -> float:
    """Python parity implementation of the exact Fortran mixture correction."""

    legacy_fraction = 1.0 - direct_fraction
    if not 0.0 <= direct_fraction < 1.0:
        raise ValueError("direct_fraction must lie in [0,1)")
    inside = (
        q2 > 0.0
        and xb > 0.0
        and t_jacobian > 0.0
        and xb_bounds[0] <= xb <= xb_bounds[1]
        and t_bounds[0] <= minus_t <= t_bounds[1]
        and phi_bounds[0] <= phi_deg % 360.0 <= phi_bounds[1]
    )
    if not inside:
        return 1.0 / legacy_fraction
    direct_to_legacy = (
        ep_range
        * 2.0
        * target_mass
        * xb**2
        / (q2 * (xb_bounds[1] - xb_bounds[0]))
        * 2.0
        * t_jacobian
        / (t_bounds[1] - t_bounds[0])
        * 360.0
        / (phi_bounds[1] - phi_bounds[0])
    )
    return 1.0 / (legacy_fraction + direct_fraction * direct_to_legacy)


def _canonical_legacy(
    records: list[str], settings: dict[str, object], events: int, envelope: float
) -> str:
    theory = int(records[0].split()[0])
    original_fmcall = float(records[16].split()[0])
    position = 17
    if original_fmcall == 0.0:
        position += 1
    w_record = records[position] if theory > 10 else None
    base = records[:17]
    bounds = settings["bounds"]  # type: ignore[assignment]
    base[3] = "4"
    base[4] = "1"
    base[11] = f"{float(settings['beam_energy']):.17g}"
    base[12] = " ".join(f"{value:.17g}" for value in bounds["Q2"])  # type: ignore[index]
    base[13] = (
        f"{float(settings['electron_p_min']):.17g} "
        f"{float(settings['beam_energy']):.17g}"
    )
    base[15] = str(events)
    base[16] = "0"
    base.append(f"{envelope:.17g}")
    if w_record is not None:
        base.append(w_record)
    return "\n".join(base) + "\n"


def _legacy_physics_settings(records: list[str]) -> dict[str, object]:
    """Return the non-campaign physics settings inherited from the template."""

    try:
        theory = int(records[0].split()[0])
        helicity = int(records[1].split()[0])
        regions = [float(value) for value in records[2].split()]
        missing_mass_cut = float(records[5].split()[0])
        target_length = float(records[6].split()[0])
        target_radius = float(records[7].split()[0])
        vertex = [float(records[index].split()[0]) for index in range(8, 11)]
        minimum_photon_energy = float(records[14].split()[0])
        original_fmcall = float(records[16].split()[0])
    except (ValueError, IndexError) as error:
        raise Mode3Error("legacy input contains malformed physics settings") from error
    if len(regions) != 4 or any(value <= 0.0 for value in regions):
        raise Mode3Error("legacy integration-region record must contain 4 positives")
    if minimum_photon_energy <= 0.0:
        raise Mode3Error("legacy minimum photon energy must be positive")
    position = 17 + int(original_fmcall == 0.0)
    theory_record = records[position].split() if theory > 10 else []
    return {
        "physics_model": theory,
        "electron_helicity_flag": helicity,
        "integration_regions": regions,
        "missing_mass_squared_cut": missing_mass_cut,
        "target_length_cm": target_length,
        "target_radius_cm": target_radius,
        "vertex_cm": vertex,
        "minimum_photon_energy_gev": minimum_photon_energy,
        "optional_theory_record": [float(value) for value in theory_record],
    }


def _trailer(
    settings: dict[str, object], *, seed: int, replica: int,
    direct_fraction: float, operation: int, trials: int, heartbeat: int
) -> str:
    bounds = settings["bounds"]  # type: ignore[assignment]
    pair = lambda name: " ".join(  # noqa: E731
        f"{float(value):.17g}" for value in bounds[name]  # type: ignore[index]
    )
    return "\n".join(
        [
            "3",
            str(seed),
            str(replica),
            f"{direct_fraction:.17g}",
            pair("xB"),
            pair("minus_t"),
            pair("phi_deg"),
            pair("Q2"),
            pair("xB"),
            pair("minus_t"),
            pair("phi_deg"),
            f"{float(settings['W_min']):.17g}",
            f"{int(bool(settings['apply_y_max']))} {float(settings['y_max']):.17g}",
            str(operation),
            str(trials),
            str(heartbeat),
        ]
    ) + "\n"


def _prepare(args: argparse.Namespace, *, operation: str) -> Path:
    output = args.output.expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    config_path = args.config.expanduser().resolve()
    input_path = args.input.expanduser().resolve()
    if not config_path.is_file() or not input_path.is_file():
        raise FileNotFoundError(config_path if not config_path.is_file() else input_path)
    settings = _configuration(config_path, args.padding_fraction)
    records = _legacy_records(input_path)
    minimum_photon_energy = getattr(args, "minimum_photon_energy", None)
    if minimum_photon_energy is not None:
        if not 0.0 < minimum_photon_energy < float(settings["beam_energy"]):
            raise ValueError(
                "--minimum-photon-energy must lie between zero and the beam energy"
            )
        records = list(records)
        records[14] = f"{minimum_photon_energy:.17g}"
    settings["generator_physics"] = _legacy_physics_settings(records)
    # The legacy input is a physics-model template.  Beam energy, Q2 bounds,
    # and electron momentum bounds are intentionally replaced by the frozen
    # analysis configuration below; requiring the template's old beam value
    # to match would reject the historical 10.6-GeV input for RGA 10.604 GeV.
    if not 0.0 <= args.direct_fraction < 1.0:
        raise ValueError("--direct-fraction must lie in [0,1)")
    if args.replicas <= 0 or args.heartbeat_interval <= 0:
        raise ValueError("replicas and heartbeat interval must be positive")
    if args.seed_base == 0:
        raise ValueError("--seed-base must be nonzero")
    if operation == "calibration":
        if args.trials <= 0:
            raise ValueError("--trials must be positive")
        requested = [0] * args.replicas
        total_events = 0
        envelope = 1.0
    else:
        calibration_path = args.calibration.expanduser().resolve()
        calibration = _load_json(calibration_path)
        if calibration.get("schema") != CALIBRATION_SCHEMA:
            raise Mode3Error(f"{calibration_path}: wrong calibration schema")
        if calibration.get("settings") != settings:
            raise Mode3Error("production settings differ from calibration")
        if not math.isclose(
            float(calibration["direct_fraction"]), args.direct_fraction,
            rel_tol=0.0, abs_tol=1.0e-12,
        ):
            raise Mode3Error("production direct fraction differs from calibration")
        envelope = float(calibration["recommended_sigr_max"])
        if envelope <= 0.0:
            raise Mode3Error("calibration has no positive envelope")
        if args.total_events <= 0 or args.events_per_job <= 0:
            raise ValueError("total events and events per job must be positive")
        if args.lund_files_per_directory <= 0:
            raise ValueError("--lund-files-per-directory must be positive")
        jobs = math.ceil(args.total_events / args.events_per_job)
        requested = [args.events_per_job] * jobs
        requested[-1] = args.total_events - args.events_per_job * (jobs - 1)
        total_events = args.total_events
    frozen = output / "inputs"
    frozen.mkdir()
    config_snapshot = frozen / "analysis_config.json"
    input_snapshot = frozen / "legacy_input.inp"
    shutil.copyfile(config_path, config_snapshot)
    shutil.copyfile(input_path, input_snapshot)
    tag = args.tag.strip()
    if not tag or any(character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-" for character in tag):
        raise ValueError("--tag contains unsupported characters")
    runs: list[dict[str, object]] = []
    for replica, events in enumerate(requested):
        seed = args.seed_base + replica
        input_text = _canonical_legacy(records, settings, max(1, events), envelope)
        input_text += _trailer(
            settings,
            seed=seed,
            replica=replica,
            direct_fraction=args.direct_fraction,
            operation=int(operation == "calibration"),
            trials=args.trials if operation == "calibration" else 0,
            heartbeat=args.heartbeat_interval,
        )
        stem = f"{tag}__g{replica:08d}"
        file_chunk_index = replica // 5_000
        input_chunk = frozen / f"chunk_{file_chunk_index:04d}"
        input_chunk.mkdir(exist_ok=True)
        generated_input = input_chunk / f"{stem}.inp"
        generated_input.write_text(input_text, encoding="utf-8")
        task_path = (
            output
            / "tasks"
            / f"chunk_{file_chunk_index:04d}"
            / f"task_{replica:08d}.json"
        )
        lund_chunk = (
            output
            / "lund"
            / f"chunk_{replica // args.lund_files_per_directory:04d}"
            if operation == "production"
            else None
        )
        runs.append(
            {
                "task_id": replica,
                "replica_index": replica,
                "seed": seed,
                "requested_trials": args.trials if operation == "calibration" else 0,
                "requested_events": events,
                "input": str(generated_input),
                "input_sha256": _sha256(generated_input),
                "run_json": str(output / "runs" / f"{stem}.json"),
                "lund": str(lund_chunk / f"{stem}.lund") if lund_chunk else None,
                "lund_chunk": str(lund_chunk) if lund_chunk else None,
                "task": str(task_path),
                "stem": stem,
            }
        )
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "created_utc": _now(),
        "operation": operation,
        "root": str(output),
        "tag": tag,
        "generator_revision": args.generator_revision,
        "direct_fraction": args.direct_fraction,
        "legacy_fraction": 1.0 - args.direct_fraction,
        "settings": settings,
        "source_config": str(config_path),
        "source_config_sha256": _sha256(config_path),
        "source_input": str(input_path),
        "source_input_sha256": _sha256(input_path),
        "frozen_config": str(config_snapshot),
        "frozen_config_sha256": _sha256(config_snapshot),
        "frozen_input": str(input_snapshot),
        "frozen_input_sha256": _sha256(input_snapshot),
        "calibration": str(args.calibration.expanduser().resolve()) if operation == "production" else None,
        "calibration_sha256": _sha256(args.calibration.expanduser().resolve()) if operation == "production" else None,
        "recommended_sigr_max": envelope,
        "trials_per_replica": args.trials if operation == "calibration" else 0,
        "total_events_requested": total_events,
        "events_per_job": args.events_per_job if operation == "production" else 0,
        "lund_files_per_directory": (
            args.lund_files_per_directory if operation == "production" else 0
        ),
        "task_count": len(runs),
        "runs": runs,
    }
    manifest_path = output / "manifest.json"
    _write_json(manifest_path, manifest)
    manifest_hash = _sha256(manifest_path)
    for run_record in runs:
        task_path = Path(str(run_record["task"]))
        _write_json(
            task_path,
            {
                "schema": TASK_SCHEMA,
                "source_manifest": str(manifest_path),
                "source_manifest_sha256": manifest_hash,
                "operation": operation,
                "root": str(output),
                "run": run_record,
            },
        )
    with (output / "tasks.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("task_id", "replica_index", "seed", "requested_trials", "requested_events", "expected_run", "expected_lund"))
        for run in runs:
            writer.writerow((run["task_id"], run["replica_index"], run["seed"], run["requested_trials"], run["requested_events"], run["run_json"], run["lund"] or ""))
    return manifest_path


def prepare_calibration(args: argparse.Namespace) -> Path:
    return _prepare(args, operation="calibration")


def prepare_production(args: argparse.Namespace) -> Path:
    return _prepare(args, operation="production")


def _parse_norm(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        if "=" not in line:
            raise Mode3Error(f"{path}:{number}: expected key=value")
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def _number(values: dict[str, str], key: str) -> float:
    try:
        return float(values[key].replace("D", "E").replace("d", "e"))
    except (KeyError, ValueError) as error:
        raise Mode3Error(f"normalization sidecar lacks numeric {key}") from error


def _run_record(manifest: dict[str, object], replica: int) -> dict[str, object]:
    matches = [run for run in manifest["runs"] if int(run["replica_index"]) == replica]  # type: ignore[index,union-attr]
    if len(matches) != 1:
        raise Mode3Error(f"manifest has {len(matches)} records for replica {replica}")
    return matches[0]


def _count_lines(path: Path) -> int:
    count = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            count += block.count(b"\n")
    return count


def _calibration_moments(path: Path) -> tuple[int, float, float, float, float]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != "# schema=aao-rad-mode3-calibration-v1":
        raise Mode3Error(f"{path}: wrong calibration CSV schema")
    reader = csv.DictReader(lines[1:])
    if reader.fieldnames != CALIBRATION_COLUMNS:
        raise Mode3Error(f"{path}: unexpected calibration columns")
    rows = 0
    total = square = maximum = maximum_ratio = 0.0
    for raw in reader:
        try:
            value = float(raw["integrand_corrected"])
            ratio = float(raw["proposal_density_ratio"])
            component = int(raw["proposal_component"])
            inside = int(raw["inside_direct"])
        except (TypeError, ValueError) as error:
            raise Mode3Error(f"{path}: malformed calibration row") from error
        if value <= 0.0 or ratio <= 0.0 or component not in (0, 1) or inside not in (0, 1):
            raise Mode3Error(f"{path}: invalid calibration row domain")
        rows += 1
        total += value
        square += value * value
        maximum = max(maximum, value)
        maximum_ratio = max(maximum_ratio, ratio)
    return rows, total, square, maximum, maximum_ratio


def run(args: argparse.Namespace) -> Path:
    campaign_path = args.manifest.expanduser().resolve()
    campaign = _load_json(campaign_path)
    if campaign.get("schema") == TASK_SCHEMA:
        record = campaign["run"]
        if int(record["replica_index"]) != args.replica_index:  # type: ignore[index]
            raise Mode3Error("task record and requested replica differ")
        manifest_path = Path(str(campaign["source_manifest"]))
        manifest_hash = str(campaign["source_manifest_sha256"])
        manifest = {
            "operation": campaign["operation"],
            "root": campaign["root"],
        }
    elif campaign.get("schema") == MANIFEST_SCHEMA:
        manifest_path = campaign_path
        manifest_hash = _sha256(manifest_path)
        manifest = campaign
        record = _run_record(manifest, args.replica_index)
    else:
        raise Mode3Error(f"{campaign_path}: wrong campaign/task schema")
    run_path = Path(str(record["run_json"]))
    if run_path.exists() and not args.overwrite:
        existing = _load_json(run_path)
        if existing.get("source_manifest_sha256") == manifest_hash:
            return run_path
        raise FileExistsError(run_path)
    executable = args.executable.expanduser().resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise FileNotFoundError(executable)
    input_path = Path(str(record["input"]))
    if _sha256(input_path) != record["input_sha256"]:
        raise Mode3Error(f"input hash changed: {input_path}")
    root = Path(str(manifest["root"]))
    scratch_parent = args.scratch_root.expanduser().resolve() if args.scratch_root else None
    if scratch_parent is not None:
        scratch_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"{record['stem']}_", dir=scratch_parent) as temporary:
        work = Path(temporary)
        completed = subprocess.run(
            [str(executable)],
            input=input_path.read_text(encoding="utf-8"),
            text=True,
            cwd=work,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0:
            raise Mode3Error(
                f"AAO failed with exit code {completed.returncode}:\n"
                f"{completed.stderr[-4000:]}"
            )
        norm_path = work / "aao_rad.norm"
        if not norm_path.is_file():
            raise Mode3Error("AAO produced no normalization sidecar")
        norm = _parse_norm(norm_path)
        if norm.get("mode3_schema") != MODE3_SCHEMA or int(_number(norm, "sampling_mode")) != 3:
            raise Mode3Error("AAO output is not mode 3")
        if int(_number(norm, "mode3_replica")) != args.replica_index:
            raise Mode3Error("AAO replica differs from manifest")
        proposals = int(_number(norm, "ntries"))
        events = int(_number(norm, "events"))
        phase_volume = _number(norm, "mode3_phase_volume")
        operation = str(manifest["operation"])
        payload: dict[str, object] = {
            "schema": RUN_SCHEMA,
            "created_utc": _now(),
            "host": socket.gethostname(),
            "operation": operation,
            "replica_index": args.replica_index,
            "seed": int(record["seed"]),
            "source_manifest": str(manifest_path),
            "source_manifest_sha256": manifest_hash,
            "input": str(input_path),
            "input_sha256": _sha256(input_path),
            "generator_executable": str(executable),
            "generator_executable_sha256": _sha256(executable),
            "ntries": proposals,
            "events": events,
            "sig_int_microbarn": _number(norm, "sig_int"),
            "sig_sum_microbarn": _number(norm, "sig_sum"),
            "sigr_max": _number(norm, "sigr_max"),
            "mcall_max": int(_number(norm, "mcall_max")),
            "direct_trials": int(_number(norm, "mode3_direct_trials")),
            "legacy_trials": int(_number(norm, "mode3_legacy_trials")),
            "internal_rows": int(_number(norm, "mode3_internal_rows")),
            "final_candidates": int(_number(norm, "mode3_final_candidates")),
            "target_candidates": int(_number(norm, "mode3_target_candidates")),
            "direct_events": int(_number(norm, "mode3_direct_events")),
            "legacy_events": int(_number(norm, "mode3_legacy_events")),
            "emitting_candidates": int(_number(norm, "mode3_emitting_candidates")),
            "duplicate_events": int(_number(norm, "mode3_duplicate_events")),
            "phase_volume": phase_volume,
        }
        run_path.parent.mkdir(parents=True, exist_ok=True)
        diagnostics = run_path.with_suffix("")
        (diagnostics.parent / f"{diagnostics.name}.stdout").write_text(completed.stdout, encoding="utf-8")
        (diagnostics.parent / f"{diagnostics.name}.stderr").write_text(completed.stderr, encoding="utf-8")
        for source_name, suffix in (
            ("aao_rad.norm", ".norm"),
            ("aao_rad.sum", ".sum"),
            ("aao_rad.out", ".out"),
            ("aao_rad.mode3.heartbeat.csv", ".heartbeat.csv"),
        ):
            source = work / source_name
            if source.exists():
                shutil.copyfile(source, diagnostics.parent / f"{diagnostics.name}{suffix}")
        if operation == "calibration":
            expected = int(record["requested_trials"])
            if proposals != expected or events != 0:
                raise Mode3Error("calibration proposal/event totals disagree")
            csv_path = work / "aao_rad.mode3.csv"
            rows, total, square, maximum, maximum_ratio = _calibration_moments(csv_path)
            estimate = phase_volume * total / proposals
            variance_sum = max(0.0, square - total * total / proposals)
            sem = phase_volume * math.sqrt(variance_sum / (proposals - 1) / proposals) if proposals > 1 else 0.0
            if not math.isclose(estimate, float(payload["sig_sum_microbarn"]), rel_tol=4.0e-6, abs_tol=1.0e-18):
                raise Mode3Error("calibration CSV and normalization disagree")
            if (work / "aao_rad.lund").stat().st_size != 0:
                raise Mode3Error("mode-3 calibration emitted LUND data")
            shutil.copyfile(csv_path, diagnostics.parent / f"{diagnostics.name}.calibration.csv")
            payload.update(
                {
                    "calibration_rows": rows,
                    "integrand_sum": total,
                    "integrand_square_sum": square,
                    "maximum_integrand": maximum,
                    "maximum_proposal_density_ratio": maximum_ratio,
                    "cross_section_sem_microbarn": sem,
                    "effective_sample_size": total * total / square if square > 0.0 else 0.0,
                }
            )
        else:
            requested = int(record["requested_events"])
            lund_source = work / "aao_rad.lund"
            lines = _count_lines(lund_source)
            if events != requested or lines != 5 * events:
                raise Mode3Error(
                    f"production expected {requested} events and {5 * events} "
                    f"LUND lines; found {events} events and {lines} lines"
                )
            lund_target = Path(str(record["lund"]))
            lund_target.parent.mkdir(parents=True, exist_ok=True)
            if lund_target.exists() and not args.overwrite:
                raise FileExistsError(lund_target)
            shutil.move(str(lund_source), lund_target)
            payload.update(
                {
                    "requested_events": requested,
                    "event_overshoot": 0,
                    "event_yield_per_proposal": events / proposals,
                    "duplicate_event_fraction": int(payload["duplicate_events"]) / events,
                    "lund": str(lund_target),
                    "lund_bytes": lund_target.stat().st_size,
                    "lund_lines": lines,
                    "event_weight_microbarn": float(payload["sig_sum_microbarn"]) / events,
                }
            )
        _write_json(run_path, payload)
    return run_path


def status(args: argparse.Namespace) -> Path:
    manifest_path = args.manifest.expanduser().resolve()
    manifest = _load_json(manifest_path)
    counts = {"complete": 0, "failed": 0, "pending": 0}
    rows = []
    expected_hash = _sha256(manifest_path)
    for record in manifest["runs"]:  # type: ignore[union-attr]
        path = Path(str(record["run_json"]))
        state = "pending"
        message = ""
        if path.exists():
            try:
                run_record = _load_json(path)
                if run_record.get("schema") != RUN_SCHEMA or run_record.get("source_manifest_sha256") != expected_hash:
                    raise Mode3Error("run provenance differs")
                state = "complete"
            except Exception as error:  # status must report, not abort
                state, message = "failed", str(error)
        counts[state] += 1
        rows.append({"task_id": record["task_id"], "replica_index": record["replica_index"], "state": state, "message": message})
    output = args.output.expanduser().resolve() if args.output else Path(str(manifest["root"])) / "status.json"
    _write_json(output, {"schema": "aao-rad-mode3-status-v1", "created_utc": _now(), "manifest": str(manifest_path), "status_counts": counts, "complete": counts["complete"] == len(rows), "runs": rows})
    return output


def _complete_runs(manifest_path: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    manifest = _load_json(manifest_path)
    expected_hash = _sha256(manifest_path)
    runs = []
    for record in manifest["runs"]:  # type: ignore[union-attr]
        path = Path(str(record["run_json"]))
        if not path.is_file():
            raise Mode3Error(f"run is incomplete: {path}")
        run_record = _load_json(path)
        if run_record.get("schema") != RUN_SCHEMA or run_record.get("source_manifest_sha256") != expected_hash:
            raise Mode3Error(f"run provenance differs: {path}")
        runs.append(run_record)
    return manifest, runs


def finalize_calibration(args: argparse.Namespace) -> Path:
    manifest_path = args.manifest.expanduser().resolve()
    manifest, runs = _complete_runs(manifest_path)
    if manifest["operation"] != "calibration":
        raise Mode3Error("not a calibration manifest")
    if args.envelope_safety_factor <= 1.0:
        raise ValueError("--envelope-safety-factor must exceed one")
    proposals = sum(int(run["ntries"]) for run in runs)
    total = sum(float(run["integrand_sum"]) for run in runs)
    square = sum(float(run["integrand_square_sum"]) for run in runs)
    maximum = max(float(run["maximum_integrand"]) for run in runs)
    volume = float(runs[0]["phase_volume"])
    cross_section = volume * total / proposals
    variance_sum = max(0.0, square - total * total / proposals)
    sem = volume * math.sqrt(variance_sum / (proposals - 1) / proposals)
    report = {
        "schema": CALIBRATION_SCHEMA,
        "created_utc": _now(),
        "source_manifest": str(manifest_path),
        "source_manifest_sha256": _sha256(manifest_path),
        "generator_revision": manifest["generator_revision"],
        "direct_fraction": manifest["direct_fraction"],
        "legacy_fraction": manifest["legacy_fraction"],
        "settings": manifest["settings"],
        "replicas": len(runs),
        "proposals": proposals,
        "target_candidates": sum(int(run["target_candidates"]) for run in runs),
        "integrated_cross_section_microbarn": cross_section,
        "integrated_cross_section_sem_microbarn": sem,
        "effective_sample_size": total * total / square if square > 0.0 else 0.0,
        "observed_maximum_integrand": maximum,
        "envelope_safety_factor": args.envelope_safety_factor,
        "recommended_sigr_max": args.envelope_safety_factor * maximum,
        "maximum_proposal_density_ratio": max(float(run["maximum_proposal_density_ratio"]) for run in runs),
        "runs": runs,
    }
    output = args.output.expanduser().resolve() if args.output else Path(str(manifest["root"])) / "envelope_calibration.json"
    _write_json(output, report)
    return output


def compare_calibrations(args: argparse.Namespace) -> Path:
    candidate_path = args.candidate.expanduser().resolve()
    control_path = args.control.expanduser().resolve()
    candidate = _load_json(candidate_path)
    control = _load_json(control_path)
    for path, report in ((candidate_path, candidate), (control_path, control)):
        if report.get("schema") != CALIBRATION_SCHEMA:
            raise Mode3Error(f"{path}: wrong calibration schema")
    if candidate["settings"] != control["settings"]:
        raise Mode3Error("candidate and control settings differ")
    c = float(candidate["integrated_cross_section_microbarn"])
    r = float(control["integrated_cross_section_microbarn"])
    combined = math.hypot(float(candidate["integrated_cross_section_sem_microbarn"]), float(control["integrated_cross_section_sem_microbarn"]))
    relative = (c - r) / r
    z_score = (c - r) / combined if combined > 0.0 else math.inf
    passed = abs(relative) <= args.maximum_relative_difference and abs(z_score) <= args.maximum_z_score
    output = args.output.expanduser().resolve()
    _write_json(output, {"schema": "aao-rad-mode3-validation-v1", "created_utc": _now(), "passed": passed, "candidate": str(candidate_path), "control": str(control_path), "candidate_cross_section_microbarn": c, "control_cross_section_microbarn": r, "relative_difference": relative, "difference_z_score": z_score, "maximum_relative_difference": args.maximum_relative_difference, "maximum_z_score": args.maximum_z_score})
    return output


def finalize_production(args: argparse.Namespace) -> Path:
    manifest_path = args.manifest.expanduser().resolve()
    manifest, runs = _complete_runs(manifest_path)
    if manifest["operation"] != "production":
        raise Mode3Error("not a production manifest")
    proposals = sum(int(run["ntries"]) for run in runs)
    events = sum(int(run["events"]) for run in runs)
    combined = sum(int(run["ntries"]) * float(run["sig_sum_microbarn"]) for run in runs) / proposals
    report = {
        "schema": WEIGHTS_SCHEMA,
        "created_utc": _now(),
        "source_manifest": str(manifest_path),
        "source_manifest_sha256": _sha256(manifest_path),
        "settings": manifest["settings"],
        "runs": len(runs),
        "requested_events": manifest["total_events_requested"],
        "total_events": events,
        "event_overshoot": events - int(manifest["total_events_requested"]),
        "total_proposals": proposals,
        "combined_sig_sum_microbarn": combined,
        "pooled_event_weight_microbarn": combined / events,
        "duplicate_events": sum(int(run["duplicate_events"]) for run in runs),
        "maximum_mcall": max(int(run["mcall_max"]) for run in runs),
        "lund_directory": str(Path(str(manifest["root"])) / "lund"),
        "lund_files": [run["lund"] for run in runs],
    }
    output = args.output.expanduser().resolve() if args.output else Path(str(manifest["root"])) / "campaign_weights.json"
    _write_json(output, report)
    return output


def emit_swif(args: argparse.Namespace) -> Path:
    manifest_path = args.manifest.expanduser().resolve()
    manifest = _load_json(manifest_path)
    executable = args.executable.expanduser().resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise FileNotFoundError(executable)
    start = args.task_start
    stop = manifest["task_count"] if args.task_stop is None else args.task_stop
    if not 0 <= start < stop <= int(manifest["task_count"]):
        raise ValueError("task interval lies outside the campaign")
    workflow = args.workflow.strip()
    if not workflow or any(character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-" for character in workflow):
        raise ValueError("--workflow contains unsupported characters")
    root = Path(str(manifest["root"]))
    driver = Path(__file__).resolve()
    wrapper = root / f"run_swif_{start:06d}_{stop:06d}.sh"
    executable_hash = _sha256(executable)
    scratch_argument = (
        f" --scratch-root {shlex.quote(str(args.scratch_root.expanduser().resolve()))}"
        if args.scratch_root
        else ""
    )
    wrapper.write_text(
        "#!/bin/bash\nset -euo pipefail\n"
        "if [ -f /etc/profile.d/modules.sh ]; then source /etc/profile.d/modules.sh; fi\n"
        "module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles 2>/dev/null || true\n"
        "module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles 2>/dev/null || true\n"
        "module load clas12/5.4 2>/dev/null || true\n"
        f"test \"$(sha256sum \"$3\" | awk '{{print $1}}')\" = {shlex.quote(executable_hash)}\n"
        f"python3 {shlex.quote(str(driver))} run \"$1\" --replica-index \"$2\" "
        f"--executable \"$3\"{scratch_argument}\n",
        encoding="utf-8",
    )
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    output = args.output.expanduser().resolve() if args.output else root / f"submit_swif_{start:06d}_{stop:06d}.sh"
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        f"WORKFLOW={shlex.quote(workflow)}",
        f"MANIFEST={shlex.quote(str(manifest_path))}",
        f"EXECUTABLE={shlex.quote(str(executable))}",
        f"WRAPPER={shlex.quote(str(wrapper))}",
        'swif2 create -workflow "$WORKFLOW"',
    ]
    for task in range(start, stop):
        lines.extend(
            [
                "swif2 add-job \\",
                '  -workflow "$WORKFLOW" \\',
                f"  -name {shlex.quote(workflow + f'_g{task:08d}')} \\",
                f"  -cores {args.cores} -ram {shlex.quote(args.ram)} -disk {shlex.quote(args.disk)} -time {shlex.quote(args.walltime)} -os el9 \\",
                f'  -- /bin/bash "$WRAPPER" {shlex.quote(str(manifest["runs"][task]["task"]))} {task} "$EXECUTABLE"',
            ]
        )
    lines.extend(['swif2 run "$WORKFLOW"', 'echo "Submitted $WORKFLOW"'])
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    output.chmod(output.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    return output


def _common_prepare(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tag", default="aao_rad_mode3")
    parser.add_argument("--padding-fraction", type=float, default=0.035)
    parser.add_argument("--direct-fraction", type=float, default=0.75)
    parser.add_argument(
        "--minimum-photon-energy",
        type=float,
        help=(
            "override the legacy input's minimum radiated-photon energy in GeV; "
            "the value is frozen into calibration/production compatibility metadata"
        ),
    )
    parser.add_argument("--seed-base", type=int, required=True)
    parser.add_argument("--heartbeat-interval", type=int, default=100_000)
    parser.add_argument("--generator-revision", required=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    calibration = commands.add_parser("prepare-calibration")
    _common_prepare(calibration)
    calibration.add_argument("--trials", type=int, required=True)
    calibration.add_argument("--replicas", type=int, default=4)
    calibration.set_defaults(
        events_per_job=0,
        total_events=0,
        calibration=None,
        lund_files_per_directory=0,
    )
    production = commands.add_parser("prepare-production")
    _common_prepare(production)
    production.add_argument("--calibration", type=Path, required=True)
    production.add_argument("--total-events", type=int, required=True)
    production.add_argument("--events-per-job", type=int, default=5000)
    production.add_argument("--lund-files-per-directory", type=int, default=5000)
    production.set_defaults(trials=0, replicas=1)
    runner = commands.add_parser("run")
    runner.add_argument("manifest", type=Path)
    runner.add_argument("--replica-index", type=int, required=True)
    runner.add_argument("--executable", type=Path, required=True)
    runner.add_argument("--scratch-root", type=Path)
    runner.add_argument("--overwrite", action="store_true")
    state = commands.add_parser("status")
    state.add_argument("manifest", type=Path)
    state.add_argument("--output", type=Path)
    final_cal = commands.add_parser("finalize-calibration")
    final_cal.add_argument("manifest", type=Path)
    final_cal.add_argument("--envelope-safety-factor", type=float, default=1.5)
    final_cal.add_argument("--output", type=Path)
    compare = commands.add_parser("compare-calibrations")
    compare.add_argument("--candidate", type=Path, required=True)
    compare.add_argument("--control", type=Path, required=True)
    compare.add_argument("--maximum-relative-difference", type=float, default=0.10)
    compare.add_argument("--maximum-z-score", type=float, default=3.0)
    compare.add_argument("--output", type=Path, required=True)
    final_prod = commands.add_parser("finalize-production")
    final_prod.add_argument("manifest", type=Path)
    final_prod.add_argument("--output", type=Path)
    swif = commands.add_parser("emit-swif")
    swif.add_argument("manifest", type=Path)
    swif.add_argument("--workflow", required=True)
    swif.add_argument("--executable", type=Path, required=True)
    swif.add_argument("--task-start", type=int, default=0)
    swif.add_argument("--task-stop", type=int)
    swif.add_argument("--cores", type=int, default=1)
    swif.add_argument("--ram", default="2gb")
    swif.add_argument("--disk", default="2gb")
    swif.add_argument("--walltime", default="24hr")
    swif.add_argument("--scratch-root", type=Path)
    swif.add_argument("--output", type=Path)
    return parser


def main() -> int:
    args = _parser().parse_args()
    functions = {
        "prepare-calibration": prepare_calibration,
        "prepare-production": prepare_production,
        "run": run,
        "status": status,
        "finalize-calibration": finalize_calibration,
        "compare-calibrations": compare_calibrations,
        "finalize-production": finalize_production,
        "emit-swif": emit_swif,
    }
    result = functions[args.command](args)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
