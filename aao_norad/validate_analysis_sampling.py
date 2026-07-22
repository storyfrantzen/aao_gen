#!/usr/bin/env python3
"""Validate all bounded AAO sampling modes over the same physical domain."""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.stats import ks_2samp


MODE_NAMES = {
    1: "direct_q2",
    2: "bounded_legacy",
    3: "optimal_inverse_q2",
}
COMPARISON_MODES = ((3, 1), (3, 2), (1, 2))
LUND_Q2_TOLERANCE = 1.0e-5
LUND_XB_TOLERANCE = 1.0e-5
LUND_MINUS_T_TOLERANCE = 1.0e-3
# List-directed single-precision LUND momentum output limits the reconstructed
# hadron-plane angle, especially close to forward pion kinematics.
LUND_PHI_TOLERANCE_DEG = 5.0e-3


@dataclass(frozen=True)
class RunResult:
    mode: int
    seed: int
    sig_sum: float
    sig_int: float
    events: int
    ntries: int
    mcall_max: int
    kinematics: np.ndarray
    max_lund_q2_error: float
    max_lund_xb_error: float
    max_lund_minus_t_error: float
    max_lund_phi_error_deg: float


def _arguments() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description=(
            "Run repeated AAO samples in direct Q2, bounded legacy, and optimal "
            "inverse-Q2 analysis coordinates, then compare normalization, "
            "unweighted shapes, and proposal efficiency."
        )
    )
    parser.add_argument("--executable", type=Path, default=here / "build" / "aao_norad")
    parser.add_argument("--output", type=Path, default=here / "validation_sampling")
    parser.add_argument("--events", type=int, default=20_000)
    parser.add_argument("--replicas", type=int, default=4)
    parser.add_argument("--seed", type=int, default=371_001)
    parser.add_argument("--physics-model", type=int, default=1)
    parser.add_argument("--beam-energy", type=float, default=6.535)
    parser.add_argument("--q2", type=float, nargs=2, default=(1.0, 3.0), metavar=("MIN", "MAX"))
    parser.add_argument("--xb", type=float, nargs=2, default=(0.15, 0.50), metavar=("MIN", "MAX"))
    parser.add_argument("--minus-t", type=float, nargs=2, default=(0.09, 1.0), metavar=("MIN", "MAX"))
    parser.add_argument("--phi", type=float, nargs=2, default=(0.0, 360.0), metavar=("MIN", "MAX"))
    parser.add_argument("--ep-min", type=float, default=0.2)
    parser.add_argument("--fmcall", type=float, default=2.0)
    parser.add_argument("--ks-threshold", type=float, default=0.01)
    parser.add_argument("--xsec-z-threshold", type=float, default=3.0)
    return parser.parse_args()


def _input_text(args: argparse.Namespace, mode: int, seed: int) -> str:
    values = (
        f"{args.physics_model}\n0\n3\n1\n{args.beam_energy}\n"
        f"{args.q2[0]} {args.q2[1]}\n"
        f"{args.ep_min} {args.beam_energy}\n"
        f"{args.events}\n{args.fmcall}\n0\n{-abs(seed)}\n{mode}\n"
        f"{args.xb[0]} {args.xb[1]}\n"
        f"{args.minus_t[0]} {args.minus_t[1]}\n"
        f"{args.phi[0]} {args.phi[1]}\n"
    )
    return values


def _norm_value(text: str, key: str, cast: type[float] | type[int]) -> float | int:
    match = re.search(rf"(?m)^\s*{re.escape(key)}\s*=\s*([^\s]+)", text)
    if not match:
        raise RuntimeError(f"Missing {key} in aao_norad.norm")
    value = float(match.group(1).replace("D", "E").replace("d", "e"))
    return int(value) if cast is int else value


def _run(args: argparse.Namespace, mode: int, seed: int, directory: Path) -> RunResult:
    directory.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        [str(args.executable.resolve())],
        input=_input_text(args, mode, seed),
        text=True,
        cwd=directory,
        capture_output=True,
        check=False,
    )
    (directory / "stdout.log").write_text(completed.stdout, encoding="utf-8")
    (directory / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"AAO mode {mode} seed {seed} failed with {completed.returncode}; "
            f"see {directory / 'stderr.log'}"
        )
    norm_path = directory / "aao_norad.norm"
    if not norm_path.is_file():
        raise RuntimeError(
            f"AAO mode {mode} seed {seed} produced no normalization file.\n"
            f"stdout:\n{completed.stdout[-4000:]}\n"
            f"stderr:\n{completed.stderr[-4000:]}"
        )
    norm_text = norm_path.read_text(encoding="utf-8")
    kinematics = np.loadtxt(directory / "aao_norad.kin", comments="#", ndmin=2)
    if kinematics.shape[1] != 7:
        raise RuntimeError(f"Unexpected kinematics sidecar shape {kinematics.shape}")
    lund_errors = _lund_consistency(directory / "aao_norad.lund", kinematics, args.beam_energy)
    return RunResult(
        mode=mode,
        seed=seed,
        sig_sum=float(_norm_value(norm_text, "sig_sum", float)),
        sig_int=float(_norm_value(norm_text, "sig_int", float)),
        events=int(_norm_value(norm_text, "events", int)),
        ntries=int(_norm_value(norm_text, "ntries", int)),
        mcall_max=int(_norm_value(norm_text, "mcall_max", int)),
        kinematics=kinematics,
        max_lund_q2_error=lund_errors[0],
        max_lund_xb_error=lund_errors[1],
        max_lund_minus_t_error=lund_errors[2],
        max_lund_phi_error_deg=lund_errors[3],
    )


def _lund_consistency(path: Path, kin: np.ndarray, beam_energy: float) -> tuple[float, ...]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) != 5 * len(kin):
        raise RuntimeError(f"Expected {5 * len(kin)} LUND lines in {path}, found {len(lines)}")
    q2 = np.empty(len(kin))
    xb = np.empty(len(kin))
    minus_t = np.empty(len(kin))
    phi = np.empty(len(kin))
    beam = np.array([0.0, 0.0, beam_energy])
    proton_mass = 0.93828
    for event in range(len(kin)):
        header = np.fromstring(lines[5 * event], sep=" ")
        electron = np.fromstring(lines[5 * event + 1], sep=" ")
        proton = np.fromstring(lines[5 * event + 2], sep=" ")
        q2[event] = header[8]
        xb[event] = header[5]
        electron_p = electron[6:9]
        proton_p = proton[6:9]
        proton_energy = math.sqrt(float(proton_p @ proton_p) + proton_mass**2)
        minus_t[event] = 2.0 * proton_mass * (proton_energy - proton_mass)
        q_vector = beam - electron_p
        lepton_normal = np.cross(beam, electron_p)
        hadron_normal = np.cross(proton_p, q_vector)
        lepton_normal /= np.linalg.norm(lepton_normal)
        hadron_normal /= np.linalg.norm(hadron_normal)
        q_hat = q_vector / np.linalg.norm(q_vector)
        phi[event] = math.degrees(
            math.atan2(
                float(q_hat @ np.cross(lepton_normal, hadron_normal)),
                float(lepton_normal @ hadron_normal),
            )
        ) % 360.0
    phi_delta = (phi - kin[:, 4] + 180.0) % 360.0 - 180.0
    return (
        float(np.max(np.abs(q2 - kin[:, 1]))),
        float(np.max(np.abs(xb - kin[:, 2]))),
        float(np.max(np.abs(minus_t - kin[:, 3]))),
        float(np.max(np.abs(phi_delta))),
    )


def _mean_sem(values: np.ndarray) -> tuple[float, float]:
    mean = float(np.mean(values))
    sem = float(np.std(values, ddof=1) / math.sqrt(values.size)) if values.size > 1 else math.nan
    return mean, sem


def _serializable_run(result: RunResult) -> dict[str, float | int]:
    return {
        "mode": result.mode,
        "seed": result.seed,
        "sig_sum": result.sig_sum,
        "sig_int": result.sig_int,
        "events": result.events,
        "ntries": result.ntries,
        "mcall_max": result.mcall_max,
        "proposal_efficiency": result.events / result.ntries,
        "max_lund_q2_error": result.max_lund_q2_error,
        "max_lund_xb_error": result.max_lund_xb_error,
        "max_lund_minus_t_error": result.max_lund_minus_t_error,
        "max_lund_phi_error_deg": result.max_lund_phi_error_deg,
    }


def _plot(output: Path, samples: dict[str, np.ndarray]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    columns = ((1, r"$Q^2$ (GeV$^2$)"), (2, r"$x_B$"), (3, r"$-t$ (GeV$^2$)"), (4, r"$\phi$ (deg)"))
    figure, axes = plt.subplots(2, 2, figsize=(10, 8))
    for axis, (column, label) in zip(axes.ravel(), columns):
        combined = np.concatenate([sample[:, column] for sample in samples.values()])
        edges = np.linspace(float(np.min(combined)), float(np.max(combined)), 31)
        for name, sample in samples.items():
            axis.hist(
                sample[:, column],
                bins=edges,
                density=True,
                histtype="step",
                linewidth=1.5,
                label=name.replace("_", " "),
            )
        axis.set_xlabel(label)
        axis.set_ylabel("normalized density")
        axis.legend()
    figure.tight_layout()
    figure.savefig(output)
    plt.close(figure)


def _cross_section_summary(runs: list[RunResult]) -> dict[str, float]:
    values = np.array([run.sig_sum for run in runs])
    mean, sem = _mean_sem(values)
    return {"mean": mean, "sem": sem}


def _cross_section_comparison(
    candidate: dict[str, float], reference: dict[str, float]
) -> dict[str, float]:
    combined_sem = math.hypot(candidate["sem"], reference["sem"])
    difference = candidate["mean"] - reference["mean"]
    if combined_sem > 0.0:
        z_score = abs(difference) / combined_sem
    else:
        z_score = 0.0 if difference == 0.0 else math.inf
    return {
        "relative_difference": difference / reference["mean"],
        "difference_z_score": z_score,
    }


def _ks_comparison(candidate: np.ndarray, reference: np.ndarray) -> dict[str, dict[str, float]]:
    names = {1: "Q2", 2: "xB", 3: "minus_t", 4: "phi_deg"}
    results = {}
    for column, name in names.items():
        statistic, pvalue = ks_2samp(candidate[:, column], reference[:, column])
        results[name] = {"statistic": float(statistic), "pvalue": float(pvalue)}
    return results


def main() -> int:
    args = _arguments()
    if args.events <= 0 or args.replicas < 2:
        raise ValueError("events must be positive and replicas must be at least two")
    if not args.executable.is_file():
        raise FileNotFoundError(f"Build the generator first: {args.executable}")
    args.output.mkdir(parents=True, exist_ok=True)

    runs: dict[int, list[RunResult]] = {mode: [] for mode in MODE_NAMES}
    with tempfile.TemporaryDirectory(prefix="aao_sampling_validation_") as temporary:
        root = Path(temporary)
        for replica in range(args.replicas):
            for offset, mode in enumerate(MODE_NAMES):
                seed = args.seed + 3 * replica + offset
                name = MODE_NAMES[mode]
                runs[mode].append(_run(args, mode, seed, root / f"{name}_{replica:02d}"))

    summaries = {MODE_NAMES[mode]: _cross_section_summary(mode_runs) for mode, mode_runs in runs.items()}
    kinematics = {
        MODE_NAMES[mode]: np.concatenate([run.kinematics for run in mode_runs])
        for mode, mode_runs in runs.items()
    }
    cross_section_comparisons = {}
    ks_comparisons = {}
    for candidate_mode, reference_mode in COMPARISON_MODES:
        candidate = MODE_NAMES[candidate_mode]
        reference = MODE_NAMES[reference_mode]
        comparison = f"{candidate}_vs_{reference}"
        cross_section_comparisons[comparison] = _cross_section_comparison(
            summaries[candidate], summaries[reference]
        )
        ks_comparisons[comparison] = _ks_comparison(
            kinematics[candidate], kinematics[reference]
        )

    all_runs = [run for mode_runs in runs.values() for run in mode_runs]
    lund_passed = all(
        run.max_lund_q2_error < LUND_Q2_TOLERANCE
        and run.max_lund_xb_error < LUND_XB_TOLERANCE
        and run.max_lund_minus_t_error < LUND_MINUS_T_TOLERANCE
        and run.max_lund_phi_error_deg < LUND_PHI_TOLERANCE_DEG
        for run in all_runs
    )
    cross_section_passed = all(
        comparison["difference_z_score"] <= args.xsec_z_threshold
        for comparison in cross_section_comparisons.values()
    )
    shapes_passed = all(
        result["pvalue"] >= args.ks_threshold
        for comparison in ks_comparisons.values()
        for result in comparison.values()
    )
    passed = lund_passed and cross_section_passed and shapes_passed
    efficiency = {
        MODE_NAMES[mode]: sum(run.events for run in mode_runs)
        / sum(run.ntries for run in mode_runs)
        for mode, mode_runs in runs.items()
    }
    legacy_efficiency = efficiency[MODE_NAMES[2]]
    report = {
        "passed": passed,
        "configuration": {
            "events_per_replica": args.events,
            "replicas": args.replicas,
            "physics_model": args.physics_model,
            "beam_energy": args.beam_energy,
            "q2": args.q2,
            "xb": args.xb,
            "minus_t": args.minus_t,
            "phi_deg": args.phi,
            "fmcall": args.fmcall,
        },
        "integrated_cross_section_microbarn": {
            "modes": summaries,
            "comparisons": cross_section_comparisons,
        },
        "ks_tests": ks_comparisons,
        "proposal_efficiency": {
            name: {
                "aggregate": value,
                "relative_to_bounded_legacy": value / legacy_efficiency,
            }
            for name, value in efficiency.items()
        },
        "lund_kinematics_consistent": lund_passed,
        "lund_tolerances": {
            "Q2": LUND_Q2_TOLERANCE,
            "xB": LUND_XB_TOLERANCE,
            "minus_t": LUND_MINUS_T_TOLERANCE,
            "phi_deg": LUND_PHI_TOLERANCE_DEG,
        },
        "all_mcall_max_le_one": all(run.mcall_max <= 1 for run in all_runs),
        "runs": {
            MODE_NAMES[mode]: [_serializable_run(run) for run in mode_runs]
            for mode, mode_runs in runs.items()
        },
    }
    report_path = args.output / "validation.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    _plot(args.output / "shape_comparison.pdf", kinematics)
    print(json.dumps(report, indent=2))
    print(f"Wrote {report_path}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
