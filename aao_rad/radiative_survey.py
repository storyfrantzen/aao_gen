#!/usr/bin/env python3
"""Run and validate deterministic fixed-trial AAO radiative surveys."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from pathlib import Path
from typing import Iterable


PROTON_MASS_GEV = 0.9382720813
SURVEY_SCHEMA = "aao-rad-survey-v1"
SURVEY_FILENAME = "aao_rad.survey.csv"
NORM_FILENAME = "aao_rad.norm"
LUND_FILENAME = "aao_rad.lund"

INTEGER_COLUMNS = {
    "replica",
    "trial",
    "final_valid",
    "candidate_status",
    "intreg",
    "helicity",
}

SURVEY_COLUMNS = [
    "replica",
    "trial",
    "final_valid",
    "candidate_status",
    "intreg",
    "helicity",
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
    "cos_theta_gamma",
    "phi_gamma_deg",
    "cos_theta_cm",
    "phi_cm_deg",
    "energy_in_vertex",
    "energy_e_pre_external",
    "energy_e_lund",
    "energy_gamma",
    "final_e_px",
    "final_e_py",
    "final_e_pz",
    "final_e_energy",
    "final_p_px",
    "final_p_py",
    "final_p_pz",
    "final_p_energy",
    "q2_leptonic",
    "xb_leptonic",
    "w_leptonic",
    "q2_hard",
    "xb_hard",
    "w_hard",
    "minus_t_hard",
    "q2_observed",
    "xb_observed",
    "minus_t_observed",
    "phi_observed_deg",
    "w_observed",
    "y_observed",
    "integrand_internal",
    "integrand_observed",
    "trial_xsec_internal_microbarn",
    "trial_xsec_observed_microbarn",
    "rotation_deg",
    "incoming_loss_fraction",
    "outgoing_loss_fraction",
]

COORDINATE_TOLERANCES = {
    "q2_observed": 3.0e-5,
    "xb_observed": 3.0e-6,
    "minus_t_observed": 3.0e-5,
    "phi_observed_deg": 3.0e-4,
    "w_observed": 3.0e-6,
    "y_observed": 3.0e-6,
}


class SurveyValidationError(RuntimeError):
    """Raised when a survey violates a required milestone-1 invariant."""


def parse_norm(path: Path) -> dict[str, str]:
    """Parse AAO's key=value normalization sidecar."""

    values: dict[str, str] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        stripped = line.strip()
        if not stripped:
            continue
        if "=" not in stripped:
            raise SurveyValidationError(
                f"{path}:{line_number}: expected key=value normalization metadata"
            )
        key, value = stripped.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def read_survey(path: Path) -> tuple[str, list[dict[str, int | float]]]:
    """Read the versioned survey CSV without requiring NumPy or pandas."""

    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith("# schema="):
        raise SurveyValidationError(f"{path}: missing '# schema=...' header")
    schema = lines[0].split("=", 1)[1].strip()
    reader = csv.DictReader(lines[1:])
    if reader.fieldnames != SURVEY_COLUMNS:
        raise SurveyValidationError(
            f"{path}: unexpected columns for {schema}; "
            f"expected {SURVEY_COLUMNS}, found {reader.fieldnames}"
        )

    rows: list[dict[str, int | float]] = []
    for line_number, raw in enumerate(reader, 3):
        try:
            converted = {
                key: int(value) if key in INTEGER_COLUMNS else float(value)
                for key, value in raw.items()
            }
        except (TypeError, ValueError) as error:
            raise SurveyValidationError(
                f"{path}:{line_number}: invalid numeric field: {error}"
            ) from error
        rows.append(converted)
    return schema, rows


def observed_coordinates(
    beam_energy: float,
    electron: tuple[float, float, float, float],
    proton: tuple[float, float, float, float],
) -> dict[str, float]:
    """Reconstruct the analysis coordinates from final-state four-vectors."""

    ex, ey, ez, electron_energy = electron
    px, py, pz, proton_energy = proton
    q_energy = beam_energy - electron_energy
    qx, qy, qz = -ex, -ey, beam_energy - ez
    qnorm2 = qx * qx + qy * qy + qz * qz
    q2 = qnorm2 - q_energy * q_energy
    if q_energy <= 0.0 or q2 <= 0.0:
        raise SurveyValidationError("final four-vectors give nonphysical DIS coordinates")

    xb = q2 / (2.0 * PROTON_MASS_GEV * q_energy)
    y = q_energy / beam_energy
    w2 = PROTON_MASS_GEV**2 + 2.0 * PROTON_MASS_GEV * q_energy - q2
    if w2 <= 0.0:
        raise SurveyValidationError("final four-vectors give nonpositive observed W^2")
    minus_t = px * px + py * py + pz * pz - (PROTON_MASS_GEV - proton_energy) ** 2

    lepton_normal = _cross((0.0, 0.0, beam_energy), (ex, ey, ez))
    hadron_normal = _cross((px, py, pz), (qx, qy, qz))
    lepton_norm = _norm(lepton_normal)
    hadron_norm = _norm(hadron_normal)
    qnorm = math.sqrt(qnorm2)
    if min(lepton_norm, hadron_norm, qnorm) <= 0.0:
        raise SurveyValidationError("final four-vectors give an undefined Trento angle")

    cosine = _dot(lepton_normal, hadron_normal) / (lepton_norm * hadron_norm)
    sine = _dot(
        (qx, qy, qz), _cross(lepton_normal, hadron_normal)
    ) / (qnorm * lepton_norm * hadron_norm)
    phi_deg = math.degrees(math.atan2(sine, cosine)) % 360.0
    return {
        "q2_observed": q2,
        "xb_observed": xb,
        "minus_t_observed": minus_t,
        "phi_observed_deg": phi_deg,
        "w_observed": math.sqrt(w2),
        "y_observed": y,
    }


def validate_survey(directory: Path) -> dict[str, object]:
    """Validate fixed-trial normalization, records, and coordinate parity."""

    directory = directory.resolve()
    norm_path = directory / NORM_FILENAME
    survey_path = directory / SURVEY_FILENAME
    lund_path = directory / LUND_FILENAME
    if not norm_path.is_file():
        raise SurveyValidationError(f"missing {norm_path}")
    if not survey_path.is_file():
        raise SurveyValidationError(f"missing {survey_path}")

    norm = parse_norm(norm_path)
    _require_norm_value(norm, "generator", "aao_rad")
    _require_norm_value(norm, "sampling_mode", "1")
    _require_norm_value(norm, "fixed_trial_survey", "1")
    _require_norm_value(norm, "survey_schema", SURVEY_SCHEMA)
    _require_norm_value(norm, "survey_emits_lund", "0")

    schema, rows = read_survey(survey_path)
    if schema != SURVEY_SCHEMA:
        raise SurveyValidationError(
            f"survey schema is {schema!r}, expected {SURVEY_SCHEMA!r}"
        )

    ntrials = _norm_int(norm, "ntries")
    requested = _norm_int(norm, "survey_ntrials_requested")
    recorded_rows = _norm_int(norm, "survey_rows")
    recorded_valid = _norm_int(norm, "survey_final_valid")
    replica = _norm_int(norm, "survey_replica")
    if ntrials != requested:
        raise SurveyValidationError(
            f"fixed-trial run stopped after {ntrials} proposals, expected {requested}"
        )
    if recorded_rows != len(rows):
        raise SurveyValidationError(
            f"survey_rows={recorded_rows}, but CSV contains {len(rows)} rows"
        )
    valid_rows = [row for row in rows if row["final_valid"] == 1]
    if recorded_valid != len(valid_rows):
        raise SurveyValidationError(
            f"survey_final_valid={recorded_valid}, but CSV contains "
            f"{len(valid_rows)} valid rows"
        )

    trial_ids = [int(row["trial"]) for row in rows]
    if trial_ids != sorted(set(trial_ids)):
        raise SurveyValidationError("CSV trial indices are not unique and increasing")
    if any(trial < 1 or trial > ntrials for trial in trial_ids):
        raise SurveyValidationError("CSV contains a trial outside the fixed trial range")
    if any(row["replica"] != replica for row in rows):
        raise SurveyValidationError("CSV replica IDs do not match survey_replica")

    phase_volume = _norm_float(norm, "survey_phase_volume")
    internal_total = 0.0
    observed_total = 0.0
    internal_square_total = 0.0
    observed_square_total = 0.0
    max_internal_contribution = 0.0
    max_observed_contribution = 0.0
    max_errors = {name: 0.0 for name in COORDINATE_TOLERANCES}
    for row in rows:
        _validate_row_domains(row)
        _validate_proposal_mapping(row, norm)
        internal_integrand = float(row["integrand_internal"])
        observed_integrand = float(row["integrand_observed"])
        internal_contribution = float(row["trial_xsec_internal_microbarn"])
        observed_contribution = float(row["trial_xsec_observed_microbarn"])
        _assert_close(
            "internal trial cross section",
            internal_contribution,
            phase_volume * internal_integrand,
            relative=2.0e-6,
            absolute=1.0e-18,
        )
        internal_total += internal_contribution
        internal_square_total += internal_contribution**2
        max_internal_contribution = max(
            max_internal_contribution, internal_contribution
        )

        if row["final_valid"] == 0:
            _assert_close(
                "invalid observed integrand", observed_integrand, 0.0, absolute=0.0
            )
            _assert_close(
                "invalid observed contribution", observed_contribution, 0.0, absolute=0.0
            )
            continue

        if row["candidate_status"] != 0:
            raise SurveyValidationError("valid row has a nonzero candidate_status")
        _assert_close(
            "valid observed integrand",
            observed_integrand,
            internal_integrand,
            relative=2.0e-6,
            absolute=1.0e-18,
        )
        _assert_close(
            "observed trial cross section",
            observed_contribution,
            phase_volume * observed_integrand,
            relative=2.0e-6,
            absolute=1.0e-18,
        )
        observed_total += observed_contribution
        observed_square_total += observed_contribution**2
        max_observed_contribution = max(
            max_observed_contribution, observed_contribution
        )

        calculated = observed_coordinates(
            _norm_float(norm, "ebeam"),
            (
                float(row["final_e_px"]),
                float(row["final_e_py"]),
                float(row["final_e_pz"]),
                float(row["final_e_energy"]),
            ),
            (
                float(row["final_p_px"]),
                float(row["final_p_py"]),
                float(row["final_p_pz"]),
                float(row["final_p_energy"]),
            ),
        )
        for name, tolerance in COORDINATE_TOLERANCES.items():
            generated = float(row[name])
            if name == "phi_observed_deg":
                error = abs((calculated[name] - generated + 180.0) % 360.0 - 180.0)
            else:
                error = abs(calculated[name] - generated)
            max_errors[name] = max(max_errors[name], error)
            if error > tolerance:
                raise SurveyValidationError(
                    f"coordinate parity failed for trial {row['trial']}: "
                    f"{name} error {error:.6g} exceeds {tolerance:.6g}"
                )

    internal_estimate = internal_total / ntrials
    observed_estimate = observed_total / ntrials
    internal_sem = _fixed_trial_sem(
        internal_total, internal_square_total, ntrials
    )
    observed_sem = _fixed_trial_sem(
        observed_total, observed_square_total, ntrials
    )
    _assert_close(
        "survey internal integral",
        internal_estimate,
        _norm_float(norm, "survey_internal_sig_sum"),
        relative=3.0e-6,
        absolute=1.0e-18,
    )
    _assert_close(
        "survey observed integral",
        observed_estimate,
        _norm_float(norm, "survey_observed_sig_sum"),
        relative=3.0e-6,
        absolute=1.0e-18,
    )
    _assert_close(
        "survey sig_sum",
        observed_estimate,
        _norm_float(norm, "sig_sum"),
        relative=3.0e-6,
        absolute=1.0e-18,
    )

    lund_bytes = lund_path.stat().st_size if lund_path.exists() else 0
    if lund_bytes != 0:
        raise SurveyValidationError(
            f"survey mode must not emit LUND events, but {lund_path} has {lund_bytes} bytes"
        )

    return {
        "passed": True,
        "schema": schema,
        "survey_seed": _norm_int(norm, "survey_seed"),
        "replica": replica,
        "proposals": ntrials,
        "internally_valid_rows": len(rows),
        "final_valid_rows": len(valid_rows),
        "final_valid_fraction": len(valid_rows) / ntrials,
        "internal_cross_section_microbarn": internal_estimate,
        "internal_cross_section_sem_microbarn": internal_sem,
        "observed_cross_section_microbarn": observed_estimate,
        "observed_cross_section_sem_microbarn": observed_sem,
        "internal_effective_sample_size": _importance_effective_size(
            internal_total, internal_square_total
        ),
        "observed_effective_sample_size": _importance_effective_size(
            observed_total, observed_square_total
        ),
        "max_internal_trial_contribution_microbarn": max_internal_contribution,
        "max_observed_trial_contribution_microbarn": max_observed_contribution,
        "max_coordinate_errors": max_errors,
        "lund_bytes": lund_bytes,
    }


def run_survey(
    executable: Path,
    input_path: Path,
    output_directory: Path,
    trials: int,
    seed: int,
    replica: int,
) -> dict[str, object]:
    """Run AAO in an empty output directory and validate its artifacts."""

    if trials <= 0:
        raise ValueError("--trials must be positive")
    if seed == 0:
        raise ValueError("--seed must be nonzero")
    if replica < 0:
        raise ValueError("--replica must be nonnegative")
    executable = executable.resolve()
    input_path = input_path.resolve()
    output_directory = output_directory.resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise FileNotFoundError(f"AAO executable is not executable: {executable}")
    if not input_path.is_file():
        raise FileNotFoundError(f"AAO input does not exist: {input_path}")

    base_input = input_path.read_text(encoding="utf-8")
    _validate_legacy_input_shape(base_input, input_path)
    output_directory.mkdir(parents=True, exist_ok=True)
    generated_names = (
        SURVEY_FILENAME,
        NORM_FILENAME,
        LUND_FILENAME,
        "aao_rad.out",
        "aao_rad.sum",
        "generator.stdout.txt",
        "generator.stderr.txt",
        "survey_input.inp",
        "validation.json",
    )
    collisions = [name for name in generated_names if (output_directory / name).exists()]
    if collisions:
        raise FileExistsError(
            f"{output_directory} already contains survey artifacts: {collisions}"
        )

    survey_input = (
        base_input.rstrip()
        + f"\n1 ! fixed-trial survey mode\n{trials}\n{seed}\n{replica}\n"
    )
    (output_directory / "survey_input.inp").write_text(survey_input, encoding="utf-8")
    completed = subprocess.run(
        [str(executable)],
        input=survey_input,
        text=True,
        cwd=output_directory,
        capture_output=True,
        check=False,
    )
    (output_directory / "generator.stdout.txt").write_text(
        completed.stdout, encoding="utf-8"
    )
    (output_directory / "generator.stderr.txt").write_text(
        completed.stderr, encoding="utf-8"
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"AAO failed with exit code {completed.returncode}; see "
            f"{output_directory / 'generator.stderr.txt'}"
        )

    result = validate_survey(output_directory)
    (output_directory / "validation.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def _validate_row_domains(row: dict[str, int | float]) -> None:
    if row["final_valid"] not in (0, 1):
        raise SurveyValidationError("final_valid must be 0 or 1")
    if row["candidate_status"] not in (0, 1, 2, 3):
        raise SurveyValidationError("candidate_status must be in [0, 3]")
    if row["intreg"] not in range(1, 7):
        raise SurveyValidationError("intreg must be in [1, 6]")
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
        if not 0.0 <= value <= 1.0:
            raise SurveyValidationError(f"{name}={value} is outside [0, 1]")
    if not -1.0 <= float(row["cos_theta_gamma"]) <= 1.0:
        raise SurveyValidationError("cos_theta_gamma is outside [-1, 1]")
    if not -1.0 <= float(row["cos_theta_cm"]) <= 1.0:
        raise SurveyValidationError("cos_theta_cm is outside [-1, 1]")
    if float(row["integrand_internal"]) <= 0.0:
        raise SurveyValidationError("recorded internal integrand must be positive")
    if not 0.0 <= float(row["incoming_loss_fraction"]) < 1.0:
        raise SurveyValidationError("incoming loss fraction is outside [0, 1)")
    if not 0.0 <= float(row["outgoing_loss_fraction"]) <= 1.0:
        raise SurveyValidationError("outgoing loss fraction is outside [0, 1]")


def _validate_proposal_mapping(
    row: dict[str, int | float], norm: dict[str, str]
) -> None:
    q2_min = _norm_float(norm, "q2_min")
    q2_max = _norm_float(norm, "q2_max")
    inverse_q2 = 1.0 / q2_max + float(row["r_u"]) * (
        1.0 / q2_min - 1.0 / q2_max
    )
    _assert_close(
        "inverse-Q2 proposal mapping",
        float(row["q2_leptonic"]),
        1.0 / inverse_q2,
        relative=3.0e-6,
    )
    ep_min = _norm_float(norm, "ep_min")
    ep_max = _norm_float(norm, "ep_max_effective")
    expected_ep = ep_max - float(row["r_ep"]) * (ep_max - ep_min)
    _assert_close(
        "electron-energy proposal mapping",
        float(row["energy_e_pre_external"]),
        expected_ep,
        relative=3.0e-6,
        absolute=1.0e-7,
    )
    _assert_close(
        "photon-energy proposal mapping",
        float(row["energy_gamma"]),
        -math.log(float(row["u_gamma"])) / 5.0,
        relative=3.0e-6,
        absolute=1.0e-8,
    )
    _assert_close(
        "hadron cosine proposal mapping",
        float(row["cos_theta_cm"]),
        -1.0 + 2.0 * float(row["hadron_cosine_base"]),
        relative=3.0e-6,
        absolute=1.0e-7,
    )
    _assert_close(
        "hadron phi proposal mapping",
        float(row["phi_cm_deg"]),
        360.0 * float(row["hadron_phi_base"]),
        relative=3.0e-6,
        absolute=1.0e-5,
    )
    beam_energy = _norm_float(norm, "ebeam")
    _assert_close(
        "incoming-loss mapping",
        float(row["energy_in_vertex"]),
        beam_energy * (1.0 - float(row["incoming_loss_fraction"])),
        relative=3.0e-6,
        absolute=1.0e-7,
    )
    _assert_close(
        "outgoing-loss mapping",
        float(row["energy_e_lund"]),
        float(row["energy_e_pre_external"])
        * (1.0 - float(row["outgoing_loss_fraction"])),
        relative=3.0e-6,
        absolute=1.0e-7,
    )
    if row["final_valid"] == 1:
        _assert_close(
            "final electron energy",
            float(row["final_e_energy"]),
            float(row["energy_e_lund"]),
            relative=3.0e-6,
            absolute=1.0e-7,
        )


def _validate_legacy_input_shape(text: str, path: Path) -> None:
    records = [
        line.split("!", 1)[0].strip()
        for line in text.splitlines()
        if line.split("!", 1)[0].strip()
    ]
    if len(records) < 17:
        raise ValueError(f"{path}: expected at least 17 legacy AAO input records")
    try:
        theory = int(records[0].split()[0])
        fmcall = float(records[16].split()[0])
    except ValueError as error:
        raise ValueError(f"{path}: cannot parse theory or fmcall") from error
    expected = 17 + (1 if fmcall == 0.0 else 0) + (1 if theory > 10 else 0)
    if len(records) != expected:
        raise ValueError(
            f"{path}: expected {expected} legacy records, found {len(records)}; "
            "pass a legacy input without an optional sampling-mode trailer"
        )


def _require_norm_value(norm: dict[str, str], key: str, expected: str) -> None:
    value = norm.get(key)
    if value != expected:
        raise SurveyValidationError(f"{key}={value!r}, expected {expected!r}")


def _norm_int(norm: dict[str, str], key: str) -> int:
    try:
        return int(norm[key])
    except (KeyError, ValueError) as error:
        raise SurveyValidationError(f"normalization metadata lacks integer {key}") from error


def _norm_float(norm: dict[str, str], key: str) -> float:
    try:
        return float(norm[key])
    except (KeyError, ValueError) as error:
        raise SurveyValidationError(f"normalization metadata lacks float {key}") from error


def _assert_close(
    label: str,
    actual: float,
    expected: float,
    *,
    relative: float = 0.0,
    absolute: float = 0.0,
) -> None:
    if not math.isclose(actual, expected, rel_tol=relative, abs_tol=absolute):
        raise SurveyValidationError(
            f"{label}: {actual:.17g} does not match {expected:.17g}"
        )


def _fixed_trial_sem(total: float, square_total: float, count: int) -> float:
    if count <= 1:
        return math.nan
    sample_variance = max(0.0, square_total - total * total / count) / (count - 1)
    return math.sqrt(sample_variance / count)


def _importance_effective_size(total: float, square_total: float) -> float:
    if square_total <= 0.0:
        return 0.0
    return total * total / square_total


def _dot(first: Iterable[float], second: Iterable[float]) -> float:
    return sum(left * right for left, right in zip(first, second))


def _cross(
    first: tuple[float, float, float], second: tuple[float, float, float]
) -> tuple[float, float, float]:
    ax, ay, az = first
    bx, by, bz = second
    return ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx


def _norm(vector: tuple[float, float, float]) -> float:
    return math.sqrt(_dot(vector, vector))


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    run = subparsers.add_parser("run", help="run and validate a fixed-trial survey")
    run.add_argument("--executable", type=Path, default=Path("build/aao_rad"))
    run.add_argument("--input", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--trials", type=int, required=True)
    run.add_argument("--seed", type=int, required=True)
    run.add_argument("--replica", type=int, required=True)

    validate = subparsers.add_parser("validate", help="validate an existing survey")
    validate.add_argument("--directory", type=Path, required=True)
    return parser


def main() -> int:
    args = _build_parser().parse_args()
    try:
        if args.command == "run":
            result = run_survey(
                args.executable,
                args.input,
                args.output,
                args.trials,
                args.seed,
                args.replica,
            )
        else:
            result = validate_survey(args.directory)
    except (FileNotFoundError, FileExistsError, ValueError, RuntimeError) as error:
        print(f"ERROR: {error}")
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
