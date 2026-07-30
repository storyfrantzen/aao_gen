#!/usr/bin/env python3
"""Run and validate deterministic fixed-trial AAO radiative surveys."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import subprocess
from pathlib import Path
from typing import Iterable


PROTON_MASS_GEV = 0.9382720813
GENERATOR_PROTON_MASS_GEV = 0.938
PI0_MASS_GEV = 0.1349
LEGACY_SURVEY_SCHEMA = "aao-rad-survey-v1"
SURVEY_SCHEMA = "aao-rad-survey-v2"
SURVEY_FILENAME = "aao_rad.survey.csv"
NORM_FILENAME = "aao_rad.norm"
LUND_FILENAME = "aao_rad.lund"
ANALYSIS_CONFIG_FILENAME = "survey_analysis_config.json"

INTEGER_COLUMNS = {
    "replica",
    "trial",
    "final_valid",
    "candidate_status",
    "intreg",
    "helicity",
    "proposal_component",
    "proposal_q2_bin",
    "proposal_xb_bin",
    "proposal_t_bin",
    "proposal_phi_bin",
}

SURVEY_COLUMNS_V1 = [
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

SURVEY_COLUMNS = (
    SURVEY_COLUMNS_V1[:6]
    + [
        "proposal_component",
        "proposal_q2_bin",
        "proposal_xb_bin",
        "proposal_t_bin",
        "proposal_phi_bin",
        "proposal_density_ratio",
    ]
    + SURVEY_COLUMNS_V1[6:]
)

SURVEY_COLUMNS_BY_SCHEMA = {
    LEGACY_SURVEY_SCHEMA: SURVEY_COLUMNS_V1,
    SURVEY_SCHEMA: SURVEY_COLUMNS,
}

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


def _records(text: str) -> list[str]:
    return [
        line.split("!", 1)[0].strip()
        for line in text.splitlines()
        if line.split("!", 1)[0].strip()
    ]


def _legacy_record_count(records: list[str], path: Path) -> int:
    if len(records) < 17:
        raise ValueError(f"{path}: expected at least 17 legacy AAO input records")
    try:
        theory = int(records[0].split()[0])
        fmcall = float(records[16].split()[0])
    except ValueError as error:
        raise ValueError(f"{path}: cannot parse theory or fmcall") from error
    return 17 + (1 if fmcall == 0.0 else 0) + (1 if theory > 10 else 0)


def _strict_edges(values: object, name: str) -> list[float]:
    if not isinstance(values, list):
        raise ValueError(f"analysis config binning.{name} must be a list")
    try:
        edges = [float(value) for value in values]
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"analysis config binning.{name} contains a nonnumeric edge"
        ) from error
    if (
        len(edges) < 2
        or any(not math.isfinite(value) for value in edges)
        or any(right <= left for left, right in zip(edges, edges[1:]))
    ):
        raise ValueError(
            f"analysis config binning.{name} must have increasing finite edges"
        )
    if len(edges) - 1 > 64:
        raise ValueError(f"analysis config binning.{name} exceeds 64 bins")
    return edges


def _load_balanced_config(
    path: Path,
    *,
    legacy_input: str | None = None,
    legacy_path: Path | None = None,
) -> tuple[dict[str, object], bytes, str]:
    raw = path.read_bytes()
    try:
        config = json.loads(raw)
    except json.JSONDecodeError as error:
        raise ValueError(f"{path}: invalid analysis JSON: {error}") from error
    try:
        target_mass = float(config["target_mass"])
        binning = config["binning"]
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            f"{path}: balanced survey requires target_mass and binning"
        ) from error
    if not math.isclose(
        target_mass, PROTON_MASS_GEV, rel_tol=0.0, abs_tol=5.0e-7
    ):
        raise ValueError(
            f"{path}: target_mass={target_mass} differs from "
            f"{PROTON_MASS_GEV}"
        )
    edges = {
        name: _strict_edges(binning.get(name), name)
        for name in ("Q2", "xB", "minus_t", "phi_deg")
    }
    if edges["xB"][0] <= 0.0:
        raise ValueError("balanced survey xB edges must be positive")
    if edges["minus_t"][0] < 0.0:
        raise ValueError("balanced survey -t edges must be nonnegative")
    if not (
        math.isclose(edges["phi_deg"][0], 0.0, abs_tol=1.0e-9)
        and math.isclose(edges["phi_deg"][-1], 360.0, abs_tol=1.0e-9)
    ):
        raise ValueError("balanced survey phi_deg edges must span [0,360]")
    if legacy_input is not None:
        input_path = legacy_path or Path("<legacy-input>")
        records = _records(legacy_input)
        legacy_count = _legacy_record_count(records, input_path)
        if len(records) != legacy_count:
            raise ValueError(
                f"{input_path}: pass a legacy input without an optional "
                "sampling-mode trailer"
            )
        try:
            q2_min, q2_max = (
                float(value) for value in records[12].split()[:2]
            )
            beam_energy = float(records[11].split()[0])
            epirea = int(records[4].split()[0])
            npart = int(records[3].split()[0])
        except (ValueError, IndexError) as error:
            raise ValueError(
                f"{input_path}: cannot parse balanced-survey settings"
            ) from error
        if edges["Q2"][0] < q2_min or edges["Q2"][-1] > q2_max:
            raise ValueError(
                "analysis Q2 edges must lie inside the generator Q2 range"
            )
        try:
            analysis_beam_energy = float(config["beam_energy"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(
                f"{path}: balanced survey requires beam_energy"
            ) from error
        if not math.isclose(
            analysis_beam_energy, beam_energy, rel_tol=2.0e-7
        ):
            raise ValueError(
                "analysis and generator beam energies differ"
            )
        if epirea != 1 or npart != 4:
            raise ValueError(
                "balanced survey currently requires epirea=1 and npart=4"
            )
    return (
        {"target_mass": target_mass, "binning": edges},
        raw,
        hashlib.sha256(raw).hexdigest(),
    )


def _survey_trailer(
    trials: int,
    seed: int,
    replica: int,
    *,
    proposal: str,
    legacy_fraction: float,
    balanced_config: dict[str, object] | None,
) -> str:
    records = ["1", str(trials), str(seed), str(replica)]
    if proposal == "legacy":
        records.append("0")
    else:
        if balanced_config is None:
            raise ValueError("balanced proposal requires an analysis config")
        records.extend(["1", f"{legacy_fraction:.17g}"])
        for name in ("Q2", "xB", "minus_t", "phi_deg"):
            edges = balanced_config["binning"][name]
            records.extend(
                [
                    str(len(edges) - 1),
                    " ".join(f"{float(value):.17g}" for value in edges),
                ]
            )
    return "\n".join(records) + "\n"


def _parse_survey_input(path: Path) -> dict[str, object]:
    records = _records(path.read_text(encoding="utf-8"))
    legacy_count = _legacy_record_count(records, path)
    trailer = records[legacy_count:]
    if len(trailer) < 4 or trailer[0].split()[0] != "1":
        raise SurveyValidationError(f"{path}: not a fixed-trial survey input")
    try:
        spec: dict[str, object] = {
            "trials": int(trailer[1].split()[0]),
            "seed": int(trailer[2].split()[0]),
            "replica": int(trailer[3].split()[0]),
            "mode": int(trailer[4].split()[0]) if len(trailer) >= 5 else 0,
            "legacy_records": records[:legacy_count],
        }
    except (ValueError, IndexError) as error:
        raise SurveyValidationError(
            f"{path}: malformed fixed-trial survey trailer"
        ) from error
    if spec["mode"] == 0:
        if len(trailer) not in (4, 5):
            raise SurveyValidationError(
                f"{path}: unexpected records after legacy survey proposal"
            )
        spec["legacy_fraction"] = 1.0
        return spec
    if spec["mode"] != 1 or len(trailer) != 14:
        raise SurveyValidationError(
            f"{path}: malformed balanced survey proposal trailer"
        )
    try:
        spec["legacy_fraction"] = float(trailer[5].split()[0])
        binning: dict[str, list[float]] = {}
        position = 6
        for name in ("Q2", "xB", "minus_t", "phi_deg"):
            count = int(trailer[position].split()[0])
            edges = [float(value) for value in trailer[position + 1].split()]
            if len(edges) != count + 1:
                raise SurveyValidationError(
                    f"{path}: {name} trailer has {len(edges)} edges "
                    f"for {count} bins"
                )
            binning[name] = _strict_edges(edges, name)
            position += 2
        spec["binning"] = binning
    except (ValueError, IndexError) as error:
        if isinstance(error, SurveyValidationError):
            raise
        raise SurveyValidationError(
            f"{path}: malformed balanced survey proposal settings"
        ) from error
    return spec


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
    try:
        expected_columns = SURVEY_COLUMNS_BY_SCHEMA[schema]
    except KeyError as error:
        raise SurveyValidationError(
            f"{path}: unsupported survey schema {schema!r}"
        ) from error
    reader = csv.DictReader(lines[1:])
    if reader.fieldnames != expected_columns:
        raise SurveyValidationError(
            f"{path}: unexpected columns for {schema}; "
            f"expected {expected_columns}, found {reader.fieldnames}"
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
    _require_norm_value(norm, "survey_emits_lund", "0")
    norm_schema = norm.get("survey_schema")
    if norm_schema not in SURVEY_COLUMNS_BY_SCHEMA:
        raise SurveyValidationError(
            f"survey_schema={norm_schema!r} is unsupported"
        )

    schema, rows = read_survey(survey_path)
    if schema != norm_schema:
        raise SurveyValidationError(
            f"CSV schema is {schema!r}, norm records {norm_schema!r}"
        )
    proposal_spec = _parse_survey_input(directory / "survey_input.inp")
    proposal_mode = int(proposal_spec["mode"])
    if schema == LEGACY_SURVEY_SCHEMA and proposal_mode != 0:
        raise SurveyValidationError(
            "v1 survey schema cannot describe a balanced proposal"
        )
    if schema == SURVEY_SCHEMA:
        _require_norm_value(
            norm, "survey_proposal_mode", str(proposal_mode)
        )
        expected_name = (
            "analysis_balanced_mixture" if proposal_mode == 1 else "legacy"
        )
        _require_norm_value(norm, "survey_proposal", expected_name)
        _assert_close(
            "survey legacy fraction",
            _norm_float(norm, "survey_legacy_fraction"),
            float(proposal_spec["legacy_fraction"]),
            relative=2.0e-7,
            absolute=1.0e-8,
        )
        if proposal_mode == 1:
            config_path = directory / ANALYSIS_CONFIG_FILENAME
            if not config_path.is_file():
                raise SurveyValidationError(
                    f"balanced survey is missing {config_path}"
                )
            frozen_config, raw_config, config_sha256 = _load_balanced_config(
                config_path
            )
            if frozen_config["binning"] != proposal_spec["binning"]:
                raise SurveyValidationError(
                    "balanced survey input edges differ from frozen analysis config"
                )
            _require_norm_value(
                norm, "survey_analysis_config_sha256", config_sha256
            )
            proposal_spec["analysis_config_sha256"] = config_sha256
            proposal_spec["analysis_config_bytes"] = len(raw_config)

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
    if (
        int(proposal_spec["trials"]) != ntrials
        or int(proposal_spec["seed"]) != _norm_int(norm, "survey_seed")
        or int(proposal_spec["replica"]) != replica
    ):
        raise SurveyValidationError(
            "survey input trailer disagrees with normalization metadata"
        )
    legacy_trials = (
        _norm_int(norm, "survey_legacy_trials")
        if schema == SURVEY_SCHEMA
        else ntrials
    )
    balanced_trials = (
        _norm_int(norm, "survey_balanced_trials")
        if schema == SURVEY_SCHEMA
        else 0
    )
    if legacy_trials + balanced_trials != ntrials:
        raise SurveyValidationError(
            "legacy and balanced trial counts do not sum to ntries"
        )

    phase_volume = _norm_float(norm, "survey_phase_volume")
    internal_total = 0.0
    observed_total = 0.0
    internal_square_total = 0.0
    observed_square_total = 0.0
    max_internal_contribution = 0.0
    max_observed_contribution = 0.0
    max_proposal_density_ratio = 0.0
    max_errors = {name: 0.0 for name in COORDINATE_TOLERANCES}
    for row in rows:
        _validate_row_domains(row, schema=schema, proposal_spec=proposal_spec)
        _validate_proposal_mapping(
            row, norm, schema=schema, proposal_spec=proposal_spec
        )
        if schema == SURVEY_SCHEMA:
            max_proposal_density_ratio = max(
                max_proposal_density_ratio,
                float(row["proposal_density_ratio"]),
            )
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

    result: dict[str, object] = {
        "passed": True,
        "schema": schema,
        "survey_seed": _norm_int(norm, "survey_seed"),
        "replica": replica,
        "proposal": (
            "legacy" if proposal_mode == 0 else "analysis_balanced_mixture"
        ),
        "legacy_fraction": float(proposal_spec["legacy_fraction"]),
        "legacy_trials": legacy_trials,
        "balanced_trials": balanced_trials,
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
        "max_proposal_density_ratio": max_proposal_density_ratio,
        "max_coordinate_errors": max_errors,
        "lund_bytes": lund_bytes,
    }
    if proposal_mode == 1:
        result["allocation"] = _allocation_summary(rows, proposal_spec)
    return result


def run_survey(
    executable: Path,
    input_path: Path,
    output_directory: Path,
    trials: int,
    seed: int,
    replica: int,
    *,
    proposal: str = "legacy",
    config_path: Path | None = None,
    legacy_fraction: float = 0.25,
) -> dict[str, object]:
    """Run AAO in an empty output directory and validate its artifacts."""

    if trials <= 0:
        raise ValueError("--trials must be positive")
    if seed == 0:
        raise ValueError("--seed must be nonzero")
    if replica < 0:
        raise ValueError("--replica must be nonnegative")
    if proposal not in ("legacy", "balanced"):
        raise ValueError("--proposal must be legacy or balanced")
    if not 0.0 < legacy_fraction < 1.0:
        raise ValueError("--legacy-fraction must lie strictly between zero and one")
    executable = executable.resolve()
    input_path = input_path.resolve()
    output_directory = output_directory.resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise FileNotFoundError(f"AAO executable is not executable: {executable}")
    if not input_path.is_file():
        raise FileNotFoundError(f"AAO input does not exist: {input_path}")

    base_input = input_path.read_text(encoding="utf-8")
    _validate_legacy_input_shape(base_input, input_path)
    balanced_config: dict[str, object] | None = None
    config_raw: bytes | None = None
    config_sha256: str | None = None
    if proposal == "balanced":
        if config_path is None:
            raise ValueError("--config is required with --proposal balanced")
        config_path = config_path.resolve()
        if not config_path.is_file():
            raise FileNotFoundError(
                f"analysis config does not exist: {config_path}"
            )
        balanced_config, config_raw, config_sha256 = _load_balanced_config(
            config_path,
            legacy_input=base_input,
            legacy_path=input_path,
        )
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
        ANALYSIS_CONFIG_FILENAME,
        "validation.json",
    )
    collisions = [name for name in generated_names if (output_directory / name).exists()]
    if collisions:
        raise FileExistsError(
            f"{output_directory} already contains survey artifacts: {collisions}"
        )

    survey_input = base_input.rstrip() + "\n" + _survey_trailer(
        trials,
        seed,
        replica,
        proposal=proposal,
        legacy_fraction=legacy_fraction,
        balanced_config=balanced_config,
    )
    (output_directory / "survey_input.inp").write_text(survey_input, encoding="utf-8")
    if config_raw is not None:
        (output_directory / ANALYSIS_CONFIG_FILENAME).write_bytes(config_raw)
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
    if config_sha256 is not None:
        with (output_directory / NORM_FILENAME).open(
            "a", encoding="utf-8"
        ) as norm_output:
            norm_output.write(
                f"survey_analysis_config_sha256={config_sha256}\n"
            )

    result = validate_survey(output_directory)
    (output_directory / "validation.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def _validate_row_domains(
    row: dict[str, int | float],
    *,
    schema: str,
    proposal_spec: dict[str, object],
) -> None:
    if row["final_valid"] not in (0, 1):
        raise SurveyValidationError("final_valid must be 0 or 1")
    if row["candidate_status"] not in (0, 1, 2, 3):
        raise SurveyValidationError("candidate_status must be in [0, 3]")
    if row["intreg"] not in range(1, 7):
        raise SurveyValidationError("intreg must be in [1, 6]")
    if schema == SURVEY_SCHEMA:
        component = int(row["proposal_component"])
        mode = int(proposal_spec["mode"])
        if component not in (0, 1):
            raise SurveyValidationError("proposal_component must be zero or one")
        if mode == 0 and component != 0:
            raise SurveyValidationError(
                "legacy survey row has a balanced proposal component"
            )
        bin_names = (
            "proposal_q2_bin",
            "proposal_xb_bin",
            "proposal_t_bin",
            "proposal_phi_bin",
        )
        if component == 0 and any(int(row[name]) != -1 for name in bin_names):
            raise SurveyValidationError(
                "legacy proposal component must use -1 target-bin indices"
            )
        if component == 1 and any(int(row[name]) < 0 for name in bin_names):
            raise SurveyValidationError(
                "balanced proposal component lacks target-bin indices"
            )
        ratio = float(row["proposal_density_ratio"])
        maximum = 1.0 / float(proposal_spec["legacy_fraction"])
        if not math.isfinite(ratio) or not 0.0 < ratio <= maximum * (1.0 + 1.0e-5):
            raise SurveyValidationError(
                "proposal_density_ratio is outside its full-support bound"
            )
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
    row: dict[str, int | float],
    norm: dict[str, str],
    *,
    schema: str,
    proposal_spec: dict[str, object],
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
    if schema == SURVEY_SCHEMA:
        expected_ratios = _proposal_density_ratio_candidates(
            row,
            norm,
            proposal_spec,
        )
        actual_ratio = float(row["proposal_density_ratio"])
        if not any(
            math.isclose(
                actual_ratio,
                expected_ratio,
                rel_tol=2.0e-5,
                abs_tol=2.0e-7,
            )
            for expected_ratio in expected_ratios
        ):
            raise SurveyValidationError(
                f"trial {row['trial']} survey proposal density ratio "
                f"{actual_ratio:.17g} does not match boundary-aware "
                f"candidates {expected_ratios}"
            )
        if int(row["proposal_component"]) == 1:
            binning = proposal_spec["binning"]
            values = (
                (
                    "proposal_q2_bin",
                    float(row["q2_leptonic"]),
                    binning["Q2"],
                    COORDINATE_TOLERANCES["q2_observed"],
                    False,
                ),
                (
                    "proposal_xb_bin",
                    float(row["xb_leptonic"]),
                    binning["xB"],
                    COORDINATE_TOLERANCES["xb_observed"],
                    False,
                ),
                (
                    "proposal_t_bin",
                    float(row["minus_t_hard"]),
                    binning["minus_t"],
                    COORDINATE_TOLERANCES["minus_t_observed"],
                    False,
                ),
                (
                    "proposal_phi_bin",
                    float(row["phi_cm_deg"]),
                    binning["phi_deg"],
                    COORDINATE_TOLERANCES["phi_observed_deg"],
                    True,
                ),
            )
            for column, value, edges, tolerance, periodic in values:
                declared_bin = int(row[column])
                if not _declared_bin_contains(
                    value,
                    edges,
                    declared_bin,
                    tolerance=tolerance,
                    periodic=periodic,
                ):
                    lower = (
                        edges[declared_bin]
                        if 0 <= declared_bin < len(edges) - 1
                        else math.nan
                    )
                    upper = (
                        edges[declared_bin + 1]
                        if 0 <= declared_bin < len(edges) - 1
                        else math.nan
                    )
                    raise SurveyValidationError(
                        f"{column}={declared_bin} declares [{lower}, {upper}], "
                        f"but the generated value is {value:.17g}"
                    )
    if row["final_valid"] == 1:
        _assert_close(
            "final electron energy",
            float(row["final_e_energy"]),
            float(row["energy_e_lund"]),
            relative=3.0e-6,
            absolute=1.0e-7,
        )


def _bin_index(value: float, edges: list[float]) -> int | None:
    for index, (lower, upper) in enumerate(zip(edges, edges[1:])):
        if lower <= value < upper:
            return index
    return None


def _declared_bin_contains(
    value: float,
    edges: list[float],
    index: int,
    *,
    tolerance: float,
    periodic: bool = False,
) -> bool:
    """Allow only numerical boundary spillover around a declared proposal bin."""

    if index < 0 or index >= len(edges) - 1:
        return False
    candidates = [value]
    if periodic:
        period = edges[-1] - edges[0]
        wrapped = (value - edges[0]) % period + edges[0]
        candidates = [wrapped, wrapped - period, wrapped + period]
    lower, upper = edges[index], edges[index + 1]
    return any(
        lower - tolerance <= candidate <= upper + tolerance
        for candidate in candidates
    )


def _bin_candidates(
    value: float,
    edges: list[float],
    *,
    tolerance: float,
    periodic: bool = False,
) -> list[int | None]:
    """Return only bin classifications ambiguous at numerical boundaries."""

    candidates: list[int | None] = [
        index
        for index in range(len(edges) - 1)
        if _declared_bin_contains(
            value,
            edges,
            index,
            tolerance=tolerance,
            periodic=periodic,
        )
    ]
    strict_value = value
    if periodic:
        strict_value = (value - edges[0]) % (edges[-1] - edges[0]) + edges[0]
    strict = _bin_index(strict_value, edges)
    if strict is not None and strict not in candidates:
        candidates.append(strict)
    if not candidates:
        return [None]
    if not periodic and (
        abs(value - edges[0]) <= tolerance
        or abs(value - edges[-1]) <= tolerance
    ):
        candidates.append(None)
    return candidates


def _allocation_summary(
    rows: list[dict[str, int | float]],
    proposal_spec: dict[str, object],
) -> dict[str, object]:
    """Summarize raw target and final-bin allocation without using weights."""

    binning = proposal_spec["binning"]
    axes = (
        ("Q2", "proposal_q2_bin", "q2_observed"),
        ("xB", "proposal_xb_bin", "xb_observed"),
        ("minus_t", "proposal_t_bin", "minus_t_observed"),
        ("phi_deg", "proposal_phi_bin", "phi_observed_deg"),
    )
    target_counts = {
        axis: [0] * (len(binning[axis]) - 1)
        for axis, _, _ in axes
    }
    observed_counts = {
        axis: [0] * (len(binning[axis]) - 1)
        for axis, _, _ in axes
    }
    target_joint: set[tuple[int, ...]] = set()
    observed_joint: set[tuple[int, ...]] = set()
    balanced_recorded_rows = 0
    final_rows_inside = 0

    for row in rows:
        if int(row["proposal_component"]) == 1:
            balanced_recorded_rows += 1
            target_indices = tuple(
                int(row[target_column])
                for _, target_column, _ in axes
            )
            for (axis, _, _), index in zip(axes, target_indices):
                target_counts[axis][index] += 1
            target_joint.add(target_indices)

        if int(row["final_valid"]) != 1:
            continue
        observed_indices = tuple(
            _bin_index(
                float(row[observed_column]) % 360.0
                if axis == "phi_deg"
                else float(row[observed_column]),
                binning[axis],
            )
            for axis, _, observed_column in axes
        )
        if any(index is None for index in observed_indices):
            continue
        final_rows_inside += 1
        integer_indices = tuple(int(index) for index in observed_indices)
        for (axis, _, _), index in zip(axes, integer_indices):
            observed_counts[axis][index] += 1
        observed_joint.add(integer_indices)

    total_joint = math.prod(len(binning[axis]) - 1 for axis, _, _ in axes)
    return {
        "note": (
            "Raw row counts diagnose allocation only; cross-section "
            "estimators and guard coverage still use corrected trial weights."
        ),
        "balanced_component_recorded_internal_rows": balanced_recorded_rows,
        "balanced_target_axis_counts_in_recorded_internal_rows": target_counts,
        "balanced_target_joint_strata_occupied_in_recorded_internal_rows": len(
            target_joint
        ),
        "final_valid_rows_inside_analysis_binning": final_rows_inside,
        "final_observed_axis_counts": observed_counts,
        "final_observed_joint_strata_occupied": len(observed_joint),
        "analysis_joint_strata_total": total_joint,
    }


def _radiative_t_jacobian(row: dict[str, int | float]) -> float:
    es = float(row["energy_in_vertex"])
    ep = float(row["energy_e_pre_external"])
    q2 = float(row["q2_leptonic"])
    photon_energy = float(row["energy_gamma"])
    photon_cosine = float(row["cos_theta_gamma"])
    nu = es - ep
    if min(es, ep, q2) <= 0.0:
        return 0.0
    q_magnitude = math.sqrt(q2 + nu * nu)
    hadron_energy = nu + GENERATOR_PROTON_MASS_GEV - photon_energy
    hadron_momentum_squared = (
        q_magnitude * q_magnitude
        + photon_energy * photon_energy
        - 2.0 * q_magnitude * photon_energy * photon_cosine
    )
    if hadron_energy <= 0.0 or hadron_momentum_squared <= 0.0:
        return 0.0
    w_squared = hadron_energy * hadron_energy - hadron_momentum_squared
    threshold = GENERATOR_PROTON_MASS_GEV + PI0_MASS_GEV
    if w_squared <= threshold * threshold:
        return 0.0
    w = math.sqrt(w_squared)
    hadron_momentum = math.sqrt(hadron_momentum_squared)
    beta = hadron_momentum / hadron_energy
    gamma = hadron_energy / w
    pstar_squared = (
        (
            w * w
            - GENERATOR_PROTON_MASS_GEV**2
            - PI0_MASS_GEV**2
        )
        ** 2
        / 4.0
        - (GENERATOR_PROTON_MASS_GEV * PI0_MASS_GEV) ** 2
    ) / (w * w)
    if pstar_squared <= 0.0:
        return 0.0
    return (
        2.0
        * PROTON_MASS_GEV
        * gamma
        * beta
        * math.sqrt(pstar_squared)
    )


def _proposal_density_ratio(
    row: dict[str, int | float],
    norm: dict[str, str],
    proposal_spec: dict[str, object],
) -> float:
    legacy_fraction = float(proposal_spec["legacy_fraction"])
    if int(proposal_spec["mode"]) == 0:
        return 1.0
    coordinates = _proposal_density_coordinates(row)
    if coordinates is None:
        return 1.0 / legacy_fraction
    q2, xb, minus_t, phi = coordinates
    binning = proposal_spec["binning"]
    indices = (
        _bin_index(q2, binning["Q2"]),
        _bin_index(xb, binning["xB"]),
        _bin_index(minus_t, binning["minus_t"]),
        _bin_index(phi % 360.0, binning["phi_deg"]),
    )
    return _proposal_density_ratio_for_indices(
        row, norm, proposal_spec, coordinates, indices
    )


def _proposal_density_ratio_candidates(
    row: dict[str, int | float],
    norm: dict[str, str],
    proposal_spec: dict[str, object],
) -> list[float]:
    legacy_fraction = float(proposal_spec["legacy_fraction"])
    if int(proposal_spec["mode"]) == 0:
        return [1.0]
    coordinates = _proposal_density_coordinates(row)
    if coordinates is None:
        return [1.0 / legacy_fraction]
    q2, xb, minus_t, phi = coordinates
    binning = proposal_spec["binning"]
    choices = (
        _bin_candidates(
            q2,
            binning["Q2"],
            tolerance=COORDINATE_TOLERANCES["q2_observed"],
        ),
        _bin_candidates(
            xb,
            binning["xB"],
            tolerance=COORDINATE_TOLERANCES["xb_observed"],
        ),
        _bin_candidates(
            minus_t,
            binning["minus_t"],
            tolerance=COORDINATE_TOLERANCES["minus_t_observed"],
        ),
        _bin_candidates(
            phi,
            binning["phi_deg"],
            tolerance=COORDINATE_TOLERANCES["phi_observed_deg"],
            periodic=True,
        ),
    )
    candidates = {
        _proposal_density_ratio_for_indices(
            row, norm, proposal_spec, coordinates, indices
        )
        for indices in itertools.product(*choices)
    }
    return sorted(candidates)


def _proposal_density_coordinates(
    row: dict[str, int | float],
) -> tuple[float, float, float, float] | None:
    q2 = float(row["q2_leptonic"])
    energy_transfer = (
        float(row["energy_in_vertex"])
        - float(row["energy_e_pre_external"])
    )
    if energy_transfer <= 0.0:
        return None
    xb = q2 / (2.0 * PROTON_MASS_GEV * energy_transfer)
    minus_t = float(row["minus_t_hard"])
    phi = float(row["phi_cm_deg"])
    return q2, xb, minus_t, phi


def _proposal_density_ratio_for_indices(
    row: dict[str, int | float],
    norm: dict[str, str],
    proposal_spec: dict[str, object],
    coordinates: tuple[float, float, float, float],
    indices: tuple[int | None, int | None, int | None, int | None],
) -> float:
    legacy_fraction = float(proposal_spec["legacy_fraction"])
    q2, xb, _, _ = coordinates
    binning = proposal_spec["binning"]
    if any(index is None for index in indices):
        return 1.0 / legacy_fraction
    iq2, ixb, it, iphi = (int(index) for index in indices)
    q2_width = binning["Q2"][iq2 + 1] - binning["Q2"][iq2]
    xb_width = binning["xB"][ixb + 1] - binning["xB"][ixb]
    t_width = binning["minus_t"][it + 1] - binning["minus_t"][it]
    phi_width = binning["phi_deg"][iphi + 1] - binning["phi_deg"][iphi]
    inverse_q2_range = (
        1.0 / _norm_float(norm, "q2_min")
        - 1.0 / _norm_float(norm, "q2_max")
    )
    electron_energy_range = (
        _norm_float(norm, "ep_max_effective")
        - _norm_float(norm, "ep_min")
    )
    t_jacobian = _radiative_t_jacobian(row)
    if t_jacobian <= 0.0:
        return 1.0 / legacy_fraction
    balanced_to_legacy = (
        inverse_q2_range
        * q2**2
        / ((len(binning["Q2"]) - 1) * q2_width)
        * electron_energy_range
        * 2.0
        * PROTON_MASS_GEV
        * xb**2
        / (q2 * (len(binning["xB"]) - 1) * xb_width)
        * 2.0
        * t_jacobian
        / ((len(binning["minus_t"]) - 1) * t_width)
        * 360.0
        / ((len(binning["phi_deg"]) - 1) * phi_width)
    )
    return 1.0 / (
        legacy_fraction + (1.0 - legacy_fraction) * balanced_to_legacy
    )


def _validate_legacy_input_shape(text: str, path: Path) -> None:
    records = _records(text)
    expected = _legacy_record_count(records, path)
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
    run.add_argument(
        "--proposal",
        choices=("legacy", "balanced"),
        default="legacy",
        help=(
            "legacy proposal or a full-support mixture whose balanced "
            "component allocates analysis-coordinate bins uniformly"
        ),
    )
    run.add_argument(
        "--config",
        type=Path,
        help="analysis config required by --proposal balanced",
    )
    run.add_argument(
        "--legacy-fraction",
        type=float,
        default=0.25,
        help=(
            "full-support legacy-mixture probability for a balanced survey "
            "(default: 0.25)"
        ),
    )

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
                proposal=args.proposal,
                config_path=args.config,
                legacy_fraction=args.legacy_fraction,
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
