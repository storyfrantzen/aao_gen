#!/usr/bin/env python3
"""Unit tests for the milestone-2 radiative guard learner."""

from __future__ import annotations

import argparse
import csv
import json
import math
import tempfile
import unittest
from pathlib import Path

import radiative_guards
import radiative_survey


def _config() -> dict:
    return {
        "beam_energy": 6.535,
        "phase_space": {"Q2_min": 1.0, "W_min": 2.0, "y_max": 0.8},
        "binning": {
            "Q2": [1.0, 2.0, 3.0],
            "xB": [0.1, 0.3],
            "minus_t": [0.09, 0.5],
            "phi_deg": [0.0, 180.0, 360.0],
        },
    }


def _row(
    *,
    replica: int,
    trial: int,
    weight: float,
    r_u: float,
    hadron_phi_base: float = 0.1,
    phi_observed_deg: float = 0.0,
    final_valid: int = 1,
    w_observed: float = 2.5,
    y_observed: float = 0.3,
) -> dict[str, str | int | float]:
    beam_energy = 6.535
    electron_momentum = 3.0
    electron_theta = math.radians(15.0)
    electron = (
        electron_momentum * math.sin(electron_theta),
        0.0,
        electron_momentum * math.cos(electron_theta),
        electron_momentum,
    )
    q = (-electron[0], 0.0, beam_energy - electron[2])
    proton_momentum = (-q[2] * 0.1, 0.0, q[0] * 0.1)
    proton_energy = math.sqrt(
        sum(component * component for component in proton_momentum)
        + radiative_survey.PROTON_MASS_GEV**2
    )
    proton = (*proton_momentum, proton_energy)
    observed = radiative_survey.observed_coordinates(
        beam_energy, electron, proton
    )
    row: dict[str, str | int | float] = {
        name: 0 for name in radiative_survey.SURVEY_COLUMNS
    }
    row.update(
        {
            "replica": replica,
            "trial": trial,
            "final_valid": final_valid,
            "candidate_status": 0 if final_valid else 1,
            "intreg": 1,
            "r_u": r_u,
            "r_ep": 0.2,
            "u_gamma": 0.5,
            "photon_cosine_base": 0.2,
            "photon_phi_base": 0.2,
            "hadron_cosine_base": 0.2,
            "hadron_phi_base": hadron_phi_base,
            "final_e_px": electron[0],
            "final_e_py": electron[1],
            "final_e_pz": electron[2],
            "final_e_energy": electron[3],
            "final_p_px": proton[0],
            "final_p_py": proton[1],
            "final_p_pz": proton[2],
            "final_p_energy": proton[3],
            "q2_observed": observed["q2_observed"],
            "xb_observed": observed["xb_observed"],
            "minus_t_observed": observed["minus_t_observed"],
            "phi_observed_deg": phi_observed_deg,
            "w_observed": (
                w_observed if w_observed != 2.5 else observed["w_observed"]
            ),
            "y_observed": (
                y_observed if y_observed != 0.3 else observed["y_observed"]
            ),
            "integrand_internal": max(weight, 1.0),
            "integrand_observed": weight if final_valid else 0.0,
            "trial_xsec_internal_microbarn": max(weight, 1.0),
            "trial_xsec_observed_microbarn": weight if final_valid else 0.0,
        }
    )
    return row


def _write_survey(
    directory: Path,
    *,
    replica: int,
    rows: list[dict[str, str | int | float]],
    ntrials: int = 100,
) -> None:
    directory.mkdir()
    total = sum(float(row["trial_xsec_observed_microbarn"]) for row in rows)
    valid = sum(int(row["final_valid"]) for row in rows)
    norm = {
        "generator": "aao_rad",
        "sampling_mode": "1",
        "fixed_trial_survey": "1",
        "survey_schema": radiative_survey.SURVEY_SCHEMA,
        "survey_emits_lund": "0",
        "ntries": str(ntrials),
        "survey_ntrials_requested": str(ntrials),
        "survey_rows": str(len(rows)),
        "survey_final_valid": str(valid),
        "survey_replica": str(replica),
        "survey_seed": str(371001 + 2 * replica),
        "survey_observed_sig_sum": repr(total / ntrials),
        "ebeam": "6.535",
        "q2_min": "1.0",
        "q2_max": "3.0",
        "ep_min": "0.2",
        "ep_max_effective": "6.0",
        "delta": "0.005",
        "epirea": "1",
        "th_opt": "5",
        "res_opt": "0",
        "survey_phase_volume": "1.0",
    }
    (directory / radiative_survey.NORM_FILENAME).write_text(
        "".join(f"{key}={value}\n" for key, value in norm.items()),
        encoding="utf-8",
    )
    legacy = [
        "5",
        "0",
        ".20 .12 .20 .20",
        "4",
        "1",
        ".2",
        "5",
        ".43",
        "0",
        "0",
        "0",
        "6.535",
        "1 3",
        ".2 6.535",
        ".005",
        "100",
        "2",
    ]
    trailer = ["1", str(ntrials), str(371001 + 2 * replica), str(replica)]
    (directory / "survey_input.inp").write_text(
        "\n".join(legacy + trailer) + "\n", encoding="utf-8"
    )
    with (directory / radiative_survey.SURVEY_FILENAME).open(
        "w", encoding="utf-8", newline=""
    ) as output:
        output.write(f"# schema={radiative_survey.SURVEY_SCHEMA}\n")
        writer = csv.DictWriter(output, fieldnames=radiative_survey.SURVEY_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


class StratumTests(unittest.TestCase):
    def test_phi_wraps_and_upper_nonperiodic_edges_are_outside(self) -> None:
        config = _config()
        wrapped = _row(
            replica=0,
            trial=1,
            weight=1.0,
            r_u=0.1,
            phi_observed_deg=360.0,
        )
        identifier, reason = radiative_guards.assign_stratum(wrapped, config)
        self.assertEqual(identifier, "s00000")
        self.assertIsNone(reason)

        wrapped["q2_observed"] = 3.0
        identifier, reason = radiative_guards.assign_stratum(wrapped, config)
        self.assertIsNone(identifier)
        self.assertEqual(reason, "outside_analysis_binning")

    def test_y_is_diagnostic_by_default_and_has_an_explicit_opt_in_cut(self) -> None:
        config = _config()
        row = _row(replica=0, trial=1, weight=1.0, r_u=0.1, y_observed=-0.1)
        self.assertEqual(
            radiative_guards.assign_stratum(row, config), ("s00000", None)
        )
        row["y_observed"] = 0.9
        self.assertEqual(
            radiative_guards.assign_stratum(row, config), ("s00000", None)
        )
        self.assertEqual(
            radiative_guards.assign_stratum(
                row, config, apply_y_max=True
            ),
            (None, "failed_analysis_phase_space"),
        )
        row["w_observed"] = 1.9
        self.assertEqual(
            radiative_guards.assign_stratum(row, config),
            (None, "failed_analysis_phase_space"),
        )


class CellTests(unittest.TestCase):
    def test_periodic_axis_dilation_wraps(self) -> None:
        partition = radiative_guards.parse_partition(
            ["r_u=4", "hadron_phi_base=4"]
        )
        seed = {(1, radiative_guards.flatten_cell((3, 0), partition))}
        expanded = radiative_guards.dilate_cells(seed, partition, 1)
        self.assertIn(
            (1, radiative_guards.flatten_cell((3, 3), partition)), expanded
        )
        self.assertIn(
            (1, radiative_guards.flatten_cell((3, 1), partition)), expanded
        )
        self.assertNotIn(
            (1, radiative_guards.flatten_cell((0, 0), partition)), expanded
        )


class WorkflowTests(unittest.TestCase):
    def test_pooling_uses_all_fixed_trials_as_the_denominator(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first"
            second = root / "second"
            _write_survey(
                first,
                replica=0,
                rows=[_row(replica=0, trial=1, weight=10.0, r_u=0.5)],
            )
            _write_survey(
                second,
                replica=1,
                rows=[_row(replica=1, trial=1, weight=10.0, r_u=0.5)],
            )
            campaign = radiative_guards.aggregate_surveys(
                [first, second],
                _config(),
                radiative_guards.parse_partition(["r_u=4"]),
            )
            self.assertEqual(campaign.proposals, 200)
            self.assertEqual(campaign.global_observed.count, 2)
            self.assertAlmostEqual(
                campaign.global_observed.total / campaign.proposals, 0.1
            )
            self.assertAlmostEqual(
                radiative_guards._ess(campaign.global_observed), 2.0
            )

    def test_weighted_learning_and_independent_validation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "config.json"
            config_path.write_text(json.dumps(_config()), encoding="utf-8")
            training = root / "training"
            training_rows = [
                _row(replica=0, trial=index, weight=1.0, r_u=0.1)
                for index in range(1, 5)
            ]
            training_rows.append(
                _row(replica=0, trial=5, weight=10.0, r_u=0.9)
            )
            _write_survey(training, replica=0, rows=training_rows)

            output = root / "guard"
            result = radiative_guards.learn_guards(
                argparse.Namespace(
                    config=config_path,
                    survey=[training],
                    output=output,
                    target_core_fraction=0.7,
                    core_probability=0.98,
                    uniform_cell_floor=0.05,
                    dilation=1,
                    iteration=0,
                    minimum_training_rows=1,
                    minimum_training_ess=0.0,
                    apply_y_max=False,
                    partition=["r_u=4", "hadron_phi_base=4"],
                    generator_revision="test-revision",
                )
            )
            self.assertTrue(result["passed"])
            manifest_path = output / "guard_manifest.json"
            self.assertTrue((output / "training_cells.csv").is_file())
            self.assertTrue((output / "training_cells.csv.sha256").is_file())
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            record = manifest["strata"]["s00000"]
            partition = radiative_guards._partition_from_manifest(manifest)
            heavy_cell = radiative_guards.flatten_cell((3, 0), partition)
            light_cell = radiative_guards.flatten_cell((0, 0), partition)
            self.assertIn(
                heavy_cell,
                record["core_cells"]["seed_cell_ids_by_channel"]["1"],
            )
            self.assertNotIn(
                light_cell,
                record["core_cells"]["seed_cell_ids_by_channel"]["1"],
            )
            self.assertGreater(record["estimated_core_fraction"], 0.7)
            self.assertGreater(
                manifest["proposal_mixture"]["tail_probability"], 0.0
            )

            validation = root / "validation"
            _write_survey(
                validation,
                replica=1,
                rows=[
                    _row(replica=1, trial=1, weight=8.0, r_u=0.9),
                    _row(replica=1, trial=2, weight=2.0, r_u=0.3),
                ],
            )
            validation_output = root / "validation_output"
            checked = radiative_guards.validate_guards(
                argparse.Namespace(
                    manifest=manifest_path,
                    survey=[validation],
                    output=validation_output,
                    minimum_core_fraction=0.75,
                )
            )
            self.assertTrue(checked["passed"])
            report = json.loads(
                (validation_output / "guard_validation.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertAlmostEqual(
                report["strata"]["s00000"]["holdout_core_fraction"], 0.8
            )
            self.assertEqual(report["validation"]["replica_ids"], [1])
            self.assertTrue(
                (validation_output / "validation_cells.csv.sha256").is_file()
            )
            summary = radiative_guards.summarize_coverage(
                argparse.Namespace(
                    manifest=manifest_path,
                    validation=validation_output / "guard_validation.json",
                    limit=1,
                )
            )
            self.assertTrue(summary["validation_passed"])
            self.assertEqual(
                summary["worst_holdout_core_fractions"][0]["stratum_id"],
                "s00000",
            )

    def test_training_replica_cannot_be_reused_for_validation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "config.json"
            config_path.write_text(json.dumps(_config()), encoding="utf-8")
            survey = root / "survey"
            _write_survey(
                survey,
                replica=0,
                rows=[_row(replica=0, trial=1, weight=1.0, r_u=0.5)],
            )
            output = root / "guard"
            radiative_guards.learn_guards(
                argparse.Namespace(
                    config=config_path,
                    survey=[survey],
                    output=output,
                    target_core_fraction=0.9,
                    core_probability=0.98,
                    uniform_cell_floor=0.05,
                    dilation=1,
                    iteration=0,
                    minimum_training_rows=1,
                    minimum_training_ess=0.0,
                    apply_y_max=False,
                    partition=["r_u=4"],
                    generator_revision="test",
                )
            )
            manifest = json.loads(
                (output / "guard_manifest.json").read_text(encoding="utf-8")
            )
            campaign = radiative_guards.aggregate_surveys(
                [survey], _config(), radiative_guards.parse_partition(["r_u=4"])
            )
            with self.assertRaisesRegex(
                radiative_guards.GuardLearningError, "overlap"
            ):
                radiative_guards.validate_manifest(
                    manifest,
                    output / "guard_manifest.json",
                    campaign,
                    radiative_guards.parse_partition(["r_u=4"]),
                    0.5,
                )


if __name__ == "__main__":
    unittest.main()
