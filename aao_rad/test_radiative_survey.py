#!/usr/bin/env python3
"""Unit tests for the milestone-1 radiative survey utilities."""

from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

import radiative_survey


class ObservedCoordinateTests(unittest.TestCase):
    def test_coordinates_match_direct_invariants_and_zero_trento_angle(self) -> None:
        beam_energy = 6.535
        electron_momentum = 3.0
        electron_theta = math.radians(20.0)
        electron = (
            electron_momentum * math.sin(electron_theta),
            0.0,
            electron_momentum * math.cos(electron_theta),
            electron_momentum,
        )
        q = (-electron[0], 0.0, beam_energy - electron[2])
        proton_momentum = (-q[2] * 0.2, 0.0, q[0] * 0.2)
        proton_energy = math.sqrt(
            sum(component * component for component in proton_momentum)
            + radiative_survey.PROTON_MASS_GEV**2
        )
        proton = (*proton_momentum, proton_energy)

        result = radiative_survey.observed_coordinates(
            beam_energy, electron, proton
        )

        q_energy = beam_energy - electron[3]
        expected_q2 = sum(component * component for component in q) - q_energy**2
        expected_t = sum(component * component for component in proton_momentum) - (
            radiative_survey.PROTON_MASS_GEV - proton_energy
        ) ** 2
        self.assertAlmostEqual(result["q2_observed"], expected_q2)
        self.assertAlmostEqual(
            result["xb_observed"],
            expected_q2
            / (2.0 * radiative_survey.PROTON_MASS_GEV * q_energy),
        )
        self.assertAlmostEqual(result["minus_t_observed"], expected_t)
        self.assertAlmostEqual(result["phi_observed_deg"], 0.0)


class InputValidationTests(unittest.TestCase):
    def test_accepts_legacy_fmcall_zero_input(self) -> None:
        text = "\n".join(
            [
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
                "0",
                ".005",
            ]
        )
        radiative_survey._validate_legacy_input_shape(text, Path("input.inp"))

    def test_rejects_an_existing_survey_trailer(self) -> None:
        text = "\n".join(
            [
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
                "1",
                "1",
                "1000",
                "371001",
            ]
        )
        with self.assertRaisesRegex(ValueError, "without an optional"):
            radiative_survey._validate_legacy_input_shape(text, Path("input.inp"))

    def test_balanced_trailer_round_trips_analysis_edges(self) -> None:
        config = {
            "target_mass": radiative_survey.PROTON_MASS_GEV,
            "binning": {
                "Q2": [1.0, 2.0, 3.0],
                "xB": [0.1, 0.2, 0.4],
                "minus_t": [0.09, 0.3, 1.0],
                "phi_deg": [0.0, 180.0, 360.0],
            },
        }
        trailer = radiative_survey._survey_trailer(
            1000,
            371001,
            0,
            proposal="balanced",
            legacy_fraction=0.25,
            balanced_config=config,
        )
        legacy = "\n".join(
            [
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
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "survey_input.inp"
            path.write_text(legacy + "\n" + trailer, encoding="utf-8")
            parsed = radiative_survey._parse_survey_input(path)

        self.assertEqual(parsed["mode"], 1)
        self.assertEqual(parsed["trials"], 1000)
        self.assertEqual(parsed["binning"], config["binning"])
        self.assertAlmostEqual(parsed["legacy_fraction"], 0.25)


class ProposalDensityTests(unittest.TestCase):
    def test_declared_bin_allows_only_small_boundary_roundoff(self) -> None:
        edges = [0.09, 1.0, 1.5, 2.0]

        self.assertTrue(
            radiative_survey._declared_bin_contains(
                1.5 - 2.0e-6,
                edges,
                2,
                tolerance=3.0e-5,
            )
        )
        self.assertTrue(
            radiative_survey._declared_bin_contains(
                2.0 + 2.0e-6,
                edges,
                2,
                tolerance=3.0e-5,
            )
        )
        self.assertFalse(
            radiative_survey._declared_bin_contains(
                1.49,
                edges,
                2,
                tolerance=3.0e-5,
            )
        )
        self.assertTrue(
            radiative_survey._declared_bin_contains(
                360.0,
                [0.0, 180.0, 360.0],
                1,
                tolerance=3.0e-4,
                periodic=True,
            )
        )

    def test_global_tail_bounds_density_ratio_outside_balanced_support(self) -> None:
        row = {
            "q2_leptonic": 4.0,
            "xb_leptonic": 0.25,
            "minus_t_hard": 0.2,
            "phi_cm_deg": 90.0,
            "energy_in_vertex": 6.0,
            "energy_e_pre_external": 3.0,
            "energy_gamma": 0.05,
            "cos_theta_gamma": 0.5,
        }
        norm = {
            "q2_min": "1.0",
            "q2_max": "5.0",
            "ep_min": "0.2",
            "ep_max_effective": "6.0",
        }
        spec = {
            "mode": 1,
            "legacy_fraction": 0.25,
            "binning": {
                "Q2": [1.0, 2.0, 3.0],
                "xB": [0.1, 0.2, 0.4],
                "minus_t": [0.09, 0.3, 1.0],
                "phi_deg": [0.0, 180.0, 360.0],
            },
        }

        self.assertAlmostEqual(
            radiative_survey._proposal_density_ratio(row, norm, spec),
            4.0,
        )

    def test_balanced_config_is_frozen_without_reformatting(self) -> None:
        payload = {
            "beam_energy": 6.535,
            "target_mass": radiative_survey.PROTON_MASS_GEV,
            "binning": {
                "Q2": [1.0, 2.0, 3.0],
                "xB": [0.1, 0.2, 0.4],
                "minus_t": [0.09, 0.3, 1.0],
                "phi_deg": [0.0, 180.0, 360.0],
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.json"
            path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
            _, raw, digest = radiative_survey._load_balanced_config(path)

        self.assertEqual(raw, (json.dumps(payload, indent=2) + "\n").encode())
        self.assertEqual(len(digest), 64)

    def test_allocation_summary_separates_target_and_observed_bins(self) -> None:
        spec = {
            "binning": {
                "Q2": [1.0, 2.0, 3.0],
                "xB": [0.1, 0.2, 0.4],
                "minus_t": [0.09, 0.3, 1.0],
                "phi_deg": [0.0, 180.0, 360.0],
            }
        }
        rows = [
            {
                "proposal_component": 1,
                "proposal_q2_bin": 1,
                "proposal_xb_bin": 0,
                "proposal_t_bin": 1,
                "proposal_phi_bin": 0,
                "final_valid": 1,
                "q2_observed": 1.5,
                "xb_observed": 0.3,
                "minus_t_observed": 0.2,
                "phi_observed_deg": 270.0,
            },
            {
                "proposal_component": 0,
                "proposal_q2_bin": -1,
                "proposal_xb_bin": -1,
                "proposal_t_bin": -1,
                "proposal_phi_bin": -1,
                "final_valid": 0,
                "q2_observed": 0.0,
                "xb_observed": 0.0,
                "minus_t_observed": 0.0,
                "phi_observed_deg": 0.0,
            },
        ]

        summary = radiative_survey._allocation_summary(rows, spec)

        self.assertEqual(
            summary["balanced_target_axis_counts_in_recorded_internal_rows"][
                "Q2"
            ],
            [0, 1],
        )
        self.assertEqual(summary["final_observed_axis_counts"]["Q2"], [1, 0])
        self.assertEqual(summary["final_observed_joint_strata_occupied"], 1)
        self.assertEqual(summary["analysis_joint_strata_total"], 16)


class SchemaTests(unittest.TestCase):
    def test_rejects_unknown_schema_columns(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / radiative_survey.SURVEY_FILENAME
            path.write_text(
                "# schema=aao-rad-survey-v1\ntrial,not_a_real_column\n1,2\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                radiative_survey.SurveyValidationError, "unexpected columns"
            ):
                radiative_survey.read_survey(path)


if __name__ == "__main__":
    unittest.main()
