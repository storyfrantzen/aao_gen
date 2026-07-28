#!/usr/bin/env python3
"""Unit tests for the milestone-1 radiative survey utilities."""

from __future__ import annotations

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
