#!/usr/bin/env python3
"""Tests for the final-electron momentum retention audit."""

from __future__ import annotations

import argparse
import json
import math
import tempfile
import unittest
from pathlib import Path

import radiative_momentum_audit as audit
from test_radiative_guards import _row, _write_survey


class MomentumAuditTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.config = self.root / "config.json"
        self.config.write_text(
            json.dumps(
                {
                    "beam_energy": 6.535,
                    "phase_space": {
                        "Q2_min": 1.0,
                        "W_min": 2.0,
                        "y_max": 0.95,
                        "electron_p_min": 1.0,
                    },
                    "binning": {
                        "Q2": [1.0, 3.0],
                        "xB": [0.1, 0.3],
                        "minus_t": [0.09, 0.5],
                        "phi_deg": [0.0, 360.0],
                    },
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _arguments(
        self,
        *,
        surveys: list[Path] | None = None,
        training: list[Path] | None = None,
        validation: list[Path] | None = None,
        output: str = "audit",
    ) -> argparse.Namespace:
        return argparse.Namespace(
            config=self.config,
            surveys=surveys,
            training_surveys=training,
            validation_surveys=validation,
            electron_p_min=None,
            apply_y_max=True,
            flat_index_file=None,
            low_nominal_fraction=0.5,
            output=self.root / output,
        )

    def test_fraction_uses_cross_section_contributions_not_row_counts(self) -> None:
        survey = self.root / "survey"
        _write_survey(
            survey,
            replica=0,
            rows=[
                _row(
                    replica=0,
                    trial=1,
                    weight=3.0,
                    r_u=0.1,
                    electron_momentum=1.2,
                    electron_theta_deg=35.0,
                ),
                _row(
                    replica=0,
                    trial=2,
                    weight=1.0,
                    r_u=0.2,
                    electron_momentum=0.8,
                    electron_theta_deg=35.0,
                ),
            ],
        )

        result = audit.audit(self._arguments(surveys=[survey]))
        payload = json.loads(result.read_text(encoding="utf-8"))
        record = payload["strata"][0]
        metrics = record["pooled"]
        self.assertEqual(record["stratum_id"], "s00000")
        self.assertAlmostEqual(
            metrics["nominal_cross_section_fraction"], 0.75
        )
        self.assertAlmostEqual(
            metrics["raw_nominal_row_fraction_diagnostic_only"], 0.5
        )
        expected_sem = math.sqrt(100.0 / 99.0 * 1.125 / 16.0)
        self.assertAlmostEqual(
            metrics["nominal_fraction_sem_delta_method"], expected_sem
        )
        self.assertAlmostEqual(
            metrics["loose_events_per_expected_nominal_event"], 4.0 / 3.0
        )
        self.assertAlmostEqual(
            metrics["denominator"]["cross_section_microbarn"], 0.04
        )
        self.assertAlmostEqual(
            metrics["nominal"]["cross_section_microbarn"], 0.03
        )
        self.assertTrue(
            payload["interpretation"][
                "constant_per_stratum_weights_remain_valid"
            ]
        )

    def test_training_validation_split_is_pooled_by_fixed_trials(self) -> None:
        training = self.root / "training"
        validation = self.root / "validation"
        _write_survey(
            training,
            replica=0,
            rows=[
                _row(
                    replica=0,
                    trial=1,
                    weight=3.0,
                    r_u=0.1,
                    electron_momentum=1.2,
                    electron_theta_deg=35.0,
                ),
                _row(
                    replica=0,
                    trial=2,
                    weight=1.0,
                    r_u=0.2,
                    electron_momentum=0.8,
                    electron_theta_deg=35.0,
                ),
            ],
        )
        _write_survey(
            validation,
            replica=1,
            rows=[
                _row(
                    replica=1,
                    trial=1,
                    weight=1.0,
                    r_u=0.3,
                    electron_momentum=1.2,
                    electron_theta_deg=35.0,
                ),
                _row(
                    replica=1,
                    trial=2,
                    weight=1.0,
                    r_u=0.4,
                    electron_momentum=0.8,
                    electron_theta_deg=35.0,
                ),
            ],
        )

        result = audit.audit(
            self._arguments(
                training=[training], validation=[validation], output="split"
            )
        )
        payload = json.loads(result.read_text(encoding="utf-8"))
        record = payload["strata"][0]
        self.assertAlmostEqual(
            record["training"]["nominal_cross_section_fraction"], 0.75
        )
        self.assertAlmostEqual(
            record["validation"]["nominal_cross_section_fraction"], 0.5
        )
        self.assertAlmostEqual(
            record["pooled"]["nominal_cross_section_fraction"], 4.0 / 6.0
        )
        self.assertEqual(payload["pooling"]["proposals"], 200)
        self.assertIsNotNone(
            record["training_validation_fraction_difference_z_score"]
        )


if __name__ == "__main__":
    unittest.main()
