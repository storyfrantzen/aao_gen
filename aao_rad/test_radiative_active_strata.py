#!/usr/bin/env python3
"""Tests for the evidence-audited active-stratum workflow."""

from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

import radiative_active_strata as active
import radiative_mode4


class ActiveStratumWorkflowTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.config = self.root / "config.json"
        self.config.write_text(
            json.dumps(
                {
                    "beam_energy": 6.535,
                    "phase_space": {"W_min": 2.0},
                    "binning": {
                        "Q2": [1.0, 2.0],
                        "xB": [0.2, 0.3],
                        "minus_t": [0.1, 0.2],
                        "phi_deg": [0.0, 360.0],
                    },
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        self.config_sha256 = hashlib.sha256(
            self.config.read_bytes()
        ).hexdigest()
        self.source = self.root / "analysis_source.json"
        self.source.write_text('{"evidence": true}\n', encoding="utf-8")

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _template(self, name: str = "relevance.json") -> Path:
        output = self.root / name
        active.create_template(
            argparse.Namespace(
                config=self.config,
                output=output,
                flat_indices=[0],
                evidence_sources=[f"analysis={self.source}"],
            )
        )
        return output

    def _claim(self, path: Path, **updates: object) -> None:
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload["strata"][0].update(
            {
                "rationale": "Test evidence for the selected stratum.",
                "source_ids": ["analysis"],
                **updates,
            }
        )
        path.write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )

    def _calibration(
        self, *, positive: bool, ready: bool, name: str = "calibration.json"
    ) -> Path:
        path = self.root / name
        sigma = 2.0e-8 if positive else 0.0
        targets = 12 if positive else 0
        payload = {
            "schema": radiative_mode4.CALIBRATION_SCHEMA,
            "analysis_config_sha256": self.config_sha256,
            "strata": [
                {
                    "stratum_id": "s00000",
                    "flat_index": 0,
                    "indices": {
                        "iq2": 0,
                        "ixb": 0,
                        "it": 0,
                        "iphi": 0,
                    },
                    "integrated_cross_section_microbarn": sigma,
                    "integrated_cross_section_sem_microbarn": 1.0e-9,
                    "inside_guard": {"target_candidates": targets},
                    "guard_complement": {"target_candidates": 0},
                    "recommendation_status": (
                        "provisional_zero_complement"
                        if ready
                        else "insufficient_inside_guard_target_support"
                    ),
                    "pilot_readiness": (
                        "ready_provisional_zero_complement"
                        if ready
                        else "not_ready"
                    ),
                    "recommended_envelope": (
                        {"sigr_max": 4.0e-8} if ready else None
                    ),
                }
            ],
        }
        path.write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )
        return path

    def _validation(self, calibration: Path) -> Path:
        path = self.root / "validation.json"
        payload = {
            "schema": radiative_mode4.PILOT_VALIDATION_SCHEMA,
            "calibration_report": str(calibration.resolve()),
            "calibration_report_sha256": hashlib.sha256(
                calibration.read_bytes()
            ).hexdigest(),
            "strata": [
                {
                    "stratum_id": "s00000",
                    "flat_index": 0,
                    "passed": True,
                    "recommendation": "ready_for_multi_stratum_pilot",
                    "relative_cross_section_difference": 0.02,
                    "cross_section_difference_z_score": 0.4,
                }
            ],
        }
        path.write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )
        return path

    def _classify(
        self,
        relevance: Path,
        *,
        calibrations: list[Path] | None = None,
        validations: list[Path] | None = None,
        output_name: str = "classification",
    ) -> dict:
        result = active.classify(
            argparse.Namespace(
                config=self.config,
                relevance=relevance,
                calibrations=calibrations,
                pilot_validations=validations,
                output=self.root / output_name,
                flat_indices=[0],
                minimum_data_events=1,
                maximum_model_fraction=0.001,
                maximum_feed_in_fraction=0.001,
                maximum_closure_fraction=0.001,
                reference_cross_section_microbarn=None,
            )
        )
        return json.loads(result.read_text(encoding="utf-8"))

    def _write_data_census_inputs(self) -> tuple[Path, Path, Path]:
        config = self.root / "data_config.json"
        config.write_text(
            json.dumps(
                {
                    "beam_energy": 6.535,
                    "target_mass": 0.9382720813,
                    "phase_space": {
                        "Q2_min": 1.0,
                        "W_min": 2.0,
                        "y_max": 0.95,
                    },
                    "binning": {
                        "Q2": [1.0, 2.0, 3.0],
                        "xB": [0.2, 0.4],
                        "minus_t": [0.1, 0.2],
                        "phi_deg": [0.0, 180.0, 360.0],
                    },
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        data = self.root / "data_events.npz"
        np.savez_compressed(
            data,
            run=np.ones(9, dtype=np.int32),
            event=np.arange(100, 109, dtype=np.int64),
            rec_Q2=np.asarray(
                [1.5, 1.5, 2.5, 2.5, 2.5, 1.1, 0.9, 1.5, np.nan]
            ),
            rec_xB=np.asarray(
                [0.3, 0.3, 0.3, 0.3, 0.2, 0.39, 0.25, 0.3, 0.3]
            ),
            rec_minus_t=np.asarray(
                [0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.25, 0.15]
            ),
            rec_trento_phi=np.asarray(
                [0.1, 3.5, 0.1, 3.5, 0.1, 0.1, 0.1, 0.1, 0.1]
            ),
        )
        mask = self.root / "selection.npy"
        values = np.ones(9, dtype=bool)
        values[3] = False
        np.save(mask, values)
        return config, data, mask

    def test_template_hashes_external_sources(self) -> None:
        template = self._template()
        payload = json.loads(template.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema"], active.EVIDENCE_SCHEMA)
        self.assertEqual(payload["strata"][0]["physical_status"], "unknown")
        self.assertEqual(
            payload["sources"][0]["sha256"],
            hashlib.sha256(self.source.read_bytes()).hexdigest(),
        )

    def test_data_census_builds_occupancy_and_relevance_evidence(self) -> None:
        config, data, mask = self._write_data_census_inputs()
        output = self.root / "data_census"
        relevance_path = active.build_data_evidence(
            argparse.Namespace(
                config=config,
                data_events=data,
                selection_mask=mask,
                output=output,
                allow_duplicate_event_keys=False,
            )
        )
        occupancy = json.loads(
            (output / "data_occupancy.json").read_text(encoding="utf-8")
        )
        self.assertEqual(occupancy["schema"], active.DATA_OCCUPANCY_SCHEMA)
        self.assertEqual(
            occupancy["catalog_summary"],
            {
                "strata_total": 4,
                "selected_data_events": 3,
                "occupied_strata": 3,
                "strata_with_at_least_5_events": 0,
                "strata_with_at_least_10_events": 0,
                "strata_with_at_least_50_events": 0,
            },
        )
        self.assertEqual(
            [item["rows_remaining"] for item in occupancy["cut_flow"]],
            [9, 8, 7, 7, 6, 5, 4, 3],
        )
        self.assertEqual(
            [record["data_events"] for record in occupancy["strata"]],
            [1, 1, 1, 0],
        )
        self.assertEqual(
            (output / "occupied_flat_indices.txt").read_text(
                encoding="utf-8"
            ),
            "0\n1\n2\n",
        )
        self.assertEqual(
            (output / "at_least_5_events_flat_indices.txt").read_text(
                encoding="utf-8"
            ),
            "",
        )
        relevance = json.loads(relevance_path.read_text(encoding="utf-8"))
        self.assertEqual(relevance["schema"], active.EVIDENCE_SCHEMA)
        self.assertEqual(
            [record["data_events"] for record in relevance["strata"]],
            [1, 1, 1, 0],
        )
        self.assertTrue(
            all(
                record["physical_status"] == "unknown"
                and record["analysis_included"] is None
                for record in relevance["strata"]
            )
        )
        classification = active.classify(
            argparse.Namespace(
                config=config,
                relevance=relevance_path,
                calibrations=None,
                pilot_validations=None,
                output=self.root / "data_classification",
                flat_indices=None,
                minimum_data_events=1,
                maximum_model_fraction=0.001,
                maximum_feed_in_fraction=0.001,
                maximum_closure_fraction=0.001,
                reference_cross_section_microbarn=None,
            )
        )
        classified = json.loads(classification.read_text(encoding="utf-8"))
        self.assertEqual(
            [record["status"] for record in classified["strata"]],
            [
                "active_needs_optimization",
                "active_needs_optimization",
                "active_needs_optimization",
                "needs_relevance_assessment",
            ],
        )

    def test_data_census_rejects_duplicate_selected_event_keys(self) -> None:
        config, data, mask = self._write_data_census_inputs()
        with np.load(data, allow_pickle=False) as sample:
            arrays = {name: sample[name] for name in sample.files}
        arrays["event"][1] = arrays["event"][0]
        np.savez_compressed(data, **arrays)
        with self.assertRaisesRegex(
            active.ActiveStratumError, "duplicate \\(run,event\\) keys"
        ):
            active.build_data_evidence(
                argparse.Namespace(
                    config=config,
                    data_events=data,
                    selection_mask=mask,
                    output=self.root / "duplicate_census",
                    allow_duplicate_event_keys=False,
                )
            )

    def test_zero_calibration_never_implies_structural_emptiness(self) -> None:
        template = self._template()
        calibration = self._calibration(positive=False, ready=False)
        payload = self._classify(template, calibrations=[calibration])
        record = payload["strata"][0]
        self.assertEqual(record["status"], "needs_relevance_assessment")
        self.assertEqual(record["physical_status"], "unknown")
        self.assertFalse(record["active"])

    def test_changed_external_evidence_is_rejected(self) -> None:
        template = self._template()
        self.source.write_text('{"evidence": false}\n', encoding="utf-8")
        with self.assertRaisesRegex(
            active.ActiveStratumError, "evidence source analysis changed"
        ):
            self._classify(template)

    def test_complete_negligible_evidence_becomes_closure_only(self) -> None:
        template = self._template()
        self._claim(
            template,
            physical_status="partially_accessible",
            analysis_included=False,
            data_events=0,
            model_cross_section_fraction=1.0e-5,
            maximum_feed_in_fraction=2.0e-5,
            global_closure_impact_fraction=3.0e-5,
        )
        calibration = self._calibration(positive=True, ready=False)
        payload = self._classify(template, calibrations=[calibration])
        record = payload["strata"][0]
        self.assertEqual(record["status"], "closure_only")
        self.assertFalse(record["active"])

    def test_relevant_ready_stratum_becomes_active_ready(self) -> None:
        template = self._template()
        self._claim(
            template,
            physical_status="nonempty",
            analysis_included=True,
        )
        calibration = self._calibration(positive=True, ready=True)
        validation = self._validation(calibration)
        payload = self._classify(
            template,
            calibrations=[calibration],
            validations=[validation],
        )
        record = payload["strata"][0]
        self.assertEqual(record["status"], "active_ready")
        self.assertTrue(record["active"])
        self.assertTrue(record["production_ready"])
        ready = self.root / "classification/production_ready_flat_indices.txt"
        self.assertEqual(ready.read_text(encoding="utf-8"), "0\n")

    def test_relevant_stratum_preserves_readiness_state(self) -> None:
        template = self._template()
        self._claim(
            template,
            physical_status="nonempty",
            analysis_included=True,
        )
        unready = self._calibration(positive=True, ready=False)
        first = self._classify(
            template,
            calibrations=[unready],
            output_name="unready",
        )["strata"][0]
        self.assertEqual(first["status"], "active_needs_optimization")

        ready = self._calibration(
            positive=True, ready=True, name="ready_calibration.json"
        )
        second = self._classify(
            template,
            calibrations=[ready],
            output_name="needs_pilot",
        )["strata"][0]
        self.assertEqual(
            second["status"], "active_needs_pilot_validation"
        )

    def test_structural_empty_claim_rejects_positive_calibration(self) -> None:
        template = self._template()
        self._claim(
            template,
            physical_status="structurally_empty",
            analysis_included=False,
            data_events=0,
            model_cross_section_fraction=0.0,
            maximum_feed_in_fraction=0.0,
            global_closure_impact_fraction=0.0,
        )
        calibration = self._calibration(positive=True, ready=False)
        with self.assertRaisesRegex(
            active.ActiveStratumError, "contradicts a positive calibration"
        ):
            self._classify(template, calibrations=[calibration])


if __name__ == "__main__":
    unittest.main()
