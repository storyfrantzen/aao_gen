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

    def _write_survey_evidence_inputs(
        self, config: Path
    ) -> tuple[Path, Path]:
        config_payload, config_sha256 = active.radiative_guards.load_analysis_config(
            config
        )
        catalog = active.radiative_guards.enumerate_strata(config_payload)
        proposals = 100

        def metrics(
            count: int, total: float, square: float, maximum: float
        ) -> dict:
            return active.radiative_guards._metrics(
                active.radiative_guards.Moment(
                    count=count,
                    total=total,
                    square_total=square,
                    maximum=maximum,
                ),
                proposals,
            )

        training_values = (
            metrics(1, 3.0, 9.0, 3.0),
            metrics(1, 2.0, 4.0, 2.0),
            metrics(0, 0.0, 0.0, 0.0),
            metrics(0, 0.0, 0.0, 0.0),
        )
        manifest = self.root / "migration_manifest.json"
        manifest_payload = {
            "schema": active.radiative_migrations.MANIFEST_SCHEMA,
            "analysis_config_sha256": config_sha256,
            "analysis_selection": {
                "q2_minimum": 1.0,
                "w_minimum": 2.0,
                "apply_y_max": True,
                "y_maximum": 0.95,
            },
            "generator_revision": "test-generator-revision",
            "generator_revision_source": "test",
            "training": {
                "total_proposals": proposals,
                "inside_analysis_partition": metrics(2, 5.0, 13.0, 3.0),
            },
            "strata": {
                stratum.identifier: {
                    "stratum_id": stratum.identifier,
                    "flat_index": stratum.flat_index,
                    "indices": active._indices(stratum),
                    "bounds": active._bounds(stratum),
                    "status": (
                        "learned" if index < 2 else "no_training_contribution"
                    ),
                    "training_total": training_values[index],
                }
                for index, stratum in enumerate(catalog)
            },
        }
        manifest.write_text(
            json.dumps(manifest_payload, indent=2) + "\n", encoding="utf-8"
        )

        validation = self.root / "migration_validation.json"
        validation_payload = {
            "schema": active.radiative_migrations.VALIDATION_SCHEMA,
            "manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
            "passed": False,
            "global_training_holdout_difference_z_score": 0.25,
            "validation": {
                "total_proposals": proposals,
                "inside_analysis_partition": metrics(2, 5.0, 17.0, 4.0),
            },
            "strata": {
                catalog[0].identifier: {
                    "holdout_total": metrics(1, 1.0, 1.0, 1.0),
                    "coverage_passed": True,
                    "training_holdout_difference_z_score": 0.5,
                },
                catalog[1].identifier: {
                    "holdout_total": metrics(0, 0.0, 0.0, 0.0),
                    "coverage_passed": None,
                    "training_holdout_difference_z_score": None,
                },
                catalog[2].identifier: {
                    "holdout_total": metrics(1, 4.0, 16.0, 4.0),
                    "coverage_passed": False,
                    "training_holdout_difference_z_score": None,
                },
            },
        }
        validation.write_text(
            json.dumps(validation_payload, indent=2) + "\n",
            encoding="utf-8",
        )
        return manifest, validation

    def _write_stratified_queue(self) -> Path:
        config = self.root / "stratified_config.json"
        config.write_text(
            json.dumps(
                {
                    "beam_energy": 6.535,
                    "phase_space": {"W_min": 2.0},
                    "binning": {
                        "Q2": [1.0, 2.0, 3.0],
                        "xB": [0.2, 0.3, 0.4],
                        "minus_t": [0.1, 0.2, 0.3],
                        "phi_deg": [0.0, 180.0, 360.0],
                    },
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        config_payload, config_sha256 = active.radiative_guards.load_analysis_config(
            config
        )
        catalog = active.radiative_guards.enumerate_strata(config_payload)
        total_weight = sum(range(1, len(catalog) + 1))
        records = []
        for weight, stratum in enumerate(catalog, start=1):
            basis = (
                "data_occupancy"
                if stratum.ixb == 0
                else "cumulative_model_tail"
            )
            records.append(
                {
                    "stratum_id": stratum.identifier,
                    "flat_index": stratum.flat_index,
                    "indices": active._indices(stratum),
                    "bounds": active._bounds(stratum),
                    "data_events": (
                        100 - stratum.flat_index
                        if basis == "data_occupancy"
                        else 0
                    ),
                    "model_cross_section_fraction": weight / total_weight,
                    "survey_model_evidence": {
                        "support_status": "independent_support",
                        "parent_coverage_status": "passed",
                        "pooled_cross_section_microbarn": weight * 1.0e-9,
                        "pooled_sem_microbarn": 1.0e-10,
                        "pooled_ess": 25.0 + weight,
                    },
                    "selected": True,
                    "selection_basis": basis,
                    "calibration_priority_rank": stratum.flat_index + 1,
                    "zero_data_model_rank": (
                        stratum.flat_index + 1
                        if basis == "cumulative_model_tail"
                        else None
                    ),
                    "zero_data_cumulative_selected_model_fraction": None,
                    "work_category": "supported_calibration",
                }
            )
        path = self.root / "stratified_queue.json"
        path.write_text(
            json.dumps(
                {
                    "schema": active.CUMULATIVE_QUEUE_SCHEMA,
                    "analysis_config": str(config.resolve()),
                    "analysis_config_sha256": config_sha256,
                    "summary": {
                        "catalog_strata": len(catalog),
                        "selected_strata": len(catalog),
                    },
                    "strata": records,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return path

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

    def test_survey_augmentation_pools_and_expands_active_set(self) -> None:
        config, data, mask = self._write_data_census_inputs()
        selected = np.asarray(np.load(mask), dtype=bool)
        selected[1:4] = False
        np.save(mask, selected)
        data_output = self.root / "data_census_for_survey"
        base_relevance = active.build_data_evidence(
            argparse.Namespace(
                config=config,
                data_events=data,
                selection_mask=mask,
                output=data_output,
                allow_duplicate_event_keys=False,
            )
        )
        manifest, validation = self._write_survey_evidence_inputs(config)
        output = self.root / "survey_evidence"
        augmented_path = active.augment_survey_evidence(
            argparse.Namespace(
                config=config,
                base_relevance=base_relevance,
                migration_manifest=manifest,
                migration_validation=validation,
                output=output,
            )
        )
        survey = json.loads(
            (output / "survey_model_evidence.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(survey["schema"], active.SURVEY_MODEL_SCHEMA)
        self.assertEqual(
            survey["support_counts"],
            {
                "independent_support": 1,
                "training_only": 1,
                "holdout_only": 1,
                "no_survey_contribution": 1,
            },
        )
        self.assertEqual(
            survey["parent_coverage_counts"],
            {
                "passed": 1,
                "failed": 1,
                "no_holdout_contribution": 1,
                "not_assessed": 1,
            },
        )
        self.assertAlmostEqual(
            survey["pooling"]["pooled_inside_analysis"][
                "cross_section_microbarn"
            ],
            0.05,
        )
        self.assertEqual(
            [
                record["model_cross_section_fraction"]
                for record in survey["strata"]
            ],
            [0.4, 0.2, 0.4, 0.0],
        )
        self.assertEqual(
            (output / "survey_nonzero_flat_indices.txt").read_text(
                encoding="utf-8"
            ),
            "0\n1\n2\n",
        )
        self.assertEqual(
            (output / "independent_survey_support_flat_indices.txt").read_text(
                encoding="utf-8"
            ),
            "0\n",
        )
        self.assertEqual(
            (output / "parent_coverage_failed_flat_indices.txt").read_text(
                encoding="utf-8"
            ),
            "2\n",
        )

        augmented = json.loads(augmented_path.read_text(encoding="utf-8"))
        self.assertEqual(len(augmented["sources"]), 2)
        self.assertEqual(
            [
                record["model_cross_section_fraction"]
                for record in augmented["strata"]
            ],
            [0.4, 0.2, 0.4, 0.0],
        )
        classification = active.classify(
            argparse.Namespace(
                config=config,
                relevance=augmented_path,
                calibrations=None,
                pilot_validations=None,
                output=self.root / "survey_classification",
                flat_indices=None,
                minimum_data_events=1,
                maximum_model_fraction=0.3,
                maximum_feed_in_fraction=0.001,
                maximum_closure_fraction=0.001,
                reference_cross_section_microbarn=None,
            )
        )
        classified = json.loads(classification.read_text(encoding="utf-8"))
        self.assertEqual(
            [record["active"] for record in classified["strata"]],
            [True, False, True, False],
        )
        self.assertEqual(
            (self.root / "survey_classification/active_flat_indices.txt").read_text(
                encoding="utf-8"
            ),
            "0\n2\n",
        )

    def test_survey_augmentation_rejects_wrong_validation_manifest(self) -> None:
        config, data, mask = self._write_data_census_inputs()
        base_relevance = active.build_data_evidence(
            argparse.Namespace(
                config=config,
                data_events=data,
                selection_mask=mask,
                output=self.root / "hash_data_census",
                allow_duplicate_event_keys=False,
            )
        )
        manifest, validation = self._write_survey_evidence_inputs(config)
        payload = json.loads(validation.read_text(encoding="utf-8"))
        payload["manifest_sha256"] = "0" * 64
        validation.write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(
            active.ActiveStratumError,
            "does not reference the supplied manifest",
        ):
            active.augment_survey_evidence(
                argparse.Namespace(
                    config=config,
                    base_relevance=base_relevance,
                    migration_manifest=manifest,
                    migration_validation=validation,
                    output=self.root / "wrong_hash_survey",
                )
            )

    def test_cumulative_queue_selects_data_and_global_model_tail(self) -> None:
        config, data, mask = self._write_data_census_inputs()
        selected = np.asarray(np.load(mask), dtype=bool)
        selected[1:4] = False
        np.save(mask, selected)
        base_relevance = active.build_data_evidence(
            argparse.Namespace(
                config=config,
                data_events=data,
                selection_mask=mask,
                output=self.root / "queue_data_census",
                allow_duplicate_event_keys=False,
            )
        )
        manifest, validation = self._write_survey_evidence_inputs(config)
        augmented = active.augment_survey_evidence(
            argparse.Namespace(
                config=config,
                base_relevance=base_relevance,
                migration_manifest=manifest,
                migration_validation=validation,
                output=self.root / "queue_survey_evidence",
            )
        )
        output = self.root / "cumulative_queue"
        result = active.build_cumulative_queue(
            argparse.Namespace(
                config=config,
                relevance=augmented,
                output=output,
                minimum_data_events=1,
                maximum_global_model_residual_fraction=0.25,
            )
        )
        payload = json.loads(result.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema"], active.CUMULATIVE_QUEUE_SCHEMA)
        self.assertEqual(
            payload["summary"],
            {
                "catalog_strata": 4,
                "data_occupied_strata": 1,
                "model_required_zero_data_strata": 1,
                "selected_strata": 2,
                "omitted_strata": 2,
                "data_occupied_model_fraction": 0.4,
                "initial_zero_data_model_fraction": 0.6000000000000001,
                "selected_zero_data_model_fraction": 0.4,
                "selected_total_model_fraction": 0.8,
                "actual_global_model_residual_fraction": 0.2,
                "work_category_counts": {
                    "supported_calibration": 1,
                    "guard_refinement": 1,
                    "targeted_discovery": 0,
                    "omitted": 2,
                },
            },
        )
        self.assertEqual(
            (output / "selected_flat_indices.txt").read_text(
                encoding="utf-8"
            ),
            "0\n2\n",
        )
        self.assertEqual(
            (output / "data_occupied_flat_indices.txt").read_text(
                encoding="utf-8"
            ),
            "0\n",
        )
        self.assertEqual(
            (output / "model_required_zero_data_flat_indices.txt").read_text(
                encoding="utf-8"
            ),
            "2\n",
        )
        self.assertEqual(
            (output / "supported_calibration_flat_indices.txt").read_text(
                encoding="utf-8"
            ),
            "0\n",
        )
        self.assertEqual(
            (output / "guard_refinement_flat_indices.txt").read_text(
                encoding="utf-8"
            ),
            "2\n",
        )
        records = {
            record["flat_index"]: record for record in payload["strata"]
        }
        self.assertEqual(records[0]["selection_basis"], "data_occupancy")
        self.assertEqual(
            records[2]["selection_basis"], "cumulative_model_tail"
        )
        self.assertEqual(
            records[1]["selection_basis"], "omitted_global_residual"
        )
        self.assertEqual(records[2]["zero_data_model_rank"], 1)
        self.assertLessEqual(
            payload["summary"]["actual_global_model_residual_fraction"],
            0.25,
        )

    def test_stratified_batch_selects_two_bases_across_q2(self) -> None:
        queue = self._write_stratified_queue()
        output = self.root / "stratified_batch"
        result = active.select_stratified_batch(
            argparse.Namespace(
                queue=queue,
                output=output,
                work_category="supported_calibration",
                selection_bases=None,
                representatives_per_q2_basis=2,
            )
        )
        payload = json.loads(result.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema"], active.STRATIFIED_BATCH_SCHEMA)
        self.assertEqual(payload["summary"]["selected_strata"], 8)
        self.assertEqual(
            payload["summary"]["basis_counts"],
            {"data_occupancy": 4, "cumulative_model_tail": 4},
        )
        self.assertEqual(
            payload["summary"]["q2_index_counts"], {"0": 4, "1": 4}
        )
        self.assertEqual(len(payload["groups"]), 4)
        self.assertTrue(
            all(group["eligible_strata"] == 4 for group in payload["groups"])
        )
        for group in payload["groups"]:
            self.assertEqual(len(group["selected_flat_indices"]), 2)
        roles = [record["selection_role"] for record in payload["strata"]]
        self.assertEqual(roles.count("largest_model_contribution_anchor"), 4)
        self.assertEqual(roles.count("normalized_index_maximin"), 4)
        self.assertTrue(
            all(
                record["minimum_normalized_distance_squared"] > 0.0
                for record in payload["strata"]
                if record["selection_role"] == "normalized_index_maximin"
            )
        )
        selected_lines = (
            output / "selected_flat_indices.txt"
        ).read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(selected_lines), 8)
        self.assertEqual(len(set(selected_lines)), 8)

        repeat_output = self.root / "stratified_batch_repeat"
        active.select_stratified_batch(
            argparse.Namespace(
                queue=queue,
                output=repeat_output,
                work_category="supported_calibration",
                selection_bases=None,
                representatives_per_q2_basis=2,
            )
        )
        self.assertEqual(
            (output / "selected_flat_indices.txt").read_text(encoding="utf-8"),
            (repeat_output / "selected_flat_indices.txt").read_text(
                encoding="utf-8"
            ),
        )

    def test_stratified_batch_requires_enough_strata_in_every_group(self) -> None:
        queue = self._write_stratified_queue()
        with self.assertRaisesRegex(
            active.ActiveStratumError, "only 4 eligible"
        ):
            active.select_stratified_batch(
                argparse.Namespace(
                    queue=queue,
                    output=self.root / "oversized_stratified_batch",
                    work_category="supported_calibration",
                    selection_bases=None,
                    representatives_per_q2_basis=5,
                )
            )

    def test_stratified_batch_audits_same_q2_basis_fallback(self) -> None:
        queue = self._write_stratified_queue()
        payload = json.loads(queue.read_text(encoding="utf-8"))
        for record in payload["strata"]:
            indices = active._normalized_indices(record["indices"])
            if indices[0] == 1 and record["selection_basis"] == "data_occupancy":
                record["work_category"] = "guard_refinement"
        queue.write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )

        with self.assertRaisesRegex(
            active.ActiveStratumError,
            "data_occupancy, Q2 index 1: only 0 eligible",
        ):
            active.select_stratified_batch(
                argparse.Namespace(
                    queue=queue,
                    output=self.root / "strict_sparse_basis_batch",
                    work_category="supported_calibration",
                    selection_bases=None,
                    representatives_per_q2_basis=2,
                    allow_basis_fallback=False,
                )
            )

        output = self.root / "fallback_sparse_basis_batch"
        result = active.select_stratified_batch(
            argparse.Namespace(
                queue=queue,
                output=output,
                work_category="supported_calibration",
                selection_bases=None,
                representatives_per_q2_basis=2,
                allow_basis_fallback=True,
            )
        )
        batch = json.loads(result.read_text(encoding="utf-8"))
        self.assertEqual(batch["summary"]["selected_strata"], 8)
        self.assertEqual(batch["summary"]["basis_fallback_strata"], 2)
        self.assertEqual(
            batch["summary"]["basis_counts"],
            {"data_occupancy": 2, "cumulative_model_tail": 6},
        )
        self.assertEqual(
            batch["summary"]["q2_index_counts"], {"0": 4, "1": 4}
        )
        sparse_data = next(
            group
            for group in batch["groups"]
            if group["selection_basis"] == "data_occupancy"
            and group["q2_index"] == 1
        )
        expanded_model = next(
            group
            for group in batch["groups"]
            if group["selection_basis"] == "cumulative_model_tail"
            and group["q2_index"] == 1
        )
        self.assertEqual(sparse_data["unfilled_primary_quota"], 2)
        self.assertEqual(sparse_data["selected_flat_indices"], [])
        self.assertEqual(len(expanded_model["primary_selected_flat_indices"]), 2)
        self.assertEqual(
            len(expanded_model["basis_fallback_selected_flat_indices"]), 2
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
