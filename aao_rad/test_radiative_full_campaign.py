#!/usr/bin/env python3
"""Tests for the scheduler-neutral full mode-4 campaign driver."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import radiative_full_campaign as full
import radiative_mode4


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FullCampaignTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.config = self.root / "analysis.json"
        self.recipes = self.root / "recipes.json"
        self.legacy = self.root / "legacy.inp"
        self.config.write_text("{}\n", encoding="utf-8")
        self.recipes.write_text("{}\n", encoding="utf-8")
        self.legacy.write_text("legacy\n", encoding="utf-8")
        self.queue = self.root / "queue.json"
        records = []
        categories = (
            "supported_calibration",
            "guard_refinement",
            "targeted_discovery",
            "supported_calibration",
            "supported_calibration",
        )
        for flat_index, category in enumerate(categories):
            records.append(
                {
                    "flat_index": flat_index,
                    "stratum_id": f"s{flat_index:05d}",
                    "selected": flat_index < 4,
                    "work_category": category,
                    "data_events": 1 if flat_index in (0, 1) else 0,
                }
            )
        _write_json(
            self.queue,
            {
                "schema": full.QUEUE_SCHEMA,
                "analysis_config_sha256": _sha256(self.config),
                "summary": {"selected_strata": 4},
                "strata": records,
            },
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def _fake_prepare(args: argparse.Namespace) -> Path:
        # Match the eager numeric conversions performed by the real
        # radiative_mode4.prepare().  This guards against constructing an
        # internal Namespace with None where argparse supplies zero/defaults.
        int(args.events_per_stratum)
        int(args.trials)
        float(args.calibration_inside_guard_fraction)
        args.output.mkdir(parents=True)
        (args.output / "analysis_config.json").write_bytes(args.config.read_bytes())
        (args.output / "continuous_guard_recipes.json").write_bytes(
            args.recipes.read_bytes()
        )
        (args.output / "legacy_input.inp").write_bytes(args.input.read_bytes())
        if args.refinements is not None:
            (args.output / "guard_refinements.json").write_bytes(
                args.refinements.read_bytes()
            )
        indices = [
            int(line)
            for line in args.flat_index_file.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        operation = (
            "calibration"
            if args.command == "prepare-calibration"
            else "generation"
        )
        runs = []
        for flat_index in indices:
            for replica in range(args.replicas):
                identifier = f"s{flat_index:05d}"
                runs.append(
                    {
                        "flat_index": flat_index,
                        "stratum_id": identifier,
                        "replica_index": replica,
                        "seed": args.seed_base + 1000 * flat_index + replica,
                        "trials_requested": args.trials or 0,
                        "events_requested": args.events_per_stratum or 0,
                        "output_stem": (
                            f"runs/{identifier}/{identifier}__g{replica:04d}"
                        ),
                    }
                )
        manifest = {
            "schema": radiative_mode4.MANIFEST_SCHEMA,
            "operation": operation,
            "generator_revision": args.generator_revision,
            "guard_candidate": args.candidate,
            "core_fraction": args.core_fraction,
            "analysis_selection": {
                "apply_y_max": bool(args.apply_y_max),
            },
            "runs": runs,
        }
        path = args.output / "manifest.json"
        _write_json(path, manifest)
        return path

    def _plan(self) -> Path:
        output = self.root / "initial"
        args = argparse.Namespace(
            queue=self.queue,
            config=self.config,
            recipes=self.recipes,
            refinements=None,
            input=self.legacy,
            output=output,
            selection="selected",
            work_category=None,
            candidate="padding_0p035",
            core_fraction=0.9,
            inside_guard_trial_fraction=0.5,
            trials=7_000_000,
            replicas=1,
            seed_base=907_001,
            heartbeat_interval=100_000,
            apply_y_max=True,
            generator_revision="physics-revision",
        )
        with mock.patch.object(
            radiative_mode4, "prepare", side_effect=self._fake_prepare
        ):
            return full.plan_calibration(args)

    def test_plan_status_and_idempotent_task_execution(self) -> None:
        campaign_path = self._plan()
        campaign = json.loads(campaign_path.read_text(encoding="utf-8"))
        self.assertEqual(campaign["task_count"], 4)
        self.assertEqual(
            campaign["selection"]["category_counts"],
            {
                "guard_refinement": 1,
                "supported_calibration": 2,
                "targeted_discovery": 1,
            },
        )
        executable = self.root / "aao_rad"
        executable.write_text("#!/bin/sh\n", encoding="utf-8")
        executable.chmod(0o755)

        def fake_run(args: argparse.Namespace) -> Path:
            manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
            record = next(
                item
                for item in manifest["runs"]
                if item["flat_index"] == args.flat_index
                and item["replica_index"] == args.replica_index
            )
            path = args.manifest.parent / (record["output_stem"] + ".json")
            _write_json(
                path,
                {
                    "schema": radiative_mode4.RUN_SCHEMA,
                    "flat_index": args.flat_index,
                    "replica_index": args.replica_index,
                    "source_manifest_sha256": _sha256(args.manifest),
                },
            )
            return path

        run_args = argparse.Namespace(
            campaign=campaign_path,
            task_id=1,
            executable=executable,
            overwrite=False,
        )
        with mock.patch.object(
            radiative_mode4, "run", side_effect=fake_run
        ) as runner:
            first = full.run_task(run_args)
            second = full.run_task(run_args)
        self.assertEqual(first, second)
        self.assertEqual(runner.call_count, 1)

        status_path = full.status(
            argparse.Namespace(campaign=campaign_path, output=None)
        )
        status = json.loads(status_path.read_text(encoding="utf-8"))
        self.assertEqual(
            status["status_counts"],
            {"complete": 1, "failed": 0, "pending": 3},
        )

    def _calibration_stratum(
        self,
        flat_index: int,
        status: str,
        readiness: str,
        inside_targets: int,
        complement_targets: int,
        inside_trials: int = 3_500_000,
        complement_trials: int = 3_500_000,
    ) -> dict:
        recommendation = (
            {"sigr_max": 1.0e-6}
            if readiness in full.READY_STATES
            else None
        )
        return {
            "flat_index": flat_index,
            "stratum_id": f"s{flat_index:05d}",
            "recommendation_status": status,
            "pilot_readiness": readiness,
            "recommended_envelope": recommendation,
            "inside_guard": {
                "target_candidates": inside_targets,
                "trials": inside_trials,
            },
            "guard_complement": {
                "target_candidates": complement_targets,
                "trials": complement_trials,
            },
        }

    def test_followup_is_grouped_and_automatically_sized(self) -> None:
        campaign_path = self._plan()
        calibration_path = self.root / "calibration.json"
        _write_json(
            calibration_path,
            {
                "schema": radiative_mode4.CALIBRATION_SCHEMA,
                "strata": [
                    self._calibration_stratum(
                        0, "provisional_zero_complement",
                        "ready_provisional_zero_complement", 1200, 0
                    ),
                    self._calibration_stratum(
                        1, "insufficient_provisional_inside_envelope_support",
                        "not_ready", 500, 0
                    ),
                    self._calibration_stratum(
                        2, "insufficient_guard_complement_target_support",
                        "not_ready", 1200, 2
                    ),
                    self._calibration_stratum(
                        3, "insufficient_zero_complement_exposure",
                        "not_ready", 1200, 0,
                        complement_trials=500_000
                    ),
                    self._calibration_stratum(
                        4, "no_candidate_meets_duplicate_limit",
                        "not_ready", 1200, 20
                    ),
                ],
            },
        )
        output = self.root / "followup"
        args = argparse.Namespace(
            campaign=campaign_path,
            calibration=calibration_path,
            output=output,
            minimum_component_targets=20,
            minimum_provisional_inside_targets=1000,
            zero_complement_confidence=0.95,
            maximum_zero_complement_target_rate=1.0e-6,
            followup_target_safety_factor=1.5,
            trial_quantum=1_000_000,
            maximum_followup_trials=100_000_000,
            discovery_trials=20_000_000,
            heartbeat_interval=100_000,
            seed_base=107_000_001,
        )
        with mock.patch.object(
            radiative_mode4, "prepare", side_effect=self._fake_prepare
        ):
            result = full.plan_followup(args)
        payload = json.loads(result.read_text(encoding="utf-8"))
        self.assertEqual(payload["task_count"], 3)
        self.assertEqual(payload["followup"]["groups"], 3)
        self.assertEqual(payload["followup"]["ready_strata"], 1)
        self.assertEqual(payload["followup"]["manual_review_strata"], 1)
        self.assertEqual(len(payload["pool_manifests"]), 4)
        text = (output / "followup_plan.tsv").read_text(encoding="utf-8")
        self.assertIn("inside_support", text)
        self.assertIn("complement_support", text)
        self.assertIn("zero_complement_exposure", text)
        self.assertIn("manual_envelope_review", text)

    def test_followup_can_exclude_capped_estimates(self) -> None:
        campaign_path = self._plan()
        calibration_path = self.root / "capped_calibration.json"
        _write_json(
            calibration_path,
            {
                "schema": radiative_mode4.CALIBRATION_SCHEMA,
                "strata": [
                    self._calibration_stratum(
                        0,
                        "provisional_zero_complement",
                        "ready_provisional_zero_complement",
                        1200,
                        0,
                    ),
                    self._calibration_stratum(
                        1,
                        "insufficient_guard_complement_target_support",
                        "not_ready",
                        1200,
                        19,
                        complement_trials=3_500_000,
                    ),
                    self._calibration_stratum(
                        2,
                        "insufficient_guard_complement_target_support",
                        "not_ready",
                        1200,
                        1,
                        complement_trials=100_000_000,
                    ),
                ],
            },
        )
        output = self.root / "uncapped_followup"
        args = argparse.Namespace(
            campaign=campaign_path,
            calibration=calibration_path,
            output=output,
            minimum_component_targets=20,
            minimum_provisional_inside_targets=1000,
            zero_complement_confidence=0.95,
            maximum_zero_complement_target_rate=1.0e-6,
            followup_target_safety_factor=1.5,
            trial_quantum=1_000_000,
            maximum_followup_trials=100_000_000,
            discovery_trials=20_000_000,
            heartbeat_interval=100_000,
            seed_base=207_000_001,
            capped_policy="exclude",
        )
        with mock.patch.object(
            radiative_mode4, "prepare", side_effect=self._fake_prepare
        ):
            result = full.plan_followup(args)
        payload = json.loads(result.read_text(encoding="utf-8"))
        self.assertEqual(payload["task_count"], 1)
        self.assertEqual(payload["followup"]["scheduled_strata"], 1)
        self.assertEqual(
            payload["followup"]["deferred_by_capped_policy"], 1
        )
        with (output / "followup_plan.tsv").open() as source:
            rows = list(csv.DictReader(source, delimiter="\t"))
        capped = next(row for row in rows if row["stratum_id"] == "s00002")
        self.assertEqual(capped["capped_at_maximum"], "True")
        self.assertEqual(capped["selected_by_capped_policy"], "False")

    @staticmethod
    def _guard() -> dict:
        box = radiative_mode4.GuardBox(
            nonperiodic={
                "r_u": (0.2, 0.6),
                "r_ep": (0.2, 0.6),
                "u_gamma": (0.2, 1.0),
                "hadron_cosine_base": (0.2, 0.6),
            },
            phi_origin=0.5,
            phi_relative=(-0.1, 0.1),
        )
        return box.manifest_record()

    def _complete_calibration_stratum(
        self, *, guard: dict, refined: bool
    ) -> dict:
        complement_targets = 0 if refined else 1
        recommendation = (
            {
                "sigr_max": 2.0e-6,
                "expected_events_per_proposal": 1.0e-4,
                "expected_duplicate_event_fraction": 0.0,
            }
            if refined
            else None
        )
        return {
            "flat_index": 0,
            "stratum_id": "s00000",
            "indices": {"iq2": 0, "ixb": 0, "it": 0, "iphi": 0},
            "bounds": {
                "Q2": [1.0, 1.5],
                "xB": [0.15, 0.2],
                "minus_t": [0.09, 0.15],
                "phi_deg": [0.0, 18.0],
            },
            "guard": guard,
            "inside_guard": {
                "target_candidates": 2000,
                "trials": 3_500_000,
            },
            "guard_complement": {
                "target_candidates": complement_targets,
                "trials": 100_000_000 if not refined else 3_500_000,
            },
            "integrated_cross_section_microbarn": 1.0e-7,
            "integrated_cross_section_sem_microbarn": 1.0e-9,
            "estimated_guard_complement_cross_section_fraction": (
                0.1 if not refined else 0.0
            ),
            "provisional_inside_envelope_support": {
                "minimum_required_targets": 1000,
                "empirical_next_target_rank_resolution": 1.0 / 2001.0,
            },
            "zero_complement_stopping_test": {
                "one_sided_upper_target_rate": (
                    None if not refined else 8.0e-7
                ),
                "maximum_allowed_target_rate": 1.0e-6,
            },
            "recommendation_status": (
                "provisional_zero_complement"
                if refined
                else "insufficient_guard_complement_target_support"
            ),
            "recommendation_basis": (
                "safety_scaled_inside_guard_observed_maximum_with_"
                "zero_complement_rate_bound"
                if refined
                else None
            ),
            "pilot_readiness": (
                "ready_provisional_zero_complement"
                if refined
                else "not_ready"
            ),
            "recommended_envelope": recommendation,
        }

    def test_capped_guard_batch_refinement_and_composite_merge(self) -> None:
        guard = self._guard()
        source_root = self.root / "calibration_source"
        csv_path = source_root / "runs/s00000/s00000__g0000.calibration.csv"
        csv_path.parent.mkdir(parents=True)
        with csv_path.open("w", encoding="utf-8", newline="") as destination:
            destination.write(
                f"# schema={radiative_mode4.MODE4_CALIBRATION_EVENT_SCHEMA}\n"
            )
            writer = csv.DictWriter(
                destination,
                fieldnames=radiative_mode4.MODE4_CALIBRATION_COLUMNS,
            )
            writer.writeheader()
            writer.writerow(
                {
                    "trial": 1,
                    "proposal_component": 0,
                    "inside_core": 0,
                    "proposal_density_ratio": 1.0,
                    "integrand_corrected": 1.0e-6,
                    "component_importance_weight": 1.0,
                    "r_u": 0.65,
                    "r_ep": 0.4,
                    "u_gamma": 0.9,
                    "hadron_cosine_base": 0.4,
                    "hadron_phi_base": 0.5,
                    "q2_observed": 1.2,
                    "xb_observed": 0.17,
                    "minus_t_observed": 0.1,
                    "phi_observed_deg": 9.0,
                    "w_observed": 2.1,
                    "y_observed": 0.5,
                }
            )
        source_manifest = source_root / "manifest.json"
        source_record = {
            "flat_index": 0,
            "stratum_id": "s00000",
            "replica_index": 0,
            "seed": 1,
            "trials_requested": 100_000_000,
            "events_requested": 0,
            "output_stem": "runs/s00000/s00000__g0000",
            "guard": guard,
        }
        _write_json(
            source_manifest,
            {
                "schema": radiative_mode4.MANIFEST_SCHEMA,
                "operation": "calibration",
                "generator_revision": "physics-revision",
                "guard_candidate": "padding_0p035",
                "core_fraction": 0.9,
                "analysis_selection": {"apply_y_max": True},
                "runs": [source_record],
            },
        )
        base_stratum = self._complete_calibration_stratum(
            guard=guard, refined=False
        )
        common = {
            "schema": radiative_mode4.CALIBRATION_SCHEMA,
            "analysis_config_sha256": _sha256(self.config),
            "guard_recipes_sha256": _sha256(self.recipes),
            "guard_refinements_sha256": None,
            "guard_candidate": "padding_0p035",
            "core_fraction": 0.9,
            "analysis_selection": {"apply_y_max": True},
            "calibration_proposal": "guard_partition",
            "envelope_safety_factor": 1.2,
            "maximum_duplicate_fraction": 0.05,
            "minimum_component_targets": 20,
            "minimum_provisional_inside_targets": 1000,
            "zero_complement_policy": {
                "enabled": True,
                "confidence_level": 0.95,
                "maximum_target_rate": 1.0e-6,
            },
            "generator_revision": "physics-revision",
            "generator_revisions": ["physics-revision"],
            "stratum_count": 1,
        }
        calibration_path = self.root / "base_calibration.json"
        _write_json(
            calibration_path,
            {
                **common,
                "source_manifests": [
                    {
                        "path": str(source_manifest),
                        "sha256": _sha256(source_manifest),
                    }
                ],
                "strata": [base_stratum],
            },
        )
        selection = self.root / "base_campaign/selection.txt"
        selection.parent.mkdir(parents=True)
        selection.write_text("0\n", encoding="utf-8")
        parent_root = self.root / "base_campaign"
        parent_payload = full._campaign_payload(
            kind="calibration_followup",
            root=parent_root,
            tasks=[],
            manifests=[source_manifest],
            pool_manifests=[source_manifest],
            frozen_inputs={
                "config": str(self.config),
                "config_sha256": _sha256(self.config),
                "recipes": str(self.recipes),
                "recipes_sha256": _sha256(self.recipes),
                "legacy_input": str(self.legacy),
                "legacy_input_sha256": _sha256(self.legacy),
                "refinements": None,
                "refinements_sha256": None,
            },
            selection={
                "flat_indices": str(selection),
                "flat_indices_sha256": _sha256(selection),
                "selected_strata": 1,
            },
        )
        parent_path = full._finish_campaign(parent_root, parent_payload, [])
        output = self.root / "refinement_campaign"
        args = argparse.Namespace(
            campaign=parent_path,
            calibration=calibration_path,
            output=output,
            minimum_component_targets=20,
            minimum_provisional_inside_targets=1000,
            zero_complement_confidence=0.95,
            maximum_zero_complement_target_rate=1.0e-6,
            followup_target_safety_factor=1.5,
            trial_quantum=1_000_000,
            maximum_followup_trials=100_000_000,
            discovery_trials=20_000_000,
            minimum_face_margin=0.005,
            excursion_margin_fraction=0.25,
            maximum_volume_ratio=4.0,
            trials=7_000_000,
            inside_guard_trial_fraction=0.5,
            replicas=1,
            seed_base=507_000_001,
            heartbeat_interval=100_000,
        )
        original_box = radiative_mode4._guard_box_from_manifest(
            {"guard": guard}
        )
        recipes = {"strata": {"s00000": {}}}
        with mock.patch.object(
                radiative_mode4,
                "_load_config_and_recipes",
                return_value=(
                    {},
                    recipes,
                    _sha256(self.config),
                    _sha256(self.recipes),
                ),
            ), mock.patch.object(
                radiative_mode4, "_candidate_padding", return_value=0.035
            ), mock.patch.object(
                radiative_mode4,
                "reconstruct_guard_box",
                return_value=original_box,
            ), mock.patch.object(
                radiative_mode4, "prepare", side_effect=self._fake_prepare
            ):
            campaign_path = full.plan_refinement(args)
        campaign = json.loads(campaign_path.read_text(encoding="utf-8"))
        self.assertEqual(campaign["kind"], "calibration_refinement")
        self.assertEqual(campaign["task_count"], 1)
        artifact_path = output / "guard_refinements.json"
        artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
        refined_guard = artifact["preview"]["s00000"]["refined_guard"]
        self.assertGreater(refined_guard["axes"]["r_u"][1], 0.65)
        self.assertEqual(
            artifact["batch_refinement"]["selected_capped_strata"], 1
        )

        refined_report_path = output / "refined_report.json"
        refined_stratum = self._complete_calibration_stratum(
            guard=refined_guard, refined=True
        )
        _write_json(
            refined_report_path,
            {
                **common,
                "guard_refinements_sha256": _sha256(artifact_path),
                "source_manifests": [],
                "strata": [refined_stratum],
            },
        )
        composite = output / "composite.json"
        merged_path = full._merge_refinement_calibrations(
            campaign_path=campaign_path,
            campaign=campaign,
            refined_path=refined_report_path,
            output=composite,
        )
        merged = json.loads(merged_path.read_text(encoding="utf-8"))
        self.assertEqual(
            merged["strata"][0]["pilot_readiness"],
            "ready_provisional_zero_complement",
        )
        self.assertTrue(
            merged["composite_calibration"][
                "pre_refinement_trials_excluded_for_changed_guards"
            ]
        )

    def test_followup_finalize_recomputes_only_touched_strata(self) -> None:
        guard = self._guard()
        first = self._complete_calibration_stratum(
            guard=guard, refined=True
        )
        second = json.loads(json.dumps(first))
        second["flat_index"] = 1
        second["stratum_id"] = "s00001"
        second["indices"]["iphi"] = 1
        common = {
            "schema": radiative_mode4.CALIBRATION_SCHEMA,
            "analysis_config_sha256": _sha256(self.config),
            "guard_recipes_sha256": _sha256(self.recipes),
            "guard_refinements_sha256": None,
            "guard_candidate": "padding_0p035",
            "core_fraction": 0.9,
            "analysis_selection": {"apply_y_max": True},
            "calibration_proposal": "guard_partition",
            "envelope_safety_factor": 1.2,
            "maximum_duplicate_fraction": 0.05,
            "minimum_component_targets": 20,
            "minimum_provisional_inside_targets": 1000,
            "zero_complement_policy": {
                "enabled": True,
                "confidence_level": 0.95,
                "maximum_target_rate": 1.0e-6,
            },
            "generator_revision": "physics-revision",
            "generator_revisions": ["physics-revision"],
            "revision_compatibility_override": {"enabled": False},
        }
        parent_manifest = self.root / "parent_stage/manifest.json"
        followup_manifest = self.root / "followup_stage/manifest.json"
        parent_run = {
            "flat_index": 0,
            "stratum_id": "s00000",
            "replica_index": 0,
            "seed": 1,
            "output_stem": "runs/s00000/s00000__g0000",
        }
        followup_run = {
            "flat_index": 1,
            "stratum_id": "s00001",
            "replica_index": 0,
            "seed": 2,
            "output_stem": "runs/s00001/s00001__g0000",
        }
        _write_json(
            parent_manifest,
            {"operation": "calibration", "runs": [parent_run]},
        )
        _write_json(
            followup_manifest,
            {"operation": "calibration", "runs": [followup_run]},
        )
        sources = [
            {"path": str(path), "sha256": _sha256(path)}
            for path in (parent_manifest, followup_manifest)
        ]
        base_path = self.root / "parent_calibration.json"
        _write_json(
            base_path,
            {
                **common,
                "stratum_count": 2,
                "source_manifests": [sources[0]],
                "strata": [first, second],
            },
        )
        campaign_root = self.root / "incremental_followup"
        campaign_root.mkdir()
        task = {
            "task_id": 1,
            "stage": "followup",
            "operation": "calibration",
            "manifest": str(followup_manifest),
            "manifest_sha256": _sha256(followup_manifest),
            "flat_index": 1,
            "stratum_id": "s00001",
            "replica_index": 0,
            "seed": 2,
            "requested_trials": 1_000_000,
            "requested_events": 0,
            "expected_run": str(self.root / "unused_run.json"),
        }
        payload = full._campaign_payload(
            kind="calibration_followup",
            root=campaign_root,
            tasks=[task],
            manifests=[followup_manifest],
            pool_manifests=[parent_manifest, followup_manifest],
            frozen_inputs={
                "config": str(self.config),
                "config_sha256": _sha256(self.config),
                "recipes": str(self.recipes),
                "recipes_sha256": _sha256(self.recipes),
                "legacy_input": str(self.legacy),
                "legacy_input_sha256": _sha256(self.legacy),
                "refinements": None,
                "refinements_sha256": None,
            },
            selection={"selected_strata": 2},
            calibration_report=base_path,
        )
        campaign_path = full._finish_campaign(
            campaign_root, payload, [task]
        )
        updated = json.loads(json.dumps(second))
        updated["inside_guard"]["target_candidates"] = 2500

        def fake_finalize(args: argparse.Namespace) -> Path:
            self.assertEqual(args.stratum_ids, {"s00001"})
            self.assertEqual(
                args.manifests, [parent_manifest, followup_manifest]
            )
            _write_json(
                args.output,
                {
                    **common,
                    "stratum_count": 1,
                    "source_manifests": sources,
                    "strata": [updated],
                },
            )
            return args.output

        output = campaign_root / "envelope_calibration.json"
        args = argparse.Namespace(
            campaign=campaign_path,
            output=output,
            envelope_safety_factor=1.2,
            maximum_duplicate_fraction=0.05,
            minimum_component_targets=20,
            minimum_provisional_inside_targets=1000,
            allow_zero_complement=True,
            zero_complement_confidence=0.95,
            maximum_zero_complement_target_rate=1.0e-6,
            full_recompute=False,
        )
        with mock.patch.object(full, "_assert_complete"), mock.patch.object(
            radiative_mode4, "finalize", side_effect=fake_finalize
        ):
            result = full.finalize(args)
        merged = json.loads(result.read_text(encoding="utf-8"))
        by_id = {item["stratum_id"]: item for item in merged["strata"]}
        self.assertEqual(
            by_id["s00000"]["inside_guard"]["target_candidates"], 2000
        )
        self.assertEqual(
            by_id["s00001"]["inside_guard"]["target_candidates"], 2500
        )
        audit = merged["incremental_calibration"]
        self.assertEqual(audit["recomputed_strata"], 1)
        self.assertEqual(audit["unchanged_strata_copied_from_parent"], 1)
        self.assertTrue(
            audit["statistically_equivalent_to_full_recomputation"]
        )

    def test_production_plan_and_swif_script_cover_every_replica(self) -> None:
        campaign_path = self._plan()
        calibration_path = self.root / "ready.json"
        strata = [
            self._calibration_stratum(
                index,
                "provisional_zero_complement",
                "ready_provisional_zero_complement",
                1200,
                0,
            )
            for index in range(4)
        ]
        _write_json(
            calibration_path,
            {"schema": radiative_mode4.CALIBRATION_SCHEMA, "strata": strata},
        )
        output = self.root / "production"
        args = argparse.Namespace(
            campaign=campaign_path,
            calibration=calibration_path,
            output=output,
            events_per_stratum=200,
            replicas=2,
            seed_base=307_000_001,
            heartbeat_interval=100_000,
            allow_incomplete=False,
        )
        with mock.patch.object(
            radiative_mode4, "prepare", side_effect=self._fake_prepare
        ):
            result = full.plan_production(args)
        payload = json.loads(result.read_text(encoding="utf-8"))
        self.assertEqual(payload["task_count"], 8)
        self.assertEqual(payload["selection"]["ready_strata"], 4)

        executable = self.root / "aao_rad"
        executable.write_text("#!/bin/sh\n", encoding="utf-8")
        executable.chmod(0o755)
        script = full.emit_swif(
            argparse.Namespace(
                campaign=result,
                workflow="rgk_mode4_test",
                executable=executable,
                output=None,
                cores=1,
                disk="2gb",
                ram="1gb",
                walltime="8hr",
            )
        )
        content = script.read_text(encoding="utf-8")
        self.assertIn("swif2 create -workflow", content)
        self.assertIn("run-task", (output / "run_swif_task.sh").read_text())
        self.assertIn("tail -n +2", content)
        subprocess.run(["bash", "-n", str(script)], check=True)
        subprocess.run(
            ["bash", "-n", str(output / "run_swif_task.sh")], check=True
        )


if __name__ == "__main__":
    unittest.main()
