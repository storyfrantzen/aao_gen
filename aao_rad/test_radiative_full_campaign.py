#!/usr/bin/env python3
"""Tests for the scheduler-neutral full mode-4 campaign driver."""

from __future__ import annotations

import argparse
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
