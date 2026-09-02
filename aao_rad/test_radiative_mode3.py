#!/usr/bin/env python3
"""Tests for the global radiative direct-mixture mode-3 workflow."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import tempfile
import unittest
from pathlib import Path

import radiative_mode3 as mode3


def _fixtures(root: Path) -> tuple[Path, Path]:
    config = {
        "beam_energy": 6.535,
        "target_mass": 0.9382720813,
        "phase_space": {
            "Q2_min": 0.2,
            "W_min": 1.08,
            "electron_p_min": 0.2,
            "y_max": None,
        },
        "binning": {
            "Q2": [0.2, 5.0],
            "xB": [0.02, 0.99],
            "minus_t": [0.0, 10.0],
            "phi_deg": [0.0, 360.0],
        },
    }
    config_path = root / "analysis.json"
    config_path.write_text(json.dumps(config) + "\n", encoding="utf-8")
    legacy = "\n".join(
        [
            "1",
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
            "6.5",
            ".2 5",
            ".2 6.5",
            ".005",
            "2",
            "0",
            "1",
        ]
    )
    input_path = root / "legacy.inp"
    input_path.write_text(legacy + "\n", encoding="utf-8")
    return config_path, input_path


def _prepare_args(
    root: Path,
    config: Path,
    legacy: Path,
    *,
    operation: str,
    direct_fraction: float = 0.75,
) -> argparse.Namespace:
    values = dict(
        config=config,
        input=legacy,
        output=root / operation,
        tag="test_mode3",
        padding_fraction=0.0,
        direct_fraction=direct_fraction,
        minimum_photon_energy=None,
        seed_base=731001,
        heartbeat_interval=100,
        generator_revision="generator-test",
        trials=2_000,
        replicas=1,
        calibration=None,
        total_events=0,
        events_per_job=0,
        lund_files_per_directory=0,
    )
    return argparse.Namespace(**values)


class Mode3WorkflowTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.config, self.legacy = _fixtures(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_rga_padding_and_proposal_density_correction(self) -> None:
        rga = {
            "beam_energy": 10.604,
            "target_mass": 0.9382720813,
            "phase_space": {
                "Q2_min": 1.0,
                "W_min": 2.0,
                "electron_p_min": 2.0,
            },
            "binning": {
                "Q2": [1.0, 10.5],
                "xB": [0.05, 0.7],
                "minus_t": [0.09, 2.0],
                "phi_deg": [0.0, 360.0],
            },
        }
        path = self.root / "rga.json"
        path.write_text(json.dumps(rga), encoding="utf-8")
        settings = mode3._configuration(path, 0.035)
        self.assertAlmostEqual(settings["bounds"]["Q2"][0], 0.6675)
        self.assertAlmostEqual(settings["bounds"]["Q2"][1], 10.8325)
        self.assertAlmostEqual(settings["bounds"]["xB"][0], 0.02725)
        self.assertAlmostEqual(settings["bounds"]["xB"][1], 0.72275)
        self.assertAlmostEqual(settings["bounds"]["minus_t"][0], 0.02315)
        self.assertAlmostEqual(settings["bounds"]["minus_t"][1], 2.06685)
        arguments = dict(
            q2=2.0,
            xb=0.3,
            minus_t=0.5,
            phi_deg=90.0,
            t_jacobian=0.8,
            ep_range=8.0,
            direct_fraction=0.75,
            xb_bounds=(0.1, 0.7),
            t_bounds=(0.09, 2.0),
            phi_bounds=(0.0, 360.0),
            target_mass=0.9382720813,
        )
        direct_to_legacy = (
            8.0 * 2.0 * 0.9382720813 * 0.3**2 / (2.0 * 0.6)
            * 2.0 * 0.8 / 1.91
        )
        expected = 1.0 / (0.25 + 0.75 * direct_to_legacy)
        self.assertAlmostEqual(mode3.proposal_density_ratio(**arguments), expected)
        arguments["minus_t"] = 3.0
        self.assertEqual(mode3.proposal_density_ratio(**arguments), 4.0)
        arguments["direct_fraction"] = 0.0
        self.assertEqual(mode3.proposal_density_ratio(**arguments), 1.0)
        arguments["direct_fraction"] = mode3.MAXIMUM_DIRECT_FRACTION
        self.assertAlmostEqual(mode3.proposal_density_ratio(**arguments), 20.0)
        arguments["minus_t"] = 0.5
        expected_bounded = 1.0 / (0.05 + 0.95 * direct_to_legacy)
        self.assertAlmostEqual(
            mode3.proposal_density_ratio(**arguments), expected_bounded
        )
        self.assertLessEqual(expected_bounded, 20.0)

    def test_prepare_snapshots_mode3_input_and_splits_production(self) -> None:
        args = _prepare_args(
            self.root, self.config, self.legacy, operation="calibration"
        )
        manifest_path = mode3.prepare_calibration(args)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["operation"], "calibration")
        task = json.loads(
            Path(manifest["runs"][0]["task"]).read_text(encoding="utf-8")
        )
        self.assertEqual(task["schema"], mode3.TASK_SCHEMA)
        self.assertEqual(task["run"]["replica_index"], 0)
        generated = Path(manifest["runs"][0]["input"]).read_text(
            encoding="utf-8"
        )
        records = mode3._records(generated)
        self.assertAlmostEqual(float(records[11]), 6.535)
        self.assertEqual(records[18], "3")
        self.assertEqual(records[21], "0.75")
        self.assertEqual(records[-3:], ["1", "2000", "100"])
        self.assertEqual(
            manifest["settings"]["generator_physics"][
                "minimum_photon_energy_gev"
            ],
            0.005,
        )

        calibration = {
            "schema": mode3.CALIBRATION_SCHEMA,
            "settings": manifest["settings"],
            "direct_fraction": 0.75,
            "recommended_sigr_max": 0.01,
        }
        calibration_path = self.root / "envelope.json"
        calibration_path.write_text(json.dumps(calibration), encoding="utf-8")
        production_args = _prepare_args(
            self.root, self.config, self.legacy, operation="production"
        )
        production_args.output = self.root / "production"
        production_args.calibration = calibration_path
        production_args.total_events = 10_001
        production_args.events_per_job = 5_000
        production_args.lund_files_per_directory = 2
        production_path = mode3.prepare_production(production_args)
        production = json.loads(production_path.read_text(encoding="utf-8"))
        self.assertEqual(production["task_count"], 3)
        self.assertEqual(
            [run["requested_events"] for run in production["runs"]],
            [5_000, 5_000, 1],
        )
        self.assertEqual(
            Path(production["runs"][0]["lund"]).name,
            "test_mode3__g00000000.lund",
        )
        self.assertEqual(
            Path(production["runs"][2]["lund"]).parent.name,
            "chunk_0001",
        )

    def test_prepare_rejects_unbounded_direct_only_proposal(self) -> None:
        args = _prepare_args(
            self.root,
            self.config,
            self.legacy,
            operation="calibration",
            direct_fraction=1.0,
        )
        with self.assertRaisesRegex(ValueError, "5% legacy support"):
            mode3.prepare_calibration(args)

    def test_production_validation_accepts_overshoot_and_rejects_truncation(self) -> None:
        valid = {
            "mode3_schema": mode3.MODE3_SCHEMA,
            "requested_events": 5_000,
            "events": 5_002,
            "mcall_max": 4,
            "event_overshoot": 2,
            "emitting_candidates": 4_990,
            "duplicate_events": 12,
        }
        mode3._validate_production_run(valid, self.root / "valid.json")
        unsafe = dict(valid)
        unsafe.update(
            events=5_000,
            mcall_max=40_930_660,
            event_overshoot=0,
            emitting_candidates=4_600,
            duplicate_events=400,
        )
        with self.assertRaisesRegex(mode3.Mode3Error, "exceeds"):
            mode3._validate_production_run(unsafe, self.root / "unsafe.json")

    def test_photon_threshold_is_frozen_across_calibration_and_production(self) -> None:
        args = _prepare_args(
            self.root, self.config, self.legacy, operation="calibration"
        )
        args.minimum_photon_energy = 0.010
        manifest_path = mode3.prepare_calibration(args)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        generated = mode3._records(
            Path(manifest["runs"][0]["input"]).read_text(encoding="utf-8")
        )
        self.assertAlmostEqual(float(generated[14]), 0.010)
        self.assertAlmostEqual(
            manifest["settings"]["generator_physics"][
                "minimum_photon_energy_gev"
            ],
            0.010,
        )

        calibration = {
            "schema": mode3.CALIBRATION_SCHEMA,
            "settings": manifest["settings"],
            "direct_fraction": 0.75,
            "recommended_sigr_max": 0.01,
        }
        calibration_path = self.root / "threshold-envelope.json"
        calibration_path.write_text(json.dumps(calibration), encoding="utf-8")
        production = _prepare_args(
            self.root, self.config, self.legacy, operation="production"
        )
        production.output = self.root / "threshold-production"
        production.calibration = calibration_path
        production.total_events = 10
        production.events_per_job = 5
        production.lund_files_per_directory = 5_000
        with self.assertRaisesRegex(
            mode3.Mode3Error, "production settings differ from calibration"
        ):
            mode3.prepare_production(production)

    def test_executable_calibration_and_production_smoke(self) -> None:
        executable = Path(__file__).resolve().parent / "build" / "aao_rad"
        if not executable.is_file():
            self.skipTest("build/aao_rad is created by the Makefile test target")
        calibration_args = _prepare_args(
            self.root, self.config, self.legacy, operation="calibration"
        )
        manifest_path = mode3.prepare_calibration(calibration_args)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        task_path = Path(manifest["runs"][0]["task"])
        run_path = mode3.run(
            argparse.Namespace(
                manifest=task_path,
                replica_index=0,
                executable=executable,
                scratch_root=None,
                overwrite=False,
            )
        )
        run = json.loads(run_path.read_text(encoding="utf-8"))
        self.assertEqual(run["ntries"], 2_000)
        self.assertEqual(run["events"], 0)
        self.assertGreater(run["target_candidates"], 0)
        self.assertGreater(run["integrand_sum"], 0.0)
        envelope_path = mode3.finalize_calibration(
            argparse.Namespace(
                manifest=manifest_path,
                envelope_safety_factor=2.0,
                output=None,
            )
        )
        envelope = json.loads(envelope_path.read_text(encoding="utf-8"))
        self.assertGreater(
            envelope["recommended_sigr_max"],
            envelope["observed_maximum_integrand"],
        )

        production_args = _prepare_args(
            self.root, self.config, self.legacy, operation="production"
        )
        production_args.output = self.root / "production-smoke"
        production_args.calibration = envelope_path
        production_args.total_events = 2
        production_args.events_per_job = 2
        production_args.lund_files_per_directory = 5_000
        production_manifest = mode3.prepare_production(production_args)
        production_run = mode3.run(
            argparse.Namespace(
                manifest=production_manifest,
                replica_index=0,
                executable=executable,
                scratch_root=None,
                overwrite=False,
            )
        )
        generated = json.loads(production_run.read_text(encoding="utf-8"))
        self.assertGreaterEqual(generated["events"], 2)
        self.assertLess(generated["event_overshoot"], generated["mcall_max"])
        self.assertEqual(generated["lund_lines"], 5 * generated["events"])
        self.assertTrue(Path(generated["lund"]).is_file())

    def test_executable_mixture_automatic_envelope_smoke(self) -> None:
        executable = Path(__file__).resolve().parent / "build" / "aao_rad"
        if not executable.is_file():
            self.skipTest("build/aao_rad is created by the Makefile test target")
        settings = mode3._configuration(self.config, 0.0)
        records = mode3._legacy_records(self.legacy)
        base = records[:17]
        base[15] = "2"
        base[16] = "2.0"
        generated = "\n".join(base) + "\n" + mode3._trailer(
            settings,
            seed=791001,
            replica=0,
            direct_fraction=0.75,
            operation=2,
            trials=2_000,
            heartbeat=100,
        )
        work = self.root / "automatic-envelope"
        work.mkdir()
        completed = subprocess.run(
            [str(executable)],
            input=generated,
            text=True,
            cwd=work,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("Mode-3 automatic sigr_max", completed.stdout)
        norm = mode3._parse_norm(work / "aao_rad.norm")
        self.assertEqual(norm["dvmp_boundary_version"], "1")
        self.assertEqual(int(mode3._number(norm, "mode3_operation")), 0)
        events = int(mode3._number(norm, "events"))
        maximum_multiplicity = int(mode3._number(norm, "mcall_max"))
        self.assertGreaterEqual(events, 2)
        self.assertLess(events - 2, maximum_multiplicity)
        self.assertGreater(mode3._number(norm, "sigr_max"), 0.0)
        self.assertEqual(
            len((work / "aao_rad.lund").read_text().splitlines()),
            5 * events,
        )

    def test_executable_nonfinite_ratio_has_distinct_diagnostic(self) -> None:
        executable = Path(__file__).resolve().parent / "build" / "aao_rad"
        if not executable.is_file():
            self.skipTest("build/aao_rad is created by the Makefile test target")
        settings = mode3._configuration(self.config, 0.0)
        base = mode3._legacy_records(self.legacy)[:17]
        base[15] = "2"
        base[16] = "NaN"  # Invalid scan multiplier must never yield usable output.
        generated = "\n".join(base) + "\n" + mode3._trailer(
            settings, seed=791001, replica=0, direct_fraction=0.75,
            operation=2, trials=2_000, heartbeat=100,
        )
        work = self.root / "nonfinite-envelope"
        work.mkdir()
        completed = subprocess.run(
            [str(executable)], input=generated, text=True, cwd=work,
            capture_output=True, timeout=60,
        )
        output = completed.stdout + completed.stderr
        self.assertNotEqual(completed.returncode, 0, output)
        self.assertIn("FATAL mode-3 non-finite cross section", output)
        self.assertNotIn("multiplicity exceeds job capacity", output)
        self.assertFalse((work / "aao_rad.lund").exists())
        self.assertFalse((work / "aao_rad.norm").exists())

    def test_executable_capacity_guard_fails_before_staging_lund(self) -> None:
        executable = Path(__file__).resolve().parent / "build" / "aao_rad"
        if not executable.is_file():
            self.skipTest("build/aao_rad is created by the Makefile test target")
        calibration_args = _prepare_args(
            self.root, self.config, self.legacy, operation="calibration"
        )
        calibration_manifest = mode3.prepare_calibration(calibration_args)
        calibration = json.loads(
            calibration_manifest.read_text(encoding="utf-8")
        )
        envelope_path = self.root / "deliberately-unsafe-envelope.json"
        envelope_path.write_text(
            json.dumps(
                {
                    "schema": mode3.CALIBRATION_SCHEMA,
                    "settings": calibration["settings"],
                    "direct_fraction": 0.75,
                    "recommended_sigr_max": 1.0e-30,
                }
            ),
            encoding="utf-8",
        )
        production_args = _prepare_args(
            self.root, self.config, self.legacy, operation="production"
        )
        production_args.output = self.root / "capacity-guard"
        production_args.calibration = envelope_path
        production_args.total_events = 2
        production_args.events_per_job = 2
        production_args.lund_files_per_directory = 5_000
        manifest_path = mode3.prepare_production(production_args)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        with self.assertRaisesRegex(
            mode3.Mode3Error, "FATAL mode-3 multiplicity exceeds job capacity"
        ):
            mode3.run(
                argparse.Namespace(
                    manifest=manifest_path,
                    replica_index=0,
                    executable=executable,
                    scratch_root=None,
                    overwrite=False,
                )
            )
        self.assertFalse(Path(manifest["runs"][0]["lund"]).exists())

    def test_swif_emission_is_chunkable_and_shell_valid(self) -> None:
        args = _prepare_args(
            self.root, self.config, self.legacy, operation="calibration"
        )
        args.replicas = 4
        manifest_path = mode3.prepare_calibration(args)
        executable = self.root / "aao_rad"
        executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        executable.chmod(0o755)
        script = mode3.emit_swif(
            argparse.Namespace(
                manifest=manifest_path,
                workflow="mode3_test",
                executable=executable,
                task_start=1,
                task_stop=3,
                cores=1,
                ram="2gb",
                disk="2gb",
                walltime="24hr",
                scratch_root=self.root / "scratch",
                output=None,
            )
        )
        text = script.read_text(encoding="utf-8")
        self.assertEqual(text.count("swif2 add-job"), 2)
        wrapper_text = (args.output / "run_swif_000001_000003.sh").read_text(
            encoding="utf-8"
        )
        self.assertIn("--scratch-root", wrapper_text)
        self.assertIn("SWIF_JOB_WORK_DIR", wrapper_text)
        self.assertIn(str(self.root / "scratch"), wrapper_text)
        self.assertLess(
            wrapper_text.index("SWIF_JOB_WORK_DIR"),
            wrapper_text.index(str(self.root / "scratch")),
        )
        completed = subprocess.run(
            ["bash", "-n", str(script)], capture_output=True, text=True
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)


if __name__ == "__main__":
    unittest.main()
