#!/usr/bin/env python3
"""Tests for the radiative bin-conditional core-plus-tail workflow."""

from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import radiative_mode4
import radiative_survey


def _raw_axis(
    name: str,
    *,
    periodic: bool = False,
    lower: float = 0.0,
    upper: float = 1.0,
) -> dict:
    if not periodic:
        return {
            "name": name,
            "periodic": False,
            "lower": lower,
            "upper": upper,
            "width": upper - lower,
        }
    return {
        "name": name,
        "periodic": True,
        "origin": 0.0,
        "lower_relative_to_origin": -0.5,
        "upper_relative_to_origin": 0.5,
        "width": 1.0,
        "full_period": True,
        "interval_start": 0.5,
        "interval_end": 0.5,
        "wraps": False,
    }


def _fixtures(root: Path) -> tuple[Path, Path, Path]:
    config = {
        "beam_energy": 6.535,
        "target_mass": radiative_survey.PROTON_MASS_GEV,
        "phase_space": {
            "Q2_min": 0.2,
            "W_min": 1.08,
            "y_max": 0.95,
        },
        "binning": {
            "Q2": [0.2, 5.0],
            "xB": [0.02, 0.99],
            "minus_t": [0.0, 10.0],
            "phi_deg": [0.0, 360.0],
        },
    }
    config_path = root / "config.json"
    config_path.write_text(
        json.dumps(config, indent=2) + "\n", encoding="utf-8"
    )
    recipes = {
        "schema": radiative_mode4.RECIPE_SCHEMA,
        "analysis_config_sha256": hashlib.sha256(
            config_path.read_bytes()
        ).hexdigest(),
        "continuous_guard_learner_revision": "learner-test",
        "generator_revision": "survey-test",
        "padding_candidates": [
            {"identifier": "padding_0p035", "base_padding": 0.035}
        ],
        "strata": {
            "s00000": {
                "bounds": config["binning"],
                "padding_scale": 1.0,
                "training_status": "learned",
                "fit_source": "local_only",
                "raw_box": {
                    "axes": [
                        _raw_axis("r_u", lower=0.2, upper=0.8),
                        _raw_axis("r_ep"),
                        _raw_axis(
                            "u_gamma", lower=0.1, upper=0.8
                        ),
                        _raw_axis("hadron_cosine_base"),
                        _raw_axis("hadron_phi_base", periodic=True),
                    ]
                },
            }
        },
    }
    recipes_path = root / "recipes.json"
    recipes_path.write_text(
        json.dumps(recipes, indent=2) + "\n", encoding="utf-8"
    )
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
            "6.535",
            ".2 5",
            ".2 6.535",
            ".005",
            "2",
            "0",
            ".005",
        ]
    )
    input_path = root / "legacy.inp"
    input_path.write_text(legacy + "\n", encoding="utf-8")
    return config_path, recipes_path, input_path


def _prepare_args(
    root: Path, config: Path, recipes: Path, legacy: Path
) -> argparse.Namespace:
    return argparse.Namespace(
        command="prepare",
        config=config,
        recipes=recipes,
        input=legacy,
        output=root / "campaign",
        candidate="padding_0p035",
        core_fraction=0.9,
        sigr_max=0.005,
        events_per_stratum=2,
        replicas=1,
        seed_base=481001,
        bin_start=0,
        bin_stop=None,
        apply_y_max=False,
        generator_revision="generator-test",
        heartbeat_interval=100,
        overwrite=False,
    )


def _calibration_args(
    root: Path, config: Path, recipes: Path, legacy: Path
) -> argparse.Namespace:
    return argparse.Namespace(
        command="prepare-calibration",
        config=config,
        recipes=recipes,
        input=legacy,
        output=root / "calibration",
        candidate="padding_0p035",
        core_fraction=0.9,
        trials=1000,
        calibration_core_fraction=0.5,
        replicas=1,
        seed_base=581001,
        bin_start=0,
        bin_stop=None,
        apply_y_max=False,
        generator_revision="generator-test",
        heartbeat_interval=100,
        overwrite=False,
    )


class ProposalTests(unittest.TestCase):
    def test_periodic_guard_wrap_and_exact_mixture_identity(self) -> None:
        box = radiative_mode4.GuardBox(
            nonperiodic={
                name: (0.2, 0.6) for name in radiative_mode4.AXES[:-1]
            },
            phi_origin=0.0,
            phi_relative=(-0.1, 0.1),
        )
        inside = {
            name: 0.4 for name in radiative_mode4.AXES[:-1]
        }
        inside["hadron_phi_base"] = 0.99
        outside = dict(inside)
        outside["hadron_phi_base"] = 0.5
        self.assertTrue(box.contains(inside))
        self.assertFalse(box.contains(outside))
        core_fraction = 0.9
        tail_fraction = 1.0 - core_fraction
        inside_ratio = radiative_mode4.proposal_density_ratio(
            inside, box, core_fraction
        )
        outside_ratio = radiative_mode4.proposal_density_ratio(
            outside, box, core_fraction
        )
        mixture_inside_mass = core_fraction + tail_fraction * box.volume
        mixture_outside_mass = tail_fraction * (1.0 - box.volume)
        self.assertAlmostEqual(
            mixture_inside_mass * inside_ratio
            + mixture_outside_mass * outside_ratio,
            1.0,
        )
        self.assertAlmostEqual(outside_ratio, 1.0 / tail_fraction)

    def test_recipe_padding_is_clipped_to_native_domain(self) -> None:
        recipe = {
            "padding_scale": 2.0,
            "raw_box": {
                "axes": [
                    {
                        "name": name,
                        "periodic": False,
                        "lower": 0.01,
                        "upper": 0.99,
                    }
                    for name in radiative_mode4.AXES[:-1]
                ]
                + [
                    {
                        "name": "hadron_phi_base",
                        "periodic": True,
                        "origin": 0.95,
                        "lower_relative_to_origin": -0.49,
                        "upper_relative_to_origin": 0.49,
                    }
                ]
            },
        }
        box = radiative_mode4.reconstruct_guard_box(recipe, 0.035)
        self.assertEqual(box.nonperiodic["r_u"], (0.0, 1.0))
        self.assertEqual(box.nonperiodic["u_gamma"], (0.0, 1.0))
        self.assertEqual(box.phi_relative, (-0.5, 0.5))
        self.assertAlmostEqual(box.volume, 1.0)

    def test_u_gamma_guard_is_anchored_at_soft_endpoint(self) -> None:
        recipe = {
            "padding_scale": 1.0,
            "raw_box": {
                "axes": [
                    _raw_axis("r_u", lower=0.2, upper=0.8),
                    _raw_axis("r_ep"),
                    _raw_axis("u_gamma", lower=0.3, upper=0.7),
                    _raw_axis("hadron_cosine_base"),
                    _raw_axis("hadron_phi_base", periodic=True),
                ]
            },
        }
        box = radiative_mode4.reconstruct_guard_box(recipe, 0.05)
        self.assertEqual(box.nonperiodic["u_gamma"], (0.25, 1.0))
        self.assertAlmostEqual(box.volume, 0.7 * 0.75)


class WorkflowTests(unittest.TestCase):
    def test_prepare_snapshots_provenance_and_writes_mode4_trailer(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config, recipes, legacy = _fixtures(root)
            manifest_path = radiative_mode4.prepare(
                _prepare_args(root, config, recipes, legacy)
            )
            manifest = json.loads(
                manifest_path.read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["schema"], radiative_mode4.MANIFEST_SCHEMA)
            self.assertEqual(manifest["sampling_mode"], 4)
            self.assertTrue(manifest["full_support_guaranteed"])
            self.assertFalse(manifest["production_ready"])
            self.assertEqual(manifest["guard_candidate"], "padding_0p035")
            self.assertAlmostEqual(manifest["legacy_tail_fraction"], 0.1)
            record = manifest["runs"][0]
            self.assertEqual(record["stratum_id"], "s00000")
            self.assertEqual(record["generation_id"], "g0000")
            self.assertIn("s00000__g0000", record["output_stem"])
            self.assertAlmostEqual(
                record["guard"]["normalized_volume"], 0.67 * 0.935
            )
            self.assertTrue(
                record["guard"]["u_gamma_soft_endpoint_anchored"]
            )
            prepared = (
                manifest_path.parent / record["input_file"]
            ).read_text(encoding="utf-8")
            parsed = radiative_survey._records(prepared)
            legacy_count = radiative_survey._legacy_record_count(
                parsed, Path("prepared.inp")
            )
            self.assertEqual(parsed[legacy_count], "4")
            self.assertEqual(parsed[legacy_count + 1], "481001")
            self.assertEqual(
                parsed[legacy_count + 3].split(),
                ["0", "0", "0", "0", "0"],
            )
            self.assertEqual(parsed[-5:-2], ["0", "0", "100"])
            self.assertAlmostEqual(float(parsed[-2]), 0.9)
            self.assertEqual(parsed[-1], "0")
            self.assertEqual(
                hashlib.sha256(
                    (manifest_path.parent / "analysis_config.json").read_bytes()
                ).hexdigest(),
                manifest["analysis_config_sha256"],
            )

    @unittest.skipUnless(
        Path(__file__).with_name("build").joinpath("aao_rad").is_file(),
        "build/aao_rad is required for the end-to-end smoke test",
    )
    def test_end_to_end_mode4_lund_and_weights(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config, recipes, legacy = _fixtures(root)
            manifest_path = radiative_mode4.prepare(
                _prepare_args(root, config, recipes, legacy)
            )
            executable = (
                Path(__file__).with_name("build").joinpath("aao_rad")
            )
            run_path = radiative_mode4.run(
                argparse.Namespace(
                    manifest=manifest_path,
                    flat_index=0,
                    replica_index=0,
                    executable=executable,
                    overwrite=False,
                )
            )
            completed = json.loads(run_path.read_text(encoding="utf-8"))
            self.assertEqual(completed["schema"], radiative_mode4.RUN_SCHEMA)
            self.assertGreaterEqual(completed["events"], 2)
            self.assertGreater(completed["ntries"], 0)
            self.assertGreater(completed["sig_sum_microbarn"], 0.0)
            self.assertEqual(
                completed["core_trials"] + completed["legacy_trials"],
                completed["ntries"],
            )
            stem = run_path.with_suffix("")
            lund_lines = [
                line
                for line in Path(str(stem) + ".lund")
                .read_text(encoding="utf-8")
                .splitlines()
                if line.strip()
            ]
            self.assertEqual(len(lund_lines), 5 * completed["events"])
            self.assertTrue(all(
                lund_lines[offset].split()[0] == "4"
                for offset in range(0, len(lund_lines), 5)
            ))
            norm = radiative_survey.parse_norm(
                Path(str(stem) + ".norm")
            )
            self.assertEqual(norm["sampling_mode"], "4")
            self.assertEqual(
                norm["mode4_guard_candidate"], "padding_0p035"
            )
            weights_path = radiative_mode4.finalize(
                argparse.Namespace(manifest=manifest_path, output=None)
            )
            weights = json.loads(weights_path.read_text(encoding="utf-8"))
            self.assertEqual(weights["schema"], radiative_mode4.WEIGHTS_SCHEMA)
            self.assertEqual(weights["stratum_count"], 1)
            stratum = weights["strata"][0]
            self.assertEqual(stratum["total_events"], completed["events"])
            self.assertAlmostEqual(
                stratum["pooled_event_weight_microbarn"]
                * stratum["total_events"],
                stratum["combined_sig_sum_microbarn"],
            )

    @unittest.skipUnless(
        Path(__file__).with_name("build").joinpath("aao_rad").is_file(),
        "build/aao_rad is required for the end-to-end calibration test",
    )
    def test_fixed_trial_calibration_heartbeat_and_envelope(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config, recipes, legacy = _fixtures(root)
            manifest_path = radiative_mode4.prepare(
                _calibration_args(root, config, recipes, legacy)
            )
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["operation"], "calibration")
            self.assertEqual(manifest["calibration_trials_per_replica"], 1000)
            self.assertAlmostEqual(
                manifest["component_core_fraction"], 0.5
            )
            self.assertEqual(
                manifest["calibration_proposal"], "guard_partition"
            )
            run_path = radiative_mode4.run(
                argparse.Namespace(
                    manifest=manifest_path,
                    flat_index=0,
                    replica_index=0,
                    executable=Path(__file__).with_name("build")
                    / "aao_rad",
                    overwrite=False,
                )
            )
            completed = json.loads(run_path.read_text(encoding="utf-8"))
            self.assertEqual(completed["events"], 0)
            self.assertEqual(completed["ntries"], 1000)
            self.assertEqual(
                completed["final_heartbeat"]["proposals"], 1000
            )
            self.assertGreater(completed["calibration_target_rows"], 0)
            stem = run_path.with_suffix("")
            self.assertTrue(
                Path(str(stem) + ".calibration.csv").is_file()
            )
            heartbeat = Path(str(stem) + ".heartbeat.csv").read_text(
                encoding="utf-8"
            )
            self.assertIn("1000,0,", heartbeat)
            second_args = _calibration_args(
                root, config, recipes, legacy
            )
            second_args.output = root / "calibration_complement_heavy"
            second_args.calibration_core_fraction = 0.2
            second_args.seed_base = 681001
            second_manifest_path = radiative_mode4.prepare(second_args)
            second_run_path = radiative_mode4.run(
                argparse.Namespace(
                    manifest=second_manifest_path,
                    flat_index=0,
                    replica_index=0,
                    executable=Path(__file__).with_name("build")
                    / "aao_rad",
                    overwrite=False,
                )
            )
            second_completed = json.loads(
                second_run_path.read_text(encoding="utf-8")
            )
            report_path = radiative_mode4.finalize(
                argparse.Namespace(
                    manifests=[manifest_path, second_manifest_path],
                    output=None,
                    envelope_safety_factor=1.2,
                    maximum_duplicate_fraction=0.05,
                    minimum_component_targets=5,
                )
            )
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(
                report["schema"], radiative_mode4.CALIBRATION_SCHEMA
            )
            self.assertEqual(report["stratum_count"], 1)
            self.assertEqual(len(report["source_manifests"]), 2)
            result = report["strata"][0]
            self.assertGreater(
                result["inside_guard"]["target_candidates"], 0
            )
            self.assertGreater(
                result["guard_complement"]["target_candidates"], 0
            )
            self.assertEqual(
                result["inside_guard"]["trials"],
                completed["core_trials"] + second_completed["core_trials"],
            )
            self.assertEqual(
                result["guard_complement"]["trials"],
                completed["noncore_trials"]
                + second_completed["noncore_trials"],
            )
            self.assertTrue(result["envelope_candidates"])
            self.assertEqual(
                result["recommendation_status"], "recommended"
            )
            self.assertGreater(
                result["recommended_envelope"]["sigr_max"], 0.0
            )
            self.assertLessEqual(
                result["recommended_envelope"][
                    "expected_duplicate_event_fraction"
                ],
                0.05,
            )


if __name__ == "__main__":
    unittest.main()
