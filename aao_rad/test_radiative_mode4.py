#!/usr/bin/env python3
"""Tests for the radiative bin-conditional core-plus-tail workflow."""

from __future__ import annotations

import argparse
import csv
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
        refinements=None,
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
        refinements=None,
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

    def test_refinement_expands_only_requested_final_faces(self) -> None:
        box = radiative_mode4.GuardBox(
            nonperiodic={
                "r_u": (0.2, 0.6),
                "r_ep": (0.3, 0.7),
                "u_gamma": (0.1, 1.0),
                "hadron_cosine_base": (0.4, 0.8),
            },
            phi_origin=0.9,
            phi_relative=(-0.1, 0.1),
        )
        specification = {
            "rationale": "Independent calibration crossed two nearby faces.",
            "evidence": [{"path": "evidence.json", "sha256": "a" * 64}],
            "faces": {
                "r_u": {"upper": 0.67},
                "r_ep": {"lower": 0.275},
            },
        }
        refined, record = radiative_mode4.apply_guard_refinement(
            box, specification, stratum_id="s04468"
        )
        self.assertEqual(refined.nonperiodic["r_u"], (0.2, 0.67))
        self.assertEqual(refined.nonperiodic["r_ep"], (0.275, 0.7))
        self.assertEqual(refined.nonperiodic["u_gamma"], (0.1, 1.0))
        self.assertEqual(box.nonperiodic["r_u"], (0.2, 0.6))
        self.assertGreater(refined.volume, box.volume)
        self.assertEqual(len(record["applied_face_changes"]), 2)
        self.assertEqual(
            record["coordinate_space"],
            radiative_mode4.REFINEMENT_COORDINATE_SPACE,
        )

    def test_refinement_rejects_guard_contraction(self) -> None:
        box = radiative_mode4.GuardBox(
            nonperiodic={
                name: (0.2, 0.8)
                for name in radiative_mode4.AXES[:-1]
            },
            phi_origin=0.0,
            phi_relative=(-0.2, 0.2),
        )
        with self.assertRaisesRegex(ValueError, "would contract"):
            radiative_mode4.apply_guard_refinement(
                box,
                {
                    "rationale": "Invalid inward movement.",
                    "evidence": [
                        {"path": "evidence.json", "sha256": "b" * 64}
                    ],
                    "faces": {"r_u": {"upper": 0.7}},
                },
                stratum_id="s00000",
            )

    def test_zero_success_upper_rate_is_exact_and_stable(self) -> None:
        trials = 5_750_000
        confidence = 0.95
        expected = 1.0 - (1.0 - confidence) ** (1.0 / trials)
        observed = radiative_mode4._zero_success_upper_rate(
            trials, confidence
        )
        self.assertAlmostEqual(observed, expected, places=15)
        self.assertLess(observed, 1.0e-6)


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

    def test_create_and_prepare_refinement_freezes_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config, recipes, legacy = _fixtures(root)
            evidence = root / "calibration.json"
            evidence.write_text(
                '{"diagnosis":"upper r_u escape"}\n', encoding="utf-8"
            )
            refinement_path = radiative_mode4.create_refinement(
                argparse.Namespace(
                    config=config,
                    recipes=recipes,
                    candidate="padding_0p035",
                    output=root / "refinement.json",
                    stratum="s00000",
                    face=["r_u:upper:0.87"],
                    rationale=(
                        "Independent complement calibration crossed r_u."
                    ),
                    evidence=[evidence],
                    overwrite=False,
                )
            )
            refinement = json.loads(
                refinement_path.read_text(encoding="utf-8")
            )
            self.assertEqual(
                refinement["schema"], radiative_mode4.REFINEMENT_SCHEMA
            )
            self.assertEqual(
                refinement["strata"]["s00000"]["evidence"][0]["sha256"],
                hashlib.sha256(evidence.read_bytes()).hexdigest(),
            )
            args = _prepare_args(root, config, recipes, legacy)
            args.output = root / "refined_campaign"
            args.refinements = refinement_path
            manifest_path = radiative_mode4.prepare(args)
            manifest = json.loads(
                manifest_path.read_text(encoding="utf-8")
            )
            self.assertEqual(
                manifest["guard_refinements_sha256"],
                hashlib.sha256(refinement_path.read_bytes()).hexdigest(),
            )
            self.assertEqual(manifest["guard_refined_strata"], ["s00000"])
            self.assertEqual(
                hashlib.sha256(
                    (
                        manifest_path.parent / "guard_refinements.json"
                    ).read_bytes()
                ).hexdigest(),
                manifest["guard_refinements_sha256"],
            )
            record = manifest["runs"][0]
            self.assertAlmostEqual(
                record["guard_original"]["axes"]["r_u"][0], 0.165
            )
            self.assertAlmostEqual(
                record["guard_original"]["axes"]["r_u"][1], 0.835
            )
            self.assertAlmostEqual(
                record["guard"]["axes"]["r_u"][0], 0.165
            )
            self.assertAlmostEqual(
                record["guard"]["axes"]["r_u"][1], 0.87
            )
            self.assertEqual(
                record["guard_refinement"]["applied_face_changes"][0][
                    "axis"
                ],
                "r_u",
            )
            self.assertAlmostEqual(
                record["guard_refinement"]["volume_ratio"],
                0.705 / 0.67,
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
            self.assertEqual(
                completed["emitting_candidates"]
                + completed["duplicate_events"],
                completed["events"],
            )
            self.assertAlmostEqual(
                completed["duplicate_event_fraction"],
                completed["duplicate_events"] / completed["events"],
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
            self.assertEqual(
                stratum["emitting_candidates"] + stratum["duplicate_events"],
                stratum["total_events"],
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
            self.assertEqual(completed["emitting_candidates"], 0)
            self.assertEqual(completed["duplicate_events"], 0)
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
            mismatch_args = _calibration_args(
                root, config, recipes, legacy
            )
            mismatch_args.output = root / "calibration_new_revision"
            mismatch_args.generator_revision = "generator-diagnostics-test"
            mismatch_args.seed_base = 781001
            mismatch_manifest_path = radiative_mode4.prepare(mismatch_args)
            radiative_mode4.run(
                argparse.Namespace(
                    manifest=mismatch_manifest_path,
                    flat_index=0,
                    replica_index=0,
                    executable=Path(__file__).with_name("build")
                    / "aao_rad",
                    overwrite=False,
                )
            )
            with self.assertRaisesRegex(
                radiative_mode4.Mode4Error,
                "generator_revision is incompatible",
            ):
                radiative_mode4.finalize(
                    argparse.Namespace(
                        manifests=[
                            manifest_path,
                            mismatch_manifest_path,
                        ],
                        output=root / "revision_mismatch_rejected.json",
                    )
                )
            override_path = radiative_mode4.finalize(
                argparse.Namespace(
                    manifests=[manifest_path, mismatch_manifest_path],
                    output=root / "revision_mismatch_audited.json",
                    envelope_safety_factor=1.2,
                    maximum_duplicate_fraction=0.05,
                    minimum_component_targets=5,
                    minimum_provisional_inside_targets=5,
                    allow_calibration_revision_mismatch=True,
                    revision_compatibility_rationale=(
                        "test-only diagnostic revision; calibration "
                        "physics and proposal are unchanged"
                    ),
                )
            )
            override = json.loads(
                override_path.read_text(encoding="utf-8")
            )
            self.assertIsNone(override["generator_revision"])
            self.assertEqual(len(override["generator_revisions"]), 2)
            self.assertTrue(
                override["revision_compatibility_override"]["enabled"]
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
            self.assertEqual(len(report["finalizer_source_sha256"]), 64)
            self.assertTrue(report["finalizer_revision"])
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
            for source_manifest in (
                manifest_path,
                second_manifest_path,
            ):
                source = json.loads(
                    source_manifest.read_text(encoding="utf-8")
                )
                for record in source["runs"]:
                    csv_path = source_manifest.parent / (
                        record["output_stem"] + ".calibration.csv"
                    )
                    lines = csv_path.read_text(encoding="utf-8").splitlines()
                    inside_only = lines[:2] + [
                        line
                        for line in lines[2:]
                        if line.split(",")[1] == "1"
                    ]
                    csv_path.write_text(
                        "\n".join(inside_only) + "\n", encoding="utf-8"
                    )
            strict_zero_path = radiative_mode4.finalize(
                argparse.Namespace(
                    manifests=[manifest_path, second_manifest_path],
                    output=root / "strict_zero.json",
                    envelope_safety_factor=1.2,
                    maximum_duplicate_fraction=0.05,
                    minimum_component_targets=5,
                    allow_zero_complement=False,
                    zero_complement_confidence=0.95,
                    maximum_zero_complement_target_rate=0.01,
                )
            )
            strict_zero = json.loads(
                strict_zero_path.read_text(encoding="utf-8")
            )["strata"][0]
            self.assertEqual(
                strict_zero["recommendation_status"],
                "insufficient_guard_complement_target_support",
            )
            insufficient_path = radiative_mode4.finalize(
                argparse.Namespace(
                    manifests=[manifest_path, second_manifest_path],
                    output=root / "insufficient_zero_exposure.json",
                    envelope_safety_factor=1.2,
                    maximum_duplicate_fraction=0.05,
                    minimum_component_targets=5,
                    allow_zero_complement=True,
                    zero_complement_confidence=0.95,
                    maximum_zero_complement_target_rate=1.0e-12,
                )
            )
            insufficient = json.loads(
                insufficient_path.read_text(encoding="utf-8")
            )["strata"][0]
            self.assertEqual(
                insufficient["recommendation_status"],
                "insufficient_zero_complement_exposure",
            )
            insufficient_inside_path = radiative_mode4.finalize(
                argparse.Namespace(
                    manifests=[manifest_path, second_manifest_path],
                    output=root / "insufficient_inside_envelope.json",
                    envelope_safety_factor=1.2,
                    maximum_duplicate_fraction=0.05,
                    minimum_component_targets=5,
                    minimum_provisional_inside_targets=1000000,
                    allow_zero_complement=True,
                    zero_complement_confidence=0.95,
                    maximum_zero_complement_target_rate=0.01,
                )
            )
            insufficient_inside = json.loads(
                insufficient_inside_path.read_text(encoding="utf-8")
            )["strata"][0]
            self.assertEqual(
                insufficient_inside["recommendation_status"],
                "insufficient_provisional_inside_envelope_support",
            )
            self.assertFalse(
                insufficient_inside[
                    "provisional_inside_envelope_support"
                ]["passes_target_threshold"]
            )
            provisional_path = radiative_mode4.finalize(
                argparse.Namespace(
                    manifests=[manifest_path, second_manifest_path],
                    output=root / "provisional.json",
                    envelope_safety_factor=1.2,
                    maximum_duplicate_fraction=0.05,
                    minimum_component_targets=5,
                    minimum_provisional_inside_targets=5,
                    allow_zero_complement=True,
                    zero_complement_confidence=0.95,
                    maximum_zero_complement_target_rate=0.01,
                )
            )
            provisional_report = json.loads(
                provisional_path.read_text(encoding="utf-8")
            )
            provisional = provisional_report["strata"][0]
            self.assertEqual(
                provisional["recommendation_status"],
                "provisional_zero_complement",
            )
            self.assertEqual(
                provisional["pilot_readiness"],
                "ready_provisional_zero_complement",
            )
            self.assertTrue(
                provisional["recommended_envelope"]["provisional"]
            )
            self.assertEqual(
                provisional["recommended_envelope"]["source"],
                "inside_guard_maximum",
            )
            self.assertEqual(
                provisional["recommended_envelope"][
                    "expected_duplicate_event_fraction"
                ],
                0.0,
            )
            self.assertEqual(
                provisional["recommended_envelope"]["sigr_max"],
                provisional["provisional_inside_envelope_support"][
                    "safety_scaled_observed_maximum_sigr_max"
                ],
            )
            source_record = json.loads(
                manifest_path.read_text(encoding="utf-8")
            )["runs"][0]
            pilot_run_path = root / "inside_only_pilot.json"
            pilot_run_path.write_text(
                json.dumps(
                    {
                        **{
                            name: source_record[name]
                            for name in (
                                "stratum_id",
                                "flat_index",
                                "indices",
                                "bounds",
                                "guard",
                            )
                        },
                        "schema": radiative_mode4.RUN_SCHEMA,
                        "operation": "generation",
                        "events": 1,
                        "noncore_events": 0,
                    },
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
            pilot_event_path = pilot_run_path.with_suffix(".mode4.csv")
            with pilot_event_path.open(
                "w", encoding="utf-8", newline=""
            ) as destination:
                destination.write(
                    "# schema="
                    f"{radiative_mode4.MODE4_KINEMATICS_SCHEMA}\n"
                )
                writer = csv.DictWriter(
                    destination,
                    fieldnames=radiative_mode4.MODE4_KINEMATICS_COLUMNS,
                )
                writer.writeheader()
                writer.writerow(
                    {
                        name: (
                            1
                            if name in ("event", "proposal_component")
                            else 9.0
                            if name == "integrand_corrected"
                            else 0.5
                        )
                        for name in (
                            radiative_mode4.MODE4_KINEMATICS_COLUMNS
                        )
                    }
                )
            pilot_floor_path = radiative_mode4.finalize(
                argparse.Namespace(
                    manifests=[manifest_path, second_manifest_path],
                    output=root / "provisional_with_pilot_floor.json",
                    envelope_safety_factor=1.2,
                    maximum_duplicate_fraction=0.05,
                    minimum_component_targets=5,
                    minimum_provisional_inside_targets=5,
                    allow_zero_complement=True,
                    zero_complement_confidence=0.95,
                    maximum_zero_complement_target_rate=0.01,
                    additional_inside_pilot_run=[pilot_run_path],
                )
            )
            pilot_floor = json.loads(
                pilot_floor_path.read_text(encoding="utf-8")
            )["strata"][0]
            self.assertEqual(
                pilot_floor["recommended_envelope"]["source"],
                "additional_inside_pilot_observed_maximum",
            )
            self.assertAlmostEqual(
                pilot_floor["recommended_envelope"]["sigr_max"], 10.8
            )
            self.assertEqual(
                len(pilot_floor["additional_inside_pilot_observations"]),
                1,
            )
            self.assertTrue(
                pilot_floor["additional_inside_pilot_observations"][0][
                    "excluded_from_fixed_trial_calibration_statistics"
                ]
            )
            self.assertTrue(
                provisional["zero_complement_stopping_test"][
                    "passes_rate_threshold"
                ]
            )
            self.assertTrue(
                provisional_report["zero_complement_policy"]["enabled"]
            )
            self.assertIn(
                "zero_complement_upper_target_rate",
                provisional_path.with_suffix(".tsv").read_text(
                    encoding="utf-8"
                ).splitlines()[0],
            )


if __name__ == "__main__":
    unittest.main()
