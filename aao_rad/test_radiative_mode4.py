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
from unittest import mock

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


def _multistratum_fixtures(
    root: Path,
) -> tuple[Path, Path, Path, Path]:
    config_path, recipes_path, input_path = _fixtures(root)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    config["binning"]["minus_t"] = [0.0, 1.0, 2.0]
    config_path.write_text(
        json.dumps(config, indent=2) + "\n", encoding="utf-8"
    )
    recipes = json.loads(recipes_path.read_text(encoding="utf-8"))
    recipes["analysis_config_sha256"] = hashlib.sha256(
        config_path.read_bytes()
    ).hexdigest()
    prototype = recipes["strata"]["s00000"]
    recipes["strata"] = {}
    guards: dict[str, dict] = {}
    for index, bounds in enumerate(((0.0, 1.0), (1.0, 2.0))):
        stratum_id = f"s{index:05d}"
        record = json.loads(json.dumps(prototype))
        record["bounds"] = {
            **config["binning"],
            "minus_t": list(bounds),
        }
        recipes["strata"][stratum_id] = record
    recipes_path.write_text(
        json.dumps(recipes, indent=2) + "\n", encoding="utf-8"
    )
    for stratum_id, recipe in recipes["strata"].items():
        guards[stratum_id] = radiative_mode4.reconstruct_guard_box(
            recipe, 0.035
        ).manifest_record()
    report = {
        "schema": radiative_mode4.CALIBRATION_SCHEMA,
        "analysis_config_sha256": hashlib.sha256(
            config_path.read_bytes()
        ).hexdigest(),
        "guard_recipes_sha256": hashlib.sha256(
            recipes_path.read_bytes()
        ).hexdigest(),
        "guard_refinements_sha256": None,
        "guard_candidate": "padding_0p035",
        "core_fraction": 0.9,
        "generator_revision": "generator-test",
        "generator_revisions": ["generator-test"],
        "analysis_selection": {
            "coordinate_definition": "final_lund_analysis",
            "w_minimum": 1.08,
            "apply_y_max": False,
            "y_maximum": None,
            "no_implicit_y_minimum": True,
        },
        "strata": [
            {
                "stratum_id": f"s{index:05d}",
                "flat_index": index,
                "indices": {
                    "iq2": 0,
                    "ixb": 0,
                    "it": index,
                    "iphi": 0,
                },
                "bounds": {
                    "Q2": [0.2, 5.0],
                    "xB": [0.02, 0.99],
                    "minus_t": list(bounds),
                    "phi_deg": [0.0, 360.0],
                },
                "guard": guards[f"s{index:05d}"],
                "recommendation_status": "recommended",
                "recommendation_basis": (
                    "both_calibration_components_observed"
                ),
                "pilot_readiness": "ready",
                "recommended_envelope": {
                    "sigr_max": 0.004 + 0.002 * index
                },
            }
            for index, bounds in enumerate(((0.0, 1.0), (1.0, 2.0)))
        ],
    }
    report_path = root / "multi_envelopes.json"
    report_path.write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return config_path, recipes_path, input_path, report_path


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


def _write_pilot_artifacts(
    root: Path,
    *,
    name: str,
    seed: int,
    sigr_max: float,
    sigma: float,
    legacy_outside_duplicate: bool,
) -> Path:
    campaign = root / name
    run_directory = campaign / "runs" / "s00000"
    run_directory.mkdir(parents=True)
    guard = {
        "axes": {
            "r_u": [0.2, 0.6],
            "r_ep": [0.2, 0.8],
            "u_gamma": [0.0, 1.0],
            "hadron_cosine_base": [0.2, 0.8],
        },
        "hadron_phi_base": {
            "origin": 0.5,
            "lower_relative_to_origin": -0.2,
            "upper_relative_to_origin": 0.2,
            "width": 0.4,
        },
        "normalized_volume": 0.0576,
        "u_gamma_soft_endpoint_anchored": True,
    }
    record = {
        "stratum_id": "s00000",
        "flat_index": 0,
        "replica_index": 0,
        "indices": {"iq2": 0, "ixb": 0, "it": 0, "iphi": 0},
        "bounds": {
            "Q2": [1.0, 2.0],
            "xB": [0.2, 0.3],
            "minus_t": [0.1, 0.5],
            "phi_deg": [0.0, 360.0],
        },
        "guard": guard,
        "seed": seed,
        "sigr_max": sigr_max,
    }
    manifest = {
        "schema": radiative_mode4.MANIFEST_SCHEMA,
        "operation": "generation",
        "generator_revision": "pilot-validator-test",
        "analysis_config_sha256": "a" * 64,
        "guard_recipes_sha256": "b" * 64,
        "guard_refinements_sha256": "c" * 64,
        "guard_candidate": "padding_0p035",
        "core_fraction": 0.9,
        "analysis_selection": {
            "w_minimum": 1.08,
            "y_maximum": 0.95,
            "apply_y_max": False,
        },
        "sigr_max": sigr_max,
        "runs": [record],
    }
    manifest_path = campaign / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    rows: list[dict[str, object]] = []
    core_events = 198 if legacy_outside_duplicate else 200
    for index in range(core_events):
        rows.append(
            {
                "event": index + 1,
                "proposal_component": 1,
                "proposal_density_ratio": 1.0 / (0.1 + 0.9 / 0.0576),
                "r_u": 0.3 + index * 1.0e-5,
                "r_ep": 0.4,
                "u_gamma": 0.99,
                "hadron_cosine_base": 0.5,
                "hadron_phi_base": 0.5,
                "q2_observed": 1.5,
                "xb_observed": 0.25,
                "minus_t_observed": 0.3,
                "phi_observed_deg": 180.0,
                "w_observed": 2.0,
                "y_observed": 0.4,
                "integrand_corrected": 1.0 + index * 1.0e-4,
            }
        )
    if legacy_outside_duplicate:
        outside = {
            "proposal_component": 0,
            "proposal_density_ratio": 10.0,
            "r_u": 0.6077,
            "r_ep": 0.4,
            "u_gamma": 0.99,
            "hadron_cosine_base": 0.5,
            "hadron_phi_base": 0.5,
            "q2_observed": 1.5,
            "xb_observed": 0.25,
            "minus_t_observed": 0.3,
            "phi_observed_deg": 180.0,
            "w_observed": 2.0,
            "y_observed": 0.4,
            "integrand_corrected": 7.2,
        }
        rows.extend(
            [
                {"event": 199, **outside},
                {"event": 200, **outside},
            ]
        )
    event_path = run_directory / "s00000__g0000.mode4.csv"
    with event_path.open("w", encoding="utf-8", newline="") as destination:
        destination.write(
            f"# schema={radiative_mode4.MODE4_KINEMATICS_SCHEMA}\n"
        )
        writer = csv.DictWriter(
            destination,
            fieldnames=radiative_mode4.MODE4_KINEMATICS_COLUMNS,
        )
        writer.writeheader()
        writer.writerows(rows)
    run_path = run_directory / "s00000__g0000.json"
    run = {
        **record,
        "schema": radiative_mode4.RUN_SCHEMA,
        "operation": "generation",
        "source_manifest": str(manifest_path),
        "source_manifest_sha256": hashlib.sha256(
            manifest_path.read_bytes()
        ).hexdigest(),
        "events": 200,
        "ntries": 1000000 + seed,
        "sig_sum_microbarn": sigma,
        "event_overshoot": 0,
        "emitting_candidates": (
            199 if legacy_outside_duplicate else 200
        ),
        "duplicate_events": 1 if legacy_outside_duplicate else 0,
        "mcall_max": 2 if legacy_outside_duplicate else 1,
        "core_events": core_events,
        "legacy_events": 2 if legacy_outside_duplicate else 0,
        "noncore_events": 2 if legacy_outside_duplicate else 0,
    }
    run_path.write_text(
        json.dumps(run, indent=2) + "\n", encoding="utf-8"
    )
    return run_path


def _write_calibration_report(
    root: Path, *, legacy_metadata: bool = False
) -> Path:
    path = root / (
        "calibration_legacy.json" if legacy_metadata else "calibration.json"
    )
    payload = {
        "schema": radiative_mode4.CALIBRATION_SCHEMA,
        "analysis_config_sha256": "a" * 64,
        "guard_recipes_sha256": "b" * 64,
        "guard_refinements_sha256": "c" * 64,
        "guard_candidate": "padding_0p035",
        "core_fraction": 0.9,
        "generator_revision": "pilot-validator-test",
        "generator_revisions": ["pilot-validator-test"],
        "analysis_selection": {
            "w_minimum": 1.08,
            "y_maximum": 0.95,
            "apply_y_max": False,
        },
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
                "bounds": {
                    "Q2": [1.0, 2.0],
                    "xB": [0.2, 0.3],
                    "minus_t": [0.1, 0.5],
                    "phi_deg": [0.0, 360.0],
                },
                "guard": {
                    "axes": {
                        "r_u": [0.2, 0.6],
                        "r_ep": [0.2, 0.8],
                        "u_gamma": [0.0, 1.0],
                        "hadron_cosine_base": [0.2, 0.8],
                    },
                    "hadron_phi_base": {
                        "origin": 0.5,
                        "lower_relative_to_origin": -0.2,
                        "upper_relative_to_origin": 0.2,
                        "width": 0.4,
                    },
                    "normalized_volume": 0.0576,
                    "u_gamma_soft_endpoint_anchored": True,
                },
                "recommendation_status": "provisional_zero_complement",
                "pilot_readiness": "ready_provisional_zero_complement",
                "recommended_envelope": {"sigr_max": 4.0},
                "integrated_cross_section_microbarn": 0.3,
                "integrated_cross_section_sem_microbarn": 0.02,
            }
        ],
    }
    if legacy_metadata:
        source_directory = root / "calibration_source"
        source_directory.mkdir()
        source_path = source_directory / "manifest.json"
        source_manifest = {
            "schema": radiative_mode4.MANIFEST_SCHEMA,
            "operation": "calibration",
            "generator_revision": "pilot-validator-test",
            "analysis_selection": payload["analysis_selection"],
            "runs": [
                {
                    "stratum_id": "s00000",
                    "guard": payload["strata"][0]["guard"],
                }
            ],
        }
        source_path.write_text(
            json.dumps(source_manifest, indent=2) + "\n",
            encoding="utf-8",
        )
        payload["source_manifests"] = [
            {
                "path": str(source_path),
                "sha256": hashlib.sha256(
                    source_path.read_bytes()
                ).hexdigest(),
            }
        ]
        payload.pop("analysis_selection")
        payload["strata"][0].pop("guard")
    path.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    return path


class ProposalTests(unittest.TestCase):
    def test_legacy_shared_envelope_manifest_remains_supported(self) -> None:
        manifest = {
            "schema": radiative_mode4.LEGACY_MANIFEST_SCHEMAS[0],
            "sigr_max": 0.005,
        }
        self.assertTrue(radiative_mode4._supported_manifest_schema(manifest))
        self.assertAlmostEqual(
            radiative_mode4._run_sigr_max(manifest, {}), 0.005
        )

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

    def test_fortran_single_precision_guard_boundary_is_reproduced(self) -> None:
        box = radiative_mode4.GuardBox(
            nonperiodic={
                "r_u": (0.7, 0.9),
                "r_ep": (0.2, 0.8),
                "u_gamma": (0.2, 1.0),
                "hadron_cosine_base": (0.2, 0.8),
            },
            phi_origin=0.5,
            phi_relative=(-0.2, 0.2),
        )
        coordinates = {
            "r_u": radiative_mode4._fortran_real32(0.7),
            "r_ep": 0.5,
            "u_gamma": 0.5,
            "hadron_cosine_base": 0.5,
            "hadron_phi_base": 0.5,
        }
        self.assertFalse(box.contains(coordinates))
        self.assertTrue(
            radiative_mode4._fortran_guard_contains(box, coordinates)
        )
        expected = radiative_mode4._proposal_density_ratio_for_membership(
            box, 0.9, True
        )
        self.assertNotAlmostEqual(
            radiative_mode4.proposal_density_ratio(coordinates, box, 0.9),
            expected,
        )

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
    def test_pilot_validation_pools_stress_runs_and_separates_geometry(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            calibration = _write_calibration_report(root)
            first = _write_pilot_artifacts(
                root,
                name="pilot_1",
                seed=1001,
                sigr_max=3.5,
                sigma=0.31,
                legacy_outside_duplicate=False,
            )
            second = _write_pilot_artifacts(
                root,
                name="pilot_2",
                seed=1003,
                sigr_max=4.0,
                sigma=0.29,
                legacy_outside_duplicate=True,
            )
            output = radiative_mode4.validate_pilots(
                argparse.Namespace(
                    calibration=calibration,
                    runs=[first, second],
                    output=root / "pilot_validation.json",
                    minimum_runs=2,
                    minimum_events=400,
                    maximum_duplicate_fraction=0.05,
                    maximum_guard_complement_fraction=0.02,
                    maximum_relative_cross_section_difference=0.10,
                    maximum_cross_section_z_score=3.0,
                    confidence=0.95,
                )
            )
            report = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(
                report["schema"],
                radiative_mode4.PILOT_VALIDATION_SCHEMA,
            )
            self.assertTrue(report["passed"])
            result = report["strata"][0]
            self.assertEqual(
                result["recommendation"],
                "ready_for_multi_stratum_pilot",
            )
            self.assertEqual(result["envelope_decision"], "hold")
            self.assertEqual(result["guard_decision"], "hold")
            self.assertEqual(result["pilot_runs"], 2)
            self.assertEqual(result["qualifying_conservative_runs"], 2)
            self.assertEqual(result["total_events"], 400)
            self.assertEqual(result["emitting_candidates"], 399)
            self.assertEqual(result["duplicate_events"], 1)
            self.assertAlmostEqual(
                result["qualifying_duplicate_event_fraction"], 0.0025
            )
            self.assertLess(
                result["duplicate_candidate_wilson_upper_fraction"],
                0.02,
            )
            self.assertEqual(result["guard_complement_events"], 2)
            self.assertAlmostEqual(
                result["guard_complement_event_fraction"], 0.005
            )
            self.assertEqual(
                result["event_classification"]["legacy_outside"]["events"],
                2,
            )
            self.assertEqual(
                result["event_classification"][
                    "guard_focused_outside"
                ]["events"],
                0,
            )
            self.assertEqual(
                result["guard_complement_faces"][0]["axis"], "r_u"
            )
            self.assertEqual(
                result["guard_complement_faces"][0]["face"], "upper"
            )
            self.assertAlmostEqual(
                result["guard_complement_faces"][0][
                    "maximum_excursion"
                ],
                0.0077,
            )
            maximum = result["maximum_integrands"]["guard_complement"]
            self.assertEqual(maximum["proposal_component_name"], "legacy")
            self.assertFalse(maximum["inside_geometric_guard"])
            self.assertAlmostEqual(maximum["ratio_to_run_sigr_max"], 1.8)
            self.assertIn(
                "qualifying_duplicate_event_fraction",
                output.with_suffix(".tsv").read_text(
                    encoding="utf-8"
                ).splitlines()[0],
            )
            legacy_calibration = _write_calibration_report(
                root, legacy_metadata=True
            )
            legacy_output = radiative_mode4.validate_pilots(
                argparse.Namespace(
                    calibration=legacy_calibration,
                    runs=[first, second],
                    output=root / "pilot_validation_legacy.json",
                    minimum_runs=2,
                    minimum_events=400,
                    maximum_duplicate_fraction=0.05,
                    maximum_guard_complement_fraction=0.02,
                    maximum_relative_cross_section_difference=0.10,
                    maximum_cross_section_z_score=3.0,
                    confidence=0.95,
                )
            )
            self.assertTrue(
                json.loads(legacy_output.read_text(encoding="utf-8"))[
                    "passed"
                ]
            )
            failed_output = radiative_mode4.validate_pilots(
                argparse.Namespace(
                    calibration=calibration,
                    runs=[first, second],
                    output=root / "pilot_validation_strict.json",
                    minimum_runs=2,
                    minimum_events=400,
                    maximum_duplicate_fraction=0.001,
                    maximum_guard_complement_fraction=0.02,
                    maximum_relative_cross_section_difference=0.10,
                    maximum_cross_section_z_score=3.0,
                    confidence=0.95,
                )
            )
            failed = json.loads(
                failed_output.read_text(encoding="utf-8")
            )
            self.assertFalse(failed["passed"])
            self.assertEqual(
                failed["strata"][0]["recommendation"],
                "increase_envelope_and_revalidate",
            )

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

    def test_prepare_uses_sparse_per_stratum_envelopes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config, recipes, legacy, report = _multistratum_fixtures(root)
            args = _prepare_args(root, config, recipes, legacy)
            args.sigr_max = None
            args.envelope_report = report
            args.flat_indices = [1, 0]
            args.allow_envelope_revision_mismatch = False
            args.envelope_revision_compatibility_rationale = None
            manifest_path = radiative_mode4.prepare(args)
            manifest = json.loads(
                manifest_path.read_text(encoding="utf-8")
            )
            self.assertEqual(
                manifest["envelope_mode"], "per_stratum_calibration"
            )
            self.assertIsNone(manifest["sigr_max"])
            self.assertEqual(
                manifest["stratum_selection"],
                {
                    "mode": "sparse_flat_indices",
                    "flat_indices": [0, 1],
                },
            )
            self.assertEqual(manifest["envelope_strata_used"], [
                "s00000",
                "s00001",
            ])
            self.assertEqual(
                hashlib.sha256(
                    (
                        manifest_path.parent / "envelope_calibration.json"
                    ).read_bytes()
                ).hexdigest(),
                manifest["envelope_calibration_sha256"],
            )
            records = {
                record["stratum_id"]: record for record in manifest["runs"]
            }
            self.assertAlmostEqual(records["s00000"]["sigr_max"], 0.004)
            self.assertAlmostEqual(records["s00001"]["sigr_max"], 0.006)
            for stratum_id, expected in (
                ("s00000", 0.004),
                ("s00001", 0.006),
            ):
                prepared = (
                    manifest_path.parent / records[stratum_id]["input_file"]
                ).read_text(encoding="utf-8")
                parsed = radiative_survey._records(prepared)
                self.assertAlmostEqual(float(parsed[17]), expected)
                self.assertEqual(
                    records[stratum_id]["envelope_pilot_readiness"],
                    "ready",
                )
            snapshot = manifest_path.parent / "envelope_calibration.json"
            snapshot.write_text("{}\n", encoding="utf-8")
            with self.assertRaisesRegex(
                radiative_mode4.Mode4Error,
                "envelope-calibration snapshot",
            ):
                radiative_mode4.run(
                    argparse.Namespace(
                        manifest=manifest_path,
                        flat_index=0,
                        replica_index=0,
                        executable=root / "missing-generator",
                        overwrite=False,
                    )
                )

    def test_envelope_report_requires_readiness_and_audited_revision(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config, recipes, legacy, report = _multistratum_fixtures(root)
            payload = json.loads(report.read_text(encoding="utf-8"))
            payload["generator_revision"] = "older-generator"
            payload["generator_revisions"] = ["older-generator"]
            report.write_text(
                json.dumps(payload, indent=2) + "\n", encoding="utf-8"
            )
            args = _prepare_args(root, config, recipes, legacy)
            args.sigr_max = None
            args.envelope_report = report
            args.flat_indices = [0]
            args.allow_envelope_revision_mismatch = False
            args.envelope_revision_compatibility_rationale = None
            with self.assertRaisesRegex(ValueError, "explicit audited override"):
                radiative_mode4.prepare(args)
            args.output = root / "audited_campaign"
            args.allow_envelope_revision_mismatch = True
            args.envelope_revision_compatibility_rationale = (
                "Wrapper-only revision; generator physics and proposal "
                "are unchanged."
            )
            manifest_path = radiative_mode4.prepare(args)
            manifest = json.loads(
                manifest_path.read_text(encoding="utf-8")
            )
            self.assertTrue(
                manifest["envelope_revision_compatibility_override"][
                    "enabled"
                ]
            )
            payload["strata"][1]["pilot_readiness"] = "not_ready"
            report.write_text(
                json.dumps(payload, indent=2) + "\n", encoding="utf-8"
            )
            args.output = root / "unready_campaign"
            args.flat_indices = [1]
            with self.assertRaisesRegex(ValueError, "not pilot-ready"):
                radiative_mode4.prepare(args)

    def test_sparse_selection_rejects_range_and_duplicates(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config, recipes, legacy, _report = _multistratum_fixtures(root)
            args = _prepare_args(root, config, recipes, legacy)
            args.flat_indices = [0]
            args.bin_stop = 1
            with self.assertRaisesRegex(ValueError, "cannot be combined"):
                radiative_mode4.prepare(args)
            args.output = root / "duplicate_selection"
            args.bin_stop = None
            args.flat_indices = [0, 0]
            with self.assertRaisesRegex(ValueError, "must be unique"):
                radiative_mode4.prepare(args)

    def test_flat_index_file_is_snapshotted_and_verified(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config, recipes, legacy, report = _multistratum_fixtures(root)
            selection = root / "production_ready_flat_indices.txt"
            selection.write_text("1\n# audited active mask\n0\n", encoding="utf-8")
            args = _prepare_args(root, config, recipes, legacy)
            args.sigr_max = None
            args.envelope_report = report
            args.flat_indices = None
            args.flat_index_file = selection
            args.allow_envelope_revision_mismatch = False
            args.envelope_revision_compatibility_rationale = None
            manifest_path = radiative_mode4.prepare(args)
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(
                manifest["stratum_selection"],
                {
                    "mode": "flat_index_file",
                    "flat_indices": [0, 1],
                    "source": str(selection.resolve()),
                    "source_sha256": hashlib.sha256(
                        selection.read_bytes()
                    ).hexdigest(),
                    "snapshot": "flat_index_selection.txt",
                },
            )
            snapshot = manifest_path.parent / "flat_index_selection.txt"
            self.assertEqual(
                snapshot.read_text(encoding="utf-8"),
                selection.read_text(encoding="utf-8"),
            )
            snapshot.write_text("0\n", encoding="utf-8")
            with self.assertRaisesRegex(
                radiative_mode4.Mode4Error, "selection snapshot"
            ):
                radiative_mode4.run(
                    argparse.Namespace(
                        manifest=manifest_path,
                        flat_index=0,
                        replica_index=0,
                        executable=root / "missing-generator",
                        overwrite=False,
                    )
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
            self.assertAlmostEqual(completed["sigr_max"], 0.005)
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
            self.assertAlmostEqual(stratum["sigr_max"], 0.005)
            self.assertAlmostEqual(
                stratum["pooled_event_weight_microbarn"]
                * stratum["total_events"],
                stratum["combined_sig_sum_microbarn"],
            )
            self.assertEqual(
                stratum["emitting_candidates"] + stratum["duplicate_events"],
                stratum["total_events"],
            )
            self.assertIn(
                "sigr_max",
                weights_path.with_suffix(".tsv")
                .read_text(encoding="utf-8")
                .splitlines()[0],
            )

    def test_calibration_subset_skips_untouched_run_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest_path = root / "manifest.json"
            records = []
            for flat_index in (0, 1):
                identifier = f"s{flat_index:05d}"
                records.append(
                    {
                        "flat_index": flat_index,
                        "stratum_id": identifier,
                        "replica_index": 0,
                        "seed": flat_index + 1,
                        "indices": {
                            "iq2": 0,
                            "ixb": 0,
                            "it": 0,
                            "iphi": flat_index,
                        },
                        "bounds": {"phi_deg": [0.0, 18.0]},
                        "guard": {"test_guard": True},
                        "output_stem": (
                            f"runs/{identifier}/{identifier}__g0000"
                        ),
                    }
                )
            manifest_path.write_text(
                json.dumps(
                    {
                        "schema": radiative_mode4.MANIFEST_SCHEMA,
                        "operation": "calibration",
                        "runs": records,
                    },
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
            manifest_hash = hashlib.sha256(
                manifest_path.read_bytes()
            ).hexdigest()
            selected_path = root / (
                records[1]["output_stem"] + ".json"
            )
            selected_path.parent.mkdir(parents=True)
            selected_path.write_text(
                json.dumps(
                    {
                        "schema": radiative_mode4.RUN_SCHEMA,
                        "source_manifest_sha256": manifest_hash,
                        "operation": "calibration",
                        "stratum_id": "s00001",
                        "flat_index": 1,
                        "replica_index": 0,
                    }
                )
                + "\n",
                encoding="utf-8",
            )
            output = root / "subset.json"

            def fake_finalize(
                args: argparse.Namespace,
                sources: list[tuple[Path, dict]],
                grouped: dict,
            ) -> Path:
                self.assertEqual(len(sources), 1)
                self.assertEqual(set(grouped), {"s00001"})
                output.write_text("{}\n", encoding="utf-8")
                return output

            with mock.patch.object(
                radiative_mode4,
                "_finalize_calibration",
                side_effect=fake_finalize,
            ):
                result = radiative_mode4.finalize(
                    argparse.Namespace(
                        manifests=[manifest_path],
                        output=output,
                        stratum_ids={"s00001"},
                    )
                )
            self.assertEqual(result, output)
            untouched_path = root / (
                records[0]["output_stem"] + ".json"
            )
            self.assertFalse(untouched_path.exists())

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
