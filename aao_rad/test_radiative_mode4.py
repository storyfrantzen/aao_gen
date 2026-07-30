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


def _raw_axis(name: str, *, periodic: bool = False) -> dict:
    if not periodic:
        return {
            "name": name,
            "periodic": False,
            "lower": 0.0,
            "upper": 1.0,
            "width": 1.0,
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
                        _raw_axis("r_u"),
                        _raw_axis("r_ep"),
                        _raw_axis("u_gamma"),
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
        self.assertEqual(box.phi_relative, (-0.5, 0.5))
        self.assertAlmostEqual(box.volume, 1.0)


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
                record["guard"]["normalized_volume"], 1.0
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
            self.assertEqual(completed["events"], 2)
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
            self.assertEqual(len(lund_lines), 10)
            self.assertTrue(all(
                lund_lines[offset].split()[0] == "4"
                for offset in (0, 5)
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
            self.assertEqual(stratum["total_events"], 2)
            self.assertAlmostEqual(
                stratum["pooled_event_weight_microbarn"]
                * stratum["total_events"],
                stratum["combined_sig_sum_microbarn"],
            )


if __name__ == "__main__":
    unittest.main()
