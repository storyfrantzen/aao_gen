#!/usr/bin/env python3
"""Unit tests for continuous native-coordinate radiative guards."""

from __future__ import annotations

import argparse
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import radiative_continuous_guards as continuous
import radiative_guards as guards
from test_radiative_guards import _config, _row, _write_survey


def _point(
    r_u: float,
    *,
    weight: float = 1.0,
    r_ep: float = 0.25,
    u_gamma: float = 0.5,
    cosine: float = 0.4,
    phi: float = 0.1,
    channel: int = 1,
) -> continuous.Point:
    return continuous.Point(
        coordinates=(r_u, r_ep, u_gamma, cosine, phi),
        weight=weight,
        channel=channel,
    )


def _moment(points: list[continuous.Point]) -> guards.Moment:
    result = guards.Moment()
    for point in points:
        result.add(point.weight)
    return result


def _campaign(
    replica: int,
    points: dict[str, list[continuous.Point]],
) -> continuous.ContinuousCampaign:
    strata = {
        identifier: _moment(values)
        for identifier, values in points.items()
    }
    global_observed = guards._combined_moment(strata.values())
    rows = sum(len(values) for values in points.values())
    base = guards.Campaign(
        proposals=100,
        rows=rows,
        final_valid_rows=rows,
        replicas=[
            {
                "directory": f"replica-{replica}",
                "replica": replica,
                "seed": 371001 + 2 * replica,
                "proposals": 100,
                "rows": rows,
                "final_valid_rows": rows,
                "observed_cross_section_microbarn": (
                    global_observed.total / 100.0
                ),
                "observed_cross_section_sem_microbarn": guards._sem(
                    global_observed, 100
                ),
            }
        ],
        global_observed=global_observed,
        outside={},
        strata=strata,
        cells={},
        channels={},
        norm_reference={
            "ebeam": "6.535",
            "q2_min": "1.0",
            "q2_max": "3.0",
            "ep_min": "0.2",
            "ep_max_effective": "6.0",
            "delta": "0.005",
            "epirea": "1",
            "th_opt": "5",
            "res_opt": "0",
            "survey_phase_volume": "1.0",
        },
        input_signature="test-input-signature",
        legacy_input_settings={"physics_model": 5},
        channel_probabilities={
            str(channel): 1.0 / 6.0 for channel in range(1, 7)
        },
    )
    return continuous.ContinuousCampaign(base=base, points=points)


class ContinuousBoxTests(unittest.TestCase):
    def test_joint_fit_excludes_a_low_weight_multiaxis_outlier(self) -> None:
        points = [
            _point(
                0.2,
                weight=49.0,
                r_ep=0.2,
                u_gamma=0.3,
                cosine=0.4,
                phi=0.1,
            ),
            _point(
                0.3,
                weight=49.0,
                r_ep=0.3,
                u_gamma=0.4,
                cosine=0.5,
                phi=0.2,
            ),
            _point(
                0.9,
                weight=2.0,
                r_ep=0.9,
                u_gamma=0.9,
                cosine=0.9,
                phi=0.6,
            ),
        ]
        box = continuous._fit_joint_box(
            points, target_fraction=0.98, minimum_axis_width=1.0e-4
        )
        self.assertTrue(box.contains(points[0].coordinates))
        self.assertTrue(box.contains(points[1].coordinates))
        self.assertFalse(box.contains(points[2].coordinates))
        self.assertLess(box.upper[0], 0.9)
        self.assertGreaterEqual(box.fit_iterations, 1)

    def test_periodic_fit_wraps_across_zero(self) -> None:
        points = [
            _point(0.3, weight=49.0, phi=0.98),
            _point(0.3, weight=49.0, phi=0.02),
            _point(0.3, weight=2.0, phi=0.5),
        ]
        box = continuous._fit_joint_box(
            points, target_fraction=0.98, minimum_axis_width=1.0e-4
        )
        self.assertTrue(box.contains(points[0].coordinates))
        self.assertTrue(box.contains(points[1].coordinates))
        self.assertFalse(box.contains(points[2].coordinates))
        phi_axis = box.metadata()["axes"][continuous.PHI_INDEX]
        self.assertTrue(phi_axis["wraps"])
        self.assertLess(phi_axis["width"], 0.1)

    def test_padding_is_continuous_and_clipped_to_unit_domain(self) -> None:
        box = continuous.Box(
            lower=(0.2, 0.0, 0.2, 0.2, -0.1),
            upper=(0.4, 0.3, 0.4, 0.4, 0.1),
            phi_origin=0.0,
        )
        padded = box.padded(0.01)
        self.assertAlmostEqual(padded.lower[0], 0.19)
        self.assertAlmostEqual(padded.upper[0], 0.41)
        self.assertEqual(padded.lower[1], 0.0)
        self.assertAlmostEqual(padded.upper[1], 0.31)
        self.assertAlmostEqual(padded.lower[continuous.PHI_INDEX], -0.11)
        self.assertAlmostEqual(padded.upper[continuous.PHI_INDEX], 0.11)


class ContinuousWorkflowTests(unittest.TestCase):
    def test_load_campaign_uses_default_partition(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            survey = root / "survey"
            _write_survey(
                survey,
                replica=0,
                rows=[_row(replica=0, trial=1, weight=1.0, r_u=0.25)],
            )

            campaign = continuous._load_campaign(
                [survey],
                _config(),
                apply_y_max=False,
            )

            self.assertEqual(campaign.proposals, 100)
            self.assertEqual(sum(len(points) for points in campaign.points.values()), 1)

    def test_empty_stratum_borrows_neighbor_and_writes_artifacts(self) -> None:
        training_points = {
            "s00000": [
                _point(
                    0.20 + 0.01 * index,
                    r_ep=0.25 + 0.005 * index,
                    phi=0.98 if index % 2 == 0 else 0.02,
                    channel=1 + index % 6,
                )
                for index in range(10)
            ]
        }
        validation_points = {
            "s00000": [
                _point(0.22, r_ep=0.26, phi=0.99, weight=5.0),
                _point(0.27, r_ep=0.28, phi=0.01, weight=5.0),
            ],
            "s00001": [
                _point(0.23, r_ep=0.27, phi=0.99, weight=4.0),
                _point(0.28, r_ep=0.29, phi=0.01, weight=4.0),
            ],
        }
        training = _campaign(0, training_points)
        validation = _campaign(1, validation_points)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "config.json"
            config_path.write_text(
                json.dumps(_config()), encoding="utf-8"
            )
            output = root / "continuous"
            with (
                mock.patch.object(
                    continuous,
                    "_load_campaign",
                    side_effect=[training, validation],
                ),
                mock.patch.object(
                    continuous.guards,
                    "_current_revision",
                    return_value="test-learner-revision",
                ),
            ):
                result = continuous.compare(
                    argparse.Namespace(
                        config=config_path,
                        training_survey=[root / "training"],
                        validation_survey=[root / "validation"],
                        output=output,
                        target_core_fraction=0.995,
                        padding=[0.0, 0.01],
                        neighbor_radius=1,
                        regularization_ess=5.0,
                        minimum_axis_width=1.0e-4,
                        minimum_training_rows=10,
                        minimum_training_ess=5.0,
                        minimum_validation_coverage=0.98,
                        iteration=0,
                        apply_y_max=False,
                        generator_revision="test-generator-revision",
                    )
                )

            self.assertTrue(result["passed"])
            comparison_path = output / "continuous_guard_comparison.json"
            comparison = json.loads(
                comparison_path.read_text(encoding="utf-8")
            )
            self.assertFalse(comparison["production_ready"])
            self.assertTrue(
                comparison["proposal_support_requirement"][
                    "future_full_support_tail_required"
                ]
            )
            recipes = json.loads(
                (
                    output / "continuous_guard_recipes.json"
                ).read_text(encoding="utf-8")
            )
            empty = recipes["strata"]["s00001"]
            self.assertEqual(
                empty["training_status"], "no_training_contribution"
            )
            self.assertEqual(empty["fit_source"], "neighbor_only")
            self.assertIn("s00000", empty["neighbor_strata"])

            with (
                output / "continuous_guard_strata.csv"
            ).open(encoding="utf-8", newline="") as source:
                rows = list(csv.DictReader(source))
            self.assertEqual(len(rows), 8)
            for filename in (
                "continuous_guard_comparison.json",
                "continuous_guard_recipes.json",
                "continuous_guard_strata.csv",
            ):
                self.assertTrue((output / f"{filename}.sha256").is_file())

            plot_output = root / "plots"
            plotted = continuous.plot(
                argparse.Namespace(
                    comparison=comparison_path,
                    output=plot_output,
                )
            )
            self.assertTrue(plotted["passed"])
            self.assertGreater(
                (
                    plot_output / "continuous_guard_comparison.pdf"
                ).stat().st_size,
                0,
            )


if __name__ == "__main__":
    unittest.main()
