#!/usr/bin/env python3
"""Unit tests for the milestone-2d radiative guard geometries."""

from __future__ import annotations

import argparse
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import radiative_guard_geometries as geometries
import radiative_guards as guards
import radiative_migrations as migrations
from test_radiative_guards import _config


def _moment(count: int, value: float) -> guards.Moment:
    result = guards.Moment()
    for _ in range(count):
        result.add(value)
    return result


def _campaign(
    *,
    replica: int,
    strata: dict[str, guards.Moment],
    migrations_by_stratum: dict[
        str, dict[tuple[str, int], guards.Moment]
    ],
) -> migrations.MigrationCampaign:
    global_observed = migrations._combined(strata.values())
    parent_totals: dict[tuple[str, int], guards.Moment] = {}
    channels: dict[int, guards.Moment] = {}
    for migration_moments in migrations_by_stratum.values():
        for component, moment in migration_moments.items():
            migrations._add_moment(
                parent_totals.setdefault(component, guards.Moment()),
                moment,
            )
            migrations._add_moment(
                channels.setdefault(component[1], guards.Moment()),
                moment,
            )
    rows = sum(moment.count for moment in strata.values())
    return migrations.MigrationCampaign(
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
        migrations=migrations_by_stratum,
        parent_totals=parent_totals,
        relationships={},
        channels=channels,
        norm_reference={"ebeam": "6.535"},
        input_signature="test-input-signature",
        legacy_input_settings={"physics_model": 5},
        channel_probabilities={
            str(channel): 1.0 / 6.0
            for channel in range(1, 7)
        },
    )


class GeometryPrimitiveTests(unittest.TestCase):
    def test_observed_neighbors_wrap_phi(self) -> None:
        config = _config()
        catalog = guards.enumerate_strata(config)
        lookup = {
            geometries._stratum_key(stratum): stratum
            for stratum in catalog
        }
        target = next(
            stratum
            for stratum in catalog
            if stratum.identifier == "s00000"
        )
        neighbors = geometries._observed_neighbor_ids(
            target,
            lookup,
            geometries._target_shape(config),
            radius=1,
        )
        self.assertIn("s00001", neighbors)
        self.assertIn("s00002", neighbors)

    def test_qx_dilation_does_not_expand_t_or_phi(self) -> None:
        config = _config()
        grid = migrations.ParentGrid.from_config(config)
        seed = migrations.ParentIndex(0, 0, 0, 0)
        seed_id = migrations.parent_identifier(seed, grid)
        expanded = geometries._dilate_parent_ids(
            {seed_id}, grid, ("Q2", "xB"), radius=1
        )
        decoded = {
            migrations.parent_from_identifier(parent_id, grid)
            for parent_id in expanded
        }
        self.assertTrue(all(parent.it == 0 for parent in decoded))
        self.assertTrue(all(parent.iphi == 0 for parent in decoded))
        self.assertIn(
            migrations.ParentIndex(-1, 0, 0, 0), decoded
        )
        self.assertIn(
            migrations.ParentIndex(1, 0, 0, 0), decoded
        )


class GeometryWorkflowTests(unittest.TestCase):
    def test_hierarchical_fallback_covers_empty_neighbor_stratum(
        self,
    ) -> None:
        config = _config()
        grid = migrations.ParentGrid.from_config(config)
        catalog = {
            stratum.identifier: stratum
            for stratum in guards.enumerate_strata(config)
        }

        def same_parent(stratum_id: str) -> str:
            target = catalog[stratum_id]
            return migrations.parent_identifier(
                migrations.ParentIndex(
                    target.iq2,
                    target.ixb,
                    target.it,
                    target.iphi,
                ),
                grid,
            )

        training = _campaign(
            replica=0,
            strata={
                "s00000": _moment(10, 1.0),
                "s00002": _moment(1, 1.0),
            },
            migrations_by_stratum={
                "s00000": {
                    (same_parent("s00000"), 1): _moment(10, 1.0)
                },
                "s00002": {
                    (same_parent("s00002"), 1): _moment(1, 1.0)
                },
            },
        )
        validation = _campaign(
            replica=1,
            strata={
                "s00000": _moment(5, 1.0),
                "s00001": _moment(5, 1.0),
                "s00002": _moment(1, 1.0),
            },
            migrations_by_stratum={
                "s00000": {
                    (same_parent("s00000"), 1): _moment(5, 1.0)
                },
                "s00001": {
                    (same_parent("s00001"), 6): _moment(5, 1.0)
                },
                "s00002": {
                    (same_parent("s00002"), 1): _moment(1, 1.0)
                },
            },
        )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "config.json"
            config_path.write_text(
                json.dumps(config), encoding="utf-8"
            )
            output = root / "geometry"
            with (
                mock.patch.object(
                    geometries.migrations,
                    "aggregate_surveys",
                    side_effect=[training, validation],
                ),
                mock.patch.object(
                    geometries.guards,
                    "_current_revision",
                    return_value="test-learner-revision",
                ),
            ):
                result = geometries.compare(
                    argparse.Namespace(
                        config=config_path,
                        training_survey=[root / "training"],
                        validation_survey=[root / "validation"],
                        output=output,
                        target_parent_fraction=0.995,
                        neighbor_radius=1,
                        iteration=0,
                        minimum_training_rows=10,
                        minimum_training_ess=5.0,
                        minimum_parent_coverage=0.98,
                        apply_y_max=False,
                        generator_revision="test-generator-revision",
                    )
                )

            self.assertTrue(result["passed"])
            comparison_path = (
                output / "guard_geometry_comparison.json"
            )
            comparison = json.loads(
                comparison_path.read_text(encoding="utf-8")
            )
            candidates = comparison["candidates"]
            self.assertAlmostEqual(
                candidates["independent_frozen"]["validation"][
                    "cross_section_weighted_core_coverage"
                ],
                6.0 / 11.0,
            )
            self.assertAlmostEqual(
                candidates["independent_all_axis_dilation"][
                    "validation"
                ]["cross_section_weighted_core_coverage"],
                6.0 / 11.0,
            )
            for identifier in (
                "hierarchical_neighbor_offsets",
                "hierarchical_neighbor_offsets_qx_dilation",
            ):
                self.assertAlmostEqual(
                    candidates[identifier]["validation"][
                        "cross_section_weighted_core_coverage"
                    ],
                    1.0,
                )
                self.assertEqual(
                    candidates[identifier]["coverage_summary"][
                        "failed_strata"
                    ],
                    0,
                )
            self.assertTrue(
                comparison["proposal_support_requirement"][
                    "future_full_support_tail_required"
                ]
            )

            recipe_path = output / "guard_geometry_recipes.json"
            recipes = json.loads(
                recipe_path.read_text(encoding="utf-8")
            )
            empty_record = recipes["strata"]["s00001"]
            self.assertEqual(
                empty_record["training_status"],
                "no_training_contribution",
            )
            self.assertEqual(empty_record["own_seed_parent_ids"], [])
            self.assertIn(
                "s00000",
                empty_record["fallback_neighbor_strata"],
            )
            self.assertGreater(
                empty_record["candidate_core_counts"][
                    "hierarchical_neighbor_offsets"
                ]["hard_cells"],
                0,
            )
            self.assertEqual(
                recipes["global_offset_template"],
                [
                    {
                        "delta_minus_t_index": 0,
                        "delta_phi_index": 0,
                        "delta_q2_index": 0,
                        "delta_xb_index": 0,
                    }
                ],
            )

            row_path = output / "guard_geometry_strata.csv"
            with row_path.open(
                newline="", encoding="utf-8"
            ) as source:
                rows = list(csv.DictReader(source))
            self.assertEqual(len(rows), 16)
            self.assertTrue(
                (
                    output / "guard_geometry_comparison.json.sha256"
                ).is_file()
            )
            self.assertTrue(
                (
                    output / "guard_geometry_recipes.json.sha256"
                ).is_file()
            )
            self.assertTrue(
                (
                    output / "guard_geometry_strata.csv.sha256"
                ).is_file()
            )

            plot_output = root / "geometry_plots"
            plotted = geometries.plot(
                argparse.Namespace(
                    comparison=comparison_path,
                    output=plot_output,
                )
            )
            self.assertTrue(plotted["passed"])
            self.assertGreater(
                (
                    plot_output / "guard_geometry_comparison.pdf"
                ).stat().st_size,
                0,
            )


if __name__ == "__main__":
    unittest.main()
