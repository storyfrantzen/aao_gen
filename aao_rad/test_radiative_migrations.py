#!/usr/bin/env python3
"""Unit tests for the milestone-2b/2c migration diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
import tempfile
import unittest
from pathlib import Path

import radiative_migrations
from test_radiative_guards import _config, _row, _write_survey


def _migration_row(
    *,
    replica: int,
    trial: int,
    weight: float,
    q2_hard: float,
    channel: int = 1,
) -> dict[str, str | int | float]:
    row = _row(
        replica=replica,
        trial=trial,
        weight=weight,
        r_u=0.2,
        hadron_phi_base=0.1,
    )
    row.update(
        {
            "intreg": channel,
            "q2_hard": q2_hard,
            "xb_hard": 0.2,
            "minus_t_hard": 0.2,
            "phi_cm_deg": 36.0,
        }
    )
    return row


class ParentIndexTests(unittest.TestCase):
    def test_parent_ids_round_trip_and_keep_feed_in(self) -> None:
        config = _config()
        inside = _migration_row(
            replica=0, trial=1, weight=1.0, q2_hard=1.5
        )
        parent = radiative_migrations.parent_from_row(inside, config)
        self.assertEqual(
            parent, radiative_migrations.ParentIndex(0, 0, 0, 0)
        )
        identifier = radiative_migrations.parent_identifier(parent, config)
        self.assertEqual(
            radiative_migrations.parent_from_identifier(identifier, config),
            parent,
        )

        feed_in = _migration_row(
            replica=0, trial=2, weight=1.0, q2_hard=0.5
        )
        outside_parent = radiative_migrations.parent_from_row(
            feed_in, config
        )
        self.assertEqual(outside_parent.iq2, -1)
        metadata = radiative_migrations.parent_metadata(
            outside_parent, config
        )
        self.assertTrue(metadata["has_underflow_or_overflow"])
        self.assertEqual(
            metadata["coordinates"]["Q2"]["region"], "underflow"
        )

    def test_parent_dilation_wraps_phi_but_not_channel(self) -> None:
        config = _config()
        parent = radiative_migrations.ParentIndex(0, 0, 0, 0)
        parent_id = radiative_migrations.parent_identifier(parent, config)
        expanded = radiative_migrations.dilate_components(
            {(parent_id, 1)}, config, 1
        )
        wrapped = radiative_migrations.ParentIndex(0, 0, 0, 1)
        self.assertIn(
            (
                radiative_migrations.parent_identifier(wrapped, config),
                1,
            ),
            expanded,
        )
        self.assertFalse(any(channel == 2 for _, channel in expanded))


class MigrationWorkflowTests(unittest.TestCase):
    def test_channel_representation_comparison_uses_frozen_footprints(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "config.json"
            config_path.write_text(json.dumps(_config()), encoding="utf-8")
            training = root / "training"
            _write_survey(
                training,
                replica=0,
                rows=[
                    _migration_row(
                        replica=0,
                        trial=1,
                        weight=8.0,
                        q2_hard=1.5,
                        channel=1,
                    ),
                    _migration_row(
                        replica=0,
                        trial=2,
                        weight=2.0,
                        q2_hard=2.5,
                        channel=2,
                    ),
                ],
            )
            validation = root / "validation"
            _write_survey(
                validation,
                replica=1,
                rows=[
                    _migration_row(
                        replica=1,
                        trial=1,
                        weight=7.0,
                        q2_hard=1.5,
                        channel=3,
                    ),
                    _migration_row(
                        replica=1,
                        trial=2,
                        weight=3.0,
                        q2_hard=2.5,
                        channel=2,
                    ),
                ],
            )
            output = root / "representation_study"
            result = radiative_migrations.compare_representations(
                argparse.Namespace(
                    config=config_path,
                    training_survey=[training],
                    validation_survey=[validation],
                    output=output,
                    target_parent_fraction=0.7,
                    parent_dilation=0,
                    iteration=0,
                    minimum_training_rows=1,
                    minimum_training_ess=0.0,
                    minimum_parent_coverage=0.65,
                    apply_y_max=False,
                    generator_revision="test-revision",
                )
            )
            self.assertTrue(result["passed"])
            comparison_path = output / "representation_comparison.json"
            comparison = json.loads(
                comparison_path.read_text(encoding="utf-8")
            )
            representations = comparison["representations"]
            self.assertAlmostEqual(
                representations["six_channel"]["validation"][
                    "cross_section_weighted_selected_parent_fraction"
                ],
                0.0,
            )
            for identifier in (
                "four_group",
                "soft_resolved",
                "channel_marginalized",
            ):
                self.assertAlmostEqual(
                    representations[identifier]["validation"][
                        "cross_section_weighted_selected_parent_fraction"
                    ],
                    0.7,
                )
            self.assertAlmostEqual(
                representations["four_group"][
                    "validation_by_native_intreg"
                ]["intreg_3"][
                    "cross_section_weighted_selected_parent_fraction"
                ],
                1.0,
            )
            self.assertEqual(
                comparison[
                    "ranking_by_heldout_coverage_then_compactness"
                ][0],
                "four_group",
            )
            footprints = json.loads(
                (
                    output / "representation_footprints.json"
                ).read_text(encoding="utf-8")
            )
            four_group_seeds = footprints["strata"]["s00000"][
                "representations"
            ]["four_group"]["seed_parent_groups"]
            self.assertEqual(
                four_group_seeds[0]["channel_group"], "peak_a"
            )
            self.assertTrue(
                (
                    output / "representation_comparison.json.sha256"
                ).is_file()
            )
            with (
                output / "representation_strata.csv"
            ).open(newline="", encoding="utf-8") as source:
                rows = list(csv.DictReader(source))
            self.assertEqual(len(rows), 16)

            plot_output = root / "representation_plots"
            plotted = radiative_migrations.plot_representations(
                argparse.Namespace(
                    comparison=comparison_path,
                    output=plot_output,
                )
            )
            self.assertTrue(plotted["passed"])
            self.assertGreater(
                (
                    plot_output / "representation_comparison.pdf"
                ).stat().st_size,
                0,
            )

    def test_weighted_parent_learning_validation_summary_and_plots(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "config.json"
            config_path.write_text(json.dumps(_config()), encoding="utf-8")

            training = root / "training"
            _write_survey(
                training,
                replica=0,
                rows=[
                    _migration_row(
                        replica=0, trial=1, weight=8.0, q2_hard=1.5
                    ),
                    _migration_row(
                        replica=0, trial=2, weight=2.0, q2_hard=2.5
                    ),
                    _migration_row(
                        replica=0, trial=3, weight=1.0, q2_hard=0.5
                    ),
                ],
            )
            output = root / "migration"
            learned = radiative_migrations.learn(
                argparse.Namespace(
                    config=config_path,
                    survey=[training],
                    output=output,
                    target_parent_fraction=0.7,
                    parent_dilation=0,
                    iteration=0,
                    minimum_training_rows=1,
                    minimum_training_ess=0.0,
                    apply_y_max=False,
                    generator_revision="test-revision",
                )
            )
            self.assertTrue(learned["passed"])
            manifest_path = output / "migration_manifest.json"
            manifest = json.loads(
                manifest_path.read_text(encoding="utf-8")
            )
            record = manifest["strata"]["s00000"]
            self.assertAlmostEqual(
                record["training_total"]["cross_section_microbarn"],
                11.0 / 100.0,
            )
            self.assertAlmostEqual(
                record["estimated_selected_parent_fraction"], 8.0 / 11.0
            )
            self.assertEqual(
                len(record["parent_footprint"]["seed_component_ids"]), 1
            )
            self.assertEqual(
                record["one_more_dilation_fraction"], 1.0
            )
            self.assertTrue(
                (output / "training_migrations.csv.sha256").is_file()
            )
            with (output / "training_migrations.csv").open(
                newline="", encoding="utf-8"
            ) as source:
                migration_rows = list(csv.DictReader(source))
            self.assertIn(
                "underflow",
                {row["hard_q2_region"] for row in migration_rows},
            )

            validation = root / "validation"
            _write_survey(
                validation,
                replica=1,
                rows=[
                    _migration_row(
                        replica=1, trial=1, weight=7.0, q2_hard=1.5
                    ),
                    _migration_row(
                        replica=1, trial=2, weight=2.0, q2_hard=2.5
                    ),
                    _migration_row(
                        replica=1, trial=3, weight=1.0, q2_hard=0.5
                    ),
                ],
            )
            validation_output = root / "migration_validation"
            checked = radiative_migrations.validate(
                argparse.Namespace(
                    manifest=manifest_path,
                    survey=[validation],
                    output=validation_output,
                    minimum_parent_coverage=0.65,
                )
            )
            self.assertTrue(checked["passed"])
            validation_path = (
                validation_output / "migration_validation.json"
            )
            report = json.loads(
                validation_path.read_text(encoding="utf-8")
            )
            holdout = report["strata"]["s00000"]
            self.assertAlmostEqual(
                report["validation"]["inside_analysis_partition"][
                    "cross_section_microbarn"
                ],
                10.0 / 100.0,
            )
            self.assertAlmostEqual(
                holdout["holdout_selected_parent_fraction"], 0.7
            )
            self.assertAlmostEqual(
                holdout["holdout_one_more_dilation_fraction"], 1.0
            )

            summary = radiative_migrations.summarize(
                argparse.Namespace(
                    manifest=manifest_path,
                    validation=validation_path,
                    limit=1,
                )
            )
            weighted = summary["weighted_parent_coverage"]
            self.assertAlmostEqual(
                weighted["validation_inside_analysis_partition"][
                    "cross_section_weighted_selected_parent_fraction"
                ],
                0.7,
            )
            self.assertEqual(
                summary["worst_holdout_parent_fractions"][0][
                    "stratum_id"
                ],
                "s00000",
            )

            plot_output = root / "plots"
            plotted = radiative_migrations.plot(
                argparse.Namespace(
                    manifest=manifest_path,
                    validation=validation_path,
                    output=plot_output,
                    top_strata=1,
                )
            )
            self.assertTrue(plotted["passed"])
            self.assertGreater(
                (plot_output / "migration_diagnostics.pdf").stat().st_size,
                0,
            )
            self.assertTrue(
                (plot_output / "migration_diagnostics.pdf.sha256").is_file()
            )

            report["manifest_sha256"] = "wrong-hash"
            mismatched = (
                validation_output / "mismatched_validation.json"
            )
            mismatched.write_text(json.dumps(report), encoding="utf-8")
            with self.assertRaisesRegex(
                radiative_migrations.MigrationError,
                "does not reference the supplied manifest",
            ):
                radiative_migrations.summarize(
                    argparse.Namespace(
                        manifest=manifest_path,
                        validation=mismatched,
                        limit=1,
                    )
                )

    def test_training_replica_cannot_be_reused_for_validation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "config.json"
            config_path.write_text(json.dumps(_config()), encoding="utf-8")
            survey = root / "survey"
            _write_survey(
                survey,
                replica=0,
                rows=[
                    _migration_row(
                        replica=0, trial=1, weight=1.0, q2_hard=1.5
                    )
                ],
            )
            output = root / "migration"
            radiative_migrations.learn(
                argparse.Namespace(
                    config=config_path,
                    survey=[survey],
                    output=output,
                    target_parent_fraction=0.9,
                    parent_dilation=0,
                    iteration=0,
                    minimum_training_rows=1,
                    minimum_training_ess=0.0,
                    apply_y_max=False,
                    generator_revision="test",
                )
            )
            with self.assertRaisesRegex(
                radiative_migrations.MigrationError, "overlap"
            ):
                radiative_migrations.validate(
                    argparse.Namespace(
                        manifest=output / "migration_manifest.json",
                        survey=[survey],
                        output=root / "invalid_validation",
                        minimum_parent_coverage=0.5,
                    )
                )


if __name__ == "__main__":
    unittest.main()
