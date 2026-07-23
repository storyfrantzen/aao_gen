#!/usr/bin/env python3
"""Unit tests for the Born bin-conditional preparation and validation helpers."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


MODULE_PATH = Path(__file__).with_name("bin_conditional.py")
SPEC = importlib.util.spec_from_file_location("bin_conditional", MODULE_PATH)
bin_conditional = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = bin_conditional
SPEC.loader.exec_module(bin_conditional)


class BinConditionalTests(unittest.TestCase):
    def test_legacy_flat_index_order(self) -> None:
        index = bin_conditional.legacy_flat_index
        dimensions = {"nq2": 2, "nt": 2, "nphi": 3}
        self.assertEqual(index(0, 0, 0, 0, **dimensions), 0)
        self.assertEqual(index(0, 0, 1, 0, **dimensions), 1)
        self.assertEqual(index(0, 0, 0, 1, **dimensions), 2)
        self.assertEqual(index(1, 0, 0, 0, **dimensions), 6)
        self.assertEqual(index(0, 1, 0, 0, **dimensions), 12)

    def test_enumeration_is_complete_and_sorted(self) -> None:
        config = self._config()
        strata = bin_conditional.enumerate_strata(config)
        self.assertEqual(len(strata), 16)
        self.assertEqual([item.flat_index for item in strata], list(range(16)))
        self.assertEqual(
            (strata[-1].iq2, strata[-1].ixb, strata[-1].it, strata[-1].iphi),
            (1, 1, 1, 1),
        )

    def test_prepare_snapshots_config_and_writes_mode4_input(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "analysis.json"
            config_path.write_text(json.dumps(self._config()), encoding="utf-8")
            output = root / "prepared"
            args = argparse.Namespace(
                config=config_path,
                output=output,
                events_per_bin=25,
                physics_model=1,
                flag_ehel=1,
                npart=3,
                epirea=1,
                fmcall=3.0,
                boso=0,
                seed_base=1000,
                bin_start=0,
                bin_stop=None,
                support_samples=8,
                include_unsupported=True,
                condition_phase_space=False,
                overwrite=False,
            )
            manifest_path = bin_conditional.prepare(args)
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["schema"], bin_conditional.MANIFEST_SCHEMA)
            self.assertEqual(manifest["sampling_mode"], 4)
            self.assertFalse(manifest["phase_space_conditioned"])
            self.assertEqual(manifest["prepared_strata"], 16)
            self.assertTrue((output / "analysis_config.json").is_file())
            self.assertEqual(len(manifest["analysis_config_sha256"]), 64)

            first = manifest["strata"][0]
            lines = (output / first["input_file"]).read_text(encoding="utf-8").splitlines()
            self.assertEqual(lines[6], "0.001 6.535")
            self.assertEqual(lines[11], "4")
            self.assertEqual(lines[15], "0 2 0.8")
            self.assertEqual(lines[16], "0 0 0 0 0")

            with self.assertRaises(FileExistsError):
                bin_conditional.prepare(args)

    def test_kinematics_validation_enforces_selection(self) -> None:
        record = {
            "events_requested": 2,
            "phase_space_conditioned": True,
            "bounds": {
                "Q2": [1.0, 1.5],
                "xB": [0.15, 0.25],
                "minus_t": [0.09, 0.3],
                "phi_deg": [0.0, 18.0],
            },
        }
        rows = np.array(
            [
                [1, 1.1, 0.18, 0.12, 4.0, 1.0, 1.0],
                [2, 1.3, 0.22, 0.20, 12.0, 1.0, 1.0],
            ]
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "events.kin"
            np.savetxt(path, rows)
            bin_conditional._validate_kinematics(
                path,
                record,
                expected_events=2,
                beam_energy=6.535,
                phase_space={"W_min": 2.0, "y_max": 0.8},
            )
            rows[1, 4] = 20.0
            np.savetxt(path, rows)
            with self.assertRaisesRegex(RuntimeError, "phi_deg"):
                bin_conditional._validate_kinematics(
                    path,
                    record,
                    expected_events=2,
                    beam_energy=6.535,
                    phase_space={"W_min": 2.0, "y_max": 0.8},
                )
            rows[1, 1:5] = [1.49, 0.15, 0.20, 12.0]
            np.savetxt(path, rows)
            with self.assertRaisesRegex(RuntimeError, "y cut"):
                bin_conditional._validate_kinematics(
                    path,
                    record,
                    expected_events=2,
                    beam_energy=6.535,
                    phase_space={"W_min": 2.0, "y_max": 0.8},
                )
            record["phase_space_conditioned"] = False
            bin_conditional._validate_kinematics(
                path,
                record,
                expected_events=2,
                beam_energy=6.535,
                phase_space={"W_min": 2.0, "y_max": 0.8},
            )

    def test_multiplicity_correction_and_event_overshoot_are_valid(self) -> None:
        record = {
            "events_requested": 25,
            "phase_space_conditioned": False,
            "flat_index": 0,
            "indices": {"iq2": 0, "ixb": 0, "it": 0, "iphi": 0},
        }
        norm = {
            "sampling_mode": "4",
            "conditional_unweighted": "1",
            "phase_space_conditioned": "0",
            "stratum_flat_index": "0",
            "stratum_iq2": "0",
            "stratum_ixb": "0",
            "stratum_it": "0",
            "stratum_iphi": "0",
            "events": "27",
            "mcall_max": "3",
        }
        bin_conditional._validate_norm_record(norm, record)

    @staticmethod
    def _config() -> dict:
        return {
            "beam_energy": 6.535,
            "phase_space": {"Q2_min": 1.0, "W_min": 2.0, "y_max": 0.8},
            "binning": {
                "Q2": [1.0, 1.5, 2.0],
                "xB": [0.1, 0.2, 0.3],
                "minus_t": [0.09, 0.2, 0.4],
                "phi_deg": [0.0, 18.0, 36.0],
            },
        }


if __name__ == "__main__":
    unittest.main()
