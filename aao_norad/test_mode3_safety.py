#!/usr/bin/env python3
"""Regression tests for fail-safe global Born mode-3 multiplicity."""

from __future__ import annotations

import re
import math
import subprocess
import tempfile
import unittest
from pathlib import Path


def _input(*, events: int, fmcall: float, seed: int) -> str:
    return "\n".join(
        [
            "1",                 # AO physics model
            "0",                 # unpolarized beam
            "3",                 # e, p, pi0 before decay
            "1",                 # pi0 channel
            "6.535",             # beam energy
            "1.0 3.0",           # Q2
            "0.2 6.535",         # scattered-electron energy
            str(events),
            f"{fmcall:.17g}",
            "0",                 # no BOS output
            str(-abs(seed)),
            "3",                 # global inverse-Q2 mode
            "0.15 0.50",         # xB
            "0.09 1.00",         # -t
            "0.0 360.0",         # phi
            "0 0.0 1.0",         # do not condition on W/y
            "",
        ]
    )


def _field(text: str, name: str) -> str:
    match = re.search(rf"(?m)^\s*{re.escape(name)}\s*=\s*([^\s]+)", text)
    if match is None:
        raise AssertionError(f"normalization sidecar lacks {name}")
    return match.group(1)


def _production_input(region: str, seed: int) -> str:
    """Exact nominal-W2 production settings that exposed the LT NaN."""
    beam, q2max, epmin = {
        "rgk": ("6.535", "6.535", "1.0"),
        "rga": ("10.604", "10.5", "2.0"),
    }[region]
    return "\n".join([
        "5", "1", "3", "1", beam, f"1.0 {q2max}", f"{epmin} {beam}",
        "5000", "2.0", "0", str(seed), "3", "0.05 0.70", "0.09 2.0",
        "0.0 360.0", "1 2.0 1.0", "",
    ])


class BornMode3SafetyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.executable = Path(__file__).resolve().parent / "build" / "aao_norad"
        if not cls.executable.is_file():
            raise unittest.SkipTest("build/aao_norad is created by the Makefile")

    def test_success_records_safe_schema_and_complete_overshoot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            completed = subprocess.run(
                [str(self.executable)],
                input=_input(events=2, fmcall=2.0, seed=812_001),
                text=True,
                cwd=work,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            norm = (work / "aao_norad.norm").read_text(encoding="utf-8")
            events = int(float(_field(norm, "events")))
            target = int(float(_field(norm, "mode3_target_events")))
            overshoot = int(float(_field(norm, "mode3_event_overshoot")))
            maximum = int(float(_field(norm, "mcall_max")))
            self.assertEqual(_field(norm, "mode3_schema"), "aao-norad-mode3-v2")
            self.assertEqual(
                int(float(_field(norm, "mode3_complete_final_multiplicity"))), 1
            )
            self.assertEqual(events - target, overshoot)
            self.assertGreaterEqual(overshoot, 0)
            self.assertLess(overshoot, maximum)
            self.assertEqual(
                len((work / "aao_norad.lund").read_text().splitlines()),
                5 * events,
            )

    def test_catastrophic_multiplicity_fails_without_event_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            completed = subprocess.run(
                [str(self.executable)],
                input=_input(events=2, fmcall=1.0e-30, seed=812_002),
                text=True,
                cwd=work,
                capture_output=True,
                check=False,
            )
            output = completed.stdout + completed.stderr
            self.assertNotEqual(completed.returncode, 0, output)
            self.assertIn(
                "FATAL Born mode-3 multiplicity exceeds job capacity", output
            )
            self.assertFalse((work / "aao_norad.lund").exists())
            self.assertFalse((work / "aao_norad.kin").exists())
            self.assertFalse((work / "aao_norad.norm").exists())

    def test_nonfinite_ratio_has_distinct_diagnostic_and_no_products(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            completed = subprocess.run(
                [str(self.executable)],
                input=_input(events=2, fmcall=float("nan"), seed=812_002),
                text=True, cwd=work, capture_output=True, timeout=60,
            )
            output = completed.stdout + completed.stderr
            self.assertNotEqual(completed.returncode, 0, output)
            self.assertIn("FATAL Born mode-3 non-finite cross section", output)
            self.assertNotIn("multiplicity exceeds job capacity", output)
            self.assertIn("Q2, xB, W, cos(theta*), phi*:", output)
            for suffix in ("lund", "kin", "norm"):
                self.assertFalse((work / f"aao_norad.{suffix}").exists())

    def test_production_nan_seeds_finish_with_finite_complete_events(self) -> None:
        # Before the fix these stopped at (ntries, events) = (95924, 866)
        # and (332374, 3503), respectively, despite a finite scan maximum.
        for region, seed in (("rgk", -2100000029), ("rga", -2110004993)):
            with self.subTest(region=region), tempfile.TemporaryDirectory() as tmp:
                work = Path(tmp)
                completed = subprocess.run(
                    [str(self.executable)], input=_production_input(region, seed),
                    text=True, cwd=work, capture_output=True, timeout=120,
                )
                output = completed.stdout + completed.stderr
                self.assertEqual(completed.returncode, 0, output)
                self.assertNotIn("IEEE_INVALID_FLAG", output)
                norm = (work / "aao_norad.norm").read_text()
                self.assertEqual(_field(norm, "dvmp_boundary_version"), "1")
                events = int(float(_field(norm, "events")))
                self.assertGreaterEqual(events, 5000)
                self.assertEqual(int(_field(norm, "mode3_target_events")), 5000)
                overshoot = int(_field(norm, "mode3_event_overshoot"))
                self.assertEqual(events - 5000, overshoot)
                self.assertLess(overshoot, int(_field(norm, "mcall_max")))
                for field in ("sig_sum", "sigr_max"):
                    value = float(_field(norm, field).replace("D", "E"))
                    self.assertTrue(math.isfinite(value), norm)
                    self.assertGreater(value, 0.0)
                lund = (work / "aao_norad.lund").read_text()
                self.assertEqual(len(lund.splitlines()), 5 * events)
                self.assertIsNone(re.search(r"(?i)\b(?:nan|inf(?:inity)?)\b", lund))


if __name__ == "__main__":
    unittest.main()
