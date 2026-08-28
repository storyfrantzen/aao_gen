#!/usr/bin/env python3
"""Regression tests for fail-safe global Born mode-3 multiplicity."""

from __future__ import annotations

import re
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


if __name__ == "__main__":
    unittest.main()
