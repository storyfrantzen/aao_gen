"""Compile the endpoint regression probe against each generator's model copy."""

from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


class DVMPBoundaryTests(unittest.TestCase):
    def test_model_copies_match(self):
        root = Path(__file__).resolve().parents[1]
        for name in ("dvmpx.F", "dvmpw.F"):
            with self.subTest(source=name):
                self.assertEqual((root / "aao_rad" / name).read_bytes(),
                                 (root / "aao_norad" / name).read_bytes())

    @unittest.skipUnless(shutil.which("gfortran"), "requires gfortran")
    def test_boundary_math_in_both_generators(self):
        root = Path(__file__).resolve().parents[1]
        probe = Path(__file__).with_suffix(".f90")
        for generator in ("aao_norad", "aao_rad"):
            with self.subTest(generator=generator), tempfile.TemporaryDirectory() as tmp:
                executable = Path(tmp) / "boundary_probe"
                build = subprocess.run([
                    "gfortran", "-fno-automatic", "-ffixed-line-length-none",
                    "-fcheck=all", str(root / generator / "dvmpx.F"),
                    str(root / generator / "dvmpw.F"), str(probe),
                    "-o", str(executable),
                ], text=True, capture_output=True, timeout=60)
                self.assertEqual(build.returncode, 0, build.stdout + build.stderr)
                result = subprocess.run([str(executable)], text=True,
                                        capture_output=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("PASS boundary grid points:", result.stdout)


if __name__ == "__main__":
    unittest.main()
