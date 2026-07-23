#!/usr/bin/env python3
"""Prepare and run Born AAO strata that are unweighted within analysis bins."""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import math
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np


PROTON_MASS_GEV = 0.93828
PI0_MASS_GEV = 0.1349
MANIFEST_SCHEMA = "aao-born-bin-conditional-v1"
SAMPLING_MODE = 4


@dataclass(frozen=True)
class Stratum:
    flat_index: int
    iq2: int
    ixb: int
    it: int
    iphi: int
    q2: tuple[float, float]
    xb: tuple[float, float]
    minus_t: tuple[float, float]
    phi_deg: tuple[float, float]

    @property
    def identifier(self) -> str:
        return (
            f"bin_{self.flat_index:05d}_q{self.iq2:02d}_x{self.ixb:02d}_"
            f"t{self.it:02d}_p{self.iphi:02d}"
        )


def legacy_flat_index(
    iq2: int,
    ixb: int,
    it: int,
    iphi: int,
    *,
    nq2: int,
    nt: int,
    nphi: int,
) -> int:
    """Match the analysis response order: xB, Q2, phi, then -t fastest."""
    return ixb * nq2 * nphi * nt + iq2 * nphi * nt + iphi * nt + it


def load_analysis_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as source:
        config = json.load(source)
    binning = config.get("binning", {})
    required_edges = ("Q2", "xB", "minus_t", "phi_deg")
    missing_edges = [name for name in required_edges if name not in binning]
    if missing_edges:
        raise ValueError(f"analysis config is missing bin edges: {', '.join(missing_edges)}")
    phase_space = config.get("phase_space", {})
    missing_cuts = [name for name in ("W_min", "y_max") if name not in phase_space]
    if missing_cuts:
        raise ValueError(f"analysis config is missing phase-space cuts: {', '.join(missing_cuts)}")
    if "beam_energy" not in config:
        raise ValueError("analysis config is missing beam_energy")
    return config


def enumerate_strata(config: dict) -> list[Stratum]:
    edges = config["binning"]
    q2_edges = _strict_edges(edges["Q2"], "Q2")
    xb_edges = _strict_edges(edges["xB"], "xB")
    t_edges = _strict_edges(edges["minus_t"], "minus_t")
    phi_edges = _strict_edges(edges["phi_deg"], "phi_deg")
    nq2 = len(q2_edges) - 1
    nt = len(t_edges) - 1
    nphi = len(phi_edges) - 1
    strata: list[Stratum] = []
    for iq2 in range(nq2):
        for ixb in range(len(xb_edges) - 1):
            for it in range(nt):
                for iphi in range(nphi):
                    strata.append(
                        Stratum(
                            flat_index=legacy_flat_index(
                                iq2,
                                ixb,
                                it,
                                iphi,
                                nq2=nq2,
                                nt=nt,
                                nphi=nphi,
                            ),
                            iq2=iq2,
                            ixb=ixb,
                            it=it,
                            iphi=iphi,
                            q2=(q2_edges[iq2], q2_edges[iq2 + 1]),
                            xb=(xb_edges[ixb], xb_edges[ixb + 1]),
                            minus_t=(t_edges[it], t_edges[it + 1]),
                            phi_deg=(phi_edges[iphi], phi_edges[iphi + 1]),
                        )
                    )
    return sorted(strata, key=lambda item: item.flat_index)


def has_selected_physical_support(
    stratum: Stratum,
    *,
    beam_energy: float,
    q2_minimum: float,
    w_minimum: float,
    y_maximum: float,
    samples: int = 32,
) -> bool:
    """Numerically preflight whether a stratum intersects selected pi0 phase space."""
    if samples <= 0:
        raise ValueError("support samples must be positive")
    q2_low = max(stratum.q2[0], q2_minimum)
    q2_high = stratum.q2[1]
    if q2_high <= q2_low:
        return False
    q2 = _midpoints(q2_low, q2_high, samples)
    xb = _midpoints(stratum.xb[0], stratum.xb[1], samples)
    q2_mesh, xb_mesh = np.meshgrid(q2, xb, indexing="ij")
    w2 = PROTON_MASS_GEV**2 + q2_mesh * (1.0 / xb_mesh - 1.0)
    y = q2_mesh / (2.0 * PROTON_MASS_GEV * beam_energy * xb_mesh)
    eprime = beam_energy * (1.0 - y)
    sin2_half = np.divide(
        q2_mesh,
        4.0 * beam_energy * eprime,
        out=np.full_like(q2_mesh, np.nan),
        where=eprime > 0.0,
    )
    selected = (
        (w2 >= w_minimum**2)
        & (y <= y_maximum)
        & (eprime > 0.0)
        & (sin2_half > 0.0)
        & (sin2_half < 1.0)
    )
    t_low, t_high = _t_limits_pi0(xb_mesh, q2_mesh)
    signed_bin_low = -stratum.minus_t[1]
    signed_bin_high = -stratum.minus_t[0]
    t_overlap = np.maximum(
        0.0,
        np.minimum(t_high, signed_bin_high) - np.maximum(t_low, signed_bin_low),
    )
    return bool(np.any(selected & np.isfinite(t_overlap) & (t_overlap > 0.0)))


def prepare(args: argparse.Namespace) -> Path:
    config_path = args.config.resolve()
    config = load_analysis_config(config_path)
    beam_energy = float(config["beam_energy"])
    phase_space = config["phase_space"]
    q2_minimum = float(phase_space.get("Q2_min", config["binning"]["Q2"][0]))
    w_minimum = float(phase_space["W_min"])
    y_maximum = float(phase_space["y_max"])
    if not (0.0 < y_maximum <= 1.0):
        raise ValueError("phase_space.y_max must satisfy 0 < y_max <= 1")
    if w_minimum <= PROTON_MASS_GEV + PI0_MASS_GEV:
        raise ValueError("phase_space.W_min must exceed the pi0 production threshold")
    if args.epirea != 1:
        raise ValueError("bin-conditional generation currently supports only epirea=1 (pi0)")
    if args.npart != 3:
        raise ValueError("bin-conditional pi0 LUND output requires npart=3")

    output = args.output.resolve()
    manifest_path = output / "manifest.json"
    if manifest_path.exists() and not args.overwrite:
        raise FileExistsError(
            f"manifest already exists: {manifest_path}; pass --overwrite to replace it"
        )
    input_dir = output / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    config_snapshot = output / "analysis_config.json"
    shutil.copy2(config_path, config_snapshot)
    all_strata = enumerate_strata(config)
    stop = len(all_strata) if args.bin_stop is None else args.bin_stop
    if args.bin_start < 0 or stop < args.bin_start or stop > len(all_strata):
        raise ValueError(f"invalid flat stratum range [{args.bin_start}, {stop})")

    records: list[dict] = []
    skipped = 0
    ep_min = max(0.001, beam_energy * (1.0 - y_maximum))
    for stratum in all_strata:
        if stratum.flat_index < args.bin_start or stratum.flat_index >= stop:
            continue
        supported = has_selected_physical_support(
            stratum,
            beam_energy=beam_energy,
            q2_minimum=q2_minimum,
            w_minimum=w_minimum,
            y_maximum=y_maximum,
            samples=args.support_samples,
        )
        if not supported and not args.include_unsupported:
            skipped += 1
            continue
        q2_bounds = (max(stratum.q2[0], q2_minimum), stratum.q2[1])
        seed = -abs(args.seed_base + stratum.flat_index)
        input_name = f"{stratum.identifier}.inp"
        input_path = input_dir / input_name
        input_path.write_text(
            generator_input(
                stratum,
                q2_bounds=q2_bounds,
                beam_energy=beam_energy,
                ep_min=ep_min,
                events=args.events_per_bin,
                physics_model=args.physics_model,
                flag_ehel=args.flag_ehel,
                npart=args.npart,
                epirea=args.epirea,
                fmcall=args.fmcall,
                boso=args.boso,
                seed=seed,
                w_minimum=w_minimum,
                y_maximum=y_maximum,
            ),
            encoding="utf-8",
        )
        records.append(
            {
                "stratum_id": stratum.identifier,
                "flat_index": stratum.flat_index,
                "indices": {
                    "iq2": stratum.iq2,
                    "ixb": stratum.ixb,
                    "it": stratum.it,
                    "iphi": stratum.iphi,
                },
                "bounds": {
                    "Q2": list(q2_bounds),
                    "xB": list(stratum.xb),
                    "minus_t": list(stratum.minus_t),
                    "phi_deg": list(stratum.phi_deg),
                },
                "selected_physical_support": supported,
                "events_requested": args.events_per_bin,
                "seed": seed,
                "input_file": str(input_path.relative_to(output)),
                "output_stem": str(Path("outputs") / stratum.identifier / stratum.identifier),
            }
        )

    manifest = {
        "schema": MANIFEST_SCHEMA,
        "generator": "aao_norad",
        "sampling_mode": SAMPLING_MODE,
        "sampling_description": (
            "inverse-Q2 unweighted sampling conditional on one selected analysis bin"
        ),
        "analysis_config": str(config_path),
        "analysis_config_snapshot": config_snapshot.name,
        "analysis_config_sha256": _sha256(config_snapshot),
        "prepared_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "generator_revision": _generator_revision(),
        "proposal_coordinates": ["1/Q2", "xB", "minus_t", "phi_deg"],
        "density_jacobian": "Q2^3/(8*M*xB^2*E*Eprime*q_cm*p_pi_cm)",
        "beam_energy": beam_energy,
        "phase_space": {
            "Q2_min": q2_minimum,
            "W_min": w_minimum,
            "y_max": y_maximum,
        },
        "flat_order": "xB, Q2, phi, then minus_t fastest",
        "events_per_bin": args.events_per_bin,
        "physics_model": args.physics_model,
        "flag_ehel": args.flag_ehel,
        "npart": args.npart,
        "epirea": args.epirea,
        "fmcall": args.fmcall,
        "boso": args.boso,
        "support_samples": args.support_samples,
        "total_analysis_bins": len(all_strata),
        "prepared_strata": len(records),
        "unsupported_strata_skipped": skipped,
        "strata": records,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    _write_tsv(output / "strata.tsv", records)
    return manifest_path


def generator_input(
    stratum: Stratum,
    *,
    q2_bounds: tuple[float, float],
    beam_energy: float,
    ep_min: float,
    events: int,
    physics_model: int,
    flag_ehel: int,
    npart: int,
    epirea: int,
    fmcall: float,
    boso: int,
    seed: int,
    w_minimum: float,
    y_maximum: float,
) -> str:
    return (
        f"{physics_model}\n{flag_ehel}\n{npart}\n{epirea}\n{beam_energy:.12g}\n"
        f"{q2_bounds[0]:.12g} {q2_bounds[1]:.12g}\n"
        f"{ep_min:.12g} {beam_energy:.12g}\n"
        f"{events}\n{fmcall:.12g}\n{boso}\n{seed}\n{SAMPLING_MODE}\n"
        f"{stratum.xb[0]:.12g} {stratum.xb[1]:.12g}\n"
        f"{stratum.minus_t[0]:.12g} {stratum.minus_t[1]:.12g}\n"
        f"{stratum.phi_deg[0]:.12g} {stratum.phi_deg[1]:.12g}\n"
        f"{w_minimum:.12g} {y_maximum:.12g}\n"
        f"{stratum.flat_index} {stratum.iq2} {stratum.ixb} "
        f"{stratum.it} {stratum.iphi}\n"
    )


def run(args: argparse.Namespace) -> Path:
    manifest_path = args.manifest.resolve()
    with manifest_path.open(encoding="utf-8") as source:
        manifest = json.load(source)
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError(f"unsupported manifest schema: {manifest.get('schema')}")
    matches = [
        record for record in manifest["strata"] if int(record["flat_index"]) == args.flat_index
    ]
    if len(matches) != 1:
        raise ValueError(f"manifest has {len(matches)} records for flat index {args.flat_index}")
    record = matches[0]
    root = manifest_path.parent
    input_path = root / record["input_file"]
    output_stem = root / record["output_stem"]
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    expected = output_stem.with_suffix(".norm")
    if expected.exists() and not args.overwrite:
        raise FileExistsError(f"output already exists: {expected}")
    executable = args.executable.resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"generator executable not found: {executable}")

    with tempfile.TemporaryDirectory(prefix=f"{record['stratum_id']}_") as temporary:
        work = Path(temporary)
        completed = subprocess.run(
            [str(executable)],
            input=input_path.read_text(encoding="utf-8"),
            text=True,
            cwd=work,
            capture_output=True,
            check=False,
        )
        (work / "stdout.log").write_text(completed.stdout, encoding="utf-8")
        (work / "stderr.log").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError(
                f"AAO failed for {record['stratum_id']} with exit code "
                f"{completed.returncode}:\n{completed.stderr[-2000:]}"
            )
        norm_path = work / "aao_norad.norm"
        lund_path = work / "aao_norad.lund"
        if not norm_path.is_file() or not lund_path.is_file():
            raise RuntimeError("AAO completed without both .norm and .lund outputs")
        norm = _norm_fields(norm_path)
        _validate_norm_record(norm, record)
        _validate_kinematics(
            work / "aao_norad.kin",
            record,
            beam_energy=float(manifest["beam_energy"]),
            phase_space=manifest["phase_space"],
        )
        _validate_lund_event_count(lund_path, _int_field(norm, "events"))

        for suffix, source_name in (
            (".lund", "aao_norad.lund"),
            (".norm", "aao_norad.norm"),
            (".kin", "aao_norad.kin"),
            (".sum", "aao_norad.sum"),
            (".out", "aao_norad.out"),
            (".stdout", "stdout.log"),
            (".stderr", "stderr.log"),
        ):
            source_path = work / source_name
            if source_path.exists():
                destination = output_stem.with_suffix(suffix)
                if destination.exists():
                    destination.unlink()
                shutil.move(str(source_path), destination)
        shutil.copy2(input_path, output_stem.with_suffix(".inp"))

    run_record = {
        **record,
        "schema": MANIFEST_SCHEMA,
        "sig_sum_microbarn": _float_field(norm, "sig_sum"),
        "events": _int_field(norm, "events"),
        "ntries": _int_field(norm, "ntries"),
        "mcall_max": _int_field(norm, "mcall_max"),
        "event_weight_microbarn": _float_field(norm, "stratum_event_weight_microbarn"),
        "proposal_efficiency": _int_field(norm, "events") / _int_field(norm, "ntries"),
    }
    record_path = output_stem.with_suffix(".json")
    record_path.write_text(json.dumps(run_record, indent=2) + "\n", encoding="utf-8")
    return record_path


def _validate_norm_record(norm: dict[str, str], record: dict) -> None:
    expected = {
        "sampling_mode": SAMPLING_MODE,
        "conditional_unweighted": 1,
        "stratum_flat_index": int(record["flat_index"]),
        "stratum_iq2": int(record["indices"]["iq2"]),
        "stratum_ixb": int(record["indices"]["ixb"]),
        "stratum_it": int(record["indices"]["it"]),
        "stratum_iphi": int(record["indices"]["iphi"]),
    }
    for key, value in expected.items():
        if _int_field(norm, key) != value:
            raise RuntimeError(f"normalization metadata mismatch for {key}")
    if _int_field(norm, "events") != int(record["events_requested"]):
        raise RuntimeError("normalization event count does not match the manifest")
    if _int_field(norm, "mcall_max") > 1:
        raise RuntimeError(
            "mcall_max exceeded 1; increase fmcall and rerun this stratum"
        )


def _validate_kinematics(
    path: Path,
    record: dict,
    *,
    beam_energy: float,
    phase_space: dict,
) -> None:
    kinematics = np.loadtxt(path, comments="#", ndmin=2)
    if kinematics.shape != (int(record["events_requested"]), 7):
        raise RuntimeError(f"unexpected kinematics sidecar shape {kinematics.shape}")
    expected_ids = np.arange(1, kinematics.shape[0] + 1)
    if not np.array_equal(kinematics[:, 0].astype(int), expected_ids):
        raise RuntimeError("kinematics sidecar has non-contiguous event identifiers")
    q2, xb, minus_t, phi_deg = (kinematics[:, index] for index in range(1, 5))
    for values, name in (
        (q2, "Q2"),
        (xb, "xB"),
        (minus_t, "minus_t"),
        (phi_deg, "phi_deg"),
    ):
        low, high = record["bounds"][name]
        tolerance = 2.0e-6 * max(1.0, abs(low), abs(high))
        if np.any(values < low - tolerance) or np.any(values > high + tolerance):
            raise RuntimeError(f"generated {name} lies outside its stratum bounds")
    y = q2 / (2.0 * PROTON_MASS_GEV * beam_energy * xb)
    w2 = PROTON_MASS_GEV**2 + q2 * (1.0 / xb - 1.0)
    if np.any(y > float(phase_space["y_max"]) + 2.0e-6):
        raise RuntimeError("generated event violates the configured y cut")
    if np.any(w2 < float(phase_space["W_min"]) ** 2 - 2.0e-6):
        raise RuntimeError("generated event violates the configured W cut")


def _validate_lund_event_count(path: Path, expected: int) -> None:
    count = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 5 and fields[:5] == ["4", "1", "1", "0", "0"]:
            count += 1
    if count != expected:
        raise RuntimeError(f"LUND output contains {count} events; expected {expected}")


def _norm_fields(path: Path) -> dict[str, str]:
    fields: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$", line)
        if match:
            fields[match.group(1)] = match.group(2)
    return fields


def _float_field(fields: dict[str, str], key: str) -> float:
    if key not in fields:
        raise RuntimeError(f"normalization output is missing {key}")
    return float(fields[key].replace("D", "E").replace("d", "e"))


def _int_field(fields: dict[str, str], key: str) -> int:
    return int(_float_field(fields, key))


def _strict_edges(values: list[float], name: str) -> tuple[float, ...]:
    edges = tuple(float(value) for value in values)
    if len(edges) < 2 or any(high <= low for low, high in zip(edges[:-1], edges[1:])):
        raise ValueError(f"{name} edges must be strictly increasing")
    return edges


def _midpoints(low: float, high: float, count: int) -> np.ndarray:
    width = (high - low) / count
    return low + (np.arange(count) + 0.5) * width


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _generator_revision() -> str:
    repository = Path(__file__).resolve().parents[1]
    completed = subprocess.run(
        ["git", "describe", "--always", "--dirty"],
        cwd=repository,
        text=True,
        capture_output=True,
        check=False,
    )
    return completed.stdout.strip() if completed.returncode == 0 else "unknown"


def _kallen(a: np.ndarray, b: float, c: float) -> np.ndarray:
    return a * a + b * b + c * c - 2.0 * a * b - 2.0 * a * c - 2.0 * b * c


def _t_limits_pi0(xb: np.ndarray, q2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xb, q2 = np.broadcast_arrays(np.asarray(xb, float), np.asarray(q2, float))
    t_low = np.full(xb.shape, np.nan)
    t_high = np.full(xb.shape, np.nan)
    valid = (xb > 0.0) & (xb < 1.0) & (q2 > 0.0)
    w2 = PROTON_MASS_GEV**2 + q2 * (1.0 / xb - 1.0)
    valid &= w2 > (PROTON_MASS_GEV + PI0_MASS_GEV) ** 2
    if not np.any(valid):
        return t_low, t_high
    w = np.sqrt(w2[valid])
    q0_cm = (w2[valid] - PROTON_MASS_GEV**2 - q2[valid]) / (2.0 * w)
    q_cm2 = q0_cm**2 + q2[valid]
    epi_cm = (w2[valid] + PI0_MASS_GEV**2 - PROTON_MASS_GEV**2) / (2.0 * w)
    ppi_cm2 = _kallen(w2[valid], PI0_MASS_GEV**2, PROTON_MASS_GEV**2) / (4.0 * w2[valid])
    subvalid = (q_cm2 > 0.0) & (ppi_cm2 > 0.0)
    valid_indices = np.flatnonzero(valid)[subvalid]
    q_cm = np.sqrt(q_cm2[subvalid])
    ppi_cm = np.sqrt(ppi_cm2[subvalid])
    t_forward = (
        PI0_MASS_GEV**2
        - q2.ravel()[valid_indices]
        - 2.0 * q0_cm[subvalid] * epi_cm[subvalid]
        + 2.0 * q_cm * ppi_cm
    )
    t_backward = (
        PI0_MASS_GEV**2
        - q2.ravel()[valid_indices]
        - 2.0 * q0_cm[subvalid] * epi_cm[subvalid]
        - 2.0 * q_cm * ppi_cm
    )
    t_low.ravel()[valid_indices] = np.minimum(t_backward, t_forward)
    t_high.ravel()[valid_indices] = np.maximum(t_backward, t_forward)
    return t_low, t_high


def _write_tsv(path: Path, records: list[dict]) -> None:
    header = (
        "flat_index\tstratum_id\tiq2\tixb\tit\tiphi\tq2_low\tq2_high\t"
        "xb_low\txb_high\tminus_t_low\tminus_t_high\tphi_low\tphi_high\t"
        "events_requested\tseed\tinput_file\toutput_stem\n"
    )
    rows = [header]
    for record in records:
        bounds = record["bounds"]
        indices = record["indices"]
        rows.append(
            "\t".join(
                str(value)
                for value in (
                    record["flat_index"],
                    record["stratum_id"],
                    indices["iq2"],
                    indices["ixb"],
                    indices["it"],
                    indices["iphi"],
                    *bounds["Q2"],
                    *bounds["xB"],
                    *bounds["minus_t"],
                    *bounds["phi_deg"],
                    record["events_requested"],
                    record["seed"],
                    record["input_file"],
                    record["output_stem"],
                )
            )
            + "\n"
        )
    path.write_text("".join(rows), encoding="utf-8")


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare", help="write one mode-4 input per stratum")
    prepare_parser.add_argument("--config", type=Path, required=True)
    prepare_parser.add_argument("--output", type=Path, required=True)
    prepare_parser.add_argument("--events-per-bin", type=int, default=10_000)
    prepare_parser.add_argument("--physics-model", type=int, default=5)
    prepare_parser.add_argument("--flag-ehel", type=int, choices=(0, 1), default=1)
    prepare_parser.add_argument("--npart", type=int, default=3)
    prepare_parser.add_argument("--epirea", type=int, default=1)
    prepare_parser.add_argument("--fmcall", type=float, default=2.0)
    prepare_parser.add_argument("--boso", type=int, choices=(0, 1), default=0)
    prepare_parser.add_argument("--seed-base", type=int, default=5_000_001)
    prepare_parser.add_argument("--bin-start", type=int, default=0)
    prepare_parser.add_argument("--bin-stop", type=int)
    prepare_parser.add_argument("--support-samples", type=int, default=32)
    prepare_parser.add_argument("--include-unsupported", action="store_true")
    prepare_parser.add_argument("--overwrite", action="store_true")

    run_parser = subparsers.add_parser("run", help="execute one prepared stratum")
    run_parser.add_argument("manifest", type=Path)
    run_parser.add_argument("--flat-index", type=int, required=True)
    run_parser.add_argument("--executable", type=Path, required=True)
    run_parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = _arguments()
    if args.command == "prepare":
        if args.events_per_bin <= 0:
            raise ValueError("--events-per-bin must be positive")
        if args.fmcall <= 0.0:
            raise ValueError("--fmcall must be positive")
        path = prepare(args)
        print(f"Wrote {path}")
    else:
        path = run(args)
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
