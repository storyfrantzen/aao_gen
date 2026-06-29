#!/usr/bin/env python3
"""
avg_d4sigma.py
--------------
Average d4sigma/dQ2 dxB dt dphi over the physical part of a
(xB, Q2, t, phi) analysis bin by calling aao_xsec on a 4D midpoint grid.

Important:
    The bin is sampled as an enclosing rectangular grid, but only
    kinematically physical points are included in the average.

    This avoids treating inaccessible regions of the rectangular
    (xB,Q2,t,phi) box as zero cross section.

Usage:
    python avg_d4sigma.py \
        --xB    0.25 0.35 \
        --Q2    1.5  2.5  \
        --t    -0.4 -0.1  \
        --phi   0   60    \
        --BeamEnergy 10.6 \
        --N 8             \
        --exe ./aao_xsec
"""

import argparse
import subprocess
import sys
import os
import itertools
import multiprocessing
from functools import partial

import numpy as np
if "MPLCONFIGDIR" not in os.environ:
    os.environ["MPLCONFIGDIR"] = os.path.join(os.getcwd(), ".matplotlib")
    os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MP = 0.938272081   # proton mass [GeV]
MPI0 = 0.1349768   # pi0 mass [GeV]
ALPHA = 1.0 / 137.036


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def kallen(a, b, c):
    return a*a + b*b + c*c - 2.0*a*b - 2.0*a*c - 2.0*b*c


def midpoint_grid(lo, hi, N):
    """Return N midpoints of equally spaced subintervals between lo and hi."""
    edges = np.linspace(lo, hi, N + 1)
    return 0.5 * (edges[:-1] + edges[1:])


def edge_pairs(edges, name):
    """Return neighboring edge pairs after checking that edges increase."""
    edges = np.asarray(edges, dtype=float)
    if edges.ndim != 1 or len(edges) < 2:
        raise ValueError(f"{name} must contain at least two edges.")
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError(f"{name} edges must be strictly increasing.")
    return list(zip(edges[:-1], edges[1:]))


def t_bin_for_aao(t_lo, t_hi, signed_t_edges):
    """
    Convert user t-bin edges to the negative-valued t convention used by aao.

    By default quilt mode treats positive t edges as -t edges, so a bin
    [0.2, 0.3] becomes actual t=[-0.3, -0.2].
    """
    if signed_t_edges:
        return [t_lo, t_hi]
    return [-t_hi, -t_lo]


# ---------------------------------------------------------------------------
# Kinematic checks
# ---------------------------------------------------------------------------

def t_limits_pi0(xB, Q2):
    """
    Return the physical t limits for gamma* p -> p pi0 at fixed xB,Q2.

    Uses t = (q - pi0)^2.

    Returns:
        (t_min, t_max) in GeV^2, where both are negative-valued in normal
        DVpi0 kinematics. If kinematics are impossible, returns (nan, nan).
    """
    W2 = MP*MP + Q2 * (1.0 / xB - 1.0)

    if W2 <= (MP + MPI0)**2:
        return float("nan"), float("nan")

    W = np.sqrt(W2)

    # gamma* energy and 3-momentum in gamma*p CM
    q0_cm = (W2 - MP*MP - Q2) / (2.0 * W)
    q_cm2 = q0_cm*q0_cm + Q2

    if q_cm2 <= 0.0:
        return float("nan"), float("nan")

    q_cm = np.sqrt(q_cm2)

    # pi0 energy and 3-momentum in gamma*p CM
    epi_cm = (W2 + MPI0*MPI0 - MP*MP) / (2.0 * W)
    ppi_cm2 = kallen(W2, MPI0*MPI0, MP*MP) / (4.0 * W2)

    if ppi_cm2 <= 0.0:
        return float("nan"), float("nan")

    ppi_cm = np.sqrt(ppi_cm2)

    # t = m_pi^2 - Q2 - 2 q0 Epi + 2 |q||pi| cos(theta*)
    t_forward = MPI0*MPI0 - Q2 - 2.0*q0_cm*epi_cm + 2.0*q_cm*ppi_cm
    t_backward = MPI0*MPI0 - Q2 - 2.0*q0_cm*epi_cm - 2.0*q_cm*ppi_cm

    return min(t_backward, t_forward), max(t_backward, t_forward)


def physical_mask(xB, Q2, t, ebeam):
    """
    True if the point is physically allowed for ep -> e p pi0.

    Checks:
      - 0 < xB < 1
      - scattered electron energy positive
      - electron scattering angle kinematics valid
      - W above p pi0 threshold
      - t lies inside physical gamma*p -> p pi0 limits
    """
    if xB <= 0.0 or xB >= 1.0 or Q2 <= 0.0:
        return False

    W2 = MP*MP + Q2 * (1.0 / xB - 1.0)
    if W2 <= (MP + MPI0)**2:
        return False

    nu = Q2 / (2.0 * MP * xB)
    eprime = ebeam - nu
    if eprime <= 0.0:
        return False

    # massless electron kinematic consistency:
    # Q2 = 4 E E' sin^2(theta/2)
    sin2_half = Q2 / (4.0 * ebeam * eprime)
    if sin2_half <= 0.0 or sin2_half >= 1.0:
        return False

    t_lo, t_hi = t_limits_pi0(xB, Q2)
    if not np.isfinite(t_lo) or not np.isfinite(t_hi):
        return False

    return (t >= t_lo) and (t <= t_hi)


# ---------------------------------------------------------------------------
# Single executable call
# ---------------------------------------------------------------------------

def call_aao_xsec(point, exe, ebeam, theory, channel, resonance, verbose_failures=False):
    """
    Call aao_xsec for one physical point.

    Returns:
        SIGU printed by executable, or NaN on failure.
    """
    xB, Q2, t, phi = point

    cmd = [
        exe,
        "-xB",          str(xB),
        "-Q2",          str(Q2),
        "-t",           str(t),
        "-phi",         str(phi),
        "-BeamEnergy",  str(ebeam),
        "-theory",      str(theory),
        "-channel",     str(channel),
        "-resonance",   str(resonance),
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        return float(result.stdout.strip())

    except subprocess.CalledProcessError as e:
        if verbose_failures:
            detail = (e.stderr or e.stdout or "").strip()
            if detail:
                detail = f"\n  executable output: {detail}"
            sys.stderr.write(
                f"WARNING: aao_xsec failed for "
                f"xB={xB:.6g} Q2={Q2:.6g} t={t:.6g} phi={phi:.6g}: {e}{detail}\n"
            )
        return float("nan")

    except ValueError as e:
        if verbose_failures:
            sys.stderr.write(
                f"WARNING: could not parse aao_xsec output for "
                f"xB={xB:.6g} Q2={Q2:.6g} t={t:.6g} phi={phi:.6g}: {e}\n"
            )
        return float("nan")


# ---------------------------------------------------------------------------
# Flux factor Gamma_v(xB, Q2, Ebeam)
# ---------------------------------------------------------------------------

def gamma_v(xB, Q2, ebeam):
    """
    Virtual-photon transverse flux factor in the Hand convention.

    Returns Gamma_v such that

        d4sigma/dQ2 dxB dt dphi = Gamma_v * SIGU

    assuming SIGU is the reduced/d2 hadronic cross section from aao_xsec.
    """
    W2 = MP*MP + Q2 * (1.0 / xB - 1.0)
    nu = Q2 / (2.0 * MP * xB)
    eprime = ebeam - nu

    if eprime <= 0.0:
        return float("nan")

    K = (W2 - MP*MP) / (2.0 * MP)

    denom_tan = 4.0 * ebeam * eprime - Q2
    if denom_tan <= 0.0:
        return float("nan")

    tan2 = Q2 / denom_tan

    eps = 1.0 / (
        1.0 + 2.0 * (1.0 + nu*nu / Q2) * tan2
    )

    if eps >= 1.0:
        return float("nan")

    gv = (
        ALPHA / (2.0 * np.pi**2)
        * (eprime / ebeam)
        * (K / Q2)
        / (1.0 - eps)
    )

    return gv


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def compute_bin_centering(args, phi_lo, phi_hi, N):
    total_grid_points = N**4

    xB_pts  = midpoint_grid(args.xB[0], args.xB[1], N)
    Q2_pts  = midpoint_grid(args.Q2[0], args.Q2[1], N)
    t_pts   = midpoint_grid(args.t[0], args.t[1], N)
    phi_pts = midpoint_grid(phi_lo, phi_hi, N)

    all_points = list(itertools.product(xB_pts, Q2_pts, t_pts, phi_pts))

    physical_points = [
        p for p in all_points
        if physical_mask(p[0], p[1], p[2], args.BeamEnergy)
    ]

    n_phys = len(physical_points)
    if n_phys == 0:
        return None

    worker = partial(
        call_aao_xsec,
        exe=args.exe,
        ebeam=args.BeamEnergy,
        theory=args.theory,
        channel=args.channel,
        resonance=args.resonance,
        verbose_failures=args.verbose_failures,
    )

    n_workers = args.workers or multiprocessing.cpu_count()

    with multiprocessing.Pool(processes=n_workers) as pool:
        sigu_values = pool.map(worker, physical_points)

    sigu_arr = np.array(sigu_values, dtype=float)
    n_xsec_failed = int(np.sum(~np.isfinite(sigu_arr)))

    d4sig_values = np.array([
        gamma_v(xB, Q2, args.BeamEnergy) * sigu
        for (xB, Q2, t, phi), sigu in zip(physical_points, sigu_arr)
    ])

    valid = np.isfinite(d4sig_values) & (d4sig_values > 0.0)

    if not np.any(valid):
        return None

    avg = np.mean(d4sig_values[valid])

    valid_points = np.array([
        p for p, is_valid in zip(physical_points, valid)
        if is_valid
    ], dtype=float)

    xB_c = np.mean(valid_points[:, 0])
    Q2_c = np.mean(valid_points[:, 1])
    t_c  = np.mean(valid_points[:, 2])

    phi_rad = np.deg2rad(valid_points[:, 3])
    phi_c = np.rad2deg(
        np.arctan2(
            np.mean(np.sin(phi_rad)),
            np.mean(np.cos(phi_rad))
        )
    ) % 360.0

    sigu_center = call_aao_xsec(
        (xB_c, Q2_c, t_c, phi_c),
        exe=args.exe,
        ebeam=args.BeamEnergy,
        theory=args.theory,
        channel=args.channel,
        resonance=args.resonance,
        verbose_failures=args.verbose_failures,
    )

    d4sig_center = gamma_v(xB_c, Q2_c, args.BeamEnergy) * sigu_center

    if not np.isfinite(d4sig_center) or d4sig_center <= 0.0:
        C_BC = float("nan")
    else:
        C_BC = avg / d4sig_center

    return {
        "N": N,
        "phi_lo": phi_lo,
        "phi_hi": phi_hi,
        "phi_c": phi_c,
        "xB_c": xB_c,
        "Q2_c": Q2_c,
        "t_c": t_c,
        "avg": avg,
        "d4sig_center": d4sig_center,
        "C_BC": C_BC,
        "n_phys": n_phys,
        "n_valid": int(np.sum(valid)),
        "n_xsec_failed": n_xsec_failed,
        "physical_fraction": n_phys / total_grid_points,
    }


def phi_scan_edges(args):
    if args.phi_edges is not None:
        return edge_pairs(args.phi_edges, "phi")
    return edge_pairs(np.linspace(args.phi[0], args.phi[1], args.phi_bins + 1), "phi")


def phi_scan_keys(include_bin_edges=False):
    keys = []
    if include_bin_edges:
        keys.extend([
            "t_bin_index", "Q2_bin_index", "xB_bin_index", "phi_bin_index",
            "xB_lo", "xB_hi", "Q2_lo", "Q2_hi",
            "t_lo_input", "t_hi_input", "t_lo_aao", "t_hi_aao",
        ])
    keys.extend([
        "N", "phi_lo", "phi_hi", "phi_c",
        "xB_c", "Q2_c", "t_c",
        "avg", "d4sig_center", "C_BC",
        "n_phys", "n_valid", "n_xsec_failed", "physical_fraction",
    ])
    return keys


def write_rows_csv(path, rows, keys):
    with open(path, "w") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in keys) + "\n")


def plot_phi_scan_rows(rows, N_list, title, plot_out):
    plt.figure(figsize=(8, 5))

    for N in N_list:
        sub = [r for r in rows if r["N"] == N]
        sub = sorted(sub, key=lambda r: r["phi_c"])

        phi_vals = [r["phi_c"] for r in sub]
        cbc_vals = [r["C_BC"] for r in sub]

        plt.scatter(phi_vals, cbc_vals, label=f"N={N}")
        plt.plot(phi_vals, cbc_vals, alpha=0.6)

    plt.xlabel(r"$\phi$ [deg]")
    plt.ylabel(r"$C_{\rm BC}$")
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_out, dpi=200)
    plt.close()


def compute_phi_scan_rows(args, phi_bins, N_list):
    rows = []

    for N in N_list:
        print(f"\nScanning phi with N={N}")

        for i, (phi_lo, phi_hi) in enumerate(phi_bins):
            print(f"  phi bin {i+1}/{len(phi_bins)}: {phi_lo:.3f} to {phi_hi:.3f}")

            result = compute_bin_centering(args, phi_lo, phi_hi, N)

            if result is not None:
                result["phi_bin_index"] = i
                rows.append(result)

    return rows


def run_phi_scan(args):
    N_list = args.N_list if args.N_list is not None else [args.N]

    phi_bins = phi_scan_edges(args)
    rows = compute_phi_scan_rows(args, phi_bins, N_list)

    if len(rows) == 0:
        sys.exit("ERROR: no valid phi-scan results.")

    keys = phi_scan_keys()
    write_rows_csv(args.csv_out, rows, keys)

    title = (
        rf"Bin-centering convergence: "
        rf"$x_B=[{args.xB[0]}, {args.xB[1]}]$, "
        rf"$Q^2=[{args.Q2[0]}, {args.Q2[1]}]$, "
        rf"$t=[{args.t[0]}, {args.t[1]}]$"
    )
    plot_phi_scan_rows(rows, N_list, title, args.plot_out)

    print(f"\nSaved plot: {args.plot_out}")
    print(f"Saved CSV:  {args.csv_out}")


def run_scan_all_bins(args):
    try:
        xB_bins = edge_pairs(args.xB_edges, "xB")
        Q2_bins = edge_pairs(args.Q2_edges, "Q2")
        t_bins = edge_pairs(args.t_edges, "t")
        phi_bins = phi_scan_edges(args)
    except ValueError as exc:
        sys.exit(f"ERROR: {exc}")

    os.makedirs(args.scan_dir, exist_ok=True)

    N_list = args.N_list if args.N_list is not None else [args.N]
    all_rows = []

    for it, (t_lo_user, t_hi_user) in enumerate(t_bins):
        t_for_aao = t_bin_for_aao(t_lo_user, t_hi_user, args.signed_t_edges)

        for iq, (Q2_lo, Q2_hi) in enumerate(Q2_bins):
            for ix, (xB_lo, xB_hi) in enumerate(xB_bins):
                scan_args = argparse.Namespace(**vars(args))
                scan_args.xB = [xB_lo, xB_hi]
                scan_args.Q2 = [Q2_lo, Q2_hi]
                scan_args.t = t_for_aao

                print(
                    "\nScanning 3D bin "
                    f"t#{it} Q2#{iq} xB#{ix}: "
                    f"xB=[{xB_lo:g}, {xB_hi:g}], "
                    f"Q2=[{Q2_lo:g}, {Q2_hi:g}], "
                    f"t=[{t_for_aao[0]:g}, {t_for_aao[1]:g}]"
                )

                rows = compute_phi_scan_rows(scan_args, phi_bins, N_list)

                for row in rows:
                    row.update({
                        "t_bin_index": it,
                        "Q2_bin_index": iq,
                        "xB_bin_index": ix,
                        "xB_lo": xB_lo,
                        "xB_hi": xB_hi,
                        "Q2_lo": Q2_lo,
                        "Q2_hi": Q2_hi,
                        "t_lo_input": t_lo_user,
                        "t_hi_input": t_hi_user,
                        "t_lo_aao": t_for_aao[0],
                        "t_hi_aao": t_for_aao[1],
                    })

                if rows:
                    plot_name = f"BC_scan_t{it:02d}_Q2{iq:02d}_xB{ix:02d}.png"
                    plot_path = os.path.join(args.scan_dir, plot_name)
                    t_label = (
                        rf"$t=[{t_lo_user:g}, {t_hi_user:g}]$"
                        if args.signed_t_edges
                        else rf"$-t=[{t_lo_user:g}, {t_hi_user:g}]$"
                    )
                    title = (
                        rf"Bin-centering convergence: "
                        rf"$x_B=[{xB_lo:g}, {xB_hi:g}]$, "
                        rf"$Q^2=[{Q2_lo:g}, {Q2_hi:g}]$, "
                        rf"{t_label}"
                    )
                    plot_phi_scan_rows(rows, N_list, title, plot_path)
                    print(f"Saved plot: {plot_path}")
                else:
                    print("  No valid phi-scan rows for this 3D bin; no plot written.")

                all_rows.extend(rows)

    if len(all_rows) == 0:
        sys.exit("ERROR: no valid scan-all-bins results.")

    keys = phi_scan_keys(include_bin_edges=True)
    csv_path = os.path.join(args.scan_dir, args.scan_csv)
    write_rows_csv(csv_path, all_rows, keys)
    print(f"\nSaved combined scan CSV: {csv_path}")


def run_quilt(args):
    """
    Build one quilt plot per t bin.

    Each tile is the existing C_BC-vs-phi scan for one 3D
    (xB, Q2, t) bin. Columns increase in xB and rows increase in Q2.
    C_BC is defined so downstream values should be divided by it.
    """
    try:
        xB_bins = edge_pairs(args.xB_edges, "xB")
        Q2_bins = edge_pairs(args.Q2_edges, "Q2")
        t_bins = edge_pairs(args.t_edges, "t")
        phi_bins = edge_pairs(args.phi_edges, "phi")
    except ValueError as exc:
        sys.exit(f"ERROR: {exc}")

    rows = []
    n_xB = len(xB_bins)
    n_Q2 = len(Q2_bins)

    for it, (t_lo_user, t_hi_user) in enumerate(t_bins):
        t_for_aao = t_bin_for_aao(t_lo_user, t_hi_user, args.signed_t_edges)
        t_label = (
            rf"$t=[{t_lo_user:g}, {t_hi_user:g}]$"
            if args.signed_t_edges
            else rf"$-t=[{t_lo_user:g}, {t_hi_user:g}]$"
        )

        fig_w = max(9.0, 1.65 * n_xB)
        fig_h = max(7.0, 1.45 * n_Q2)
        fig, axes = plt.subplots(
            n_Q2,
            n_xB,
            figsize=(fig_w, fig_h),
            sharex=True,
            sharey=True,
            squeeze=False,
        )

        print(f"\nBuilding quilt {it + 1}/{len(t_bins)} for {t_label}")

        for iq, (Q2_lo, Q2_hi) in enumerate(Q2_bins):
            row = n_Q2 - 1 - iq

            for ix, (xB_lo, xB_hi) in enumerate(xB_bins):
                ax = axes[row, ix]
                tile_rows = []

                tile_args = argparse.Namespace(**vars(args))
                tile_args.xB = [xB_lo, xB_hi]
                tile_args.Q2 = [Q2_lo, Q2_hi]
                tile_args.t = t_for_aao

                print(
                    "  tile "
                    f"Q2=[{Q2_lo:g}, {Q2_hi:g}], "
                    f"xB=[{xB_lo:g}, {xB_hi:g}]"
                )

                for iphi, (phi_lo, phi_hi) in enumerate(phi_bins):
                    result = compute_bin_centering(tile_args, phi_lo, phi_hi, args.N)
                    if result is None:
                        continue

                    result.update({
                        "t_bin_index": it,
                        "Q2_bin_index": iq,
                        "xB_bin_index": ix,
                        "phi_bin_index": iphi,
                        "xB_lo": xB_lo,
                        "xB_hi": xB_hi,
                        "Q2_lo": Q2_lo,
                        "Q2_hi": Q2_hi,
                        "t_lo_input": t_lo_user,
                        "t_hi_input": t_hi_user,
                        "t_lo_aao": t_for_aao[0],
                        "t_hi_aao": t_for_aao[1],
                    })
                    tile_rows.append(result)
                    rows.append(result)

                if tile_rows:
                    tile_rows = sorted(tile_rows, key=lambda r: r["phi_c"])
                    ax.plot(
                        [r["phi_c"] for r in tile_rows],
                        [r["C_BC"] for r in tile_rows],
                        marker="o",
                        markersize=2.5,
                        linewidth=0.9,
                    )
                else:
                    ax.text(
                        0.5,
                        0.5,
                        "no data",
                        ha="center",
                        va="center",
                        transform=ax.transAxes,
                        fontsize=8,
                        color="0.45",
                    )
                    ax.set_facecolor("0.96")

                ax.grid(True, alpha=0.25, linewidth=0.5)
                ax.tick_params(labelsize=7, length=2)

                if row == n_Q2 - 1:
                    ax.set_xlabel(rf"$x_B$ [{xB_lo:g}, {xB_hi:g}]", fontsize=8)
                if ix == 0:
                    ax.set_ylabel(rf"$Q^2$ [{Q2_lo:g}, {Q2_hi:g}]", fontsize=8)

        fig.supxlabel(r"$\phi$ [deg]", fontsize=11)
        fig.supylabel(r"$C_{\rm BC}$; rows increase in $Q^2$", fontsize=11)
        fig.suptitle(
            rf"Bin-centering quilts, {t_label}, N={args.N}",
            fontsize=13,
        )
        fig.tight_layout(rect=[0.03, 0.03, 1.0, 0.96])

        plot_out = f"{args.quilt_prefix}_{it:02d}.png"
        fig.savefig(plot_out, dpi=200)
        plt.close(fig)
        print(f"Saved quilt: {plot_out}")

    if len(rows) == 0:
        sys.exit("ERROR: no valid quilt results.")

    keys = [
        "t_bin_index", "Q2_bin_index", "xB_bin_index", "phi_bin_index",
        "xB_lo", "xB_hi", "Q2_lo", "Q2_hi",
        "t_lo_input", "t_hi_input", "t_lo_aao", "t_hi_aao",
        "N", "phi_lo", "phi_hi", "phi_c",
        "xB_c", "Q2_c", "t_c",
        "avg", "d4sig_center", "C_BC",
        "n_phys", "n_valid", "n_xsec_failed", "physical_fraction",
    ]

    with open(args.quilt_csv, "w") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in keys) + "\n")

    print(f"\nSaved quilt CSV: {args.quilt_csv}")

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Average d4sigma over the physical part of a 4D kinematic bin "
            "using a midpoint Riemann sum."
        )
    )

    parser.add_argument("--xB", nargs=2, type=float, default=None,
                        metavar=("LO", "HI"), help="xB bin edges")
    parser.add_argument("--Q2", nargs=2, type=float, default=None,
                        metavar=("LO", "HI"), help="Q2 bin edges [GeV^2]")
    parser.add_argument("--t", nargs=2, type=float, default=None,
                        metavar=("LO", "HI"), help="t bin edges [GeV^2]")
    parser.add_argument("--phi", nargs=2, type=float, default=None,
                        metavar=("LO", "HI"), help="phi bin edges [degrees]")
    parser.add_argument("--BeamEnergy", type=float, required=True,
                        metavar="E", help="Beam energy [GeV]")
    parser.add_argument("--N", type=int, default=10,
                        help="Points per dimension, so total grid is N^4")
    parser.add_argument(
        "--exe",
        type=str,
        default="/work/clas12/storyf/SF_analysis_software/aao_print/aao_gen/aao_norad/build/aao_xsec",
        help="Path to aao_xsec executable"
    )
    parser.add_argument("--theory", type=int, default=5)
    parser.add_argument("--channel", type=int, default=1)
    parser.add_argument("--resonance", type=int, default=0)
    parser.add_argument("--workers", type=int, default=None,
                        help="Number of parallel workers")
    parser.add_argument("--verbose-failures", action="store_true",
                        help="Print one warning per failed aao_xsec call.")
    
    parser.add_argument("--scan-phi", action="store_true",
                    help="Scan bin-centering correction vs phi.")
    parser.add_argument("--scan-all-bins", action="store_true",
                        help="Run --scan-phi style convergence plots for every xB/Q2/t edge bin.")
    parser.add_argument("--phi-bins", type=int, default=18,
                        help="Number of phi bins for --scan-phi.")
    parser.add_argument("--N-list", nargs="+", type=int, default=None,
                        help="List of N values for convergence scan, e.g. --N-list 4 6 8 10.")
    parser.add_argument("--plot-out", type=str, default="bc_phi_scan.png",
                        help="Output plot name for --scan-phi.")
    parser.add_argument("--csv-out", type=str, default="bc_phi_scan.csv",
                        help="Output CSV name for --scan-phi.")
    parser.add_argument("--scan-dir", type=str, default="BC_scans",
                        help="Output directory for --scan-all-bins plots and CSV.")
    parser.add_argument("--scan-csv", type=str, default="BC_scans.csv",
                        help="Combined CSV filename inside --scan-dir for --scan-all-bins.")
    parser.add_argument("--quilt", action="store_true",
                        help="Make C_BC-vs-phi quilt plots over xB, Q2, and t edge arrays.")
    parser.add_argument("--xB-edges", nargs="+", type=float, default=None,
                        help="xB edge list for --quilt.")
    parser.add_argument("--Q2-edges", nargs="+", type=float, default=None,
                        help="Q2 edge list for --quilt [GeV^2].")
    parser.add_argument("--t-edges", nargs="+", type=float, default=None,
                        help="t edge list for --quilt. By default positive values are treated as -t.")
    parser.add_argument("--phi-edges", nargs="+", type=float, default=None,
                        help="phi edge list for --quilt [degrees].")
    parser.add_argument("--signed-t-edges", action="store_true",
                        help="In --quilt mode, use --t-edges as signed physical t values instead of -t.")
    parser.add_argument("--quilt-prefix", type=str, default="bc_quilt_t",
                        help="Output PNG prefix for --quilt. One file per t bin is written.")
    parser.add_argument("--quilt-csv", type=str, default="bc_quilt.csv",
                        help="Output CSV name for --quilt.")

    args = parser.parse_args()

    if not os.path.isfile(args.exe):
        sys.exit(f"ERROR: executable not found: {args.exe}")
    if not os.access(args.exe, os.X_OK):
        sys.exit(f"ERROR: file is not executable: {args.exe}")

    if args.N <= 0:
        sys.exit("ERROR: --N must be positive")

    if args.quilt:
        missing = [
            name for name, value in [
                ("--xB-edges", args.xB_edges),
                ("--Q2-edges", args.Q2_edges),
                ("--t-edges", args.t_edges),
                ("--phi-edges", args.phi_edges),
            ]
            if value is None
        ]
        if missing:
            sys.exit(f"ERROR: --quilt requires {' '.join(missing)}")
        run_quilt(args)
        return

    if args.scan_all_bins:
        missing = [
            name for name, value in [
                ("--xB-edges", args.xB_edges),
                ("--Q2-edges", args.Q2_edges),
                ("--t-edges", args.t_edges),
            ]
            if value is None
        ]
        if args.phi_edges is None and args.phi is None:
            missing.append("--phi-edges or --phi")
        if missing:
            sys.exit(f"ERROR: --scan-all-bins requires {' '.join(missing)}")
        run_scan_all_bins(args)
        return

    missing = [
        name for name, value in [
            ("--xB", args.xB),
            ("--Q2", args.Q2),
            ("--t", args.t),
            ("--phi", args.phi),
        ]
        if value is None
    ]
    if missing:
        sys.exit(f"ERROR: single-bin and --scan-phi modes require {' '.join(missing)}")

    N = args.N
    total_grid_points = N**4

    print(f"Grid: {N}^4 = {total_grid_points} midpoint cells")

    if args.scan_phi:
        run_phi_scan(args)
        return

    # Build enclosing rectangular midpoint grid
    xB_pts = midpoint_grid(args.xB[0], args.xB[1], N)
    Q2_pts = midpoint_grid(args.Q2[0], args.Q2[1], N)
    t_pts = midpoint_grid(args.t[0], args.t[1], N)
    phi_pts = midpoint_grid(args.phi[0], args.phi[1], N)

    all_points = list(itertools.product(xB_pts, Q2_pts, t_pts, phi_pts))

    # Keep only kinematically physical points
    physical_points = [
        p for p in all_points
        if physical_mask(p[0], p[1], p[2], args.BeamEnergy)
    ]

    n_phys = len(physical_points)

    if n_phys == 0:
        sys.exit(
            "ERROR: no physical midpoint cells in this bin. "
            "Try increasing --N or check the bin edges."
        )

    physical_fraction = n_phys / total_grid_points

    print(f"Physical midpoint cells: {n_phys} / {total_grid_points}")
    print(f"Physical fraction:       {physical_fraction:.6f}")

    worker = partial(
        call_aao_xsec,
        exe=args.exe,
        ebeam=args.BeamEnergy,
        theory=args.theory,
        channel=args.channel,
        resonance=args.resonance,
        verbose_failures=args.verbose_failures,
    )

    n_workers = args.workers or multiprocessing.cpu_count()
    print(f"Dispatching {n_phys} executable calls to {n_workers} workers ...")

    with multiprocessing.Pool(processes=n_workers) as pool:
        sigu_values = pool.map(worker, physical_points)

    sigu_arr = np.array(sigu_values, dtype=float)
    n_xsec_failed = int(np.sum(~np.isfinite(sigu_arr)))

    d4sig_values = np.array([
        gamma_v(xB, Q2, args.BeamEnergy) * sigu
        for (xB, Q2, t, phi), sigu in zip(physical_points, sigu_arr)
    ])

    valid = np.isfinite(d4sig_values) & (d4sig_values > 0.0)

    n_nan = np.sum(~np.isfinite(d4sig_values))
    n_nonpos = np.sum(np.isfinite(d4sig_values) & (d4sig_values <= 0.0))

    if n_nan > 0:
        sys.stderr.write(
            f"WARNING: {n_nan}/{n_phys} physical points returned NaN/inf.\n"
        )

    if n_nonpos > 0:
        sys.stderr.write(
            f"WARNING: {n_nonpos}/{n_phys} physical points had non-positive d4sigma "
            f"and were excluded from the positive-cross-section average.\n"
        )

    if not np.any(valid):
        sys.exit("ERROR: no valid positive cross-section points.")

    # -----------------------------------------------------------------------
    # Average over the physical sampled region.
    #
    # This is the bin-centering denominator:
    #
    #   <sigma>_phys = integral_phys sigma dV / integral_phys dV
    #
    # Since all midpoint cells have equal dV, the dV cancels.
    # -----------------------------------------------------------------------

    avg = np.mean(d4sig_values[valid])

    # Enclosing rectangular volume
    delta_xB = args.xB[1] - args.xB[0]
    delta_Q2 = args.Q2[1] - args.Q2[0]
    delta_t = args.t[1] - args.t[0]
    delta_phi = args.phi[1] - args.phi[0]

    rect_vol = delta_xB * delta_Q2 * delta_t * delta_phi

    # Estimated physical volume from cell counting
    cell_vol = rect_vol / total_grid_points
    physical_vol = n_phys * cell_vol

    integral = avg * physical_vol

    # Optional center value and bin-centering factor
    # Use the geometric centroid of the physically allowed sampled bin.
    # Since the midpoint grid is uniform, every accepted cell has equal weight.

    valid_points = np.array([
        p for p, is_valid in zip(physical_points, valid)
        if is_valid
    ], dtype=float)

    xB_c = np.mean(valid_points[:, 0])
    Q2_c = np.mean(valid_points[:, 1])
    t_c  = np.mean(valid_points[:, 2])

    phi_rad = np.deg2rad(valid_points[:, 3])
    phi_c = np.rad2deg(
        np.arctan2(
            np.mean(np.sin(phi_rad)),
            np.mean(np.cos(phi_rad))
        )
    ) % 360.0

    if not physical_mask(xB_c, Q2_c, t_c, args.BeamEnergy):
        print("WARNING: physical-bin centroid is not physically allowed.")

    sigu_center = call_aao_xsec(
        (xB_c, Q2_c, t_c, phi_c),
        exe=args.exe,
        ebeam=args.BeamEnergy,
        theory=args.theory,
        channel=args.channel,
        resonance=args.resonance,
        verbose_failures=args.verbose_failures,
    )

    gamma_center = gamma_v(xB_c, Q2_c, args.BeamEnergy)
    d4sig_center = gamma_center * sigu_center

    if not np.isfinite(d4sig_center) or d4sig_center <= 0.0:
        bin_centering_factor = float("nan")
    else:
        bin_centering_factor = avg / d4sig_center

    print("\nResults:")
    print(f"  <d4sigma/dQ2 dxB dt dphi> over physical bin = {avg:.6e}")
    print(f"  Enclosing rectangular volume                 = {rect_vol:.6e}")
    print(f"  Estimated physical volume                    = {physical_vol:.6e}")
    print(f"  Physical fraction                            = {physical_fraction:.6f}")
    print(f"  Integrated cross section over physical bin   = {integral:.6e}")
    print(f"  Physical points                              = {n_phys} / {total_grid_points}")
    print(f"  Valid positive xsec points                   = {np.sum(valid)} / {n_phys}")
    print(f"  Failed aao_xsec calls                        = {n_xsec_failed} / {n_phys}")

    print("\nBin center:")
    print(f"  xB   = {xB_c:.6g}")
    print(f"  Q2   = {Q2_c:.6g}")
    print(f"  t    = {t_c:.6g}")
    print(f"  phi  = {phi_c:.6g}")

    print(f"xB midpoint  = {0.5*(args.xB[0]+args.xB[1]):.6f}")
    print(f"xB centroid  = {xB_c:.6f}")

    print(f"Q2 midpoint  = {0.5*(args.Q2[0]+args.Q2[1]):.6f}")
    print(f"Q2 centroid  = {Q2_c:.6f}")

    print(f"t midpoint   = {0.5*(args.t[0]+args.t[1]):.6f}")
    print(f"t centroid   = {t_c:.6f}")

    print("\nBin-centering:")
    print(f"  d4sigma(center)             = {d4sig_center:.6e}")
    print(f"  <d4sigma>_bin               = {avg:.6e}")
    print(f"  C_BC = <d4sigma>_bin/d4sigma(center) = {bin_centering_factor:.6e}")
    print("  Apply downstream as: centered value = bin-averaged value / C_BC")

if __name__ == "__main__":
    main()
