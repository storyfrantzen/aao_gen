# aao_print
modification of https://github.com/JeffersonLab/aao_gen by codeX, to print xsec value

## Build products are not tracked

Generated files under `aao_rad/build/` and `aao_norad/build/` are intentionally
not tracked. Build the needed executable once on the target machine before
launching jobs, then point production scripts at the executable in `build/`.
Do not have many batch jobs invoke a wrapper that runs `make` in the shared
repository; on NFS this can create many stale `.nfs*` files while build outputs
are open.

## Build and plot reduced cross section vs phi

All paths below are relative to this repository root:

```bash
cd aao_norad
make aao_xsec_bin
```

The compiled single-point cross-section executable is:

```bash
aao_norad/build/aao_xsec
```

It prints the reduced cross section `sigma0` used by the `aao_norad` MC
generator before multiplying by the virtual-photon flux, Jacobian, or phase
space normalization.

Single-point example:

```bash
aao_norad/build/aao_xsec \
  -xB 0.3 \
  -Q2 2.0 \
  -t -0.2 \
  -phi 30 \
  -BeamEnergy 10.6 \
  -phys 5 \
  -epirea 1
```

Default model options are:

```text
phys = 5
epirea = 1
resonance = 0
BeamEnergy = 10.6
```

The kinematic options `-xB`, `-Q2`, `-t`, `-phi`, and `-BeamEnergy` must be
given explicitly.

For bin-centering or other batch scans, build `aao_norad/build/aao_xsec` once
before submission and pass that path directly to the analysis driver, for
example:

```bash
python3 analysis/run_analysis.py bin-centering \
  --config configs/analysis/rgk/6.535.json \
  --exe external/aao_gen/aao_norad/build/aao_xsec \
  --output /path/to/C_BC_part.npz \
  --N 4 --bin-chunks 648 --bin-chunk-index 0
```

Avoid invoking the make-running helper wrapper inside large batch workflows.
Use the compiled `build/aao_xsec` executable directly instead.

## Build event generators

Build the non-radiative generator:

```bash
cd aao_norad
make
```

This creates:

```text
aao_norad/build/aao_norad
aao_norad/build/aao_xsec
```

Build the radiative generator:

```bash
cd aao_rad
make
```

This creates:

```text
aao_rad/build/aao_rad
```

## Generator normalization sidecars

The event generators write machine-readable normalization sidecars at the end of
each successful run:

```text
aao_gen/aao_norad/aao_norad.norm
aao_gen/aao_rad/aao_rad.norm
```

These files preserve the run-level integrated cross section needed to normalize
accepted unweighted LUND events. Downstream analyses should use `sig_sum` as the
generator integrated cross section. The older `.sum` files still contain the
same value in the line `Integrated cross section = <sig_int> <sig_sum>`.

To make the phi plot with ROOT:

```bash
cd aao_norad
root -l -q 'plot_xsec_phi.C(0.3,2.0,-0.2,10.6)'
```

The arguments are:

```text
plot_xsec_phi(xB, Q2, t, BeamEnergy)
```

Optional scan range and output prefix:

```bash
root -l -q 'plot_xsec_phi.C(0.3,2.0,-0.2,10.6,0,360,5,"my_phi_plot")'
```

This writes:

```text
aao_gen/aao_norad/my_phi_plot.dat
aao_gen/aao_norad/my_phi_plot.png
aao_gen/aao_norad/my_phi_plot.pdf
```

Clean build products before committing:

```bash
cd aao_norad
make clean
```
