# aao_print
modification of https://github.com/JeffersonLab/aao_gen by codeX, to print xsec value

## Build and plot reduced cross section vs phi

All paths below are relative to this repository root:

```bash
cd aao_gen/aao_norad
make aao_xsec_bin
```

The wrapper program is:

```bash
aao_gen/aao_norad/aao_xsec
```

It prints the reduced cross section `sigma0` used by the `aao_norad` MC
generator before multiplying by the virtual-photon flux, Jacobian, or phase
space normalization.

Single-point example:

```bash
./aao_xsec \
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
cd aao_gen/aao_norad
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
cd aao_gen/aao_norad
make clean
```
