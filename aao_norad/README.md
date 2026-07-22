# aao_norad

## Build products

Files under `build/` are generated locally and are not tracked. Build the
needed target once before running local or batch workflows.

## How to run `aao_norad`:
```
cd aao_norad
make
./build/aao_norad < aao_norad.inp
```

The single-point cross-section executable used by bin-centering scans is built
with:

```
cd aao_norad
make aao_xsec_bin
./build/aao_xsec -xB 0.3 -Q2 2.0 -t -0.2 -phi 30 -BeamEnergy 10.6
```

For large batch scans, point jobs directly at `build/aao_xsec`. Avoid invoking
the `aao_xsec` helper wrapper from every job because it can run `make` in this
shared directory.

The content of input file `aao_norad.inp` configures the generator parameters.
Refer to `input.txt` for documentation on the content:
```
The input file must have the following format. Please see the script aaorun 
and example input file aao_norad_input.inp included in this package for an 
example of how to run the program:

phys         Physics model: 1=A0, 4=MAID98, 5=MAID2000
flag_ehel    0= no polarized electron, 1=polarized electron
npart        number of particles in BOS banks: 2=(e-,h+), 3=(e-,h+,h0)
epirea       final state hadron: 1=pi0, 3=pi+
ebeam        incident electron beam energy in GeV
q2min q2max  minimum and maximum Q^2 limits in GeV^2
epmin epmax  minimum and maximum scattered electron energy limits in GeV
nmax         number of output events
fmcall       factor to adjust the maximum cross section, used in M.C. selection
boso         1=bos output, 0=no bos output
seed         0=clock/farm seed, otherwise an explicit integer seed
```

The historical input ends after `seed` and remains valid.  To select a new
sampling mode, append:

```text
sampling_mode
xbmin xbmax
minus_t_min minus_t_max
phi_min phi_max
```

The modes are:

- `0`: historical uniform sampling in `(1/Q2, E', cos(theta*), phi*)`;
- `1`: uniform sampling in `(Q2, xB, -t, phi*)` over the appended bounds;
- `2`: historical coordinates with the appended analysis bounds applied.  This
  is the control mode for validating mode 1 over an identical physical domain.

For mode 1, the generator reconstructs

```text
nu = Q2/(2 M xB)
E' = E - nu
W2 = M2 + Q2*(1/xB - 1)
```

and obtains the pion center-of-mass angle from the proposed `t`.  The density
used by acceptance-rejection is

```text
Gamma_v * sigma0 * Q2 /
    (8 M xB^2 E E' |q*| |p_pi*|)
```

where `Q2` denotes the kinematic variable Q squared, rather than `Q` squared a
second time.  The associated integration volume is

```text
2*pi * Delta(Q2) * Delta(xB) * Delta(-t) * Delta(phi radians)
```

The factor `2*pi` is the unobserved electron azimuth.  Invalid points in the
rectangular analysis-coordinate domain have zero integrand and remain counted
in `ntries`.

Every run also writes `aao_norad.kin`, with one row per LUND event:

```text
event Q2 xB minus_t phi_deg signr jacobian
```

## Sampling validation

Build the generator and compare direct sampling with bounded legacy sampling:

```bash
make
python3 validate_analysis_sampling.py \
  --physics-model 1 --events 20000 --replicas 4
```

The study compares replicated `sig_sum` estimates and performs two-sample KS
tests for the unweighted `Q2`, `xB`, `-t`, and `phi` distributions.  It writes
`validation_sampling/validation.json` and `shape_comparison.pdf`.  Model 5 also
requires the normal MAID table configuration, such as `MAID07_TBL`.

Each successful `aao_norad` run writes `aao_norad.norm` in addition to the LUND,
`.kin`, `.out`, and `.sum` outputs. Use the `sig_sum` value from
`aao_norad.norm` as the Born integrated cross section when normalizing
generated-event histograms.
