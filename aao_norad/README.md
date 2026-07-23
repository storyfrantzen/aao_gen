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
  is the control mode for validating the direct modes over an identical
  physical domain;
- `3`: uniform sampling in `(1/Q2, xB, -t, phi*)`. This is the recommended
  mode for globally unweighted samples because it retains direct analysis
  bounds while preferentially proposing the low-`Q2` region where the cross
  section is largest;
- `4`: bin-conditional uniform sampling in `(1/Q2, xB, -t, phi*)`. This is the
  high-statistics production mode for generating a fixed number of unweighted
  Born events inside each analysis bin.

Mode 4 requires two additional lines after the analysis bounds:

```text
condition_phase_space W_min y_max
flat_index iq2 ixb it iphi
```

With `condition_phase_space=0`, mode 4 generates the full rectangular bin; the
configured `Q2`, `W`, and `y` values are retained as metadata but are not
generation cuts. With `condition_phase_space=1`, it generates only the
intersection of the bin and those phase-space cuts. The choice and bin identity
are written to `aao_norad.norm`.

For modes 1, 3, and 4, the generator reconstructs

```text
nu = Q2/(2 M xB)
E' = E - nu
W2 = M2 + Q2*(1/xB - 1)
```

and obtains the pion center-of-mass angle from the proposed `t`. In mode 1, the
density used by acceptance-rejection is

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

Mode 3 uses `u = 1/Q2` in place of `Q2`. Its density includes the additional
coordinate derivative

```text
|dQ2/du| = Q2^2
```

so that the acceptance-rejection density and integration volume are

```text
Gamma_v * sigma0 * Q2^3 /
    (8 M xB^2 E E' |q*| |p_pi*|)

2*pi * Delta(1/Q2) * Delta(xB) * Delta(-t) * Delta(phi radians)
```

Here each `Q2` denotes the full kinematic variable (Q^2). Mode 3 changes only
the proposal efficiency; after acceptance-rejection, its globally unweighted
events follow the same physical distribution as modes 1 and 2.

Mode 4 uses the same proposal density and Jacobian as mode 3, but its proposal
rectangle is one analysis bin. Events are unweighted *within* a stratum.
Different strata generally have different integrated cross sections, so when
strata are combined the physical weight for every event from stratum `i` is

```text
w_i = sig_sum_i / N_i
```

This is also the weight to use when filling generated-bin and reconstructed-bin
histograms for a response matrix. Migration is retained because an event keeps
its source stratum from the manifest while its reconstructed bin is determined
after detector simulation. Mode 4 supplies the response columns for truth
events inside the generated domain. Full-bin generation preserves events on
both sides of the `W/y` boundary so downstream cuts can measure those
migrations. Phase-space-conditioned generation is more efficient when those
events are definitely unwanted, but a global mode-3 or dedicated guard-region
sample is then needed to estimate feed-in and perform an independent closure
test.

AAO's acceptance step permits `mcall > 1`. If
`r = signr/sigr_max`, it emits

```text
mcall = floor(r) + Bernoulli(r - floor(r))
```

events, so `E[mcall | r] = r` even when `r > 1`. The emitted records still have
equal unit event weight; they are multiplicity copies of the same proposed hard
kinematics, with independently generated rotation and decay variables.
`mcall_max > 1` therefore indicates that the preliminary envelope was exceeded,
not that the distribution correction failed. A large or frequent multiplicity
does increase correlations and can overshoot the requested event count on the
last proposal. Mode-4 normalization uses the actual emitted count in
`sig_sum/events`; increasing `fmcall` remains useful when fewer duplicate
kinematics are desired.

Every run also writes `aao_norad.kin`, with one row per LUND event:

```text
event Q2 xB minus_t phi_deg signr jacobian
```

## Sampling validation

Build the generator and compare uniform-`Q2` direct sampling, bounded legacy
sampling, and optimal inverse-`Q2` direct sampling:

```bash
make
python3 validate_analysis_sampling.py \
  --physics-model 1 --events 20000 --replicas 4
```

The study compares repeated `sig_sum` estimates pairwise across all three modes
and performs pairwise two-sample KS tests for the unweighted `Q2`, `xB`, `-t`,
and `phi` distributions. It also reports aggregate proposal efficiencies. It
writes `validation_sampling/validation.json` and `shape_comparison.pdf`.
Because three independent runs are made per repeat, the total event count is
`3 * events * replicas`. Model 5 also requires the normal MAID table
configuration, such as `MAID07_TBL`.

Each successful `aao_norad` run writes `aao_norad.norm` in addition to the LUND,
`.kin`, `.out`, and `.sum` outputs. Use the `sig_sum` value from
`aao_norad.norm` as the Born integrated cross section when normalizing
generated-event histograms.

## Bin-conditional Born production

`bin_conditional.py` converts an analysis JSON configuration into one or more
mode-4 generator invocations per selected physical analysis bin. Strata use the
canonical ID `sNNNNN`, and independent invocations use `gNNNN`, for example
`s01440__g0001`. The script snapshots the configuration and writes a manifest
containing the bin bounds, indices, replica number, seed, generator revision,
and output location:

```bash
python3 bin_conditional.py prepare \
  --config /path/to/configs/analysis/rgk/6.535.json \
  --output born_rgk_conditional \
  --events-per-bin 100000 \
  --replicas 2 \
  --physics-model 5
```

The default is full-bin generation. Add `--condition-phase-space` to restrict
each job to the bin's intersection with the configured `Q2`, `W`, and `y`
selection.

After building `aao_norad` and configuring the MAID tables in the usual way,
run one manifest entry with:

```bash
python3 bin_conditional.py run born_rgk_conditional/manifest.json \
  --flat-index 1440 \
  --replica-index 1 \
  --executable build/aao_norad
```

The `run` command is intended to be called once per farm array job. It checks
the `.norm`, `.kin`, and LUND products before moving them into the stratum
output directory. Each output JSON records `sig_sum`,
`event_weight_microbarn`, requested and actual event counts, any final event
overshoot, `mcall_max`, whether the multiplicity correction was used, proposal
count, and event yield per proposal.

After every invocation has completed, finalize the campaign:

```bash
python3 bin_conditional.py finalize born_rgk_conditional/manifest.json
```

This writes `campaign_weights.json` and `.tsv`. Replicas estimate the same
physical stratum cross section, so their cross sections are averaged rather
than added. The combined estimator pools their Monte Carlo integration
proposals,

```text
combined_sig_sum =
    sum(ntries_r * sig_sum_r) / sum(ntries_r)
```

and every event in the pooled stratum receives

```text
pooled_event_weight = combined_sig_sum / sum(events_r).
```

Consequently, the weights of all replicas together sum once—not once per
replica—to the combined stratum cross section. The finalized weights file is
the input to stratum-preserving OSG repacking and downstream event weighting.
