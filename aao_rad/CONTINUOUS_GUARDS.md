# Milestone 2e continuous radiative guards

`radiative_continuous_guards.py` tests a simpler mode-4 core
representation: one continuous axis-aligned box per final-LUND analysis
stratum in AAO's normalized native proposal coordinates.

This learner remains a diagnostic and does not itself write LUND events.
`radiative_mode4.py` and AAO sampling mode 4 now consume its recipes through
an exact core-plus-global-tail mixture. `production_ready` remains false until
the mode-4 implementation and tuning are frozen and tested on fresh replicas.

## Coordinate definitions

The target stratum uses final-LUND generator-truth coordinates:

```text
(Q2_observed, xB_observed, minus_t_observed, phi_observed)
```

The fitted core uses:

```text
r_u
r_ep
u_gamma
hadron_cosine_base
hadron_phi_base
```

These variables lie in `[0,1]` under the legacy proposal.
`hadron_phi_base` is periodic. Native `intreg` sampling remains unchanged,
and the first continuous core leaves photon-angle local randoms and external
loss variables unrestricted.

A box in these coordinates has a directly calculable normalized volume. That
is required for the future exact mode-4 proposal correction. The hard-parent
coordinates remain useful diagnostics but are not assumed to be directly
sampleable.

## Weighted joint-box fit

For each analysis stratum, the learner uses training trials that finish in
that stratum and their fixed-trial cross-section contributions.

The optimizer begins with the complete unit box. It repeatedly considers
moving each lower or upper face inward to a weighted empirical-CDF
breakpoint. A breakpoint is an actual recorded proposal-coordinate value;
between breakpoints the empirical coverage is unchanged.

Every candidate move is evaluated jointly with all other current bounds. It
is accepted only if the retained weighted cross section remains at or above
`--target-core-fraction`. Candidate moves are ranked by normalized volume
saved per weighted contribution removed. Starting from the full box avoids a
greedy-growth failure in which a correlated point lies outside two faces and
moving either face by itself captures nothing.

The fitted face positions are not treated as final statistical limits. The
comparison evaluates continuous outward padding values after fitting.

## Sparse-stratum regularization

Training strata are classified as:

- `learned`: enough rows and effective sample size to fit locally;
- `learned_low_support`: a local sample exists but is not well supported;
- `no_training_contribution`: no positive training contribution was seen.

Well-supported strata use their own weighted distribution. Low-support
strata blend their local distribution with populated adjacent final-LUND
strata. The local fraction is:

```text
local_fraction = ESS / (ESS + regularization_ESS)
```

Empty strata use populated adjacent strata. If none exist within the
configured Manhattan radius, the learner uses one frozen global fallback
box fitted to all training contributions inside the analysis.

The observed phi-bin index wraps when finding neighbors.

## Continuous padding scan

The default base paddings are:

```text
0
0.0025
0.005
0.01
0.02
```

All are expressed in the normalized `[0,1]` proposal coordinates. The
effective padding is support adaptive:

```text
learned:
    effective_padding = base_padding

low or empty:
    effective_padding =
        base_padding * (1 + (1 - local_fraction))
```

Thus an empty stratum receives twice the base padding, while a low-support
stratum interpolates continuously between one and two times the base
padding.

The development recommendation is the smallest total core volume among
padding candidates whose aggregate held-out weighted coverage reaches
`--minimum-validation-coverage`. This is not a final production choice.
Tail ESS, dominant tail weights, per-channel coverage, and material-stratum
coverage must also be inspected.

## Run the RGK comparison

Run from `aao_rad`:

```bash
python3 radiative_continuous_guards.py compare \
  --config ../../../configs/analysis/rgk/6.535.json \
  --training-survey \
    survey_rgk_balanced_replica000 \
    survey_rgk_balanced_replica001 \
    survey_rgk_balanced_replica002 \
  --validation-survey \
    survey_rgk_balanced_replica003 \
    survey_rgk_balanced_replica004 \
  --target-core-fraction 0.995 \
  --neighbor-radius 1 \
  --regularization-ess 5 \
  --minimum-axis-width 0.0001 \
  --minimum-training-rows 10 \
  --minimum-training-ess 5 \
  --minimum-validation-coverage 0.98 \
  --iteration 0 \
  --generator-revision GENERATOR_COMMIT_USED_FOR_SURVEYS \
  --output migration_rgk_balanced_continuous_guards_iteration000
```

Replace `GENERATOR_COMMIT_USED_FOR_SURVEYS` with the exact revision used to
build the executable. In `tcsh`, capture it before starting the replicas with:

```tcsh
set generator_revision = `git rev-parse HEAD`
```

and pass `--generator-revision $generator_revision`.

These commands intentionally use only balanced-proposal surveys in one
comparison. Do not pool legacy-proposal and balanced-proposal replicas: their
fixed-trial contributions are individually valid, but the guard tools require
one frozen proposal definition per campaign.

The default applies the configured final-LUND `Q2` and `W` requirements and
does not apply an analysis-level `y` cut. Add `--apply-y-max` only after the
analysis selection is intentionally frozen.

The output directory is immutable and must not already exist.

## Plot the comparison

```bash
python3 radiative_continuous_guards.py plot \
  --comparison \
    migration_rgk_balanced_continuous_guards_iteration000/continuous_guard_comparison.json \
  --output migration_rgk_balanced_continuous_guards_iteration000_plots
```

## Artifacts

The comparison writes:

- `continuous_guard_comparison.json`: aggregate, support-class,
  material-stratum, and native-`intreg` held-out coverage for every padding;
  normalized core volumes; tail statistics; pass/fail counts; and the
  development recommendation;
- `continuous_guard_recipes.json`: exact raw boxes, circular-phi
  representation, neighbor provenance, ESS regularization fractions, and
  padding scales for all analysis strata;
- `continuous_guard_strata.csv`: one row per
  `(final-LUND analysis stratum, padding candidate)`;
- SHA-256 sidecars for every artifact;
- `continuous_guard_comparison.pdf`: coverage, coverage-versus-volume,
  residual-tail ESS, and per-stratum pass/fail plots.

The recipe stores a raw box plus an exact padding rule rather than duplicating
all padded edges. A mode-4 implementation can reconstruct every candidate
deterministically.

## Validation discipline

Replicas 3 and 4 are development data because earlier diagnostics already
influenced this design. Use them to choose and refine the continuous
procedure. Once the algorithm, regularization, and padding are frozen, assess
the result on fresh replicas that were not used in any design decision.

The mode-4 proposal retains:

```text
g_i = alpha * g_core,i + (1 - alpha) * g_global
```

with `1 - alpha` strictly positive. The global component supplies complete
legacy support. Mode 4 evaluates the sum of both mixture densities at every
trial before unweighting. See `RADIATIVE_MODE4.md`.

## Tests

Run the focused tests with:

```bash
python3 -m unittest -v test_radiative_continuous_guards.py
```

Run every radiative workflow test with:

```bash
make test
```
