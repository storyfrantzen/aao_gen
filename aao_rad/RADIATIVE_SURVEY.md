# Milestone 1 radiative survey

Milestone 1 adds a diagnostic path to `aao_rad`. It does not yet implement
bin-conditional radiative production.

The two radiative run-control modes are:

- `0`: unchanged legacy unweighted generation;
- `1`: deterministic fixed-trial survey.

Mode `4` is now a development bin-conditional unweighted prototype. Its
core-plus-legacy-tail proposal and workflow are documented in
`RADIATIVE_MODE4.md`; fresh-replica production validation is still required.

Mode `1` now has two proposal choices:

- `legacy`: the original AAO radiative proposal;
- `balanced`: a full-support mixture of the legacy proposal and a component
  that allocates trials uniformly across the configured analysis-coordinate
  bins.

The balanced proposal is the recommended choice for learning guards. The
legacy proposal remains available for reproducibility and normalization
cross-checks.

## Balanced proposal

For each fixed trial, the balanced survey first chooses a mixture component.
The default probabilities are 25% legacy and 75% analysis-balanced.

The analysis-balanced component independently:

1. chooses a configured `Q2` bin uniformly, then samples `Q2` uniformly
   inside that bin;
2. chooses an `xB` bin uniformly, samples inside it, and maps `(Q2,xB)` to
   the pre-outgoing-loss electron energy;
3. chooses a `minus_t` bin uniformly, samples inside it, and uses the exact
   affine hard-vertex mapping from `-t` to `cos(theta*)`;
4. chooses a `phi_deg` bin uniformly and samples inside it.

All photon variables, external-loss variables, and radiative-channel choices
retain their original AAO sampling. The chosen bins therefore refer to the
pre-outgoing-loss leptonic `Q2,xB` and hard hadronic `-t,phi`, not directly to
the final observed coordinates. Radiation, physical-boundary rejection, and
the final external loss can migrate a trial to another observed bin or make
it invalid. The final observed bins will consequently be more even, but not
exactly uniform.

The mixture is corrected exactly:

```text
q_mix = alpha * q_legacy + (1 - alpha) * q_balanced

corrected integrand =
    legacy-transformed integrand * q_legacy / q_mix
```

For the default `alpha=0.25`, the correction is positive and no larger than
`1/alpha = 4`. The nonzero legacy component guarantees support anywhere the
legacy generator has support. The survey validator independently recomputes
this density ratio for every recorded proposal.

Raw occupancy counts diagnose whether trials were allocated more evenly.
Cross sections and guard coverage must still use the corrected
`trial_xsec_*` contributions; raw counts are not cross-section weights.

## Survey input

Existing legacy input files remain valid and select mode `0` at end of file.
The checked Python wrapper should normally construct the survey trailer. A
legacy-proposal survey trailer has five records after the legacy input:

```text
1          ! fixed-trial survey mode
1000000    ! exact number of unrestricted proposal trials
371001     ! explicit nonzero random seed
0          ! explicit nonnegative replica ID
0          ! legacy survey proposal
```

The balanced trailer additionally freezes the mixture fraction and all four
sets of bin edges. `radiative_survey.py` builds that trailer from the analysis
JSON, copies the exact JSON bytes into the output directory, and records its
SHA-256 digest. This avoids maintaining a long trailer by hand.

The balanced survey currently supports the neutral-pion, four-particle
configuration: `epirea=1` and `npart=4`.

The event-count and `fmcall` records remain in the legacy part of the input
for compatibility, but survey termination does not use either one. It runs
exactly the requested number of proposal trials and does not perform the
legacy maximum scan or acceptance-rejection.

## Recommended command

From `aao_rad`, first build the executable:

```bash
make
```

Then run a balanced survey through the checked wrapper:

```bash
python3 radiative_survey.py run \
  --executable build/aao_rad \
  --input aao_input.inp \
  --output survey_rgk_balanced_replica000 \
  --trials 1000000 \
  --seed 371001 \
  --replica 0 \
  --proposal balanced \
  --config ../../../configs/analysis/rgk/6.535.json \
  --legacy-fraction 0.25
```

The output directory must not already contain survey artifacts. This avoids
silently mixing or overwriting replicas.

To revalidate an existing output:

```bash
python3 radiative_survey.py validate \
  --directory survey_rgk_balanced_replica000
```

## Outputs and semantics

New output uses schema `aao-rad-survey-v2`. Readers remain compatible with
the earlier `aao-rad-survey-v1` legacy surveys. Each CSV row is an
internally valid unrestricted proposal. Proposal trials rejected before the
internal integrand exists are absent from the CSV but remain zero
contributions through the fixed `ntries` denominator in `aao_rad.norm`.

The CSV records:

- the selected proposal component and, for balanced-component trials, the
  four target-bin indices;
- the exact `q_legacy/q_mix` density correction;
- base proposal coordinates and radiative channel;
- incoming, pre-outgoing-loss, and final electron energies;
- internal-photon energy and angles;
- leptonic and hard-vertex coordinates;
- final electron and proton four-vectors;
- observed analysis coordinates derived from those final particles;
- internal and observed integrands;
- per-trial cross-section contributions;
- incoming and outgoing external-loss fractions;
- validity status.

`candidate_status` means:

- `0`: complete final candidate and valid observed coordinates;
- `1`: outgoing external loss placed the electron below `ep_min`;
- `2`: the post-loss call to `missm` did not produce a physical candidate;
- `3`: final particles did not define valid observed coordinates.

An invalid final candidate keeps its internal contribution but has zero
observed contribution. This distinction lets later guard studies quantify
the cross section that survives all steps needed to assign a final observed
analysis stratum.

The cross-section estimators are:

```text
survey_internal_sig_sum =
    sum(trial_xsec_internal_microbarn) / ntries

survey_observed_sig_sum =
    sum(trial_xsec_observed_microbarn) / ntries
```

The denominator is the exact number of unrestricted proposals, not the number
of CSV rows and not the number of valid final candidates.

Survey mode intentionally writes no LUND events. The validator requires
`aao_rad.lund` to be empty and independently recomputes observed
`Q2`, `xB`, `-t`, `phi`, `W`, and `y` from the recorded final four-vectors.
Its JSON summary also reports fixed-trial standard errors, importance-sampling
effective sample sizes, and the largest single-trial contribution. These are
important because radiative survey weights can have long tails; a large raw
trial count alone does not guarantee a precise integral. Balanced summaries
also report target-axis counts, final-observed-axis counts, and the number of
occupied joint analysis strata.

## Deterministic replicas

A replica means an independent run at identical physical settings with a
different explicit seed. Keep the trial count and all legacy physics inputs
and proposal settings fixed, change only the seed, and write each replica to
a separate directory.
Repeating a run with the same compiler, executable, input, trial count, and
seed should reproduce the survey CSV byte for byte.

## Milestone 2: learn and validate guards

`radiative_guards.py` consumes survey directories without loading their full
CSVs into memory. It uses the configured observed `(Q2, xB, -t, phi)` bins,
applies the configured observed `Q2` and `W` selection, and aggregates
**cross-section contribution**, its square, and its maximum by proposal cell
and radiative channel. By default, `y_observed` is recorded but no lower or
upper `y` cut is applied.

The recommended split for the five RGK pilot replicas is:

- replicas `0`, `1`, and `2`: training;
- replicas `3` and `4`: independent held-out validation.

From `aao_rad`, learn iteration zero with:

```bash
python3 radiative_guards.py learn-guards \
  --config ../../../configs/analysis/rgk/6.535.json \
  --survey \
    survey_rgk_balanced_replica000 \
    survey_rgk_balanced_replica001 \
    survey_rgk_balanced_replica002 \
  --generator-revision GENERATOR_COMMIT_USED_FOR_SURVEYS \
  --output guard_rgk_balanced_iteration000
```

Replace the revision placeholder with the revision actually used to build the
surveyed executable. If the option is omitted, the learner records the current
checkout but marks that provenance as an assumption.

Once the RGK upper-`y` selection is finalized, it can be enabled explicitly
with `--apply-y-max`. That reads `phase_space.y_max` from the analysis
configuration and records the active value in the frozen manifest. Held-out
validation always inherits the manifest's selection policy; it cannot
silently introduce or remove a `y` cut.

The default proposal partition is deliberately modest:

```text
r_u=8, r_ep=8, u_gamma=6,
hadron_cosine_base=6, hadron_phi_base=8
```

`intreg` is always kept separate. `hadron_phi_base` is periodic, so dilation
wraps across its zero/one boundary. Photon-angle and external-loss base
variables remain unrestricted in this first guard. This preserves their full
legacy support while avoiding an excessively sparse initial grid.

Validate the frozen manifest only on replicas that were not used to learn it:

```bash
python3 radiative_guards.py validate-guards \
  --manifest guard_rgk_balanced_iteration000/guard_manifest.json \
  --survey \
    survey_rgk_balanced_replica003 \
    survey_rgk_balanced_replica004 \
  --output guard_rgk_balanced_iteration000_validation
```

Print a compact campaign summary with:

```bash
python3 radiative_guards.py summarize-coverage \
  --manifest guard_rgk_balanced_iteration000/guard_manifest.json \
  --validation \
    guard_rgk_balanced_iteration000_validation/guard_validation.json
```

The summary retains the original per-stratum pass counts and also reports a
`weighted_coverage` section. That section uses the exact fixed-trial
cross-section contributions—not raw stratum counts—to provide:

- training and held-out core and tail fractions;
- fixed-trial cross sections, standard errors, and effective sample sizes;
- coverage after one additional dilation step;
- the inside-analysis training/validation difference and z-score;
- held-out coverage separated by training-support and pass/fail status;
- coverage among the strata carrying 50%, 90%, 95%, and 99% of the cross
  section;
- radiative-channel fractions.

It also prints the frozen analysis selection and proposal partition and refuses
to combine a validation artifact with a different manifest hash.

The learner writes:

- `guard_manifest.json` and its SHA-256 sidecar;
- `training_coverage.csv`, with one row for every configured analysis stratum;
- `training_cells.csv`, preserving `S`, `S2`, and `M` for every occupied
  stratum/channel/proposal cell.

The validator writes:

- `guard_validation.json` and its SHA-256 sidecar;
- `validation_coverage.csv`, covering strata seen in training or validation;
- `validation_cells.csv`, labeling every occupied held-out cell as seed, core,
  one-more-dilation, or tail.

Each CSV also receives its own SHA-256 sidecar.

Output directories must not already exist. This is intentional: learned
manifests are immutable campaign inputs, and a later iteration must receive a
new directory and iteration number.

Each stratum's core is stored compactly as a union of weighted seed cells plus
an explicit number of axis-neighbor dilation steps. The complement of that
derived core is the complete legacy proposal tail, and the manifest assigns it
a nonzero mixture probability. Sparse strata are labeled
`learned_low_support`; they are not silently treated as physically empty.

The validator reports training and held-out core fractions, tail fractions,
fixed-trial cross sections and standard errors, importance-sampling effective
sample sizes, radiative-channel fractions, the largest tail cell, and the
fraction obtained after one additional dilation step. It also reports the
training-versus-holdout cross-section difference in combined-standard-error
units. `summarize-coverage` lists the lowest-coverage strata first.

The validation command exits with status `2` when the configured held-out
coverage threshold is not met. That means the proposed guard needs another
learning iteration; it does not mean that the survey files or cross-section
calculation failed.

This is still a learning artifact: `production_ready` is false. The manifest
does not alter radiative sampling or LUND output until milestone 3 implements
the exact core-plus-tail proposal correction and unweighting.

## Milestones 2b–2e: migrations and guard geometry

The independent proposal-cell guard is deliberately general but can be sparse.
`radiative_migrations.py` provides a complementary, physically organized
diagnostic that maps Born-like hard-parent strata and radiative channels into
the final-LUND analysis strata. The follow-on
`radiative_guard_geometries.py` diagnostic compares support-adaptive,
channel-marginalized core guards without changing AAO sampling. See
[`RADIATIVE_MIGRATIONS.md`](RADIATIVE_MIGRATIONS.md) for the coordinate
definitions, farm commands, artifacts, and interpretation.

`radiative_continuous_guards.py` then replaces whole-cell dilation with
weighted joint boxes in the normalized native proposal coordinates and
compares continuous padding values. See
[`CONTINUOUS_GUARDS.md`](CONTINUOUS_GUARDS.md).
