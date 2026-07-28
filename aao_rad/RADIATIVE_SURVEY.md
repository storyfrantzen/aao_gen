# Milestone 1 radiative survey

Milestone 1 adds a diagnostic path to `aao_rad`. It does not yet implement
bin-conditional radiative production.

The two radiative run-control modes are:

- `0`: unchanged legacy unweighted generation;
- `1`: deterministic fixed-trial survey.

Mode `4` remains reserved for the future bin-conditional unweighted radiative
generator described in `GUARD_PROPOSAL_WORKFLOW.md`.

## Survey input

Existing legacy input files remain valid and select mode `0` at end of file.
To request a survey manually, append four records after the last legacy
record:

```text
1          ! fixed-trial survey mode
1000000    ! exact number of unrestricted proposal trials
371001     ! explicit nonzero random seed
0          ! explicit nonnegative replica ID
```

The survey currently supports the neutral-pion, four-particle configuration:
`epirea=1` and `npart=4`.

The event-count and `fmcall` records remain in the legacy part of the input
for compatibility, but survey termination does not use either one. It runs
exactly the requested number of proposal trials and does not perform the
legacy maximum scan or acceptance-rejection.

## Recommended command

From `aao_rad`, first build the executable:

```bash
make
```

Then run a survey through the checked wrapper:

```bash
python3 radiative_survey.py run \
  --executable build/aao_rad \
  --input aao_input.inp \
  --output survey_rgk_replica000 \
  --trials 1000000 \
  --seed 371001 \
  --replica 0
```

The output directory must not already contain survey artifacts. This avoids
silently mixing or overwriting replicas.

To revalidate an existing output:

```bash
python3 radiative_survey.py validate \
  --directory survey_rgk_replica000
```

## Outputs and semantics

`aao_rad.survey.csv` has schema `aao-rad-survey-v1`. Each CSV row is an
internally valid unrestricted proposal. Proposal trials rejected before the
internal integrand exists are absent from the CSV but remain zero
contributions through the fixed `ntries` denominator in `aao_rad.norm`.

The CSV records:

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
trial count alone does not guarantee a precise integral.

## Deterministic replicas

A replica means an independent run at identical physical settings with a
different explicit seed. Keep the trial count and all legacy physics inputs
fixed, change only the seed, and write each replica to a separate directory.
Repeating a run with the same compiler, executable, input, trial count, and
seed should reproduce the survey CSV byte for byte.
