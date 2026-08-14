# Global radiative mode 3

Radiative mode 3 is the global, equal-event-weight alternative to the
bin-conditional mode 4 workflow. It keeps AAO's radiative proposal variables
and replaces only two inefficient legacy coordinates:

```text
legacy:  (1/Q_l^2, E',       photon variables, cos(theta*), phi_h)
direct:  (1/Q_l^2, x_l,      photon variables, -t_h,        phi_h)
```

Here `Q_l^2` and `x_l` are formed from the incoming and outgoing electron at
the radiative vertex, while `-t_h` is the transfer to the hadronic system after
the internally radiated photon is included. The generated target is defined
from the final LUND electron and proton, after external energy loss.

## Why it is a mixture

The direct proposal is efficient inside the requested padded analysis box, but
by itself it cannot describe radiative feed-in from every legacy-coordinate
point. Mode 3 therefore draws from

```text
q_mix = f_direct q_direct + (1 - f_direct) q_legacy
```

and multiplies each trial integrand by `q_legacy/q_mix`. The legacy component
has nonzero density throughout the original AAO domain and preserves support.
The default `f_direct=0.75` spends three quarters of proposals on the efficient
coordinates and one quarter on unrestricted legacy coverage.

The exact density ratio inside direct support is

```text
q_direct/q_legacy =
  [Delta(E') 2 M x_l^2 / (Q_l^2 Delta(x_l))]
  [2 |d(-t_h)/d cos(theta*)| / Delta(-t_h)]
  [360 degrees / Delta(phi_h)]
```

Outside direct support it is zero. The acceptance step uses stochastic
multiplicity, so a missed envelope maximum does not bias the sample; it creates
diagnostic duplicates instead. Production stops at exactly the requested event
count, and the LUND record remains the legacy four-particle, five-line format.

## Required validation

Before production, run two fixed-trial calibrations over identical settings:

1. the candidate mixture, normally `--direct-fraction 0.75`;
2. a legacy-proposal control with `--direct-fraction 0`.

`compare-calibrations` checks that their integrated cross sections agree within
the requested relative and statistical tolerances. The candidate calibration
also supplies the production envelope. Run a small production pilot before a
large campaign and inspect `maximum_mcall`, `duplicate_events`, proposal yield,
and the exact LUND line count.

## Simple automatic-envelope production

The fixed-envelope workflow above is useful when one wants a reusable,
validated envelope and a quantitative proposal-closure test.  It is not a
mathematical requirement for unbiased generation.  Mode-3 operation `2`
instead performs a short prescan with the actual mode-3 proposal, sets

```text
sigr_max = fmcall * maximum_observed_integrand
```

discards the prescan statistics, and immediately starts production.  As in
legacy AAO, an underestimated maximum does not bias the generated sample:
stochastic multiplicity emits `floor(sigr/sigr_max)` events plus one with the
remaining fractional probability.  Envelope crossings therefore appear as
duplicate events and reduce efficiency rather than changing the distribution.

`submit_aao_rad_mode3_native_rga10604_200M.csh` uses this path with a pure
analysis-coordinate proposal (`direct_fraction=1`), 10,000 prescan trials per
job, and the padded RGA domain.  This deliberately forgoes the unrestricted
legacy component: the configured direct hard-coordinate box must therefore be
wide enough for the radiative feed-in that the study intends to retain.

## RGA 10.604-GeV padded domain

With `configs/analysis/rga/10.604.json` and the default padding fraction 0.035,
mode 3 targets final-LUND coordinates in

```text
Q2       0.6675  to 10.8325 GeV2
xB       0.02725 to 0.72275
-t       0.02315 to 2.06685 GeV2
phi      0       to 360 degrees
W        greater than 2 GeV
electron momentum greater than 2 GeV
```

No `y_max` is applied because the RGA analysis configuration does not define
one. Padding expands `Q2`, `xB`, and `-t`; it does not relax the configured `W`
or final-electron-momentum selections.

The minimum radiated-photon energy is inherited from the legacy input unless
`--minimum-photon-energy` is supplied.  This setting is frozen into the
calibration metadata; production with a different threshold is rejected.
The RGA submission helper uses `input/rga10604_mode3_legacy.inp`, which
preserves the polarized-electron flag, target geometry, and other non-campaign
settings from the earlier RGA legacy-production script.

## Campaign artifacts

`radiative_mode3.py` creates immutable inputs, a manifest, a task table,
per-run JSON/normalization diagnostics, and LUND outputs. For a 100-million
event campaign with 5,000 events per file it creates 20,000 tasks and four
`lund/chunk_NNNN` directories of 5,000 files each. This keeps each directory
below the 10,000-job type-2 OSG submission limit. All events in the finalized
global campaign share the `pooled_event_weight_microbarn` recorded in
`campaign_weights.json`; no per-stratum weights are needed.

Generated SWIF wrappers use `$SWIF_JOB_WORK_DIR` for temporary files, falling
back to `$TMPDIR`, an explicitly requested scratch root, and finally the job's
current directory. A login-node `/scratch/$USER` path is not assumed to exist
on worker nodes.

Use `python3 radiative_mode3.py --help` or a subcommand's `--help` for the full
command interface.
