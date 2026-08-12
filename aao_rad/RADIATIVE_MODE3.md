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

## Campaign artifacts

`radiative_mode3.py` creates immutable inputs, a manifest, a task table,
per-run JSON/normalization diagnostics, and LUND outputs. For a 100-million
event campaign with 5,000 events per file it creates 20,000 tasks and four
`lund/chunk_NNNN` directories of 5,000 files each. This keeps each directory
below the 10,000-job type-2 OSG submission limit. All events in the finalized
global campaign share the `pooled_event_weight_microbarn` recorded in
`campaign_weights.json`; no per-stratum weights are needed.

Use `python3 radiative_mode3.py --help` or a subcommand's `--help` for the full
command interface.
