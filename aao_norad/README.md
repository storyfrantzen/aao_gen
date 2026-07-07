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
```

Each successful `aao_norad` run writes `aao_norad.norm` in addition to the LUND,
`.out`, and `.sum` outputs. Use the `sig_sum` value from `aao_norad.norm` as the
Born integrated cross section when normalizing generated-event histograms.
