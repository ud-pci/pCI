# Sr I Example

Sr I using CI+all-order and CI+MBPT methods.

## Overview

This example computes energy levels and E1, E2, E3, M1, M2, M3 transition matrix elements for Sr, using configuration interaction with all-order and MBPT corrections.

`for_portal: True` is used here to generate a complete set of transition matrix elements. In ci.py, it collects the unique J values from both parities and generates CI directories for every parity/J combination (even0, even1, odd0, odd1) instead of just the single J per parity specified in the config. In dtm.py, it generates all four TM combinations (tm_even0_odd1, tm_odd0_even1, tm_even0_even1, tm_odd0_odd1) with the appropriate E1/E2/E3/M1/M2/M3 operators for each, instead of only the pair specified in the config. This is useful any time you want a complete set of transition matrix elements, not only for portal output.

The two-method setup is used to calculate uncertainties by taking the difference between the two methods for each energy level and matrix element.

## Key config.yml parameters

To adapt this example for a different atom, the main parameters to change are:

- `atom.name`, `atom.isotope` — element symbol and mass number
- `atom.code_method` — list of methods to run (e.g. `[ci+all-order]` for single method)
- `basis.orbitals.core` / `valence` — core and valence orbital list for the HFD step
- `add.ref_configs` — reference configurations for even and odd parity
- `add.basis_set` — active orbital set (e.g. `17spdfg` means orbitals up to n=17, l=4)
- `add.orbitals.active` — active orbitals and their occupation number ranges
- `add.excitations` — which excitation levels to include (single, double, triple)
- `conf.odd.J` / `conf.even.J` — total angular momentum for each parity block
- `conf.odd.num_energy_levels` / `conf.even.num_energy_levels` — number of levels to compute

## Output files

After running, the main results are:

- `pconf.csv` — energy levels for a given parity, one row per level
- `tm.csv` — transition matrix elements and rates between two parity blocks, one row per transition
- `FINAL.RES` — final energy levels from CI in an overview format

## Directory structure

```
Sr/
    config.yml          - input configuration
    ref/                - reference results to compare against
        ci+all-order/
            basis/
                HFD.INP       - HFD input (from basis.py)
                HFD.RES       - HFD output (from hfd)
                BASS.INP      - BASS input (from basis.py)
                BASS.RES      - BASS output (from bass)
                BAS_INFO.RES  - basis set summary (from bas_info)
                bas_wj.in     - bas_wj input (from basis.py)
                spl.in        - bspl input (from basis.py)
                inf.aov       - all-order input (from basis.py)
                inf.vw        - all-order input (from basis.py)
            even0/
                CONF.INP  - CONF input (from ci.py)
                ADD.INP   - ADD input (from ci.py)
                FINAL.RES - CI convergence output
                pconf.csv - energy levels
            even1/  odd0/  odd1/    - same as even0/
            tm_even0_odd1/
                dtm.in    - DTM input (from dtm.py)
                MBPT.INP  - MBPT input (from dtm.py)
                E1.RES  E3.RES  M2.RES  TM.RES  tm.csv
            tm_odd0_even1/
                dtm.in  MBPT.INP
                E1.RES  E3.RES  M2.RES  TM.RES  tm.csv
            tm_even0_even1/
                dtm.in  MBPT.INP
                E2.RES  M1.RES  M3.RES  TM.RES  tm.csv
            tm_odd0_odd1/
                dtm.in  MBPT.INP
                E2.RES  M1.RES  M3.RES  TM.RES  tm.csv
        ci+second-order/
            basis/            - same files as ci+all-order/basis/
            even0/  even1/  odd0/  odd1/    - same as ci+all-order
            tm_even0_odd1/  tm_odd0_even1/  tm_even0_even1/  tm_odd0_odd1/
                - same files as ci+all-order tm dirs
```

## Running the example

Run each script from this directory in order:

```bash
# 1. Generate basis set
python basis.py
# -> creates ci+all-order/basis/ and ci+second-order/basis/

# 2. Generate configuration lists
python ci.py
# -> creates ci+all-order/even0/, ci+all-order/even1/, ci+all-order/odd0/, ci+all-order/odd1/, etc.

# 3. Run CI (on cluster via job scripts, or interactively)

# 4. Set up DTM input files
python dtm.py
# -> creates tm_even0_odd1/, tm_odd0_even1/, etc. in each method directory

# 5. Run DTM (on cluster)
```

