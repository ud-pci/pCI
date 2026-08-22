# Fe16+ Example

Fe16+ using a pure CI method with a Dirac-Fock basis set.

## Overview

This example computes energy levels and E1 transition matrix elements for Fe16+.

## Key config.yml parameters

To adapt this example for a different ion, the main parameters to change are:

- `atom.name`, `atom.isotope` - element symbol and mass number
- `basis.method` - `dirac-fock` or `b-splines`
- `basis.orbitals.core` / `valence` - core and valence orbital list for the HFD step
- `basis.orbitals.order` - groups of orbitals for sequential HFD runs (one HFD input file per group)
- `add.ref_configs` - reference configurations for even and odd parity
- `add.basis_set` - active orbital set (e.g. `12spdfg` means orbitals up to n=12, l=4)
- `add.orbitals.active` - active orbitals and their occupation number ranges
- `add.excitations` - which excitation levels to include (single, double, triple)
- `conf.odd.J` / `conf.even.J` - total angular momentum for each parity block
- `conf.odd.J_selection` - if True, only keeps levels with the specified J value
- `conf.odd.num_energy_levels` / `conf.even.num_energy_levels` - number of levels to compute
- `dtm.TM.from` / `dtm.TM.to` - parity and level range for transition matrix elements

## Output files

After running, the main results are:

- `pconf.csv` - energy levels for a given parity, one row per level
- `tm.csv` - transition matrix elements and rates between two parity blocks, one row per transition
- `FINAL.RES` - final energy levels from CI in an overview format

## Directory structure

```
Fe16+/
    config.yml          - input configuration
    ref/                - reference results to compare against
        basis/
            HFD1.INP      - HFD input for first orbital group (from basis.py)
            HFD2.INP      - HFD input for second orbital group (from basis.py)
            HFD3.INP      - HFD input for third orbital group (from basis.py)
            BASS.INP      - BASS input (from basis.py)
            BAS_INFO.RES  - basis set summary (from bas_info)
            bass.in       - bass input (from basis.py)
        even0/
            CONF.INP  - CONF input (from ci.py)
            ADD.INP   - ADD input (from ci.py)
            FINAL.RES - final energy levels from CI in an overview format
            pconf.csv - energy levels
        odd1/         - same as even0/
        tm/
            dtm.in    - DTM input (from dtm.py)
            E1.RES    - E1 matrix elements
            tm.csv    - transition matrix elements and rates
```

## Running the example

Run each script from this directory in order:

```bash
# 1. Generate basis set
python basis.py
# -> creates basis/

# 2. Generate configuration lists
python ci.py
# -> creates even0/ and odd1/

# 3. Run CI (on cluster via job scripts, or interactively)

# 4. Set up DTM input files
python dtm.py
# -> creates tm/

# 5. Run DTM (on cluster)
```
