# Sr I Example

Sr I using CI+all-order and CI+MBPT methods.

## Overview

This example computes energy levels and E1, E2, E3, M1, M2, M3 transition matrix elements for Sr, using configuration interaction with all-order and MBPT corrections.

`for_portal: True` is used here to generate a complete set of transition matrix elements. In ci.py, it collects the unique J values from both parities and generates CI directories for every parity/J combination (even0, even1, odd0, odd1) instead of just the single J per parity specified in the config. In dtm.py, it generates all four TM combinations (tm_even0_odd1, tm_odd0_even1, tm_even0_even1, tm_odd0_odd1) with the appropriate E1/E2/E3/M1/M2/M3 operators for each, instead of only the pair specified in the config. This is useful any time you want a complete set of transition matrix elements, not only for portal output.

The two-method setup is used to calculate uncertainties by taking the difference between the two methods for each energy level and matrix element.

## Key config.yml parameters

To adapt this example for a different atom, the main parameters to change are:

- `atom.name`, `atom.isotope` - element symbol and mass number
- `atom.code_method` - list of methods to run (e.g. `[ci+all-order]` for single method)
- `basis.orbitals.core` / `valence` - core and valence orbital list for the HFD step
- `add.ref_configs` - reference configurations for even and odd parity
- `add.basis_set` - active orbital set (e.g. `17spdfg` means orbitals up to n=17, l=4)
- `add.orbitals.active` - active orbitals and their occupation number ranges
- `add.excitations` - which excitation levels to include (single, double, triple)
- `conf.odd.J` / `conf.even.J` - total angular momentum for each parity block
- `conf.odd.num_energy_levels` / `conf.even.num_energy_levels` - number of levels to compute

For the polarizability calculation (ine.py / pine program):

- `ine.parity` - parity of the level to compute the polarizability for
- `ine.N0`, `ine.N2` - record numbers of the initial (X0) and final (X2) levels in CONF0.XIJ; N0 and N2 are the same for a scalar polarizability, different for a transition polarizability
- `ine.rhs`, `ine.lhs` - operators for the right- and left-hand sides of the Sternheimer equation: 1 = H_pnc, 2 = E1 length gauge, 3 = H_am (lhs only), 4 = E1 velocity gauge (lhs only), 5 = E2
- `ine.field_type` - `static` for static polarizability (omega=0), `dynamic` for frequency-dependent, or `static, dynamic` for both
- `ine.wavelength_range`, `ine.step_size` - wavelength range in nm and step size for dynamic polarizability
- `ine.tm_dir` - (optional) override the TM directory to copy DTM.INT from; derived automatically from `for_portal` and conf J values if not set

## Output files

After running, the main results are:

- `pconf.csv` - energy levels for a given parity, one row per level
- `tm.csv` - transition matrix elements and rates between two parity blocks, one row per transition
- `pine.csv` - polarizabilities from pine, one row per wavelength; columns are wavelength, alpha_0, alpha_1, alpha_2, alpha_tot, and alpha(m) for each m from 0 to J. This example computes the static and dynamic E1 polarizability of the 5s2 1S0 ground state from wavelength=1000 to 1001 nm in steps of 0.5 nm.
- `FINAL.RES` - final energy levels from CI in an overview format

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

# 6. Set up polarizability directory (after DTM is done)
python ine.py
# -> creates ine_5s2_1S0/ with ine.in, ine.qs, CI files, and DTM_RPA.INT (or DTM.INT)

# 7. Run INE (on cluster)
```

