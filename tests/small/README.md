# Small Fe16+ regression test

A small version of `examples/Fe16+` that checks the hfd, bass, pbasc, pconf, and pdtm chain in under a minute.
It covers E1, M2, and E3 transitions, then g-factors and hyperfine constants of the odd levels.
It also repeats pconf for each parity in the Emu regime with `Nemu= 100`.

| Parity | Jm | Levels kept | Configurations | Determinants |
|--------|----|-------------|----------------|--------------|
| even   | 0  | all J       | 247            | 2614         |
| odd    | 1  | J=1 only    | 454            | 5419         |

It differs from the example in these ways, all set in `config.yml` except `Nemu`, which `run.sh` adds:

- basis of 7spdfg instead of 12spdfg
- `nmax: 7` instead of 24
- single excitations only, since keeping doubles gives 8552 even configurations
- transitions to even levels 1-5 instead of level 1 only, since M2 and E3 are forbidden from J=1 to J=0
- an added DM run for the odd levels

## Running

```bash
bash run.sh       # 1 rank
bash run.sh 4     # 4 ranks
CONF_USE_BIT_DET=FALSE bash run.sh 4     # integer determinants
```

pconf uses bitstring determinants by default here, so the last line covers the integer determinant path.

Executables are taken from `<repo>/bin`. Set `BIN` to use a different directory.
Each run works in a new temp directory and prints its path, so the output can be inspected.
The script exits with an error if a program fails, a step takes over 5 minutes, a level does not converge, or a result differs from the reference.
Energies and matrix elements must agree to 1e-6 a.u., g-factors to 2e-5, and rates and hyperfine constants to a fraction of 2e-4.
Configuration and term labels must match exactly.
Matrix elements are compared by magnitude, since their overall sign is arbitrary.

## Reference values

`ref/*/pconf.csv` holds energies from pconf at commit 7e03eb7, run on 4 ranks.
`ref/tm/tm.csv` and `ref/dm/*.RES` hold pdtm results from the same commit.
`ref/*_emu/pconf.csv` holds the Emu energies from pconf at the same commit.
Commit 7e03eb7 predates the FormH neighbor list change.
Its bass hangs, so the basis came from bass with the `conf_init.f90` fix.

Only replace these files when a change in energies is expected and understood.
To replace them, run the test, then copy `even0/pconf.csv`, `odd1/pconf.csv`, `even0_emu/pconf.csv`, `odd1_emu/pconf.csv`, `tm/tm.csv`, and `dm/*.RES` from the printed temp directory into `ref/`.

## Regenerating inputs

`inputs/` came from running `py-lib/basis.py`, `py-lib/ci.py`, and `py-lib/dtm.py` on `config.yml`, with `bin_directory` set to your executables.
