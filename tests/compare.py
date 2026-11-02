"""Compare results in two output files of the same type.

Usage: compare.py ref new
Handles pconf.csv, tm.csv, and the GF.RES, A_hf.RES, and B_hf.RES tables from pdtm.
Exits with an error message if the results do not match.
"""
import csv
import sys

# Largest allowed difference for energies and matrix elements, in a.u.
TOL = 1e-6

# g-factors are printed to 5 decimal places.
TOL_GF = 2e-5

# Largest allowed fractional difference for rates and hyperfine constants.
# These span many orders of magnitude and are printed to about 5 digits.
REL_TOL = 2e-4

# Energies given in cm^-1 are converted to a.u. to use the same tolerance.
CM_PER_AU = 219474.63


def read_levels(rows, path):
    """Return results from the rows of pconf.csv."""
    results = []
    for row in rows:
        if row['converged'].strip() != 'True':
            sys.exit(f"FAIL: level {row['n']} did not converge in {path}")
        label = f"level {row['n']} {row['configuration']} {row['term']}"
        results.append((label, float(row['energy_au']), TOL))
        results.append((f'{label} excitation energy', float(row['energy_cm']) / CM_PER_AU, TOL))
        # J=0 levels have no g-factor.
        if row['gf'].strip():
            results.append((f'{label} g-factor', float(row['gf']), TOL_GF))
    return results


def read_transitions(rows):
    """Return results from the rows of tm.csv."""
    results = []
    for row in rows:
        state_one = f"{row['state_one_configuration']} {row['state_one_term']}"
        state_two = f"{row['state_two_configuration']} {row['state_two_term']}"
        label = f"{row['n1']} -> {row['n2']} {state_one} -> {state_two} {row['operator']}"
        rate = float(row['transition_rate'])
        # The overall sign of a matrix element is arbitrary, so compare magnitudes.
        results.append((label, abs(float(row['matrix_element_value'])), TOL))
        results.append((f'{label} transition energy', float(row['transition_energy_cm']) / CM_PER_AU, TOL))
        results.append((f'{label} rate', rate, abs(rate) * REL_TOL))
    return results


def read_table(path):
    """Return results from a table of the form: n, configuration, term, value, energy."""
    with open(path) as f:
        lines = [line.split() for line in f if line.strip()]

    results = []
    is_gf = 'gfactor' in lines[0]
    for fields in lines[1:]:
        label = f"level {fields[0]} {' '.join(fields[1:-3])} {fields[-3]}"
        value = float(fields[-2])
        results.append((label, value, TOL_GF if is_gf else abs(value) * REL_TOL))
    return results


def read_results(path):
    """Return a list of (label, value, tolerance), one per compared value."""
    if path.endswith('.RES'):
        return read_table(path)

    with open(path) as f:
        rows = list(csv.DictReader(f))
    if rows and 'operator' in rows[0]:
        return read_transitions(rows)
    return read_levels(rows, path)


def main():
    ref_path, new_path = sys.argv[1], sys.argv[2]
    ref = read_results(ref_path)
    new = read_results(new_path)

    if not ref:
        sys.exit(f'FAIL: no results found in {ref_path}')
    if len(ref) != len(new):
        sys.exit(f'FAIL: {len(ref)} values in {ref_path}, {len(new)} values in {new_path}')

    # Labels include the configurations and terms, so this also catches wrong labels.
    for (ref_label, _, _), (new_label, _, _) in zip(ref, new):
        if ref_label != new_label:
            sys.exit(f"FAIL: labels differ: '{ref_label}' in reference, '{new_label}' in new results")

    # Written this way so that a NaN counts as a failure.
    num_agree = len(ref)
    for (label, a, tol), (_, b, _) in zip(ref, new):
        if not abs(a - b) <= tol:
            print(f'{label}: ref = {a:.8g}, new = {b:.8g}, diff = {abs(a - b):.2e}, allowed = {tol:.2e}')
            num_agree -= 1
    if num_agree != len(ref):
        sys.exit(f'FAIL: {num_agree}/{len(ref)} values agree in {new_path}')

    print(f'PASS: {num_agree}/{len(ref)} values agree')


if __name__ == '__main__':
    main()
