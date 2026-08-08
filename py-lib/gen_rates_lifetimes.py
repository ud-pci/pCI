"""
Script to generate transition rates and lifetimes from energy levels and matrix elements.

Multipole transition probabilities A_wv(Tk) from upper state v to lower state w:

                    2.02613e18
    A_wv(E1) =  ------------------  x S(E1)
                (2*J_v+1)*lambda^3

                    1.11995e18
    A_wv(E2) =  ------------------  x S(E2)
                (2*J_v+1)*lambda^5

                    3.14441e17
    A_wv(E3) =  ------------------  x S(E3)
                (2*J_v+1)*lambda^7

                    2.69735e13
    A_wv(M1) =  ------------------  x S(M1)
                (2*J_v+1)*lambda^3

                    1.49097e13
    A_wv(M2) =  ------------------  x S(M2)
                (2*J_v+1)*lambda^5

                    4.18610e12
    A_wv(M3) =  ------------------  x S(M3)
                (2*J_v+1)*lambda^7

where J_v is the total angular momentum of the upper state, lambda is the transition
wavelength in Angstroms, and S(Tk) = |<j_w || T^k || j_v>|^2 is the line strength.

Reference: https://www1.udel.edu/atom/about

Input files:
- {atom}_Energies.csv
- {atom}_Matrix_Elements_Theory.csv

Output files:
- {atom}_Transition_Rates.csv              (Step 1 - raw rates per matrix element row)
- {atom}_Transition_Properties.csv         (Step 2 - upper->lower with display columns)
- {atom}_Lifetimes.csv                     (Step 2 - per-state lifetimes)
- {atom}_Transition_Rates_Error_Check.csv  (Step 3 - parenthetical uncertainty display)
- {atom}_Lifetimes_Error_Check.csv         (Step 3 - parenthetical uncertainty display)
- {atom}_transitions.txt                   (Step 3 - per-state decay table:
                                            header: conf termJ - energy cm^-1 - lifetime(unc) ns ->
                                            col header + rows: Lower State, Op, Wavelength(nm), ME, Rate(s^-1), BR)
"""

import math
import decimal

import pandas as pd
import sys
import os


def _rate_formula(operator):
    """
    Return (constant, wavelength_power) for:
    
              C
    A = --------------- x S,
        (2J+1)*lambda^p

    where C = constant, p = wavelength_power, and lambda must be in Angstroms.
    """
    table = {
        'E1': (2.02613e18, 3),
        'E2': (1.11995e18, 5),
        'E3': (3.14441e17, 7),
        'M1': (2.69735e13, 3),
        'M2': (1.49097e13, 5),
        'M3': (4.18610e12, 7),
    }
    if operator not in table:
        raise ValueError('Unknown operator: ' + operator)
    return table[operator]


def _round_half_up(value, n_dec):
    """Round to n_dec decimal places using round-half-up."""
    d = decimal.Decimal(str(value))
    if n_dec >= 1:
        q = decimal.Decimal('0.' + '0' * n_dec)
    elif n_dec == 0:
        q = decimal.Decimal('1')
    else:
        q = decimal.Decimal('1E' + str(-n_dec))
    return d.quantize(q, rounding=decimal.ROUND_HALF_UP)


def _fmt5(x):
    """5 decimal places, switching to scientific notation outside [1e-4, 1e5]."""
    if x != 0.0 and (abs(x) >= 1e5 or abs(x) < 1e-4):
        return f'{x:.5e}'
    return f'{x:.5f}'


def parse_termJ(termJ):
    """
    Parse a combined termJ string into separate term and J values.

    Handles both standard terms like "3P1" -> ("3P", "1")
    and parenthetical terms like "(1/2,7/2)3" -> ("(1/2,7/2)", "3")

    Parameters:
    -----------
    termJ : str
        Combined term and J string

    Returns:
    --------
    tuple : (term, J)
    """
    if termJ.startswith('('):
        # Find the matching closing parenthesis
        close_paren = termJ.rfind(')')
        if close_paren != -1:
            term = termJ[:close_paren + 1]
            J = termJ[close_paren + 1:]
        else:
            # Fallback if no closing parenthesis found
            term = termJ[:-1]
            J = termJ[-1]
    else:
        # Standard term like "3P1" - term is everything except last char
        term = termJ[:-1]
        J = termJ[-1]

    return term, J


def format_with_uncertainty(value, unc, n_sig_figs=2):
    """
    Core parenthetical formatter used by the type-specific helpers below.

    Rounds the uncertainty to n_sig_figs significant figures, then expresses
    both the value and uncertainty to that same last decimal place.

    Examples (n_sig_figs=2):
        (1.8059, 0.0018)     -> "1.8059(18)"
        (28.64,  0.18)       -> "28.64(18)"
        (119000, 65000)      -> "119000(65000)"

    Examples (n_sig_figs=1):
        (454.9813, 0.000207) -> "454.9813(2)"
        (1050.918, 0.007)    -> "1050.918(7)"
    """
    if value == 0:
        return "0"
    if unc <= 0:
        return str(_round_half_up(value, 4))

    unc_exp = math.floor(math.log10(abs(unc)))
    n_dec   = -(unc_exp - (n_sig_figs - 1))

    if n_dec <= 0:
        scale = 10 ** (-n_dec)
        v = int(_round_half_up(value / scale, 0)) * scale
        u = int(_round_half_up(unc   / scale, 0)) * scale
        return f"{v}({u})"

    n_dec   = min(n_dec, 12)
    val_str = str(_round_half_up(value, n_dec))
    unit    = 10.0 ** (-n_dec)
    u_int   = int(_round_half_up(unc / unit, 0))

    if u_int == 0:
        return str(_round_half_up(value, n_dec))

    return f"{val_str}({u_int})"


# ---------------------------------------------------------------------------
# Private helper: format a plain number to n_dec decimal places (no unc).
# ---------------------------------------------------------------------------
def _fmt_plain(value, n_dec):
    if n_dec <= 0:
        scale = 10 ** (-n_dec)
        return str(int(_round_half_up(value / scale, 0)) * scale)
    return str(_round_half_up(value, min(n_dec, 15)))


# ---------------------------------------------------------------------------
# Type-specific formatters (following ATOM portal v3 formatting conventions)
# ---------------------------------------------------------------------------

def format_wavelength(wl, unc, show_upper_bound=True):
    """
    Wavelength display (nm).
      - 1 sig fig in uncertainty
      - Primary: decimal (e.g. 454.9813(2))
      - Exception: scientific when unc > 950 nm (e.g. 2.37(1)E+05)
      - When unc > wl and show_upper_bound=True: <(wl + unc)
      - When unc > wl and show_upper_bound=False: value(unc) as usual
    """
    if unc <= 0:
        return str(_round_half_up(wl, 4))

    unc_exp = math.floor(math.log10(abs(unc)))
    n_dec   = -(unc_exp)        # 1 sig fig: last digit at unc_exp

    if unc > wl and show_upper_bound:
        upper = wl + unc
        if unc > 950:
            exp  = math.floor(math.log10(abs(upper)))
            sign = "+" if exp >= 0 else "-"
            return f"<{upper / 10**exp:.2f}E{sign}{abs(exp):02d}"
        return f"<{_fmt_plain(upper, n_dec)}"

    # Normal: switch to scientific if unc > 950
    if unc > 950:
        exp  = math.floor(math.log10(abs(wl)))
        sign = "+" if exp >= 0 else "-"
        coefficient_str = format_with_uncertainty(wl / 10**exp, unc / 10**exp, n_sig_figs=1)
        return f"{coefficient_str}E{sign}{abs(exp):02d}"

    return format_with_uncertainty(wl, unc, n_sig_figs=1)


def format_matrix_element(me, unc, show_upper_bound=True):
    """
    Matrix element display (a.u.).
      - 2 sig figs in uncertainty
      - Primary: decimal, no exceptions
      - When unc > me and show_upper_bound=True: <(me + unc)
      - When unc > me and show_upper_bound=False: value(unc) as usual
    """
    if unc <= 0:
        return str(_round_half_up(me, 4))

    unc_exp = math.floor(math.log10(abs(unc)))
    n_dec   = -(unc_exp - 1)    # 2 sig figs

    if unc > me and show_upper_bound:
        upper = me + unc
        return f"<{_fmt_plain(upper, n_dec)}"

    return format_with_uncertainty(me, unc, n_sig_figs=2)


def format_rate_scientific(value, unc, show_upper_bound=True):
    """
    Transition rate display (s^-1).
      - 2 sig figs in uncertainty
      - Primary: scientific with parenthetical coefficient (e.g. 1.661(12)E+06)
      - Exception: E+00 is omitted for normal display (e.g. 1.5(10))
      - When unc > value and show_upper_bound=True: <(value + unc) in scientific
      - When unc > value and show_upper_bound=False: value(unc) as usual
    """
    if value == 0:
        if unc > 0 and show_upper_bound:
            upper_exp = math.floor(math.log10(abs(unc)))
            sign = "+" if upper_exp >= 0 else "-"
            return f"<{_round_half_up(unc / 10**upper_exp, 1)}E{sign}{abs(upper_exp):02d}"
        return "0"

    if unc > value and show_upper_bound:
        upper     = value + unc
        upper_exp = math.floor(math.log10(abs(upper)))
        sign      = "+" if upper_exp >= 0 else "-"
        return f"<{_round_half_up(upper / 10**upper_exp, 1)}E{sign}{abs(upper_exp):02d}"

    exp      = math.floor(math.log10(abs(value)))
    coefficient = value / (10.0 ** exp)
    unc_scaled = unc / (10.0 ** exp)
    coefficient_str    = format_with_uncertainty(coefficient, unc_scaled) if unc > 0 \
               else str(_round_half_up(coefficient, 4))

    # If rounding pushed coefficient to 10, bump the exponent and reformat
    if float(coefficient_str.split('(')[0]) >= 10.0:
        exp += 1
        coefficient_str = format_with_uncertainty(value / 10.0**exp, unc / 10.0**exp) if unc > 0 \
                else str(_round_half_up(value / 10.0**exp, 4))

    # Omit E+00 in normal display
    if exp == 0:
        return coefficient_str

    sign = "+" if exp >= 0 else "-"
    return f"{coefficient_str}E{sign}{abs(exp):02d}"


def format_branching_ratio(B, unc, show_upper_bound=True):
    """
    Branching ratio display (dimensionless, in [0, 1]).
      - 2 sig figs in uncertainty (1 sig fig when unc < 1e-10)
      - Primary: decimal, no exceptions
      - When unc > 0 but rounds to zero: shows (0) explicitly
      - When unc > B and show_upper_bound=True: <(B + unc), capped at <1
      - When unc > B and show_upper_bound=False: value(unc) as usual
    """
    if B == 0:
        if unc > 0 and show_upper_bound:
            upper = unc
            if upper >= 1.0:
                return "<1"
            n_sig_z     = 1 if unc < 1e-10 else 2
            upper_exp   = math.floor(math.log10(abs(upper)))
            upper_n_dec = -(upper_exp - (n_sig_z - 1))
            return f"<{format(_round_half_up(upper, upper_n_dec), 'f').rstrip('0').rstrip('.')}"
        return "0.0"

    if unc <= 0:
        if abs(B - 1.0) < 1e-12:
            return "1.0"
        val_exp = math.floor(math.log10(abs(B)))
        return format(_round_half_up(B, min(-(val_exp - 1), 15)), 'f')

    n_sig   = 1 if unc < 1e-10 else 2
    unc_exp = math.floor(math.log10(abs(unc)))
    n_dec   = -(unc_exp - (n_sig - 1))

    if unc > B and show_upper_bound:
        upper = B + unc
        if upper >= 1.0:
            return "<1"
        upper_exp   = math.floor(math.log10(abs(upper))) if upper > 0 else 0
        upper_n_dec = -(upper_exp - (n_sig - 1))
        return f"<{format(_round_half_up(upper, upper_n_dec), 'f').rstrip('0').rstrip('.')}"

    if n_dec <= 0:
        scale = 10 ** (-n_dec)
        v = int(_round_half_up(B   / scale, 0)) * scale
        u = int(_round_half_up(unc / scale, 0)) * scale
        return f"{v}({u})"

    n_dec   = min(n_dec, 15)
    val_str = format(_round_half_up(B, n_dec), 'f')
    unit    = 10.0 ** (-n_dec)
    u_int   = int(_round_half_up(unc / unit, 0))

    if u_int == 0:
        return f"{val_str}(0)"

    return f"{val_str}({u_int})"


def format_lifetime(lt_ns, unc, show_upper_bound=True):
    """
    Lifetime display (ns, or scientific for metastable states).
      - 2 sig figs in uncertainty
      - Primary: decimal (e.g. 1.8059(18))
      - Exception: scientific when lt_ns > 1e5 (e.g. 3.863(77)E+09)
      - When unc > lt_ns and show_upper_bound=True: <(lt_ns + unc)
      - When unc > lt_ns and show_upper_bound=False: value(unc) as usual
    """
    if unc <= 0:
        if lt_ns > 1e5:
            exp  = math.floor(math.log10(abs(lt_ns)))
            sign = "+" if exp >= 0 else "-"
            return f"{_round_half_up(lt_ns / 10**exp, 2)}E{sign}{abs(exp):02d}"
        return str(_round_half_up(lt_ns, 4))

    unc_exp = math.floor(math.log10(abs(unc)))
    n_dec   = -(unc_exp - 1)    # 2 sig figs

    if unc > lt_ns and show_upper_bound:
        upper = lt_ns + unc
        if upper > 1e5:
            exp  = math.floor(math.log10(abs(upper)))
            sign = "+" if exp >= 0 else "-"
            return f"<{_round_half_up(upper / 10**exp, 1)}E{sign}{abs(exp):02d}"
        return f"<{_fmt_plain(upper, n_dec)}"

    # Scientific for large lifetimes
    if lt_ns > 1e5:
        exp   = math.floor(math.log10(abs(lt_ns)))
        coefficient_str = format_with_uncertainty(lt_ns / 10**exp, unc / 10**exp, n_sig_figs=2)
        # If rounding pushed coefficient to 10, bump the exponent and reformat
        if float(coefficient_str.split('(')[0]) >= 10.0:
            exp += 1
            coefficient_str = format_with_uncertainty(lt_ns / 10**exp, unc / 10**exp, n_sig_figs=2)
        sign = "+" if exp >= 0 else "-"
        return f"{coefficient_str}E{sign}{abs(exp):02d}"

    return format_with_uncertainty(lt_ns, unc, n_sig_figs=2)


def calculate_transition_rates(energies_file, matrix_elements_file, atom_name):
    """
    Calculate transition rates from energies and matrix elements.

    Parameters:
    -----------
    energies_file : str
        Path to the energies CSV file
    matrix_elements_file : str
        Path to the matrix elements CSV file
    atom_name : str
        Name of the atom (e.g., 'Ba1')
    """

    # Read input files
    print(f"Reading {energies_file}...")
    energies_df = pd.read_csv(energies_file)

    print(f"Reading {matrix_elements_file}...")
    matrix_df = pd.read_csv(matrix_elements_file)

    # Create a lookup dictionary for energies
    energy_lookup = {}
    for _, row in energies_df.iterrows():
        key = (row['state_configuration'], row['state_term'], str(row['state_J']))
        energy_lookup[key] = float(row['energy'])

    # Calculate transition rates
    transition_rates = []

    for _, row in matrix_df.iterrows():
        # Get state information
        conf1 = row['state_one_configuration']
        term1 = row['state_one_term']
        J1 = str(row['state_one_J'])

        conf2 = row['state_two_configuration']
        term2 = row['state_two_term']
        J2 = str(row['state_two_J'])

        operator = str(row['operator'])
        matrix_element = row['matrix_element']
        matrix_element_unc = row['matrix_element_uncertainty']

        # Look up energies
        key1 = (conf1, term1, J1)
        key2 = (conf2, term2, J2)

        if key1 not in energy_lookup or key2 not in energy_lookup:
            print(f"Warning: Could not find energy for {key1} or {key2}, skipping...")
            continue

        energy1 = energy_lookup[key1]
        energy2 = energy_lookup[key2]

        # Calculate wavelength (nm) from energy difference (cm^-1)
        energy_diff = abs(energy1 - energy2)
        if energy_diff == 0:
            print(f"Warning: Zero energy difference for transition {key1} -> {key2}, skipping...")
            continue

        wavelength = 1e7 / energy_diff  # nm

        # J of the upper state (v) goes in the (2J_v+1) denominator
        if energy1 > energy2:
            J_upper = float(J1)
        else:
            J_upper = float(J2)

        C, p = _rate_formula(operator)
        wl_angstrom = abs(wavelength) * 10
        transition_rate = C / ((2 * J_upper + 1) * wl_angstrom ** p) * matrix_element ** 2

        # Add to list
        transition_rates.append({
            'state_one_configuration': conf1,
            'state_one_term': term1,
            'state_one_J': J1,
            'state_two_configuration': conf2,
            'state_two_term': term2,
            'state_two_J': J2,
            'operator': operator,
            'matrix_element': matrix_element,
            'matrix_element_uncertainty': matrix_element_unc,
            'energy1(cm-1)': f"{energy1:.3f}",
            'energy2(cm-1)': f"{energy2:.3f}",
            'wavelength(nm)': f"{wavelength:.2f}",
            'transition_rate(s-1)': f"{transition_rate:.4e}"
        })

    # Create dataframe and save
    tr_df = pd.DataFrame(transition_rates)
    output_file = f"{atom_name}_Transition_Rates.csv"
    tr_df.to_csv(output_file, index=False)
    print(f"{output_file} has been written with {len(tr_df)} transitions")

    return output_file


def calculate_lifetimes_and_branching_ratios(tr_file):
    """
    Calculate lifetimes and branching ratios from transition rates.

    Parameters:
    -----------
    tr_file : str
        Path to the transition rates CSV file
    """

    # Read transition rates from file
    print(f"Reading {tr_file}...")
    df = pd.read_csv(tr_file)

    atom = tr_file.split('_')[0]

    # Open files to write lifetimes and branching ratios
    filename_lifetimes = atom + '_Lifetimes.csv'
    filename_br_ratios = atom + '_Transition_Properties.csv'

    lifetimes_df = pd.DataFrame(columns=['state_configuration', 'state_term', 'state_J', 'lifetime', 'lifetime_uncertainty'])
    br_ratios_df = pd.DataFrame(columns=['state_one_configuration', 'state_one_term', 'state_one_J',
                                         'state_two_configuration', 'state_two_term', 'state_two_J',
                                         'operator', 'wavelength_display', 'matrix_element_display',
                                         'branching_ratio_display', 'transition_rate_display', 'transition_rate', 'transition_rate_uncertainty'])

    # Load energies early for uncertainty lookup and energy-based sorting
    energies_file = atom + '_Energies.csv'
    energy_unc_lookup = {}
    energies_df_for_sort = None
    try:
        energies_df = pd.read_csv(energies_file)
        energies_df['state_J'] = energies_df['state_J'].astype(str)
        for _, erow in energies_df.iterrows():
            key = (erow['state_configuration'], erow['state_term'], str(erow['state_J']))
            energy_unc_lookup[key] = float(erow.get('energy_uncertainty', 0.0))
        energies_df_for_sort = energies_df
    except FileNotFoundError:
        print(f"Warning: {energies_file} not found, skipping energy sorting and uncertainty for lifetimes")

    # Parse transition rates and sum up rates for each upper state
    # Note: Matrix_Elements_Theory.csv is already deduplicated by gen_portal_csv.py
    tr_rates = {}
    upper_energies = {}

    for _, row in df.iterrows():
        energy1 = float(row['energy1(cm-1)'])
        energy2 = float(row['energy2(cm-1)'])

        configuration1 = row['state_one_configuration']
        term1 = row['state_one_term']
        J1 = str(row['state_one_J'])
        termJ1 = term1 + J1

        configuration2 = row['state_two_configuration']
        term2 = row['state_two_term']
        J2 = str(row['state_two_J'])
        termJ2 = term2 + J2

        operator = str(row['operator'])
        _, wl_power = _rate_formula(operator)
        matrix_element = row['matrix_element']
        me_unc = float(row['matrix_element_uncertainty'])
        wavelength = float(row['wavelength(nm)'])
        tr_rate = float(row['transition_rate(s-1)'])

        config1 = configuration1 + ' ' + termJ1
        config2 = configuration2 + ' ' + termJ2

        # Determine which state is upper (higher energy)
        if energy2 < energy1:
            upper_state = config1
            upper_energies[config1] = energy1
            lower_state = config2
        else:
            upper_state = config2
            upper_energies[config2] = energy2
            lower_state = config1

        # Add transition to upper state's decay list
        if upper_state not in tr_rates:
            tr_rates[upper_state] = []

        tr_rates[upper_state].append([lower_state, tr_rate, matrix_element, wavelength, me_unc, operator, wl_power])

    # Calculate lifetimes and branching ratios
    lifetimes = []

    # Calculate lifetimes and branching ratios
    for config, rates in tr_rates.items():
        configuration = config.split(' ')[0]
        termJ = config.split(' ')[1]
        term, J = parse_termJ(termJ)

        # Calculate total decay rate
        total_rates = sum(rate[1] for rate in rates)

        # Calculate lifetime in nanoseconds
        lifetime = 1/total_rates * 1e9

        # Sum squared rate uncertainties over all decay channels for lifetime uncertainty
        upper_key = (configuration, term, J)
        e_unc_upper = energy_unc_lookup.get(upper_key, 0.0)
        sum_rate_unc_sq = 0.0
        for rate_entry in rates:
            lower_conf  = rate_entry[0].split(' ')[0]
            lower_term_p, lower_J_p = parse_termJ(rate_entry[0].split(' ')[1])
            e_unc_lower = energy_unc_lookup.get((lower_conf, lower_term_p, lower_J_p), 0.0)
            tr_rate_i  = rate_entry[1]
            me_i       = rate_entry[2]
            wl_i       = rate_entry[3]
            me_unc_i   = rate_entry[4]
            p_i        = rate_entry[6]
            wl_unc_i   = (wl_i ** 2 / 1e7) * math.sqrt(e_unc_upper ** 2 + e_unc_lower ** 2)
            if me_i != 0 and wl_i != 0:
                rate_unc_i = tr_rate_i * math.sqrt((2.0 * me_unc_i / abs(me_i)) ** 2 + (p_i * wl_unc_i / wl_i) ** 2)
            else:
                rate_unc_i = 0.0
            sum_rate_unc_sq += rate_unc_i ** 2
        lt_unc = (lifetime ** 2 / 1e9) * math.sqrt(sum_rate_unc_sq)

        lifetimes.append([configuration, term, J, lifetime, lt_unc])

        # Build one branching ratio row per decay channel
        for rate in rates:
            configuration2 = rate[0].split(' ')[0]
            termJ2 = rate[0].split(' ')[1]
            term2, J2 = parse_termJ(termJ2)
            tr_rate        = rate[1]
            matrix_element = rate[2]
            wavelength     = rate[3]
            op             = rate[5]
            branching_ratio = tr_rate / total_rates

            row_data = {
                'state_one_configuration': configuration,
                'state_one_term': term,
                'state_one_J': J,
                'state_two_configuration': configuration2,
                'state_two_term': term2,
                'state_two_J': J2,
                'operator': op,
                'wavelength_display': f"{wavelength:.2f}",
                'matrix_element_display': matrix_element,
                'branching_ratio_display': format(_round_half_up(branching_ratio, (-math.floor(math.log10(branching_ratio)) + 2) if branching_ratio > 0 else 3), 'f'),
                'transition_rate_display': '',
                'transition_rate': f"{tr_rate:.3e}",
                'transition_rate_uncertainty': ''
            }
            br_ratios_df.loc[len(br_ratios_df.index)] = row_data

    # Write lifetimes to dataframe
    for lifetime in lifetimes:
        row = {
            'state_configuration': lifetime[0],
            'state_term': lifetime[1],
            'state_J': lifetime[2],
            'lifetime': lifetime[3],
            'lifetime_uncertainty': lifetime[4]
        }
        lifetimes_df.loc[len(lifetimes_df.index)] = row

    # Sort lifetimes by energy
    if energies_df_for_sort is not None:
        lifetimes_df['state_J'] = lifetimes_df['state_J'].astype(str)
        lifetimes_df = lifetimes_df.merge(
            energies_df_for_sort[['state_configuration', 'state_term', 'state_J', 'energy']],
            on=['state_configuration', 'state_term', 'state_J'],
            how='left'
        )
        lifetimes_df['energy'] = lifetimes_df['energy'].astype(float)
        lifetimes_df = lifetimes_df.sort_values(by='energy').reset_index(drop=True)
        lifetimes_df = lifetimes_df.drop(columns=['energy'])

    # Save output files
    br_ratios_df.to_csv(filename_br_ratios, index=False)
    print(f"{filename_br_ratios} has been written with {len(br_ratios_df)} transitions")

    lifetimes_df.to_csv(filename_lifetimes, index=False)
    print(f"{filename_lifetimes} has been written with {len(lifetimes_df)} states")

    return tr_rates, upper_energies, lifetimes


def add_display_formats(atom_name, tr_rates=None, upper_energies=None, lifetimes_raw=None, br_cutoff=1e-6):
    """
    Step 3: Overwrite the display columns in the Error_Check CSVs with
    parenthetical uncertainty notation.

    Reads:
      - {atom}_Transition_Rates.csv        (matrix_element, matrix_element_uncertainty,
                                             energy1, energy2, transition_rate)
      - {atom}_Energies.csv                (energy_uncertainty per state)
      - {atom}_Transition_Properties.csv   (organised upper->lower rows from Step 2)
      - {atom}_Lifetimes.csv               (per-state lifetimes from Step 2)

    Writes:
      - {atom}_Transition_Rates_Error_Check.csv
      - {atom}_Lifetimes_Error_Check.csv

    Uncertainty formulas
    --------------------
    Wavelength:
                     lambda^2
        d_lambda =  ----------  x sqrt(dE1^2 + dE2^2)
                       1e7

    Rate:
                  d_me
        dA = 2A x ----
                  |me|

    Lifetime:
                     tau_ns^2
        d_tau_ns =  ----------  x sqrt(sum dAi^2)
                       1e9
    """

    tr_file        = f"{atom_name}_Transition_Rates.csv"
    props_file     = f"{atom_name}_Transition_Properties.csv"
    lifetimes_file = f"{atom_name}_Lifetimes.csv"
    lifetimes_exp  = f"{atom_name}_Lifetimes_Exp.csv"
    energies_file  = f"{atom_name}_Energies.csv"
    tr_check_out   = f"{atom_name}_Transition_Rates_Error_Check.csv"
    tr_check_zero_out = f"{atom_name}_Transition_Rates_Error_Check_Zero_BR.csv"
    lifetimes_out  = f"{atom_name}_Lifetimes_Error_Check.csv"

    print(f"Reading {tr_file}...")
    tr_df = pd.read_csv(tr_file)

    # Build energy lookup (full precision) and energy-uncertainty lookup
    energy_lookup     = {}   # (conf, term, J_str) -> energy (cm^-1)
    energy_unc_lookup = {}   # (conf, term, J_str) -> energy_uncertainty
    try:
        energies_df = pd.read_csv(energies_file)
        for _, row in energies_df.iterrows():
            key = (str(row['state_configuration']).strip(),
                   str(row['state_term']).strip(),
                   str(row['state_J']).strip())
            energy_lookup[key] = float(row['energy'])
            if 'energy_uncertainty' in energies_df.columns:
                energy_unc_lookup[key] = float(row['energy_uncertainty'])
    except FileNotFoundError:
        print(f"Warning: {energies_file} not found; wavelength uncertainties will be 0")

    # ----------------------------------------------------------------------- #
    # Build per-transition lookup                                             #
    # key: (upper_conf, upper_term, upper_J, lower_conf, lower_term, lower_J) #
    # ----------------------------------------------------------------------- #
    trans_data  = {}   # key -> {me, me_unc, wavelength, wl_unc, rate, rate_unc}
    upper_trans = {}   # (upper_conf, upper_term, upper_J) -> [(rate, rate_unc), ...]

    for _, row in tr_df.iterrows():
        energy1 = float(row['energy1(cm-1)'])
        energy2 = float(row['energy2(cm-1)'])

        conf1 = str(row['state_one_configuration']).strip()
        term1 = str(row['state_one_term']).strip()
        J1    = str(row['state_one_J']).strip()

        conf2 = str(row['state_two_configuration']).strip()
        term2 = str(row['state_two_term']).strip()
        J2    = str(row['state_two_J']).strip()

        operator = str(row['operator']).strip()
        C, p     = _rate_formula(operator)
        me      = float(row['matrix_element'])
        me_unc  = float(row['matrix_element_uncertainty'])
        rate    = float(row['transition_rate(s-1)'])

        # Recompute wavelength with full precision from Energies.csv
        e1 = energy_lookup.get((conf1, term1, J1), energy1)
        e2 = energy_lookup.get((conf2, term2, J2), energy2)
        wl = 1e7 / abs(e1 - e2) if abs(e1 - e2) > 0 else float(row['wavelength(nm)'])

        # Wavelength uncertainty via error propagation
        e_unc1  = energy_unc_lookup.get((conf1, term1, J1), 0.0)
        e_unc2  = energy_unc_lookup.get((conf2, term2, J2), 0.0)
        wl_unc  = (wl ** 2 / 1e7) * math.sqrt(e_unc1 ** 2 + e_unc2 ** 2)

        # Identify upper and lower states
        if energy1 > energy2:
            upper = (conf1, term1, J1)
            lower = (conf2, term2, J2)
            J_upper_str = J1
        else:
            upper = (conf2, term2, J2)
            lower = (conf1, term1, J1)
            J_upper_str = J2

        # Rate uncertainty:
        #
        #   dA        /  (2*d_me)^2   (p*d_lambda)^2  \
        #   -- = sqrt(  ---------  + ----------------  )
        #    A        \    me^2          lambda^2     /
        # 
        # where p is the wavelength power for the operator (3 for E1/M1, 5 for E2/M2, 7 for E3/M3).
        # When me == 0 but me_unc > 0, evaluate the rate formula at me = me_unc for an upper bound.
        if me != 0 and wl != 0:
            rate_unc = rate * math.sqrt((2.0 * me_unc / abs(me)) ** 2 + (p * wl_unc / wl) ** 2)
        elif me_unc > 0 and wl != 0:
            j_parts     = J_upper_str.split('/')
            j_val       = float(j_parts[0]) / float(j_parts[1]) if len(j_parts) == 2 \
                          else float(J_upper_str)
            j_degen     = 2 * j_val + 1
            wl_angstrom = wl * 10.0
            rate_unc    = C / j_degen * (me_unc ** 2) / (wl_angstrom ** p)
        else:
            rate_unc = 0.0

        trans_data[upper + lower + (operator,)] = {
            'me':      me,      'me_unc':   me_unc,
            'wl':      wl,      'wl_unc':   wl_unc,
            'rate':    rate,    'rate_unc': rate_unc,
        }

        upper_trans.setdefault(upper, []).append((rate, rate_unc))

    # ------------------------------------------------------------------ #
    # Write Transition_Rates_Error_Check.csv                             #
    # ------------------------------------------------------------------ #
    print(f"Reading {props_file}...")
    props_df = pd.read_csv(props_file)

    wl_disp, me_disp, br_disp, rate_disp, rate_unc_list = [], [], [], [], []
    br_values = []

    for _, row in props_df.iterrows():
        upper = (str(row['state_one_configuration']).strip(),
                 str(row['state_one_term']).strip(),
                 str(row['state_one_J']).strip())
        lower = (str(row['state_two_configuration']).strip(),
                 str(row['state_two_term']).strip(),
                 str(row['state_two_J']).strip())
        operator = str(row['operator']).strip()
        key = upper + lower + (operator,)

        if key in trans_data:
            t = trans_data[key]
            wl_disp.append(format_wavelength(t['wl'],      t['wl_unc']))
            me_disp.append(format_matrix_element(t['me'],  t['me_unc']))
            rate_disp.append(format_rate_scientific(t['rate'], t['rate_unc']))
            rate_unc_list.append(f"{t['rate_unc']:.3e}")

            # Branching ratio and its uncertainty via error propagation
            #
            #          A_i
            # B_i = ---------
            #        A_total
            #
            #          sqrt[ (1-B_i)^2 x dA_i^2 + B_i^2 x sum_{j!=i} dA_j^2 ]
            # dB_i = -----------------------------------------------------------
            #                               A_total
            #
            A_i  = t['rate']
            dA_i = t['rate_unc']
            trans_list = upper_trans.get(upper, [(A_i, dA_i)])
            A_total    = sum(r for r, _ in trans_list)
            B_i        = A_i / A_total if A_total > 0 else 0.0
            sum_sq_all = sum(u ** 2 for _, u in trans_list)
            dB_i = math.sqrt(
                ((1 - B_i) / A_total) ** 2 * dA_i ** 2 +
                (B_i       / A_total) ** 2 * (sum_sq_all - dA_i ** 2)
            ) if A_total > 0 else 0.0
            br_disp.append(format_branching_ratio(B_i, dB_i))
            br_values.append(B_i)
        else:
            print(f"  Warning: no {operator} transition data for {upper} -> {lower}, keeping old values")
            wl_disp.append(str(row['wavelength_display']))
            me_disp.append(str(row['matrix_element_display']))
            br_disp.append(str(row['branching_ratio_display']))
            rate_disp.append(str(row.get('transition_rate_display', '')))
            rate_unc_list.append('')
            br_values.append(None)

    # Write Transition_Properties.csv with plain rate and uncertainty columns
    props_df['transition_rate_display']     = rate_disp
    props_df['transition_rate_uncertainty'] = rate_unc_list
    props_df.to_csv(props_file, index=False)
    print(f"{props_file} updated with transition_rate_uncertainty")

    # Build Error_Check df — only display format columns
    error_check_cols = [c for c in props_df.columns
                        if c not in ('transition_rate', 'transition_rate_uncertainty')]
    tr_check_df = props_df[error_check_cols].copy()
    tr_check_df['wavelength_display']      = wl_disp
    tr_check_df['matrix_element_display']  = me_disp
    tr_check_df['branching_ratio_display'] = br_disp
    tr_check_df['transition_rate_display'] = rate_disp

    # Split: rows with branching ratio below cutoff go to a separate file
    # so the main error-check CSV only contains significant transitions.
    zero_mask = pd.Series([(b is not None and b < br_cutoff) for b in br_values],
                          index=tr_check_df.index)
    tr_check_nonzero = tr_check_df[~zero_mask]
    tr_check_zero    = tr_check_df[zero_mask]

    tr_check_nonzero.to_csv(tr_check_out, index=False)
    print(f"{tr_check_out} written with {len(tr_check_nonzero)} transitions")
    if len(tr_check_zero) > 0:
        tr_check_zero.to_csv(tr_check_zero_out, index=False)
        print(f"{tr_check_zero_out} written with {len(tr_check_zero)} near-zero BR transitions")

    # ------------------------------------------------------------------ #
    # Write Lifetimes_Error_Check.csv                                    #
    # ------------------------------------------------------------------ #
    print(f"Reading {lifetimes_file}...")
    lifetimes_df = pd.read_csv(lifetimes_file)

    # Load experimental lifetimes if the file exists
    exp_lifetime_lookup = {}   # (conf, term, J) -> {'lifetime_display': ..., 'lifetime_ref': ...}
    if os.path.exists(lifetimes_exp):
        print(f"Reading {lifetimes_exp}...")
        lifetimes_exp_df = pd.read_csv(lifetimes_exp)
        for _, erow in lifetimes_exp_df.iterrows():
            key = (str(erow['state_configuration']).strip(),
                   str(erow['state_term']).strip(),
                   str(erow['state_J']).strip())
            exp_lifetime_lookup[key] = {
                'lifetime_display': str(erow['lifetime_display']),
                'lifetime_ref': str(erow['lifetime_ref']) if pd.notna(erow['lifetime_ref']) else ''
            }
        print(f"  {len(exp_lifetime_lookup)} experimental lifetimes loaded")
    else:
        print(f"  {lifetimes_exp} not found, using theory lifetimes only")

    lt_disp = []
    lt_refs = []

    for _, row in lifetimes_df.iterrows():
        upper = (str(row['state_configuration']).strip(),
                 str(row['state_term']).strip(),
                 str(row['state_J']).strip())

        lt_ns = float(row['lifetime'])

        if upper in exp_lifetime_lookup:
            exp = exp_lifetime_lookup[upper]
            lt_disp.append(exp['lifetime_display'])
            lt_refs.append(exp['lifetime_ref'])
        elif upper in upper_trans:
            rate_uncs      = [u for _, u in upper_trans[upper]]
            total_rate_unc = math.sqrt(sum(u ** 2 for u in rate_uncs))
            lt_unc         = (lt_ns ** 2 / 1e9) * total_rate_unc
            lt_disp.append(format_lifetime(lt_ns, lt_unc))
            lt_refs.append('')
        else:
            print(f"  Warning: no transition data for lifetime state {upper}, keeping old value")
            lt_disp.append(str(lt_ns))
            lt_refs.append('')

    lt_check_df = lifetimes_df.copy()
    lt_check_df['lifetime_display'] = lt_disp
    lt_check_df = lt_check_df.drop(columns=['lifetime', 'lifetime_uncertainty'])
    lt_check_df['lifetime_ref'] = lt_refs
    lt_check_df.to_csv(lifetimes_out, index=False)
    print(f"{lifetimes_out} written with {len(lt_check_df)} states")

    if tr_rates is not None and upper_energies is not None and lifetimes_raw is not None:
        lt_lookup = {(lt[0], lt[1], lt[2]): (lt[3], lt[4]) for lt in lifetimes_raw}
        transitions_file = f"{atom_name}_transitions.txt"
        with open(transitions_file, 'w') as f:
            for config, rates in sorted(tr_rates.items(), key=lambda kv: upper_energies.get(kv[0], 0.0)):
                configuration = config.split(' ')[0]
                termJ = config.split(' ')[1]
                term, J = parse_termJ(termJ)
                upper = (configuration, term, J)
                energy_cm = upper_energies.get(config, 0.0)
                lt_ns, lt_unc = lt_lookup.get(upper, (0.0, 0.0))
                total_rates = sum(rate[1] for rate in rates)
                lt_str = format_lifetime(lt_ns, lt_unc, show_upper_bound=False)
                if lt_unc > lt_ns:
                    lt_str += '*'
                f.write(f'{config} - {energy_cm:.2f} cm^-1 - {lt_str} ns ->\n')
                W_ST = 20; W_OP = 4; W_WL = 16; W_ME = 16; W_RT = 18
                f.write(f'      {"Lower State":<{W_ST}} {"Op":<{W_OP}} {"Wavelength(nm)":<{W_WL}} {"ME":<{W_ME}} {"Rate(s^-1)":<{W_RT}} BR\n')
                for rate in rates:
                    lower_conf = rate[0].split(' ')[0]
                    lower_term, lower_J = parse_termJ(rate[0].split(' ')[1])
                    lower = (lower_conf, lower_term, lower_J)
                    operator = rate[5]
                    key = upper + lower + (operator,)
                    if key in trans_data:
                        t = trans_data[key]
                        A_i  = t['rate']
                        dA_i = t['rate_unc']
                        trans_list = upper_trans.get(upper, [(A_i, dA_i)])
                        A_total    = sum(r for r, _ in trans_list)
                        B_i        = A_i / A_total if A_total > 0 else 0.0
                        if B_i < br_cutoff:
                            continue
                        sum_sq_all = sum(u ** 2 for _, u in trans_list)
                        dB_i = math.sqrt(
                            ((1 - B_i) / A_total) ** 2 * dA_i ** 2 +
                            (B_i       / A_total) ** 2 * (sum_sq_all - dA_i ** 2)
                        ) if A_total > 0 else 0.0
                        wl_str   = format_wavelength(t['wl'], t['wl_unc'], show_upper_bound=False)
                        me_str   = format_matrix_element(t['me'], t['me_unc'], show_upper_bound=False)
                        rate_str = format_rate_scientific(A_i, dA_i, show_upper_bound=False)
                        br_str   = format_branching_ratio(B_i, dB_i, show_upper_bound=False)
                        if t['wl_unc'] > t['wl']:
                            wl_str += '*'
                        if t['me_unc'] > abs(t['me']):
                            me_str += '*'
                        if dA_i > A_i:
                            rate_str += '*'
                        if dB_i > B_i:
                            br_str += '*'
                        f.write(f'      {rate[0]:<{W_ST}} {operator:<{W_OP}} {wl_str:<{W_WL}} {me_str:<{W_ME}} {rate_str:<{W_RT}} {br_str}\n')
                    else:
                        br = rate[1] / total_rates if total_rates > 0 else 0.0
                        wl_s  = _fmt5(rate[3]); me_s = _fmt5(rate[2])
                        rt_s  = _fmt5(rate[1]); br_s = _fmt5(br)
                        f.write(f'      {rate[0]:<{W_ST}} {operator:<{W_OP}} {wl_s:<{W_WL}} {me_s:<{W_ME}} {rt_s:<{W_RT}} {br_s}\n')
        print(f'{transitions_file} has been written')


def main():
    """Main function to process command line arguments and run the script."""
    import argparse
    parser = argparse.ArgumentParser(
        description='Generate transition rates and lifetimes from energy levels and matrix elements.'
    )
    parser.add_argument('atom_name', help='Atom name (e.g. Ba1, In2)')
    parser.add_argument('--br-cutoff', type=float, default=1e-6,
                        help='Branching ratio cutoff; transitions below this are excluded from '
                             '{atom}_transitions.txt and moved to the Zero_BR CSV (default: 1e-6)')
    args = parser.parse_args()
    atom_name = args.atom_name
    br_cutoff = args.br_cutoff

    # Check if input files exist
    energies_file = f"{atom_name}_Energies.csv"
    matrix_elements_file = f"{atom_name}_Matrix_Elements_Theory.csv"

    if not os.path.exists(energies_file):
        print(f"Error: {energies_file} not found!")
        sys.exit(1)

    if not os.path.exists(matrix_elements_file):
        print(f"Error: {matrix_elements_file} not found!")
        sys.exit(1)

    print(f"Processing {atom_name}...")
    print("=" * 60)

    # Step 1: Calculate transition rates
    print("\nStep 1: Calculating transition rates from matrix elements...")
    tr_file = calculate_transition_rates(energies_file, matrix_elements_file, atom_name)

    # Step 2: Calculate lifetimes and branching ratios
    print("\nStep 2: Calculating lifetimes and branching ratios...")
    tr_rates, upper_energies, lifetimes_raw = calculate_lifetimes_and_branching_ratios(tr_file)

    # Step 3: Add parenthetical uncertainty display formats and write {atom}_transitions.txt
    print("\nStep 3: Adding parenthetical uncertainty display formats...")
    add_display_formats(atom_name, tr_rates, upper_energies, lifetimes_raw, br_cutoff)

    print("\n" + "=" * 60)
    print("Processing complete!")


if __name__ == "__main__":
    main()
