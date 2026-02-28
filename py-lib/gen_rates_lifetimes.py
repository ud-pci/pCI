"""
Script to generate transition rates and lifetimes from energy levels and matrix elements.

Input files:
- {atom}_Energies.csv
- {atom}_Matrix_Elements_Theory.csv

Output files:
- {atom}_Transition_Rates.csv              (Step 1 - raw rates per matrix element row)
- {atom}_Transition_Properties.csv         (Step 2 - upper->lower with display columns)
- {atom}_Lifetimes.csv                     (Step 2 - per-state lifetimes)
- {atom}_Transition_Rates_Error_Check.csv  (Step 3 - parenthetical uncertainty display)
- {atom}_Lifetimes_Error_Check.csv         (Step 3 - parenthetical uncertainty display)
"""

import math

import pandas as pd
import sys
import os


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
        return f"{value:.4f}"

    unc_exp = math.floor(math.log10(abs(unc)))
    n_dec   = -(unc_exp - (n_sig_figs - 1))

    if n_dec <= 0:
        scale = 10 ** (-n_dec)
        v = int(round(value / scale)) * scale
        u = int(round(unc   / scale)) * scale
        return f"{v}({u})"

    n_dec   = min(n_dec, 12)
    val_str = f"{value:.{n_dec}f}"
    unit    = 10.0 ** (-n_dec)
    u_int   = round(unc / unit)

    if u_int == 0:
        return f"{value:.{n_dec}f}"

    return f"{val_str}({u_int})"


# ---------------------------------------------------------------------------
# Private helper: format a plain number to n_dec decimal places (no unc).
# ---------------------------------------------------------------------------
def _fmt_plain(value, n_dec):
    if n_dec <= 0:
        scale = 10 ** (-n_dec)
        return str(int(round(value / scale)) * scale)
    return f"{value:.{min(n_dec, 15)}f}"


# ---------------------------------------------------------------------------
# Type-specific formatters (following ATOM portal v3 formatting conventions)
# ---------------------------------------------------------------------------

def format_wavelength(wl, unc):
    """
    Wavelength display (nm).
      - 1 sig fig in uncertainty
      - Primary: decimal (e.g. 454.9813(2))
      - Exception: scientific when unc > 950 nm (e.g. 2.37(1)E+05)
      - Version 3: unc > wl  ->  <(wl + unc)
    """
    if unc <= 0:
        return f"{wl:.4f}"

    unc_exp = math.floor(math.log10(abs(unc)))
    n_dec   = -(unc_exp)        # 1 sig fig: last digit at unc_exp

    # Version 3: uncertainty exceeds value
    if unc > wl:
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
        m_str = format_with_uncertainty(wl / 10**exp, unc / 10**exp, n_sig_figs=1)
        return f"{m_str}E{sign}{abs(exp):02d}"

    return format_with_uncertainty(wl, unc, n_sig_figs=1)


def format_matrix_element(me, unc):
    """
    Matrix element display (a.u.).
      - 2 sig figs in uncertainty
      - Primary: decimal, no exceptions
      - Version 3: unc > me  ->  <(me + unc)
    """
    if unc <= 0:
        return f"{me:.4f}"

    unc_exp = math.floor(math.log10(abs(unc)))
    n_dec   = -(unc_exp - 1)    # 2 sig figs

    # Version 3
    if unc > me:
        upper = me + unc
        return f"<{_fmt_plain(upper, n_dec)}"

    return format_with_uncertainty(me, unc, n_sig_figs=2)


def format_rate_scientific(value, unc):
    """
    Transition rate display (s^-1).
      - 2 sig figs in uncertainty
      - Primary: scientific with parenthetical mantissa (e.g. 1.661(12)E+06)
      - Exception: E+00 is omitted for normal display (e.g. 1.5(10))
      - Version 3: unc > value  ->  <(value + unc) in scientific; E+00 kept
    """
    if value == 0:
        if unc > 0:
            upper_exp = math.floor(math.log10(abs(unc)))
            sign = "+" if upper_exp >= 0 else "-"
            return f"<{unc / 10**upper_exp:.2f}E{sign}{abs(upper_exp):02d}"
        return "0"

    # Version 3: uncertainty exceeds value — always full scientific, E+00 kept
    if unc > value:
        upper     = value + unc
        upper_exp = math.floor(math.log10(abs(upper)))
        sign      = "+" if upper_exp >= 0 else "-"
        return f"<{upper / 10**upper_exp:.2f}E{sign}{abs(upper_exp):02d}"

    exp      = math.floor(math.log10(abs(value)))
    mantissa = value / (10.0 ** exp)
    m_str    = format_with_uncertainty(mantissa, unc / (10.0 ** exp)) if unc > 0 \
               else f"{mantissa:.4f}"

    # Omit E+00 in normal display
    if exp == 0:
        return m_str

    sign = "+" if exp >= 0 else "-"
    return f"{m_str}E{sign}{abs(exp):02d}"


def format_branching_ratio(B, unc):
    """
    Branching ratio display (dimensionless, in [0, 1]).
      - 2 sig figs in uncertainty (1 sig fig when unc < 1e-10)
      - Primary: decimal, no exceptions
      - When unc > 0 but rounds to zero: shows (0) explicitly
      - Version 3: unc > B  ->  <(B + unc), capped at <1
    """
    if B == 0:
        if unc > 0:
            upper = unc
            if upper >= 1.0:
                return "<1"
            n_sig_z   = 1 if unc < 1e-10 else 2
            upper_exp = math.floor(math.log10(abs(upper)))
            upper_n_dec = -(upper_exp - (n_sig_z - 1))
            return f"<{upper:.{upper_n_dec}f}"
        return "0.0"

    if unc <= 0:
        if abs(B - 1.0) < 1e-12:
            return "1.0"
        val_exp = math.floor(math.log10(abs(B)))
        return f"{B:.{min(-(val_exp - 1), 15)}f}"

    n_sig   = 1 if unc < 1e-10 else 2
    unc_exp = math.floor(math.log10(abs(unc)))
    n_dec   = -(unc_exp - (n_sig - 1))

    # Version 3
    if unc > B:
        upper = B + unc
        if upper >= 1.0:
            return "<1"
        # Compute n_dec from the upper limit so the last digit is non-zero
        upper_exp   = math.floor(math.log10(abs(upper))) if upper > 0 else 0
        upper_n_dec = -(upper_exp - (n_sig - 1))
        return f"<{upper:.{upper_n_dec}f}"

    if n_dec <= 0:
        scale = 10 ** (-n_dec)
        v = int(round(B   / scale)) * scale
        u = int(round(unc / scale)) * scale
        return f"{v}({u})"

    n_dec   = min(n_dec, 15)
    val_str = f"{B:.{n_dec}f}"
    unit    = 10.0 ** (-n_dec)
    u_int   = round(unc / unit)

    if u_int == 0:
        return f"{val_str}(0)"

    return f"{val_str}({u_int})"


def format_lifetime(lt_ns, unc):
    """
    Lifetime display (ns, or scientific for metastable states).
      - 2 sig figs in uncertainty
      - Primary: decimal (e.g. 1.8059(18))
      - Exception: scientific when lt_ns > 1e5 (e.g. 3.863(77)E+09)
      - Version 3: unc > lt_ns  ->  <(lt_ns + unc)
    """
    if unc <= 0:
        if lt_ns > 1e5:
            exp  = math.floor(math.log10(abs(lt_ns)))
            sign = "+" if exp >= 0 else "-"
            return f"{lt_ns / 10**exp:.2f}E{sign}{abs(exp):02d}"
        return f"{lt_ns:.4f}"

    unc_exp = math.floor(math.log10(abs(unc)))
    n_dec   = -(unc_exp - 1)    # 2 sig figs

    # Version 3
    if unc > lt_ns:
        upper = lt_ns + unc
        if upper > 1e5:
            exp  = math.floor(math.log10(abs(upper)))
            sign = "+" if exp >= 0 else "-"
            return f"<{upper / 10**exp:.2f}E{sign}{abs(exp):02d}"
        return f"<{_fmt_plain(upper, n_dec)}"

    # Scientific for large lifetimes
    if lt_ns > 1e5:
        exp  = math.floor(math.log10(abs(lt_ns)))
        sign = "+" if exp >= 0 else "-"
        m_str = format_with_uncertainty(lt_ns / 10**exp, unc / 10**exp, n_sig_figs=2)
        return f"{m_str}E{sign}{abs(exp):02d}"

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

        wavelength = 1e7 / energy_diff  # Convert cm^-1 to nm

        # Calculate transition rate using formula
        # trate = (2.02613*10**18)/((2*J+1)*(abs(wavelength)*10)**3)*matrix_element**2
        # The J value used depends on which state has higher energy
        if energy1 > energy2:
            # Transition from state 1 to state 2
            J_denominator = int(J1)
        else:
            # Transition from state 2 to state 1
            J_denominator = int(J2)

        transition_rate = (2.02613e18) / ((2*J_denominator + 1) * (abs(wavelength)*10)**3) * matrix_element**2

        # Add to list
        transition_rates.append({
            'state_one_configuration': conf1,
            'state_one_term': term1,
            'state_one_J': J1,
            'state_two_configuration': conf2,
            'state_two_term': term2,
            'state_two_J': J2,
            'matrix_element': matrix_element,
            'matrix_element_uncertainty': matrix_element_unc,
            'energy1(cm-1)': f"{energy1:.2f}",
            'energy2(cm-1)': f"{energy2:.2f}",
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
                                         'wavelength_display', 'matrix_element_display',
                                         'branching_ratio_display', 'transition_rate_display'])

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

        matrix_element = row['matrix_element']
        me_unc = float(row['matrix_element_uncertainty'])
        wavelength = float(row['wavelength(nm)'])
        tr_rate = float(row['transition_rate(s-1)'])

        config1 = configuration1 + ' ' + termJ1
        config2 = configuration2 + ' ' + termJ2

        # Determine which state is upper (higher energy)
        if energy2 < energy1:
            upper_state = config1
            lower_state = config2
        else:
            upper_state = config2
            lower_state = config1

        # Add transition to upper state's decay list
        if upper_state not in tr_rates:
            tr_rates[upper_state] = []

        tr_rates[upper_state].append([lower_state, tr_rate, matrix_element, wavelength, me_unc])

    # Calculate lifetimes and branching ratios
    lifetimes = []

    # Write transitions to text file for inspection
    with open('transitions.txt', 'w') as f:
        for config, rates in tr_rates.items():
            configuration = config.split(' ')[0]
            termJ = config.split(' ')[1]
            term, J = parse_termJ(termJ)

            # Calculate total decay rate
            total_rates = sum(rate[1] for rate in rates)

            # Calculate lifetime in nanoseconds
            lifetime = 1/total_rates * 1e9

            # Calculate lifetime uncertainty via error propagation
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
                wl_unc_i   = (wl_i ** 2 / 1e7) * math.sqrt(e_unc_upper ** 2 + e_unc_lower ** 2)
                if me_i != 0 and wl_i != 0:
                    rate_unc_i = tr_rate_i * math.sqrt((2.0 * me_unc_i / abs(me_i)) ** 2
                                                       + (3.0 * wl_unc_i / wl_i) ** 2)
                else:
                    rate_unc_i = 0.0
                sum_rate_unc_sq += rate_unc_i ** 2
            lt_unc = (lifetime ** 2 / 1e9) * math.sqrt(sum_rate_unc_sq)

            lifetimes.append([configuration, term, J, lifetime, lt_unc])

            # Write to transitions file
            f.write(config + ' ->\n')
            for rate in rates:
                f.write('      ' + rate[0] + ' ' + str(rate[1]) + ' ' + str(rate[2]) + ' ' + str(rate[3]) + '\n')

                # Add to branching ratios dataframe
                configuration2 = rate[0].split(' ')[0]
                termJ2 = rate[0].split(' ')[1]
                term2, J2 = parse_termJ(termJ2)
                matrix_element = rate[2]
                wavelength = rate[3]
                tr_rate = rate[1]
                branching_ratio = tr_rate / total_rates

                row_data = {
                    'state_one_configuration': configuration,
                    'state_one_term': term,
                    'state_one_J': J,
                    'state_two_configuration': configuration2,
                    'state_two_term': term2,
                    'state_two_J': J2,
                    'wavelength_display': f"{wavelength:.2f}",
                    'matrix_element_display': matrix_element,
                    'branching_ratio_display': f"{branching_ratio:.3e}",
                    'transition_rate_display': f"{tr_rate:.3e}"
                }
                br_ratios_df.loc[len(br_ratios_df.index)] = row_data

    print(f"transitions.txt has been written")

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


def add_display_formats(atom_name):
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
    Wavelength:   δλ = (λ² / 1e7) x sqrt(δE₁² + δE₂²)
    Rate:         δA = 2 x A x (δd / |d|)
    Lifetime:     δτ_ns = (τ_ns² / 1e9) x sqrt(Σ δAᵢ²)
    """

    tr_file        = f"{atom_name}_Transition_Rates.csv"
    props_file     = f"{atom_name}_Transition_Properties.csv"
    lifetimes_file = f"{atom_name}_Lifetimes.csv"
    energies_file  = f"{atom_name}_Energies.csv"
    tr_check_out   = f"{atom_name}_Transition_Rates_Error_Check.csv"
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

        # Rate uncertainty: δA/A = sqrt( (2δd/d)² + (3δλ/λ)² )
        # When me == 0 but me_unc > 0, the linear formula gives δA = 0 (not useful).
        # Instead, evaluate the rate formula at me = me_unc to get an upper bound:
        #   A_max = (C / (2J+1)) x me_unc² / (λ_Å)³
        # i.e. the largest rate consistent with |me| ≤ me_unc.
        if me != 0 and wl != 0:
            rate_unc = rate * math.sqrt((2.0 * me_unc / abs(me)) ** 2 + (3.0 * wl_unc / wl) ** 2)
        elif me_unc > 0 and wl != 0:
            j_parts     = J_upper_str.split('/')
            j_val       = float(j_parts[0]) / float(j_parts[1]) if len(j_parts) == 2 \
                          else float(J_upper_str)
            j_degen     = 2 * j_val + 1
            wl_angstrom = wl * 10.0
            rate_unc    = (2.02613e18 / j_degen) * (me_unc ** 2) / (wl_angstrom ** 3)
        else:
            rate_unc = 0.0

        trans_data[upper + lower] = {
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

    wl_disp, me_disp, br_disp, rate_disp = [], [], [], []

    for _, row in props_df.iterrows():
        upper = (str(row['state_one_configuration']).strip(),
                 str(row['state_one_term']).strip(),
                 str(row['state_one_J']).strip())
        lower = (str(row['state_two_configuration']).strip(),
                 str(row['state_two_term']).strip(),
                 str(row['state_two_J']).strip())
        key = upper + lower

        if key in trans_data:
            t = trans_data[key]
            wl_disp.append(format_wavelength(t['wl'],      t['wl_unc']))
            me_disp.append(format_matrix_element(t['me'],  t['me_unc']))
            rate_disp.append(format_rate_scientific(t['rate'], t['rate_unc']))

            # Branching ratio and its uncertainty via error propagation
            # B_i = A_i / A_total
            # δB_i = (1/A_total) x sqrt[ (1-B_i)² δA_i² + B_i² x Σ_{j≠i} δA_j² ]
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
        else:
            print(f"  Warning: no transition data for {upper} -> {lower}, keeping old values")
            wl_disp.append(str(row['wavelength_display']))
            me_disp.append(str(row['matrix_element_display']))
            br_disp.append(str(row['branching_ratio_display']))
            rate_disp.append(str(row['transition_rate_display']))

    tr_check_df = props_df.copy()
    tr_check_df['wavelength_display']       = wl_disp
    tr_check_df['matrix_element_display']   = me_disp
    tr_check_df['branching_ratio_display']  = br_disp
    tr_check_df['transition_rate_display']  = rate_disp
    tr_check_df.to_csv(tr_check_out, index=False)
    print(f"{tr_check_out} written with {len(tr_check_df)} transitions")

    # ------------------------------------------------------------------ #
    # Write Lifetimes_Error_Check.csv                                    #
    # ------------------------------------------------------------------ #
    print(f"Reading {lifetimes_file}...")
    lifetimes_df = pd.read_csv(lifetimes_file)

    lt_disp = []

    for _, row in lifetimes_df.iterrows():
        upper = (str(row['state_configuration']).strip(),
                 str(row['state_term']).strip(),
                 str(row['state_J']).strip())

        lt_ns = float(row['lifetime'])

        if upper in upper_trans:
            rate_uncs      = [u for _, u in upper_trans[upper]]
            total_rate_unc = math.sqrt(sum(u ** 2 for u in rate_uncs))
            lt_unc         = (lt_ns ** 2 / 1e9) * total_rate_unc
            lt_disp.append(format_lifetime(lt_ns, lt_unc))
        else:
            print(f"  Warning: no transition data for lifetime state {upper}, keeping old value")
            lt_disp.append(str(lt_ns))

    lt_check_df = lifetimes_df.copy()
    lt_check_df['lifetime_display'] = lt_disp
    lt_check_df = lt_check_df.drop(columns=['lifetime', 'lifetime_uncertainty'])
    lt_check_df['lifetime_ref'] = ''
    lt_check_df.to_csv(lifetimes_out, index=False)
    print(f"{lifetimes_out} written with {len(lt_check_df)} states")


def main():
    """Main function to process command line arguments and run the script."""

    if len(sys.argv) < 2:
        print("Usage: python generate_transition_rates.py <atom_name>")
        print("Example: python generate_transition_rates.py Ba1")
        print("\nThis script expects the following files to exist:")
        print("  - {atom_name}_Energies.csv")
        print("  - {atom_name}_Matrix_Elements_Theory.csv")
        print("\nIt will generate:")
        print("  - {atom_name}_Transition_Rates.csv              (Step 1)")
        print("  - {atom_name}_Transition_Properties.csv         (Step 2)")
        print("  - {atom_name}_Lifetimes.csv                     (Step 2)")
        print("  - {atom_name}_Transition_Rates_Error_Check.csv  (Step 3)")
        print("  - {atom_name}_Lifetimes_Error_Check.csv         (Step 3)")
        print("  - transitions.txt")
        sys.exit(1)

    atom_name = sys.argv[1]

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
    calculate_lifetimes_and_branching_ratios(tr_file)

    # Step 3: Add parenthetical uncertainty display formats
    print("\nStep 3: Adding parenthetical uncertainty display formats...")
    add_display_formats(atom_name)

    print("\n" + "=" * 60)
    print("Processing complete!")


if __name__ == "__main__":
    main()
