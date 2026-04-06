import pandas as pd
import numpy as np
from pathlib import Path
import re


_LEADING_ORB = re.compile(r"^(\d+[spdfg](?:2|6|10|14))\.")

def _strip_to_theory(nist_conf, theory_confs):
    """Strip leading fully-occupied orbital groups from a NIST config one at a time
    until the remainder matches a known theory config, or no more can be stripped."""
    config = str(nist_conf)
    while config not in theory_confs:
        m = _LEADING_ORB.match(config)
        if m:
            config = config[len(m.group(0)):]  # remove leading orbital + dot
        else:
            break
    return config


def load_nist_lines(filepath):
    df = pd.read_csv(filepath)
    # Strip brackets from energies, e.g. [73719.4] -> 73719.4
    for col in ['E_lower(cm-1)', 'E_upper(cm-1)']:
        df[col] = df[col].astype(str).str.replace(r'[\[\]]', '', regex=True)
        df[col] = pd.to_numeric(df[col], errors='coerce').round(3)
    df['Aki(s^-1)'] = pd.to_numeric(df['Aki(s^-1)'], errors='coerce')
    for col in ['ritz_wl_vac(nm)', 'unc_ritz_wl', 'obs_wl_vac(nm)', 'unc_obs_wl']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col].astype(str).str.rstrip('+-'), errors='coerce')
    # Normalize terms for matching: strip "*" and replace letter-prefix space with dot (e.g. "a 3D" -> "a.3D")
    for src, key in [('term_lower', 'term_lower_key'), ('term_upper', 'term_upper_key')]:
        df[key] = (df[src].str.replace('*', '', regex=False)
                          .str.replace(r'^([a-zA-Z])\s+', r'\1.', regex=True))
    return df


def load_theory_lines(filepath):
    df = pd.read_csv(filepath)
    # Determine lower/upper from energy
    swap = df['energy1(cm-1)'] < df['energy2(cm-1)']
    df['E_lower'] = np.where(swap, df['energy1(cm-1)'], df['energy2(cm-1)']).round(3)
    df['E_upper'] = np.where(swap, df['energy2(cm-1)'], df['energy1(cm-1)']).round(3)
    df['conf_lower'] = np.where(swap, df['state_one_configuration'], df['state_two_configuration'])
    df['term_lower'] = np.where(swap, df['state_one_term'], df['state_two_term'])
    df['J_lower']    = np.where(swap, df['state_one_J'], df['state_two_J'])
    df['conf_upper'] = np.where(swap, df['state_two_configuration'], df['state_one_configuration'])
    df['term_upper'] = np.where(swap, df['state_two_term'], df['state_one_term'])
    df['J_upper']    = np.where(swap, df['state_two_J'], df['state_one_J'])
    return df


def compare_lines(nist_path, theory_path, atom):
    nist   = load_nist_lines(nist_path)
    theory = load_theory_lines(theory_path)

    # Strip core from NIST configs using theory configs as the target set
    theory_confs = set(theory['conf_lower']) | set(theory['conf_upper'])
    for col in ['conf_lower', 'conf_upper']:
        nist[col] = nist[col].apply(lambda c: _strip_to_theory(c, theory_confs))

    for col in ['J_lower', 'J_upper']:
        theory[col] = theory[col].astype(str)
        nist[col]   = nist[col].astype(str)

    props_path = f"{atom}_Transition_Properties.csv"
    if Path(props_path).exists():
        props_df = pd.read_csv(props_path)
        props_df['transition_rate_uncertainty'] = pd.to_numeric(props_df['transition_rate_uncertainty'], errors='coerce')
        # state_one is always upper, state_two is always lower in {atom}_Transition_Properties.csv
        props_df = props_df.rename(columns={
            'state_one_configuration': 'conf_upper', 'state_one_term': 'term_upper', 'state_one_J': 'J_upper',
            'state_two_configuration': 'conf_lower', 'state_two_term': 'term_lower', 'state_two_J': 'J_lower',
        })
        props_df['J_upper'] = props_df['J_upper'].astype(str)
        props_df['J_lower'] = props_df['J_lower'].astype(str)
        theory = theory.merge(
            props_df[['conf_upper', 'term_upper', 'J_upper',
                       'conf_lower', 'term_lower', 'J_lower',
                       'transition_rate_uncertainty']],
            on=['conf_upper', 'term_upper', 'J_upper',
                'conf_lower', 'term_lower', 'J_lower'],
            how='left'
        )

    # Merge on conf, term (no *), J only — energy is shown for comparison only
    merged = theory.merge(
        nist,
        left_on  =['conf_lower', 'term_lower', 'J_lower',
                   'conf_upper', 'term_upper', 'J_upper'],
        right_on =['conf_lower', 'term_lower_key', 'J_lower',
                   'conf_upper', 'term_upper_key', 'J_upper'],
        how='inner',
    )

    wl_col  = 'ritz_wl_vac(nm)' if 'ritz_wl_vac(nm)' in merged.columns else 'obs_wl_vac(nm)'
    unc_col = 'unc_ritz_wl'     if 'unc_ritz_wl'     in merged.columns else 'unc_obs_wl'

    base_cols = [
        'conf_lower', 'term_lower_x', 'J_lower',
        'E_lower', 'E_lower(cm-1)',
        'conf_upper', 'term_upper_x', 'J_upper',
        'E_upper', 'E_upper(cm-1)',
        'wavelength(nm)', wl_col, 'matrix_element', 'matrix_element_uncertainty',
        'transition_rate(s-1)', 'Aki(s^-1)', 'Acc',
    ]
    if 'transition_rate_uncertainty' in merged.columns:
        base_cols.append('transition_rate_uncertainty')
    rename_map = {
        'term_lower_x':         'term_lower',
        'term_upper_x':         'term_upper',
        'E_lower':              'E_lower_theory(cm-1)',
        'E_lower(cm-1)':        'E_lower_NIST(cm-1)',
        'E_upper':              'E_upper_theory(cm-1)',
        'E_upper(cm-1)':        'E_upper_NIST(cm-1)',
        'wavelength(nm)':       'wavelength_theory(nm)',
        wl_col:                 'wavelength_NIST(nm)',
        'transition_rate(s-1)':      'Aki_theory(s^-1)',
        'transition_rate_uncertainty': 'Aki_unc_theory(s^-1)',
        'Aki(s^-1)':            'Aki_NIST(s^-1)',
        'Acc':                  'Acc_NIST',
    }
    if unc_col in merged.columns:
        base_cols.insert(base_cols.index(wl_col) + 1, unc_col)
        rename_map[unc_col] = 'wavelength_uncertainty(nm)'

    result = merged[base_cols].rename(columns=rename_map)

    result['dE_lower(cm-1)'] = result['E_lower_theory(cm-1)'] - result['E_lower_NIST(cm-1)']
    result['dE_upper(cm-1)'] = result['E_upper_theory(cm-1)'] - result['E_upper_NIST(cm-1)']
    result['%diff_Aki(theory-NIST)'] = (
        (result['Aki_theory(s^-1)'] - result['Aki_NIST(s^-1)'])
        / result['Aki_NIST(s^-1)'] * 100
    ).round(2)
    result['matrix_element_uncertainty(%)'] = (result['matrix_element_uncertainty'] / result['matrix_element'] * 100).abs().round(2)

    _ACC_UNC = {'AA': 0.01, 'A+': 0.02, 'A': 0.03,
                'B+': 0.07, 'B': 0.10,
                'C+': 0.18, 'C': 0.25,
                'D+': 0.40, 'D': 0.50}
    result['Aki_unc_NIST(s^-1)'] = (result['Acc_NIST'].map(_ACC_UNC) * result['Aki_NIST(s^-1)']).round(2)

    has_theory_unc = 'Aki_unc_theory(s^-1)' in result.columns
    if has_theory_unc:
        result['unc_%diff(theory-NIST)'] = (
            (result['Aki_unc_theory(s^-1)'] - result['Aki_unc_NIST(s^-1)'])
            / result['Aki_unc_NIST(s^-1)'] * 100
        ).round(2)
        result['unc_quadrature(s^-1)'] = np.sqrt(
            result['Aki_unc_theory(s^-1)'] ** 2 + result['Aki_unc_NIST(s^-1)'] ** 2
        ).round(2)

    # Reorder columns
    wl_cols = ['wavelength_theory(nm)', 'wavelength_NIST(nm)']
    if 'wavelength_uncertainty(nm)' in result.columns:
        wl_cols.append('wavelength_uncertainty(nm)')
    result = result[[
        'conf_lower', 'term_lower', 'J_lower',
        'E_lower_theory(cm-1)', 'E_lower_NIST(cm-1)', 'dE_lower(cm-1)',
        'conf_upper', 'term_upper', 'J_upper',
        'E_upper_theory(cm-1)', 'E_upper_NIST(cm-1)', 'dE_upper(cm-1)',
        *wl_cols,
        'matrix_element', 'matrix_element_uncertainty(%)',
        *(['Aki_theory(s^-1)', 'Aki_unc_theory(s^-1)'] if has_theory_unc else ['Aki_theory(s^-1)']),
        'Aki_NIST(s^-1)', 'Aki_unc_NIST(s^-1)',
        *(['unc_%diff(theory-NIST)', 'unc_quadrature(s^-1)'] if has_theory_unc else []),
        'Acc_NIST', '%diff_Aki(theory-NIST)',
    ]]

    output_path = f"{atom}_lines_comparison.csv"
    result.to_csv(output_path, index=False)
    print(f"Matched {len(result)} of {len(theory)} theory transitions to NIST")
    print(f"Saved to {output_path}")
    
    return result


def print_statistics(result):
    pct = result['%diff_Aki(theory-NIST)'].replace([np.inf, -np.inf], np.nan)
    acc_order = ['AA', 'A+', 'A', 'B+', 'B', 'C+', 'C', 'D+', 'D', 'E']

    print(f"\n{'Acc':<6} {'N':>4} {'Mean %diff':>12} {'Min |%diff|':>13} {'Max |%diff|':>13}")
    print("-" * 54)

    # Per accuracy class
    present = [a for a in acc_order if a in result['Acc_NIST'].values]
    for acc in present:
        r = pct[result['Acc_NIST'] == acc].dropna()
        if r.empty:
            continue
        print(f"{acc:<6} {len(r):>4} {r.mean():>12.4f} {r.abs().min():>13.4f} {r.abs().max():>13.4f}")

    # Overall
    r = pct.dropna()
    print("-" * 54)
    print(f"{'All':<6} {len(r):>4} {r.mean():>12.4f} {r.abs().min():>13.4f} {r.abs().max():>13.4f}")


_ROMAN = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X']


def parse_atom(name):
    """Accept 'Y2', 'Y_II', or 'Y II' and return (atom_dir, nist_stem).
    atom_dir  e.g. 'Y2'
    nist_stem e.g. 'Y_II'
    """
    name = name.strip()
    # numeric suffix: Y2, Ba1, Al2
    m = re.fullmatch(r'([A-Za-z]+)(\d+)', name)
    if m:
        elem, num = m.group(1), int(m.group(2))
        return name, f"{elem}_{_ROMAN[num - 1]}"
    # roman suffix with underscore or space: Y_II, Ba_I, Al II
    m = re.fullmatch(r'([A-Za-z]+)[_ ](I{1,3}|IV|VI{0,3}|IX|X)', name, re.IGNORECASE)
    if m:
        elem, roman = m.group(1), m.group(2).upper()
        num = _ROMAN.index(roman) + 1
        return f"{elem}{num}", f"{elem}_{roman}"
    raise ValueError(f"Cannot parse atom name: {name!r}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 3:
        nist_path, theory_path = sys.argv[1], sys.argv[2]
        atom = Path(theory_path).stem.replace('_Transition_Rates', '')
    else:
        name = sys.argv[1] if len(sys.argv) == 2 else input("Atom (e.g. Y2 or Y_II)? ")
        atom_dir, nist_stem = parse_atom(name)
        nist_path   = f"{nist_stem}_NIST_Lines.csv"
        theory_path = f"{atom_dir}_Transition_Rates.csv"
        atom = atom_dir
    result = compare_lines(nist_path, theory_path, atom)
    print_statistics(result)