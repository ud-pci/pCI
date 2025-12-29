"""
Script to generate transition rates and lifetimes from energy levels and matrix elements.

Input files:
- {atom}_Energies.csv
- {atom}_Matrix_Elements_Theory.csv

Output files:
- {atom}_Transition_Rates.csv
- {atom}_Transition_Rates_Error_Check.csv
- {atom}_Lifetimes_Error_Check.csv
"""

import pandas as pd
import sys
import os


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
    filename_lifetimes = atom + '_Lifetimes_Error_Check.csv'
    filename_br_ratios = atom + '_Transition_Rates_Error_Check.csv'

    lifetimes_df = pd.DataFrame(columns=['state_configuration', 'state_term', 'state_J', 'lifetime_display'])
    br_ratios_df = pd.DataFrame(columns=['state_one_configuration', 'state_one_term', 'state_one_J',
                                         'state_two_configuration', 'state_two_term', 'state_two_J',
                                         'wavelength_display', 'matrix_element_display',
                                         'branching_ratio_display', 'transition_rate_display'])

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

        tr_rates[upper_state].append([lower_state, tr_rate, matrix_element, wavelength])

    # Calculate lifetimes and branching ratios
    lifetimes = []

    # Write transitions to text file for inspection
    with open('transitions.txt', 'w') as f:
        for config, rates in tr_rates.items():
            configuration = config.split(' ')[0]
            term = config.split(' ')[1][0:2]
            J = config.split(' ')[1][-1]

            # Calculate total decay rate
            total_rates = sum(rate[1] for rate in rates)

            # Calculate lifetime in nanoseconds
            lifetime = 1/total_rates * 1e9
            lifetimes.append([configuration, term, J, round(lifetime, 3)])

            # Write to transitions file
            f.write(config + ' ->\n')
            for rate in rates:
                f.write('      ' + rate[0] + ' ' + str(rate[1]) + ' ' + str(rate[2]) + ' ' + str(rate[3]) + '\n')

                # Add to branching ratios dataframe
                configuration2 = rate[0].split(' ')[0]
                term2 = rate[0].split(' ')[1][0:2]
                J2 = rate[0].split(' ')[1][-1]
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
            'lifetime_display': lifetime[3]
        }
        lifetimes_df.loc[len(lifetimes_df.index)] = row

    # Save output files
    br_ratios_df.to_csv(filename_br_ratios, index=False)
    print(f"{filename_br_ratios} has been written with {len(br_ratios_df)} transitions")

    lifetimes_df.to_csv(filename_lifetimes, index=False)
    print(f"{filename_lifetimes} has been written with {len(lifetimes_df)} states")


def main():
    """Main function to process command line arguments and run the script."""

    if len(sys.argv) < 2:
        print("Usage: python generate_transition_rates.py <atom_name>")
        print("Example: python generate_transition_rates.py Ba1")
        print("\nThis script expects the following files to exist:")
        print("  - {atom_name}_Energies.csv")
        print("  - {atom_name}_Matrix_Elements_Theory.csv")
        print("\nIt will generate:")
        print("  - {atom_name}_Transition_Rates.csv")
        print("  - {atom_name}_Transition_Rates_Error_Check.csv")
        print("  - {atom_name}_Lifetimes_Error_Check.csv")
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

    print("\n" + "=" * 60)
    print("Processing complete!")


if __name__ == "__main__":
    main()
