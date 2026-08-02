""" Atomic Term Symbol Calculator, Charles Cheung (2024)

This script allows the user to print all possible atomic term symbols from an electron configuration.
Electron configurations should be formatted with a '.' or ' ' separating different orbitals (e.g. 2s.2p3 or 2s 2p3).

This file can also be imported as a module and contains the following functions:

    * calc_microstates - returns total number of microstates
    * calc_term_symbols - the main function of the script

"""
import re
import math
from itertools import product, combinations
import numpy as np
import pandas as pd
from fractions import Fraction

'''
This python script determines all possible terms from an electron configuration
'''


def calc_microstates(num_orbital_positions, num_electrons):
    """ Calculates the total number of microstates given the number of orbital positions and number of electrons
    
    The total number of microstates is calculated using the following formula:
    N = n!/(e!(n-e)!), where
    n = number of possible orbital positions
    e = number of electrons
    
    Parameters
    ----------
    num_orbital_positions : integer
        the number of orbital positions
    num_electrons : integer
        the number of electrons
    
    Returns
    -------
    num_microstates : integer
        the total number of microstates
    
    """
    num_microstates = int(math.factorial(num_orbital_positions) / (math.factorial(num_electrons) * math.factorial(num_orbital_positions - num_electrons)))
    
    return num_microstates

def calc_term_symbols(configuration):
    """ Calculates all possible term symbols for a given electron configuration
    
    Parameters
    ----------
    configuration : string
        the electron configuration formatted with '.' or ' ' between each orbital (e.g. 2s.2p3 or 2s 2p3)
    
    Returns
    -------
    term_symbols : list
        a list of strings representing possible term symbols
    """
    
    print('Calculating atomic term symbols for ' + configuration + '...')
    
    # Separate the configuration into an array of orbitals
    orbitals = configuration.replace(' ', '.').split('.')
        
    # Define orbital angular momentum quantum numbers and letters, and max number of orbitals
    Ls = 'SPDFGHIKLMNOQRTUVWXYZ'
    l_values = {}
    L_terms = {}
    max_occupancies = {}
    for i in range(len(Ls)):
        l_values[Ls[i].lower()] = i
        L_terms[i] = Ls[i]
        max_occupancies[Ls[i].lower()] = 2*(2*i+1)
        
    num_electrons = 0
    num_microstates = 1
    
    active_orbitals = []
    num_orbital_positions = []
    configs = []
    occupancies = {}
    
    ml = [] # all possible values of ml = [-l, -l+1, ..., l-1, l] for each orbital

    for orbital in orbitals:
        l = re.findall('[a-zA-Z]+', orbital)[0]
        nn = int(re.findall('[0-9]+', orbital)[0])
        try:
            qq = int(re.findall('[0-9]+', orbital)[1])
        except IndexError:
            qq = 1
        # Check if shell is filled and skip to next orbital if filled
        if qq == max_occupancies[l]:
            continue
        for i in range(int(qq)):
            configs.append([nn,l])
        ml += list(range(-l_values[l], l_values[l]+1))
        num_orbital_positions.append(max_occupancies[l])
        num_microstates *= calc_microstates(max_occupancies[l], qq)
        num_electrons += qq
        occupancies[l] = qq
        active_orbitals.append(l)
        
    # If no active orbitals, i.e. all orbitals filled, term would be 1S0
    if not active_orbitals:
        return ['1S0']
        
    # All possible values of m_s. We work with 2*m_s instead of m_s so that we only have to handle integers
    ms2 = list(range(-num_electrons,num_electrons+1))
    if num_electrons % 2 == 0:
        ms2 = [i for i in ms2 if i % 2 == 0]
    else:
        ms2 = [i for i in ms2 if i % 2 != 0]
    
    # Generate a list of microstates by taking combinations of electrons orbiting possible orbital positions
    num_orbital_positions = sum(num_orbital_positions)
    microstates = list(combinations(range(num_orbital_positions), num_electrons))

    # List each possible microstate
    configs = [] # full expanded microstate configuration
    reduced_configs = [] # reduced microstate configuration in orbitals
    irr_configs = [] # further reduced configuration of orbitals, where (-1/2,1/2)=(1/2,-1/2) is set to (2,2)
    
    # Create an array of combinations of electron spins
    list_electrons = list(product(*[[-1,1]] * num_electrons))
        
    # number of zeros that should be found for each microstate configuration
    num_zeros = num_orbital_positions - num_electrons
        
    # loop over microstates
    for microstate in microstates:
        # loop over list of combinations of electron spins
        for electrons in list_electrons:
            # initialize a microstate configuration, i.e. an array of orbital positions for electrons to occupy
            config = np.zeros(num_orbital_positions, dtype=int)
            
            # add electrons to each orbital position in the microstate configuration
            for i in range(len(microstate)):
                config[microstate[i]] += electrons[i]
                
            # remove instances where microstate configuration does not have the right number of electrons
            # l_array is structured in the following way: [s1/2, p-1/2, p0, p1/2, d-3/2, d-1/2, d0, d1/2, d3/2, ...]
            actual_occupancies = {}
            
            first_active_orbital = active_orbitals[0]
            actual_occupancies[first_active_orbital] = sum(abs(config[0:max_occupancies[first_active_orbital]]))
            sum_occ = max_occupancies[first_active_orbital]
            end_occ = sum_occ
            for active_index in range(1,len(active_orbitals)):
                start_occ = end_occ
                end_occ = sum_occ + max_occupancies[active_orbitals[active_index]]
                sum_occ += max_occupancies[active_orbitals[active_index]]
                actual_occupancies[active_orbitals[active_index]] = sum(abs(config[start_occ:end_occ]))

            skip_config = False
            for li in active_orbitals:
                if actual_occupancies[li] != occupancies[li]:
                    skip_config = True

            if skip_config:
                continue
                            
            # apply Pauli exclusion principle by removing instances where an orbital is occupied by two electrons with the same spin
            reduced_config = [list(config[i:i + 2]) for i in range(0,len(config), 2)]
            
            for orbital in reduced_config:
                if orbital[0] == orbital[1] and orbital[0] == -1 or orbital[0] == 1:
                    skip_config = True
            if skip_config:
                continue
            
            # remove duplicate instances of orbitals with paired electrons
            irr_config = []
            for orbital in reduced_config:
                if orbital == [-1,1] or orbital == [1,-1]:
                    irr_config.append(2)
                elif orbital == [-1,0] or orbital == [0,-1]:
                    irr_config.append(-1)
                elif orbital == [1,0] or orbital == [0,1]:
                    irr_config.append(1)
                else:
                    irr_config.append(0)
            
            if reduced_config in reduced_configs:
                skip_config = True
            if skip_config:
                continue

            # add microstate configuration to list of configurations
            if list(config) not in configs and irr_config not in irr_configs:
                configs.append(list(config))
                reduced_configs.append(reduced_config)
                irr_configs.append(irr_config)
    
    # assign max value of m_l
    ml_max = 2*max(ml)+1
    
    # create a table of all possible different microstates
    # each row represents a total orbital magnetic quantum number
    # each column represents a total spin magnetic quantum number
    df = pd.DataFrame(0, index=np.arange(-ml_max,ml_max+1), columns=ms2)

    for index in range(len(reduced_configs)):
        MS2 = sum(configs[index])
        ML = sum(abs(irr_configs[index][i])*ml[i] for i in range(len(ml)))
        df.at[ML, MS2] += 1
    
    # Now we need to start subtracting term symbols from the table
    term_descriptions = []
        
    # drop rows with all zeros
    df = df.loc[(df != 0).any(axis=1)]
    
    # extract M_L and M_S from the table
    while len(df.index) != 0:        
        # save indices of columns with no zeros
        col_indices = []  
        for col in df.columns:
            num_zeros = df[col].eq(0).sum()
            if num_zeros == 0:
                col_indices.append(col)
        
        # remove 1 from each row of columns with no zeros and record M_L, M_S
        M_S = len(col_indices)
        for col in col_indices:
            M_L = abs(df.index[0])
            df[col] = df[col] - 1
        
        # add M_L, M_S, and term to array
        term_descriptions.append([M_L, M_S, str(M_S) + L_terms[abs(M_L)]])
        
        # drop rows with all zeros
        df = df.loc[(df != 0).any(axis=1)]
        
    # find each multiplet term symbol including J
    term_symbols = []
    for terms in term_descriptions:
        L = terms[0]     # total orbital angular momentum
        M_S2 = terms[1]  # spin multiplicity 2S+1
        S = (M_S2 - 1)/2 # total spin
        term = terms[2]
        
        # total angular momentum |L-S| <= J <= |L+S|
        J_min = abs(L-S)
        J_max = abs(L+S)
        Js = np.arange(J_min, J_max+1, 1)
        for J in Js:
            term_J = term + str(Fraction(J))
            if term_J not in term_symbols:
                term_symbols.append(term_J)

    return term_symbols

if __name__ == '__main__':
    terms = calc_term_symbols(input('Enter configuration: '))
    print('# of terms: ' + str(len(terms)))
    print(terms)
    