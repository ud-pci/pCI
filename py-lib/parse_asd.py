from io import StringIO
import numpy as np
import pandas as pd
import requests
import sys
import re
from collections import Counter

pd.options.mode.chained_assignment = None

def Convert_Type(s): # detect and correct the 'type' of object to 'float', 'integer', 'string' while reading data
    s = s.replace(" ", ".")
    try:
        f = float(s)
        i = int(f)
        return i if i == f else f
    
    except ValueError:
        try: return eval(s)
        except: return s
    

def RemoveInitial(l): # find core config to remove it for minimal representation
    s =""
    n=0
    for i in l:
        s+=i
        if i==".":n+=1
        if n==1:break  
    return s

def convert_res_to_csv(filename): # Modified convert .RES to csv

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    csvfile = "DATA_Filtered/UD/"+filename.split('/')[-1].split('.')[0] + '.csv' # output file path

    f = open(csvfile, 'w')
    f.write('n, conf, , term, E_n (a.u.), DEL (cm^-1), S, L, , gf, conf%, conf2, , conf2% \n')

    for line in lines[1:]:
        newline = re.sub(r'\s{3,}', '  ', line).replace('  ', ',') + '\n'
        if newline[0] == ',':
            newline = newline[1:]
        
        string_to_list = newline.split(",")
        modified_list = [str(Convert_Type(i)) for i in string_to_list]
        newline = ",".join(modified_list)
        f.write(newline)
    f.close()

    return


def generate_asd_url(spectrum): # ex: spectrum = Sr I
    """
    This function generates the url for NIST Atomic Spectral Database levels data
    """ 
    url = "https://physics.nist.gov/cgi-bin/ASD/energy1.pl?"

    spectrum_post = 'spectrum=' + str(spectrum).replace(' ', '+') + '&submit=Retrieve+Data&'
    post_req = ('units=0' + '&' 
                + 'format=2' + '&' # read asd in csv text format
                + 'output=0' + '&' 
                + 'page_size=15' + '&' 
                + 'multiplet_ordered=0' + '&' 
                + 'conf_out=on' + '&' 
                + 'term_out=on' + '&' 
                + 'level_out=on' + '&'
                + 'unc_out=1' + '&'
                + 'j_out=on' + '&' 
                + 'lande_out=on' + '&'
                + 'perc_out=on' + '&' 
                + 'biblio=on' + '&'
                + 'temp=')

    full_url = url + spectrum_post + post_req

    return full_url

def generate_df_from_asd(url):
    """
    This function returns a dataframe containing the energies from the NIST Atomic Spectra Database
    """
    r = requests.get(url)
    asd_df = pd.read_table(StringIO(r.text), delimiter = ',')
    new_df = asd_df[['Configuration','Term','J','Level (cm-1)', 'Uncertainty (cm-1)', 'Reference']].replace({'=','"'},'', regex=True)

    # remove rows indicating ionization limit
    new_df = new_df[~new_df['Term'].str.contains("Limit", regex=True, na=False)]
    
    # find all orbitals with fully filled orbitals
    pattern = r"\d+[spdfg](?:2|6|10|14)"
    all_orbitals = new_df['Configuration'].str.findall(pattern)
    
    # Count frequencies
    counts = Counter([orb for sublist in all_orbitals for orb in sublist])
    
    # Keep only those that appear in all rows
    nrows = len(new_df)
    common_full_orbitals = {orb for orb, count in counts.items() if count == nrows}
    
    def strip_common(config):
        for orb in common_full_orbitals:
            config = re.sub(rf"{orb}\.?", "", config)
        config = re.sub(r"\.{2,}", ".", config).strip(".")
        return config
    
    new_df['Configuration'] = new_df['Configuration'].apply(strip_common)
    
    return new_df

def reformat_df_to_atomdb(asd_df, theory_J=None):
    """
    This function reformats the ASD dataframe for use in the UD Atom database
    """

    # Rename columns to AtomDB format
    asd_df.rename(columns={'Configuration':'state_configuration',
                            'Term':'state_term',
                            'J':'state_J',
                            'Level (cm-1)':'energy',
                            'Uncertainty (cm-1)':'energy_uncertainty',
                            'Reference':'is_from_theory'},
                            inplace=True)

    # Expand rows with multiple J values
    # Filter out rows with missing J first
    asd_df = asd_df[asd_df['state_J'] != ""]
    expanded_rows = []
    for idx, row in asd_df.iterrows():
        j_values = str(row['state_J']).split(',')
        for j_val in j_values:
            new_row = row.copy()
            # Keep J as string for now (theory_J contains strings)
            new_row['state_J'] = j_val.strip()
            expanded_rows.append(new_row)

    asd_df = pd.DataFrame(expanded_rows).reset_index(drop=True)

    # Only keep configurations where J is in theory results
    if theory_J:
        asd_df = asd_df[(((asd_df['state_J'].isin(theory_J['even'])) & (asd_df['state_term'].str[-1] != '*'))) |
                        (((asd_df['state_J'].isin(theory_J['odd'])) & (asd_df['state_term'].str[-1] == '*')))]

    # Fill in empty cells
    asd_df['state_configuration'] = asd_df['state_configuration'].replace(r'^\s*$', np.nan, regex=True)
    asd_df['state_configuration'] = asd_df['state_configuration'].ffill()
    asd_df['state_term'] = asd_df['state_term'].replace(r'^\s*$', np.nan, regex=True)
    asd_df['state_term'] = asd_df['state_term'].ffill()
    asd_df['energy'] = asd_df['energy'].fillna('N/A')
    asd_df['energy_uncertainty'] = asd_df['energy_uncertainty'].fillna(0)

    # Replace reference ID column to 'FALSE' as a value of is_from_theory
    asd_df['is_from_theory'] = 'FALSE'

    return asd_df

def NIST_Discrepancies(asd_df,ri=False):
    # ri : replace initials
    # All Filters and Discrepencies Here
    # if ri==True: 
    #     s = RemoveInitial(asd_df["state_configuration"].values[0])
    #     asd_df["state_configuration"] = asd_df["state_configuration"].str.replace(s, "") # replace initial from the configuration

    # missing J
    asd_df["term_original"] = asd_df["state_term"]
    asd_df["state_J"] = asd_df["state_J"].apply(lambda x: Convert_Type(x))
    asd_df["state_term"] = asd_df["state_term"].str.replace('[a-z]', '',regex=True) # replace x,y,z from the terms
    asd_df["state_term"] = asd_df["state_term"].str.replace(' ', '',regex=True) # replace spaces from the terms
    asd_df["state_term"] = asd_df["state_term"].str.replace("\\(.*?\\)","",regex=True) # replace '()' and whats inside it from the configuration
    asd_df["state_term"] = asd_df["state_term"].str.replace("\\[.*?\\]","",regex=True) # replace '[]' and whats inside it from the configuration

    # asd_df["state_configuration"] = asd_df["state_configuration"].str.replace(".", " ") # replace '.' with " " from the configuration
    asd_df["state_configuration"] = asd_df["state_configuration"].str.replace("n", "0") # replace "4s np" with "4s 0p" from the configuration
    asd_df["state_configuration"] = asd_df["state_configuration"].str.replace(".\\(.*?\\)","",regex=True) # replace '.()' and whats inside it from the configuration
    asd_df["state_configuration"] = asd_df["state_configuration"].str.replace("\\(.*?\\)","",regex=True) # replace '()' and whats inside it from the configuration
    asd_df["state_configuration"] = asd_df["state_configuration"].str.replace("\\<.*?\\>","",regex=True) # replace '<>' and whats inside it from the configuration
    asd_df["state_configuration"] = asd_df["state_configuration"].str.replace("\\.\\.",'.') # replace multiple '.' with single '.' from the configuration

    asd_df["energy"] = asd_df["energy"].replace(r'^s*$', float('NaN'), regex = True) # replace blank energies with nan from the configuration
    asd_df.dropna(inplace = True) # remove states with abscent energy values

    return asd_df

def df_to_csv(asd_df, filename, parity=None,ri=False):
    """
    This function converts a dataframe into a csv file
    """
    
    # Filtering Odd and Even data before converting it to csv
    if parity=="odd":
        asd_df = asd_df[asd_df['state_term'].str[-1] == '*']
        asd_df["state_term"] = asd_df["state_term"].str.replace("*","")
        filename = filename.replace(" ","_")+"_NIST_Odd"+".csv"

    elif parity=="even":
        asd_df = asd_df[asd_df['state_term'].str[-1] != '*']
        filename = filename.replace(" ","_")+"_NIST_Even"+".csv"

    else:
        filename = filename.replace(" ","_")+"_NIST_All"+".csv"

    asd_df = NIST_Discrepancies(asd_df,ri)

    asd_df.to_csv(filename, index=False)


def generate_lines_url(spectrum, min_accur=''):
    """
    This function generates the url for NIST Atomic Spectral Database lines data.
    min_accur: minimum accuracy class filter (e.g. 'AA', 'A', 'B', 'C', 'D', 'E'), or '' for no filter
    """
    url = "https://physics.nist.gov/cgi-bin/ASD/lines1.pl?"

    spectra_post = 'spectra=' + str(spectrum).replace(' ', '+') + '&submit=Retrieve+Data&'
    post_req = ('output_type=0' + '&'
                + 'low_w=' + '&'
                + 'upp_w=' + '&'
                + 'unit=1' + '&'
                + 'de=0' + '&'
                + 'I_scale_type=1' + '&'
                + 'format=2' + '&'  # CSV format
                + 'line_out=0' + '&'
                + 'remove_js=on' + '&'
                + 'en_unit=0' + '&'
                + 'output=0' + '&'
                + 'bibrefs=1' + '&'
                + 'page_size=15' + '&'
                + 'show_obs_wl=1' + '&'
                + 'show_calc_wl=1' + '&'
                + 'unc_out=1' + '&'
                + 'order_out=0' + '&'
                + 'show_av=3' + '&'
                + 'tsb_value=0' + '&'
                + 'A_out=0' + '&'
                + 'intens_out=on' + '&'
                + 'allowed_out=1' + '&'
                + 'forbid_out=1' + '&'
                + 'min_accur=' + str(min_accur) + '&'
                + 'conf_out=on' + '&'
                + 'term_out=on' + '&'
                + 'enrg_out=on' + '&'
                + 'J_out=on')

    return url + spectra_post + post_req


def generate_lines_df_from_asd(url):
    """
    This function returns a dataframe containing the lines data from the NIST Atomic Spectra Database.
    Columns: obs_wl_vac(nm), unc_obs_wl, ritz_wl_vac(nm), unc_ritz_wl, intens, Aki(s^-1), Acc,
             Ei(cm-1), Ek(cm-1), conf_i, term_i, J_i, conf_k, term_k, J_k, Type, tp_ref, line_ref
    """
    r = requests.get(url)
    with open('lines_raw.csv', 'w') as f:
        f.write(r.text)
    if r.text.lstrip().startswith('<'):
        print('No lines data available on ASD for this spectrum')
        return pd.DataFrame()
    lines_df = pd.read_csv(StringIO(r.text))

    # Strip leading/trailing whitespace from column names
    lines_df.columns = lines_df.columns.str.strip()

    # Drop trailing empty column produced by trailing comma in header
    lines_df = lines_df.loc[:, ~lines_df.columns.str.startswith('Unnamed')]

    # Strip the ="..." Excel formula quoting (same approach as generate_df_from_asd)
    lines_df = lines_df.replace({'=', '"'}, '', regex=True)

    return lines_df


def lines_df_to_csv(lines_df, spectrum):
    """
    This function saves a lines dataframe to a CSV file named after the spectrum,
    keeping only the observed wavelength, uncertainty, Aki, accuracy, and lower/upper
    level configuration, term, J, and energy.
    """
    if 'unc_ritz_wl' in lines_df.columns:
        wl_col, unc_wl_col = 'ritz_wl_vac(nm)', 'unc_ritz_wl'
    elif 'unc_obs_wl' in lines_df.columns:
        wl_col, unc_wl_col = 'obs_wl_vac(nm)', 'unc_obs_wl'
    else:
        wl_col = 'ritz_wl_vac(nm)' if 'ritz_wl_vac(nm)' in lines_df.columns else 'obs_wl_vac(nm)'
        unc_wl_col = None
    cols = ['conf_i', 'term_i', 'J_i', 'Ei(cm-1)',
            'conf_k', 'term_k', 'J_k', 'Ek(cm-1)',
            wl_col, 'Aki(s^-1)', 'Acc']
    new_cols = ['conf_lower', 'term_lower', 'J_lower', 'E_lower(cm-1)',
                'conf_upper', 'term_upper', 'J_upper', 'E_upper(cm-1)',
                wl_col, 'Aki(s^-1)', 'Acc']
    if unc_wl_col is not None:
        cols.insert(cols.index(wl_col) + 1, unc_wl_col)
        new_cols.insert(new_cols.index(wl_col) + 1, unc_wl_col)
    filename = str(spectrum).replace(' ', '_') + '_NIST_Lines.csv'
    lines_df[cols].rename(columns=dict(zip(cols, new_cols))).to_csv(filename, index=False)
    return filename


if __name__ == "__main__":
    spectrum = input("Name of system? ")
    data_type = input("Get energies or lines data? (1=energies, 2=lines): ")

    if data_type == '1':
        asd_url = generate_asd_url(spectrum)

        asd_df = pd.DataFrame()

        try:
            asd_df = generate_df_from_asd(asd_url)
        except:
            print(spectrum + ' does not exist on ASD')
            sys.exit()

        reformat_df_to_atomdb(asd_df)
        filename = df_to_csv(asd_df, str(spectrum).replace(' ', '_'))
        print(f'Written to {filename}')
    elif data_type == '2':
        min_accur = input("Minimum uncertainty (e.g. AA, A, B, C, D, E — press Enter for no filter): ")
        lines_url = generate_lines_url(spectrum, min_accur)
        try:
            lines_df = generate_lines_df_from_asd(lines_url)
        except:
            print(spectrum + ' does not exist on ASD')
            sys.exit()
        if not lines_df.empty:
            filename = lines_df_to_csv(lines_df, spectrum)
            print(f'Written to {filename}')