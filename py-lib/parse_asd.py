import pandas as pd
import requests
import sys
from io import StringIO

def generate_asd_url(spectrum):
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

    return new_df

def reformat_df_to_atomdb(asd_df):
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
    
    # Fill in empty cells 
    asd_df['state_configuration'].fillna(method='ffill', inplace=True)
    asd_df['state_term'].fillna(method='ffill', inplace=True)
    asd_df['energy'].fillna('N/A', inplace=True)
    asd_df['energy_uncertainty'].fillna(0, inplace=True)

    # Remove spaces in energies
    asd_df['energy'] = asd_df['energy'].str.replace(' ', '')

    # Add a '.' between orbitals in configurations
    

    # Replace reference ID column to 'FALSE' as a value of is_from_theory
    asd_df['is_from_theory'] = 'FALSE'

    return asd_df

def df_to_csv(df, filename):
    """
    This function converts a dataframe into a csv file
    """

    asd_df.to_csv(filename, index=False)

    print("Dataframe converted to csv")


if __name__ == "__main__":
    spectrum = input("Name of system? ")
    asd_url = generate_asd_url(spectrum)

    asd_df = pd.DataFrame()

    try:
        asd_df = generate_df_from_asd(asd_url)
    except:
        print(spectrum + ' does not exist on ASD')
        sys.exit()

    reformat_df_to_atomdb(asd_df)
    df_to_csv(asd_df, str(spectrum).replace(' ', '_'))
    
