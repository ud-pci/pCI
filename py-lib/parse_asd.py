from io import StringIO
import pandas as pd
import requests
import sys
import re

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
        newline = re.sub('\s{3,}', '  ', line).replace('  ', ',') + '\n'
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

    return new_df

def reformat_df_to_atomdb(asd_df, theory_J): # Modified
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

    # Only keep configurations where J is in theory results
    asd_df = asd_df[(((asd_df['state_J'].isin(theory_J['even'])) & (asd_df['state_term'].str[-1] != '*'))) |
                    (((asd_df['state_J'].isin(theory_J['odd'])) & (asd_df['state_term'].str[-1] == '*')))]

    # Fill in empty cells 
    asd_df['state_configuration'].fillna(method='ffill', inplace=True)
    asd_df['state_term'].fillna(method='ffill', inplace=True)
    asd_df['energy'].fillna('N/A', inplace=True)
    asd_df['energy_uncertainty'].fillna(0, inplace=True)

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
    asd_df = asd_df[asd_df['state_J'] != ""]
    asd_df["state_J"] = asd_df["state_J"].str.split(',').str[0]
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