from UDRead import *
from parse_asd import *
from get_atomic_term import *

parity = "even"
order = "E"
# nist_max=0

Atom = "Ca I"
ri = True
nist_max = 51 if parity=="odd" else 45

# Atom = "Sr I"
# ri = False
# nist_max = 50 if parity== "odd" else 50

# Atom = "In II"
# ri = False
# nist_max = 36 if parity== "odd" else 32

# Atom = "In I"
# ri = False
# nist_max = 34 if parity== "odd" else 27

# Atom = "Y II"
# ri = False
# nist_max = 43 if parity== "odd" else 40

# Atom = "Cs VII"
# ri = True # fac=2 for Sr, fac=4 for Ca
# nist_max = 32 if parity== "odd" else 31

Name = Atom.replace(" ","_")

# # DATA NIST
# Store Filtered data of even or odd parity in DATA_Filtered/NIST/ 
url_nist = generate_asd_url(Atom,parity)
data_nist = generate_df_from_asd(url_nist)
data_nist = reformat_df_to_atomdb(data_nist)
df_to_csv(data_nist,"DATA_Filtered/NIST/"+Atom,parity,ri=ri)

## DATA UD
# Sorted formatting, Converted .RES to .csv and stored Filtered data of even or odd parity in DATA_Filtered/UD/ 
convert_res_to_csv("DATA_RAW/UD/"+Name+"_UD_"+parity.capitalize()+".RES")

path_nist = "DATA_Filtered/NIST/"+Name+"_NIST_"+parity.capitalize()+".csv"
path_ud = "DATA_Filtered/UD/"+Name+"_UD_"+parity.capitalize()+".csv"

## Main Code
data = MainCode(path_nist, path_ud,nist_max,Ordering=order)

## Finding Missing Levels
data_final = Missing_Levels(data)


# # Export
path = "DATA_Output/"+Name+"_"+parity.capitalize()+".txt" #specify path for export
ConvertToTXT(data_final,path)
