# TO RUN

**Add .RES file of the atom naming it as it should be given to NIST input**

Import Dependencies
```python
from UDRead import *
from parse_asd import *
from get_atomic_term import *
```

**Inputs Required:**

**parity** - allows to select specific parity state from NIST database to compare UD data with.

**order** - mention order to get by energy or term ordered final output.

```python
parity = "even" # can be "even" or "odd"
order = "E"     # can be "E" for energy or "T" for Term ordered
```

**Atom** - Name of Atom same as we give in NIST website

**ri** - To replace ground state configuration which NIST gives sometime throughout data. For example for "Ca I" NIST adds "3p6" for all levels, ri==True  "3p6.4s2" --> "4s2". This needs to check manually in NIST database.

**nist_max** - maximum number of states to be considered from NIST data. Ideally it should be such that the last configuration of UD data is the nist_max<sup>th</sup>  state in NIST data as somtimes if more NIST states are provided sometime wrong configurations in higher states start matching better with UD state than the actual state it corresponds.


```python
Atom = "Ca I"
ri = True # can be True of False
nist_max = 51 if parity=="odd" else 45
```

*Some examples of inputs are provided in Example.py*

**Producing CSV files:**

**NIST** - This part extract NIST data from the website and makes csv for "odd" or "even" as provided and corrects the format of NIST data. The csv file is stored in "DATA_Filtered/NIST/"

```python
Name = Atom.replace(" ","_")

url_nist = generate_asd_url(Atom,parity)
data_nist = generate_df_from_asd(url_nist)
data_nist = reformat_df_to_atomdb(data_nist)
df_to_csv(data_nist,"DATA_Filtered/NIST/"+Atom,parity,ri=ri)
```

**UD** - This part converts .RES file to csv and store it to "DATA_Filtered/UD/"

```python
convert_res_to_csv("DATA_RAW/UD/"+Name+"_UD_"+parity.capitalize()+".RES")

# Paths of csv's stored in "DATA_Filtered/"
path_nist = "DATA_Filtered/NIST/"+Name+"_NIST_"+parity.capitalize()+".csv"
path_ud = "DATA_Filtered/UD/"+Name+"_UD_"+parity.capitalize()+".csv"
```


**Finding correspondance:**

The correspondance is find between NIST and UD data. The levels which are missing in both NIST and UD data are then added to the end of data. The final output is then stored in `"DATA_Output/"`

```python
## Main Code
data = MainCode(path_nist, path_ud,nist_max,Ordering=order)

## Finding Missing Levels
data_final = Missing_Levels(data)

## Export
path = "DATA_Output/"+Name+"_"+parity.capitalize()+".txt" #specify path for export
ConvertToTXT(data_final,path)
```