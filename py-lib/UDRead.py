import pandas as pd
import numpy as np
import os
from fractions import Fraction

pd.set_option("display.max_rows", None, "display.max_columns", None)

def LongSubString(arr): # Longest common substring from array of strings
    n,s = len(arr),arr[0]
    l,res = len(s),""
    for i in range(l):
        for j in range(i + 1, l + 1):
            stem,k = s[i:j],1
            for k in range(1, n):
                if stem not in arr[k]:break
            if (k + 1 == n and len(res) < len(stem)): res = stem
                
    return res


def Convert_Type(s): # detect and correct the 'type' of object to 'float', 'integer', 'string' while reading data
    s = s.replace(" ", "")
    try:
        f = float(s)
        i = int(f)
        return i if i == f else f
    
    except ValueError:
        return s
    
    
def correct_config(config,val): # correction to configuraion by value "val", ex: '5s4d' --> '5s5d' if val = 1
    List = list(config)
    idx = [i for i in range(len(List)) if List[i].isdigit()==True]

    if len(idx)==2:
        List[idx[-1]] = str(int(List[idx[-1]])+val)
        return ''.join(List)
    
    else: return config


def inverse_config(config_nist):
    # 5s6p <--> 6p5s as somtime configuration text in ud and nist are inverse of each other
    # so correspondence will be checked with both 5s6p and 6p5s
    n=0
    first = ""
    second = ""
    for i in(config_nist):
        if n==0: first+=i

        if n==1:
            second+=i

        if i.isalpha()==True:
            n+=1
    return second+first


def term_correct(term_ao,corr=1): # correction to term symbol by value "corr", ex : 2D --> 3D if corr = 1
    term_symbol = list(term_ao)
    term_symbol[0]=str(int(term_symbol[0])+corr)
    term_symbol="".join(term_symbol)
    
    return term_symbol


def Sep_Config(config): # sperate confguration : "5s6p" return --> 5,'s',6,'p'
    d1,a1,d2,a2 = '','','',''
    for i in config:
        if            a1=='' and i.isdigit()==True: d1+=i
        if d1!='' and d2=='' and i.isalpha()==True: a1+=i
        if a1!='' and a2=='' and i.isdigit()==True: d2+=i
        if d2!='' and            i.isalpha()==True: a2+=i

    return int(d1),a1,int(d2),a2


def Data_UD(ith,df_ao,corr_config=True): # reading filtered csv data
    
    OrgConfig1 = df_ao[ith][1] # original first possible configuration from data
    config1 = df_ao[ith][1] if corr_config==False else  df_ao[ith][-1] # corrected configuration
    config2 = df_ao[ith][10] # second possible configuration
    Term = df_ao[ith][2][0:2] # Term
    J = str(Fraction(df_ao[ith][7])) # J value
    Level = round(df_ao[ith][4],3) # Level
    Levelau = df_ao[ith][3] # Level in atomic units
    Per1 = df_ao[ith][9] # Percentage Contribution of first configuration
    Per2 = df_ao[ith][11] # Percentage Contribution of second configuration
    uncer_ud = df_ao[ith][12] # Percentage Contribution of second configuration

    return OrgConfig1,config1,config2,Term,J,Level,Levelau,Per1,Per2,uncer_ud


def Data_Nist(ith,df_nist): # Reading filtered csv data
    config_nist = df_nist[ith][0]
    term_nist = df_nist[ith][1]
    j_nist = df_nist[ith][2]
    level_nist = round(df_nist[ith][3],3)
    uncer_nist = round(df_nist[ith][4],3)
    return config_nist,term_nist,j_nist,level_nist,uncer_nist


def Dataframe(path_nist,path_ao,nist_max,parity,fac,gs_parity):
    # read csv
    df_nist = pd.read_csv(path_nist)
    df_ao = pd.read_csv(path_ao).fillna('') # fill blank for nan
    df_ao = df_ao.sort_values(by=df_ao.columns[4], ascending=True)

    # Storing as numpy arrays
    df_nist = df_nist.values
    df_ao = df_ao.values

    ref_E = df_nist[:,3][0] if parity!=gs_parity else 0
    if parity!=gs_parity: df_nist[:,3] = df_nist[:,3]- ref_E # making zero reference energy in odd parity states in NIST Data

    new_config = Correct_Config(df_ao,df_nist,nist_max,fac)
    df_ao = np.append(df_ao, np.stack([new_config],axis=1), axis=1)
    
    return df_nist,df_ao,ref_E


## Main Correspondance Part
# Update : Automated Corrections to Configuration
# Need generalisation and final version to be written again from scratch after automating manual correction
#################################################################################################################
def FindJthAll(i,df_nist,df_ud,fac,corr_config=True): # find ao config corresponding to ith nist config
    j1,j2,Eperc1,Eperc2=[],[],[],[]
    config_nist, term_nist, j_nist, level_nist, uncer_nist = Data_Nist(i,df_nist)

    for j in range(len(df_ud)):
        OrgConfig1,config1_ao, config2_ao, term_ao, j_ao, level_ao,Levelau, per1,per2,uncer_ud = Data_UD(j,df_ud,corr_config)
        Ediff = round(abs(level_nist-level_ao),1)
        Eperc = round(100*Ediff/level_nist,2) if level_nist!=0 else 0
        
        # First Pass : if J_Nist = J_AO and Percentage Difference is less than 2
        # fac = 2
        EpercLim = fac-i*(fac-0.75)/len(df_ud)
        if str(j_nist)==j_ao and Eperc<EpercLim:
            inv_config = inverse_config(config_nist)
            
            # Second Pass : if Term_NIST = Term_AO
            if term_nist==term_ao:
                # Third Pass : 
                if config_nist==config1_ao or inv_config==config1_ao:j1.append(j),Eperc1.append(Eperc)
                if config_nist==config2_ao or inv_config==config2_ao:j2.append(j),Eperc2.append(Eperc)
            
            # Fourth Pass : If Term_NIST != Term_AO --> Check using term correction by +-1
            if term_nist!=term_ao: # correction to terms in ao results
                if config_nist==config1_ao or inv_config==config1_ao:
                    if term_nist==term_correct(term_ao,1):j1.append(j),Eperc1.append(Eperc)
                    if term_nist==term_correct(term_ao,-1):j1.append(j),Eperc1.append(Eperc)
            
                if config_nist==config2_ao or inv_config==config2_ao:
                    if term_nist==term_correct(term_ao,1):j2.append(j),Eperc2.append(Eperc)
                    if term_nist==term_correct(term_ao,-1):j2.append(j),Eperc2.append(Eperc)
    
    return j1,j2,Eperc1,Eperc2 # j1,j2 is the index of ud configurations 1 and 2 corresponding to NIST


# Fifth Pass : to choose between j1 and j2 if both present to find the best correspondance
def ChooseJth(J1,J2,Eperc1,Eperc2,df_nist,df_ud,nist_max,fac,corr_config=True):

    JJ1=0
    EEEperc1=0
    j_final,jf=-1,0 # jf is "which config selected" : jf=0 --> no config, jf=1 --> 1st Config, jf=2 --> 2nd Config
    
    # Choosing minimum percentage defference index if more than one values present in jj1 and jj2
    if J1!=[]:
        idx = np.argmin(Eperc1)
        j1,Eperc1 = J1[idx],Eperc1[idx]
    if J2!=[]:
        idx = np.argmin(Eperc2)
        j2,Eperc2 = J2[idx],Eperc2[idx]
    
    # choosing between j1 and j2
    if J1!=[]: j_final,jf=j1,1
        
    if J1==[] and J2!=[]:
        for i in range(nist_max):
            jj1,jj2,EEperc1,EEperc2=FindJthAll(i,df_nist,df_ud,fac,corr_config)
            if jj1!=-1 and jj1==j2:
                JJ1,EEEperc1=jj1,EEperc1
        
        if EEEperc1<Eperc2:j_final=-1
        else: j_final,jf=j2,2
        if JJ1==0: j_final,jf=j2,2
                    
    return j_final,jf


def FindJth(i,df_nist,df_ud,nist_max,fac,which_j=False,corr_config=True): # Final J value surviving all five passes
    
    j1,j2,Eperc1,Eperc2 = FindJthAll(i,df_nist,df_ud,fac,corr_config)

    # j1=-1 if len(j1)==0 else j1[0]
    # j2=-1 if len(j2)==0 else j2[0]
    # Eperc1=-1 if len(Eperc1)==0 else Eperc1[0]
    # Eperc2=-1 if len(Eperc2)==0 else Eperc2[0]

    j_final,jf = ChooseJth(j1,j2,Eperc1,Eperc2,df_nist,df_ud,nist_max,fac,corr_config)

    if which_j==True:
        return j_final,jf
    else: return j_final


def finalsort(data):
    data_final = []
    config = data[:,0]
    term = data[:,1]
    energy = np.char.replace(data[:,11],"-","0.0")
    used = np.array([])
    for i in range(len(data)):
        if (i in used.astype(int)) == True: continue
        idx = np.where((config==config[i]) & (term==term[i]))[0] if config[i]!="-" else [i]
        used = np.append(used,idx)
        sort = energy[idx].argsort()
        data_final.append(data[idx][sort])
    data_final=np.concatenate(data_final)
    data_final = np.char.replace(data_final,"nan","")
    return data_final

#################################################################################################################


def IsBeingUsed(j,df_nist,df_ud,nist_max,fac,corr_config=False): # return the list of 
    jf=[]
    for i in range(nist_max):
        j_final,jw=FindJth(i,df_nist,df_ud,nist_max,fac,True,corr_config)
        if j_final == j: jf.append(jw)

    if len(jf)==0: jf=[0]
    return jf


def Correct_Config(df_ao,df_nist,nist_max,fac):
    # Automation needs corrections, need to done by comparing with NIST, need to combine with FindJthAll
    list_config,list_term,list_number,count=np.array([]),np.array([]),np.array([]),np.array([])

    # j=-1
    new_config=[]
    final_config = []
    for j in range(len(df_ao)): #len(df_ao)-3
        s = int(df_ao[j][5])
        mul = int(df_ao[j][2][0])

        jf = IsBeingUsed(j,df_nist,df_ao,nist_max,fac)[0]
        if jf==2:
            config = df_ao[j][1]

            d1,a1,d2,a2=Sep_Config(config)
            dummy_config = "1"+a1+"1"+a2 if a2 != '' else "1"+a1+"1"
            idx_dummy = np.where((list_config.astype(str)==dummy_config) & (list_term.astype(int)==mul))[0] # index of dummy config present in list_config    

            # Store term and config for comparison
            if len(idx_dummy)==0:
                list_config=np.append(list_config,dummy_config)
                list_term=np.append(list_term,mul)
                list_number=np.append(list_number,d2-1 if a2 != '' else d1-1)
                count = np.append(count,0)

            config = df_ao[j][10]
        
        else: config = df_ao[j][1]

        ## Manual Correction
        if config=="5s9p" and mul==1: config="4d5p"

        d1,a1,d2,a2=Sep_Config(config)

        dummy_config = "1"+a1+"1"+a2 if a2 != '' else "1"+a1+"1"
        idx_dummy = np.where((list_config.astype(str)==dummy_config) & (list_term.astype(int)==mul))[0] # index of dummy config present in list_config    

        # Store term and config for comparison
        if len(idx_dummy)==0:
            list_config=np.append(list_config,dummy_config)
            list_term=np.append(list_term,mul)
            list_number=np.append(list_number,d2 if a2 != '' else d1)
            count = np.append(count,0)


        # main correction part
        idx_dummy = np.where((list_config.astype(str)==dummy_config) & (list_term.astype(int)==mul))[0] # index of dummy config present in list_config
        if a2=='s':mul=1
        if count[idx_dummy[0]]!=mul: count[idx_dummy[0]]+=1
        if a2=='':
            d1_prev = list_number[idx_dummy[0]]
            if abs(d1_prev-d1)>1: 
                d1=d1-int(abs(d1_prev-d1)-1)

        else :
            d2_prev = list_number[idx_dummy[0]]
            if abs(d2_prev-d2)>1:
                d2=d2-int(abs(d2_prev-d2)-1)

        # re-initian for the next loop
        if count[idx_dummy[0]]==mul:
            list_number[idx_dummy[0]]=d2 if a2 != '' else d1
            count[idx_dummy[0]]=0

        new_config.append(str(d1)+a1+str(d2)+a2 if jf != 2 else df_ao[j][1])
    
    return new_config



def MainCode(path_nist, path_ud, nist_max,fac, parity, gs_parity):
    
    df_nist,df_ud,ref_E = Dataframe(path_nist,path_ud,nist_max,parity,fac,gs_parity)
    df_nist_org = np.copy(df_nist)

    primary_config_ud = np.array([df_ud[i][1] for i in range(len(df_ud))])
    
    
    
    data_csv = []
    ao_used = []
    roundd = 3

    # ref_E = df_nist_org[:,3][0] if parity=="odd" else 0
    ref_E = 0
    for i in range(nist_max):
        config_nist, term_nist, j_nist, level_nist, uncer_nist = Data_Nist(i,df_nist)
        data_nist = [config_nist, term_nist, j_nist, round(level_nist+ref_E,roundd),uncer_nist]

        nist_used = 0

        jth,jf = FindJth(i,df_nist,df_ud,nist_max,fac,which_j=True)
        if jf==2 and (1 in IsBeingUsed(jth,df_nist,df_ud,nist_max,fac,True))==True:jth=-1 # avoiding duplicates
        if jth!=-1:
            ao_used.append(jth)
            nist_used = 1
            OrgConfig1,config1_ao, config2_ao, term_ao, j_ao, level_ao,level_au, per1,per2,uncer_ud = Data_UD(jth,df_ud)
            final_config = config2_ao if jf==2 else config1_ao # final config
            Ediff = round(abs(level_nist-level_ao),1)
            Eperc = round(100*Ediff/level_nist,2) if level_nist!=0 else 0
            EpercStr = str(Eperc)+"%"
            if config2_ao=='': config2_ao='-'
            data_csv.append(data_nist+[final_config,OrgConfig1,config1_ao, config2_ao, term_ao, j_ao, round(level_ao+ref_E,roundd),uncer_ud,level_au,Ediff,EpercStr, per1,per2])


        ## remaining nist states
        if nist_used==0:
            data_csv.append([config_nist, term_nist, j_nist, round(level_nist+ref_E,roundd),uncer_nist,"-","-","-","-","-","-",round(level_nist+ref_E,roundd),"-","-","-","-","","-"])

    ## remaining ao states
    for j in range(len(df_ud)):
        OrgConfig1,config1_ao, config2_ao, term_ao, j_ao, level_ao,level_au, per1,per2,uncer_ud = Data_UD(j,df_ud)
        if config2_ao=='': config2_ao='-'
        if (j in ao_used)==False:
            data_csv.append(["-","-","-",round(level_ao+ref_E,roundd),"-","-",OrgConfig1,config1_ao, config2_ao, term_ao, j_ao, round(level_ao+ref_E,roundd),uncer_ud,level_au,"-","-", per1,per2])


    header = [" Config","  Term","  J","  Level (cm⁻¹)","Uncertainty (cm-1)","Final Config","   OrgConfig1","  CorConfig1","  Config2","  Term","  J","  Level (cm⁻¹)","  Uncertainty (cm⁻¹)","  Level (au)","    \u0394E","    \u0394E%","  I","  II"]
    data = np.array(data_csv)

    # sort data
    sort = data.T[11].astype(float).argsort()
    data = data[sort]
    
    for i in range(len(data)):
        if data[i][0]=="-": data[i][3] = "-"
        if data[i][0]=="": data[i][3] = ""

    for i in range(len(data)):# 6 and 10 neends to change if added any new column between them
        if data[i][6]=="-": data[i][11] = "-"
        if data[i][6]=="": data[i][11] = ""
    
    # Force energies in a.u. to have 8 decimal places
    for i in range(len(data)):
        if data[i][13].strip('.').isnumeric():
            data[i][13] = "{:.8f}".format(float(data[i][13]))

    data[data=='nan'] = '-'

    print(np.where(np.bincount(ao_used) > 1)[0]) # duplicates
    
    data_final = finalsort(data)
    
    return [header,data_final]




def ConvertToTXT(data_final,path):
    remove = [-1,-2]
    data_new = np.delete(data_final[1], remove, axis=1)
    header_new = np.delete(data_final[0], remove, axis=0)
    
    data_filtered = pd.DataFrame(data_new,columns=header_new)
    
    if os.path.exists(path)==True:os.remove(path)

    open(path, 'a').close()
    #export DataFrame to text file
    with open(path, 'a') as f:
        df_string = data_filtered.to_string(header=True, index=False)
        f.write(df_string)