import pandas as pd
import numpy as np
from get_atomic_term import scrape_term
import os

pd.set_option("display.max_rows", None, "display.max_columns", None)

def Convert_Type(s): # detect and correct the 'type' of object to 'float', 'integer', 'string' while reading data
    # s = s.replace(" ", "")
    try:
        f = float(s)
        i = int(f)
        return i if i == f else f
    
    except ValueError:
        try: return eval(s)
        except: return s

def Sep_Config(config):
    #  divided configuration in parts: "5s2 10g"
    #  n=0: 5s2 --> l=0:  5 and l=1: s2
    #  n=1: 10g --> l=0: 10 and l=1: g
    #  and storing in list as [5,s2,10,g]
    #  can add more parts for when required

    conf = np.full(30,"",dtype=object)
    p=0 # first part of configuration
    l=0
    for i in config:
        # if i==" ":# next part of configuration after space
        if i==".":# next part of configuration after space
            p+=2
            l=0

        if i.isalpha()==True: l=1
        conf[p]+=i if l==0 and i.isdigit()==True else ""
        conf[p+1]+=i if l==1 else ""
            
    return [i for i in conf if i!=""]



def Combine_Config(conf):
    # combine list of seprated configuration
    config=""
    for i in range(len(conf)):
        config+="."if i%2==0 and i>0 else ""
        # config+=" "if i%2==0 and i>0 else ""
        config+=conf[i]
    return config



def Inverse_Config(config):
    # 5s 6p <--> 6p 5s as somtime configuration text in ud and nist are inverse of each other
    # so correspondence will be checked with both 5s 6p and 6p 5s
    conf=Sep_Config(config)
    if len(conf)==4: return conf[2]+conf[3]+"."+conf[0]+conf[1]
    # if len(conf)==4: return conf[2]+conf[3]+" "+conf[0]+conf[1]
    else:  return config


def Correct_Config(config,val,pos):
    # correction to configuraion by value "val", ex: '5s4d' --> '5s5d' if val = 1 and pos=2
    conf=Sep_Config(config)
    conf[2*(pos-1)] = str(int(conf[2*(pos-1)])+val)
    return Combine_Config(conf)


def Term_Correct(term,corr=1):
    # correction to term symbol by value "corr", ex : 2D --> 3D if corr = 1
    term_symbol = list(term)
    for i in range(len(term_symbol)):
        if term_symbol[i].isdigit()==True: term_symbol[i]=str(int(term_symbol[i])+corr)
    term_symbol="".join(term_symbol)
    
    return term_symbol


def Extract_Multiplicity(term): # extract multiplicity from term
    mul=""
    for i in term:
        if i.isdigit():mul+=i
    return int(mul)


def Dummy_Config(config):
    conf=Sep_Config(config)
    for i in range(len(conf)):
        if len(conf)>2 and i==0:continue
        if i%2==0: conf[i]="1"
    return Combine_Config(conf)


def Mark_UD(df_ud):
    listt = np.array([])
    mark = list('abcdefghijklmnopqrstuvwxyz')
    for i in range(len(df_ud)):
        if (i in listt)==True:continue
        idx = np.where((df_ud[:,1]==df_ud[i][1]) & (df_ud[:,2]==df_ud[i][2]))[0]
        for j in range(1,len(idx)):
            df_ud[idx[j]][2] = mark[j-1]+df_ud[idx[j]][2]
    
    listt = np.concatenate([listt,idx])
    return df_ud



def Mark_NIST(df_nist):
    listt = np.array([])
    mark = list('abcdefghijklmnopqrstuvwxyz')
    for i in range(len(df_nist)):
        if (i in listt)==True:continue
        idx = np.where((df_nist[:,0]==df_nist[i][0]) & (df_nist[:,1]==df_nist[i][1]) & (df_nist[:,2]==df_nist[i][2]))[0]
        for j in range(1,len(idx)):
            df_nist[idx[j]][1] = mark[j-1]+df_nist[idx[j]][1]

        listt = np.concatenate([listt,idx])
    return df_nist


def Unmark_Term(Term): # remove marker from terms
    term,s=[],0
    for i in Term:
        if i.isdigit():s=1
        if s==1:term+=i
    return "".join(term)


def Data_UD(ith,df_ao): # reading filtered csv data
    
    config1 = df_ao[ith][1] # original first possible configuration from data
    config2 = df_ao[ith][10] # second possible configuration
    Term = df_ao[ith][2][:-1] # Term
    # J = int(round(df_ao[ith][7],0)) # J value
    J = df_ao[ith][7] # J value
    Level = round(df_ao[ith][4],3) # Level
    Levelau = df_ao[ith][3] # Level in atomic units
    Per1 = df_ao[ith][9] # Percentage Contribution of first configuration
    Per2 = df_ao[ith][11] # Percentage Contribution of second configuration
    uncer_ud = df_ao[ith][12]

    return config1,config2,Term,J,Level,Levelau,Per1,Per2,uncer_ud


def Data_Nist(ith,df_nist):
    # Reading filtered csv data
    config_nist = Convert_Type(df_nist[ith][0])
    term_nist = Convert_Type(df_nist[ith][1])
    j_nist = Convert_Type(df_nist[ith][2])
    level_nist = round(Convert_Type(df_nist[ith][3]),3)
    uncer_nist = round(Convert_Type(df_nist[ith][4]),3)
    
    return config_nist,term_nist,j_nist,level_nist,uncer_nist


def Dataframe(path_nist,path_ud,nist_max=0):
    # read csv
    df_nist = pd.read_csv(path_nist)
    df_ud = pd.read_csv(path_ud).fillna('') # fill blank for nan
    df_ud = df_ud.sort_values(by=df_ud.columns[4], ascending=True)

    # Storing as numpy arrays
    df_ud = df_ud.values
    df_nist = df_nist.values
    if nist_max==0:
        maxx = df_ud[:,1][-1]
        nist_max =  max(np.where(df_nist[:,0]==maxx)[0])-8
    
    print("nist_max : ",nist_max)

    df_nist = df_nist[:nist_max]

    df_nist=Mark_NIST(df_nist)
    df_ud=Mark_UD(df_ud)

    ref_E = df_nist[:,3][0]
    df_nist[:,3] = df_nist[:,3]- ref_E
    print("Generating Dataframes...!")
    return df_nist,df_ud,ref_E



def FindJthAll(i,df_nist,df_ud,corr_config=[]): # find ao config corresponding to ith nist config
    j1,j2,Eperc1,Eperc2=[],[],[],[]
    config_nist, term_nist, j_nist, level_nist, uncer_nist = Data_Nist(i,df_nist)
    corr=1 if len(corr_config)==len(df_ud) else 0 # corr=True --> Corrected Configurations are proveded and to use those
    for j in range(len(df_ud)):
        config1_ao, config2_ao, term_ao, j_ao, level_ao,Levelau, per1,per2,uncer_ud = Data_UD(j,df_ud)
        if corr==1: config1_ao=corr_config[j]
            
        Ediff = round(abs(level_nist-level_ao),1)
        Eperc = round(100*Ediff/level_nist,2) if level_nist!=0 else 0
        
        # First Pass : if J_Nist = J_AO
        if j_nist==j_ao:
            inv_config = Inverse_Config(config_nist)
            
            # Second Pass : if Term_NIST = Term_AO
            if term_nist==term_ao:
                # Third Pass :
                if config_nist==config1_ao or inv_config==config1_ao:j1.append(j),Eperc1.append(Eperc)
                if corr==0:
                    if config_nist==config2_ao or inv_config==config2_ao:j2.append(j),Eperc2.append(Eperc)
            
            # Fourth Pass : If Term_NIST != Term_AO --> Check using term correction by +-1
            if term_nist!=term_ao: # correction to terms in ao results
                if config_nist==config1_ao or inv_config==config1_ao:
                    if term_nist==Term_Correct(term_ao,1):j1.append(j),Eperc1.append(Eperc)
                    if term_nist==Term_Correct(term_ao,-1):j1.append(j),Eperc1.append(Eperc)

                if corr==0:
                    if config_nist==config2_ao or inv_config==config2_ao:
                        if term_nist==Term_Correct(term_ao,1):j2.append(j),Eperc2.append(Eperc)
                        if term_nist==Term_Correct(term_ao,-1):j2.append(j),Eperc2.append(Eperc)
    
    return j1,j2,Eperc1,Eperc2 # j1,j2 is the index of ud configurations 1 and 2 corresponding to NIST


# Fifth Pass : to choose between j1 and j2 if both present to find the best correspondance
def ChooseJth(ith,J1,J2,Eperc1,Eperc2,df_nist,df_ud,corr_config=[]):
    JJ1=0
    EEEperc1=0
    j_final,jf,dE=-1,0,0 # jf is "which config selected" : jf=0 --> no config, jf=1 --> 1st Config, jf=2 --> 2nd Config
    
    # Choosing minimum percentage defference index if more than one values present in jj1 and jj2
    if J1!=[]:
        if len(J1)>1:
            for i in range(len(J1)):
                config1,config2,Term,J,Level,Levelau,Per1,Per2,uncer_ud=Data_UD(J1[i],df_ud)
                config_nist,term_nist,j_nist,level_nist,uncer_nist=Data_Nist(ith,df_nist)
                if Term!=term_nist:
                    J1.remove(J1[i])
                    Eperc1.remove(Eperc1[i])
                    break

        idx = np.argmin(J1)
        j1,Eperc1 = J1[idx],Eperc1[idx]

    if J2!=[]:
        if len(J2)>1:
            for i in range(len(J2)):
                config1,config2,Term,J,Level,Levelau,Per1,Per2,uncer_ud=Data_UD(J2[i],df_ud)
                config_nist,term_nist,j_nist,level_nist,uncer_nist=Data_Nist(ith,df_nist)
                if Term!=term_nist:
                    J2.remove(J2[i])
                    Eperc2.remove(Eperc2[i])
                    break

        idx = np.argmin(J2)
        j2,Eperc2 = J2[idx],Eperc2[idx]

    
    # if J1!=[] and J2==[]: j_final,jf,dE=j1,1,Eperc1
    # elif J1==[] and J2!=[]: j_final,jf,dE=j2,2,Eperc2
    # elif J1!=[] and J2!=[]:
    #     jf=1 if j1<=j2 else 2
    #     if jf==1: j_final,jf,dE=j1,1,Eperc1
    #     if jf==2: j_final,jf,dE=j2,2,Eperc2

    # else: j_final=-1

    # choosing between j1 and j2
    if J1!=[]: j_final,jf,dE=j1,1,Eperc1
        
    if J1==[] and J2!=[]:
        for i in range(len(df_nist)):
            jj1,jj2,EEperc1,EEperc2=FindJthAll(i,df_nist,df_ud,corr_config)
            if jj1!=-1 and jj1==j2:
                JJ1,EEEperc1=jj1,EEperc1
        
        if EEEperc1<Eperc2:j_final=-1
        else: j_final,jf,dE=j2,2,Eperc2
        if JJ1==0: j_final,jf,dE=j2,2,Eperc2
        
    return j_final,jf,dE



def FindJth(i,df_nist,df_ud,corr_config=[],which_j=False): # Final J value surviving all five passes
    
    j1,j2,Eperc1,Eperc2 = FindJthAll(i,df_nist,df_ud,corr_config)

    j_final,jf,dE = ChooseJth(i,j1,j2,Eperc1,Eperc2,df_nist,df_ud,corr_config)

    if which_j==True:
        return j_final,jf,dE
    else: return j_final



def IsBeingUsed(j,df_nist,df_ud,corr_config=[],OutAll=False): # return the list of
    jf=[]
    nf=[]
    de=[]
    for i in range(len(df_nist)):
        j_final,jw,dE=FindJth(i,df_nist,df_ud,corr_config,which_j=True)
        if j_final == j:
            jf.append(jw)
            nf.append(i)
            de.append(dE)
    
    
    if OutAll==False:
        if len(jf)==0: return 0,0,0
        else:
            jth = np.where(np.array(de)==min(de))[0][0]
            return nf[jth],jf[jth],de[jth]
    
    if OutAll==True:
        if len(jf)==0:jf=[0]
        return nf,jf,de


def StartingConfigs(df_nist):
    list_config,list_term,list_number,count=np.array([]),np.array([]),np.array([]),np.array([])
    for j in range(len(df_nist)): #len(df_ao)-3
        # mul = int(df_nist[j][1][0])
        mul = Extract_Multiplicity(df_nist[j][1])
        config = df_nist[j][0]
        conf=Sep_Config(config)
        dummy_config = Dummy_Config(config)
        idx_dummy = np.where((list_config.astype(str)==dummy_config) & (list_term.astype(int)==mul))[0]
        
        n = int(conf[-2])-1
        # Store term and config for comparison
        if len(idx_dummy)==0:
            list_config=np.append(list_config,dummy_config)
            list_term=np.append(list_term,mul)
            list_number=np.append(list_number,n)
            count = np.append(count,0)
        
    return list_config,list_term,list_number,count


def Correct_Config(jth,config,mul,new_config,list_config,list_term,list_number,count,Final=False):

    conf=Sep_Config(config)
    dummy_config = Dummy_Config(config)
    idx_dummy = np.where((list_config.astype(str)==dummy_config) & (list_term.astype(int)==mul))[0]

    n = int(conf[-2])
    # Store term and config for comparison
    if len(idx_dummy)==0:
        list_config=np.append(list_config,dummy_config)
        list_term=np.append(list_term,mul)
        list_number=np.append(list_number,n)
        count = np.append(count,0)

    # main correction part
    idx_dummy = np.where((list_config.astype(str)==dummy_config) & (list_term.astype(int)==mul))[0] # index of dummy config present in list_config
    if conf[-1]=='s':mul=1
    n_prev = list_number[idx_dummy[0]]

    if abs(n_prev-n)>1:
        n=n-int(abs(n_prev-n)-1)

    conf[-2]=str(n)
    new_config[jth]=Combine_Config(conf)
    
    if Final==True:
        if count[idx_dummy[0]]!=mul: count[idx_dummy[0]]+=1
        # re-initian for the next loop
        if count[idx_dummy[0]]==mul:
            list_number[idx_dummy[0]]=n
            count[idx_dummy[0]]=0
            
        return new_config,list_config,list_term,list_number,count

    else: return new_config


def Corrected_Config(df_ud,df_nist,ManCorr=False): # ManCorr : Manual Correction
    print("Correcting ud configurations")
    new_config = np.full(len(df_ud),"",dtype=object)
    list_config,list_term,list_number,count = StartingConfigs(df_nist)

    for j in range(len(df_ud)):
        nf,jf,de = IsBeingUsed(j,df_nist,df_ud)
        config1_ao, config2_ao, term_ao, j_ao, level_ao,Levelau, per1,per2,uncer_ud = Data_UD(j,df_ud)
        # mul = int(term_ao[0]) if jf==0 else int(df_nist[nf][1][0])
        mul = Extract_Multiplicity(term_ao) if jf==0 else Extract_Multiplicity(df_nist[nf][1])

        if jf==1: config = config1_ao
        if jf==0 or jf==2:
            config = config1_ao
            new_config_check = Correct_Config(j,config,mul,new_config,list_config,list_term,list_number,count,Final=False)
            nf1,jf1,de1 = IsBeingUsed(j,df_nist,df_ud,corr_config=new_config_check)

            new_config_check[j] = config2_ao
            nf2,jf2,de2 = IsBeingUsed(j,df_nist,df_ud,corr_config=new_config_check)

            if de1!=0 and de2==0: config=config1_ao
            elif de1==0 and de2!=0: config=config2_ao
            elif de1!=0 and de2!=0: config=config1_ao if de1<=de2 else config2_ao
            else: config=config1_ao

        # Manual Corrections
        if ManCorr==True:
            if jf==1 and config=="5s.9p" and mul==1:config="4d.5p"

        new_config,list_config,list_term,list_number,count = Correct_Config(j,config,mul,new_config,list_config,list_term,list_number,count,Final=True)
    print("Correcting ud configurations : Complete..!")
    return new_config



def FinalSortE(data):
    data_final = []
    energyu = np.char.replace(data[:,10],"-","0.0").astype(float)

    # "-" in data is filled with corresponding UDel Data Values
    confign,termn,energyn = np.copy(data[:,5]),np.copy(data[:,1]),np.copy(data[:,3])
    idxx = np.where(confign=="-")
    confign[idxx],termn[idxx] = data[:,5][idxx],data[:,8][idxx]

    # seprate ao which are not corresponding to UD and order them according to UD odering
    used = np.array([])
    skip = False
    for i in range(len(data)):
        if (i in used.astype(int)) == True: continue
        idx = np.where((confign==confign[i]) & (termn==termn[i]))[0] if energyn[i]=="-" else [i]
        
        skip = False
        if len(idx)==1:
            idxcheck = np.where((confign==confign[i]) & (termn==termn[i]))[0]
            for j in idxcheck:
                if energyn[j]=="-" and len(idxcheck)>1:skip=True
            if skip==True: continue
                
        used = np.append(used,idx)
        sort = energyu[idx].argsort()
        data_final.append(data[idx][sort])
        
    sort = np.array([ESort(i,data_final) for i in range(len(data_final))]).argsort()
    # print([ESort(i,data_final) for i in range(len(data_final))])
    data_final = np.array(data_final,dtype=object)
    data_final=np.concatenate(data_final[sort])
    data_final = np.char.replace(data_final.astype(str),"nan","")
    print("Sorted by energies...!")
    return data_final


def FinalSortT(data):
    data_final = []

    energyu = np.char.replace(data[:,10],"-","0.0").astype(float)

    # "-" in data is filled with corresponding UDel Data Values
    config,term,energyn = np.copy(data[:,5]),np.copy(data[:,1]),np.copy(data[:,3])
    idxx = np.where(config=="-")
    config[idxx],term[idxx],energyn[idxx] = data[:,5][idxx],data[:,8][idxx],data[:,10][idxx]

    

    # print(energy)
    used = np.array([])
    for i in range(len(data)):
        if (i in used.astype(int)) == True: continue
        idx = np.where((config==config[i]) & (term==term[i]))[0] if config[i]!="-" else [i]
        used = np.append(used,idx)
        sort = energyu[idx].argsort()
        data_final.append(data[idx][sort])
        
    sort = np.array([ESort(i,data_final) for i in range(len(data_final))]).argsort() # According to NIST Levels
    # sort = np.array([float(data_final[i][:,10][-1]) for i in range(len(data_final))]).argsort() # according to UD Levels
    data_final = np.array(data_final,dtype=object)
    data_final=np.concatenate(data_final[sort])
    # data_final=np.concatenate(data_final) # according to UD Levels same as sorting as data is already sorted in main code
    data_final = np.char.replace(data_final,"nan","")
    print("Sorted by terms...!")
    return data_final


def ESort(i,data_final):
    idxx = np.where(data_final[i][:,3]!="-")[0]
    if len(idxx)==0: return float(data_final[i][:,10][0])
    else : return float(data_final[i][:,3][idxx][-1])



# j,df_nist,df_ud,corr_config=[],OutAll=False
def MainCode(path_nist,path_ud,nist_max,Ordering="E"):

    df_nist,df_ud,ref_E = Dataframe(path_nist,path_ud,nist_max)
    ManualCorrection=True if ("Sr_I" in path_nist.split("/")[-1])==True else False
    new_config = Corrected_Config(df_ud,df_nist,ManCorr=ManualCorrection)

    data_csv = []
    ao_used = []
    roundd = 3

    print("Finding Correspondance")
    ref_E = 0
    for i in range(len(df_nist)):
        config_nist, term_nist, j_nist, level_nist, uncer_nist = Data_Nist(i,df_nist)
        term_nist=Unmark_Term(term_nist)
        data_nist = [config_nist, term_nist, j_nist, round(level_nist+ref_E,roundd),uncer_nist]

        nist_used = 0
        jth,jf,de = FindJth(i,df_nist,df_ud,corr_config=new_config,which_j=True)
        if (jth in ao_used)==True:jth=-1 # avoiding duplicates
        if jth!=-1:
            ao_used.append(jth)
            nist_used = 1
            config1_ao, config2_ao, term_ao, j_ao, level_ao,level_au, per1,per2,uncer_ud = Data_UD(jth,df_ud)
            final_config = config2_ao if jf==2 else config1_ao # final config
            Ediff = round(abs(level_nist-level_ao),1)
            Eperc = round(100*Ediff/level_nist,2) if level_nist!=0 else 0
            EpercStr = str(Eperc)+"%"
            if config2_ao=='': config2_ao='-'
            term_ao=Unmark_Term(term_ao)
            data_csv.append(data_nist+[new_config[jth],config1_ao,config2_ao, term_ao, j_ao, round(level_ao+ref_E,roundd),uncer_ud,level_au,Ediff,EpercStr, per1,per2])


        ## remaining nist states
        if nist_used==0:
            data_csv.append([config_nist, term_nist, j_nist, round(level_nist+ref_E,roundd),uncer_nist,"-","-","-","-","-",round(level_nist+ref_E,roundd),"-","-","-","-","","-"])

    ## remaining ao states
    for j in range(len(df_ud)):
        config1_ao, config2_ao, term_ao, j_ao, level_ao,level_au, per1,per2,uncer_ud = Data_UD(j,df_ud)
        term_ao=Unmark_Term(term_ao)
        if config2_ao=='': config2_ao='-'
        if (j in ao_used)==False:
            data_csv.append(["-","-","-",round(level_ao+ref_E,roundd),"-",new_config[j],config1_ao, config2_ao, term_ao, j_ao, round(level_ao+ref_E,roundd),uncer_ud,level_au,"-","-", per1,per2])


    header = [" Config","  Term","  J","  Level (cm⁻¹)","Uncertainty (cm-1)","Final Config","  Config1","  Config2","  Term","  J","  Level (cm⁻¹)","  Uncertainty (cm⁻¹)","  Level (au)","    \u0394E","    \u0394E%","  I","  II"]
    data = np.array(data_csv)

    # sort data
    sort = data.T[3].astype(float).argsort()
    data = data[sort]

    for i in range(len(data)):
        if data[i][0]=="-": data[i][3] = "-"
        if data[i][0]=="": data[i][3] = ""

    for i in range(len(data)):# 6 and 10 neends to change if added any new column between them
        if data[i][6]=="-": data[i][10] = "-"
        if data[i][6]=="": data[i][10] = ""

    # Force energies in a.u. to have 8 decimal places
    for i in range(len(data)):
        if data[i][12].replace('.','').isnumeric():
            data[i][12] = "{:.8f}".format(float(data[i][12]))
    
    # Fill blanks with '-'
    data[data=='nan'] = '-'
    
    print("Finding Correspondance : Complete...!")
    print("Duplicates : ",np.where(np.bincount(ao_used) > 1)[0]) # duplicates

    if Ordering=="E": data_final = FinalSortE(data)
    if Ordering=="T": data_final = FinalSortT(data)
    
    return [header,data_final]


def Missing_Levels(data):
    print("Finding missing levels")
    header,data_copy=data
    data_final = np.copy(data_copy)

    config_nist,term_nist,J_nist,config_ud,term_ud,J_ud = data_final[:,0],data_final[:,1],data_final[:,2],data_final[:,5],data_final[:,8],data_final[:,9]

    idx = np.where(config_nist=="-")[0]
    config_nist[idx],term_nist[idx],J_nist[idx]=config_ud[idx],term_ud[idx],J_ud[idx]
    terms_all = np.array([term_nist[i]+J_nist[i] for i in range(len(term_nist))])

    lst = np.array([])
    config_exp,term_exp,J_exp=[],[],[]
    for i in range(len(config_nist)):
        config = config_nist[i]
        if (i in lst)==True or ("0" in config)==True: continue
        terms_expected = scrape_term(config)
        for j in range(len(terms_expected)):
            term_str = terms_expected[j][:2]+str(Convert_Type(terms_expected[j][2:]))
            idx = np.where((config_nist==config) & (terms_all==term_str))[0]
            lst = np.concatenate([lst,idx])
            if len(idx)==0:
                config_exp.append(config_nist[i])
                term_exp.append(terms_expected[j][:2])
                J_exp.append(terms_expected[j][2:])
                
    data_csv=[]
    data_csv.append(["","","","","","","","","","","","","","","","",""])
    for j in range(len(config_exp)):
        data_csv.append(["~"+config_exp[j],term_exp[j],J_exp[j],"","","","","","","","","","","","","",""])

    print("Finding missing levels : complete...!")
    return [header,np.concatenate([data_copy,data_csv])]


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

    print(f"Produced Output File at {path}")
    print("Done...!")