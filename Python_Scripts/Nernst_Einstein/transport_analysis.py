#!/usr/bin/env python3
#importing libraries 
"""
Created on Sat Apr  8 04:29:51 2023

@author: Ashutosh Kumar Verma and Amey Thorat
"""
import pandas as pd
import numpy as np
import glob
import multiprocessing
import dask.dataframe as dd
from pathlib import Path
import time
import shutil

st=time.time()

#define number of cores
n_cores=multiprocessing.cpu_count()
print("Found {} cores".format(n_cores))

#charge multiplication list:
cation_list=["EMI","BMI","OMI"]
anion_list=["NC","BF","SCN","Br","Cl","MS","NO3","NSC","PF6","TCM","TFO"]
solvent_list=["ETGLY"]

#charge_multiplier
charge_multiplier=0.8
print("extracting center of mass from com_files")
def extract_COMs(f):
    x=[]
    y=[]
    z=[]
    t=[]
    with open("com_files/com_{}.pdb".format(f)) as file:
            matrix=[line.split() for line in file]
            #print("Extracting from com_files/com_{}.pdb".format(f))
    for m in matrix:
        if m[0]=="TITLE":       #TITLE line contains information on t= 0.00
            temp_t=m[3]
        elif m[0]=="CRYST1":    #CRYST1 line contains the box_length
            l=m[1]
        elif (m[0]=="ATOM"):    #All lines containing ATOM (or HETATOM) have coordinates
            res_name=m[3]
            if ('.' in m[5]):   #resID >62
                p=4
            else:
                p=5

            k=0   #no joint strings

            if len(m[p+1]) >=10:  # x is joint string
                k=1
            elif len(m[p+2]) >=10:  # y is joint string
                k=2

            if(k==0):  #no joint strings
                temp_x=m[p+1]
                temp_y=m[p+2]
                temp_z=m[p+3]
                
            elif(k!=0):   #if joint strings present -> slice
                a=m[p+k]

                #generate list containing positions of "-" sign in joint string
                lis=[i for i, letter in enumerate(m[p+k]) if letter == "-"]
                i=len(lis) # i gives number of -ve values in original joint string

                #complete the list to contain 3 elements
                if(i==1):                 # either y or z are -ve
                    lis.insert(0,0)       # insert 0 as first element and len(a) as last element of list
                    lis.append(len(a))
                elif(i==2):               # either (x,y) or (y,z) are -ve
                    if(lis[0]==0):lis.append(len(a))    #if lis[0]=0 then (x,y) are -ve, append l(a)
                    if(lis[0]!=0):lis.insert(0,0)       #if lis[0]!=0 then (y,z) are -ve, insert 0 at lis[0]

                #default slices into temp_x,y,z
                temp_x=a[:lis[1]]
                temp_y=a[lis[1]:lis[2]]
                temp_z=a[lis[2]:len(a)]

                #check conditions if x is joint string
                if(k==1):   
                    if(i==1): # -y joint with x then z = m[p+2]
                        temp_z=m[p+2]
                    elif(i==2)&(len(m)==p+5): # (x,y) joint, then z is m[p+2]
                        temp_z=m[p+2]
                    #if x,y,z are joint default slicing works irrespective of signs

                #check conditions if y is joint string
                elif(k==2):   # (y,z) joint then first slice is y, second slice is z, x is m[p+1]
                    temp_z=temp_y
                    temp_y=temp_x
                    temp_x=m[p+1]

            #append to actual x,y,z
            t.append(temp_t)
            x.append(temp_x)
            y.append(temp_y)
            z.append(temp_z)

    #charge multiplication factor +ve for cations, -ve for anions:
    multiplier = 1
    if (res_name in cation_list):
        multiplier=charge_multiplier
    elif (res_name in anion_list):
        multiplier=(-1)*charge_multiplier
    elif (res_name in solvent_list):
        multiplier=(0)*charge_multiplier

    data=(np.vstack((t,x,y,z))).T
    cols=["t","xu","yu","zu"]
    df=pd.DataFrame(data=data,columns=cols)
    df.iloc[:,1:4]=(df.iloc[:,1:4].values).astype(float)*multiplier
    df.to_csv("data/species/species_{}.csv".format(f), index=None)

def pool_handler():
    global n_cores
    p = multiprocessing.Pool(n_cores)
    p.map(extract_COMs,list_n_steps)
    return

all_filenames = [i for i in glob.glob("com_files/*.pdb")]
n_count=len(all_filenames)
list_n_steps=[i for i in range(1,n_count+1)]

if __name__ == '__main__':
    pool_handler()
    
#define number of cores
n_cores=multiprocessing.cpu_count()
print("Found {} cores".format(n_cores))


dt = 20
 

def msd_calc_ions(ions):
    df_ions = dd.read_csv(f"data/species/species_{ions}.csv").compute()
    n_snapshots = len(df_ions)
    x, y, z = df_ions[["xu", "yu", "zu"]].values.T
    msd_particle = np.empty((1, n_snapshots-1))

    for delta in range(1, n_snapshots):
        r = np.sqrt((x[delta:] - x[:-delta]) ** 2 +
                    (y[delta:] - y[:-delta]) ** 2 +
                    (z[delta:] - z[:-delta]) ** 2)
        msd_particle[0, delta-1] = np.mean(r ** 2)

    msd_particle = pd.DataFrame(msd_particle,
                                columns=[f"MSD_{delta}" for delta in range(1, n_snapshots)])
    msd_particle.to_pickle(f"data/msd/msd_{ions}.pkl", compression='infer')

def pool_handler_ions():
    global n_cores
    p = multiprocessing.Pool(n_cores)
    p.map(msd_calc_ions, list_ions_steps)

if __name__ == '__main__':
    n_cores = multiprocessing.cpu_count()
    print("Found {} cores".format(n_cores))
    data_path = Path("data/species/")
    data_filenames = list(data_path.glob("species_*.csv"))
    n_molecules = len(data_filenames)
    list_ions_steps = list(range(1,n_molecules+1))

    pool_handler_ions()

    msd_path = Path("data/msd/")
    msd_filenames = list(msd_path.glob("msd_*.pkl"))
    combined_csv = pd.concat([pd.read_pickle(f) for f in msd_filenames])
    combined_csv.to_pickle("data/msd.pkl", compression='infer')

    msd_ions = pd.read_pickle("data/msd.pkl")
    print(msd_ions)
    ions_avg = msd_ions.iloc[:, 1:].mean(axis=0).values
    time_t = np.arange(1, len(ions_avg) + 1) * dt
    summary = pd.DataFrame(np.column_stack([time_t, ions_avg]), columns=["delta", "product"])
    summary['product'] = summary['product']*int(n_molecules)
    summary.to_csv("Nernst_Einstein.csv", index=None)


et=time.time()
print(et-st)

shutil.rmtree("data")
shutil.rmtree("com_files")
shutil.rmtree("index_files")
