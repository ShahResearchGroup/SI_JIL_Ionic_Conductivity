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

"""

In this section, we will sum coordinates for all the species

"""
print("calculating sum of coordinates")
    
all_filenames=[j for j in glob.glob("data/species/*.csv")]
species_count=int(len(all_filenames)+1)

x_sum = 0
y_sum = 0
z_sum = 0
for f in range(1,species_count,1):
    df_species = pd.read_csv("data/species/species_{}.csv".format(f))
    t_count=df_species["t"].values
    x_sum=x_sum+df_species["xu"].values
    y_sum=y_sum+df_species["yu"].values
    z_sum=z_sum+df_species["zu"].values

data=np.vstack((t_count,x_sum,y_sum,z_sum)).T
df_sum = pd.DataFrame(data=data,columns=['t', 'xu', 'yu','zu'])
del data
df_sum.to_csv("data/species.csv",index=False)


"""

In this section, we will calculate displacement and product

"""
print("calculating displacement")


def subtract(delta):
    df_species = pd.read_csv("data/species.csv")
    x = df_species["xu"].values
    y = df_species["yu"].values
    z = df_species["zu"].values
    delta_count = 0
    x_temp,y_temp,z_temp=0,0,0
    for i in range(delta,n_steps,1):
        delta_count+=1
        x_temp=x_temp+(x[i] - x[i - delta])**2
        y_temp=y_temp+(y[i] - y[i - delta])**2
        z_temp=z_temp+(z[i] - z[i - delta])**2
    x_temp=x_temp/delta_count
    y_temp=y_temp/delta_count
    z_temp=z_temp/delta_count
    sum_temp=x_temp+y_temp+z_temp

    return (delta,sum_temp)

print("calculating product")
def pool_handler_delta():
    global n_cores
    delt=[]
    disp=[]
    p = multiprocessing.Pool(n_cores)
    t_steps=p.map(subtract, list_n_steps)
    for t in t_steps:
        delt.append(t[0])
        disp.append(t[1])
        
    data = np.vstack((delt,disp)).T
    cols = ["delta", "product"]
        
    df=pd.DataFrame(data=data, columns=cols)
    df['delta'] = (df['delta'])*20
    df.to_csv("Einstein.csv",index=False)
    return

df_species = pd.read_csv("data/species.csv")
print(df_species)
n_steps = len(df_species)
list_n_steps = range(1, n_steps-1, 1)

if __name__ == '__main__':
  pool_handler_delta()
  

et=time.time()
print(et-st)

shutil.rmtree("data")
shutil.rmtree("com_files")
shutil.rmtree("index_files")
