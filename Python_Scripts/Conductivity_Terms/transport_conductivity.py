#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Apr 19 23:09:03 2023

@author: Ashutosh Verma and Amey Thorat

Dr. Jindal Shah Research Lab

"""
# In this section, We import all the relevant Libraries
import pandas as pd
import numpy as np
import glob
import multiprocessing
import shutil


#define number of cores
n_cores=multiprocessing.cpu_count()
print("Found {} cores".format(n_cores))

#charge multiplication list: (Please update list based on your .itp files)
cation_list=["EMI","BMI","OMI"]
anion_list=["NC","BF","SCN","Br","Cl","MS","NO3","NSC","PF6","TCM","TFO"]
solvent_list=["ETGLY"]

#charge_multiplier
charge_multiplier=0.8  # Please adjust this values according to force-field parameters
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


#calculating displacement by considering time origin shifting
print("calculating displacements")
def disp_delta(f):
    df_species = pd.read_csv("data/species/species_{}.csv".format(f))
    x=df_species["xu"].values
    y=df_species["yu"].values
    z=df_species["zu"].values
    n_deltas = 500
    n_steps = int(len(df_species))

    disp_x = []
    disp_y = []
    disp_z = []
    delta_count = []
    cols = []
    for delta in range(0,n_deltas):
        for i in range(delta,n_steps):
            delta_count.append(delta)
            disp_x.append(x[i]-x[i-delta])
            disp_y.append(y[i]-y[i-delta])
            disp_z.append(z[i]-z[i-delta])

    data=np.vstack((delta_count,disp_x,disp_y,disp_z)).T
    cols=["delta","xu","yu","zu"]
    df=pd.DataFrame(data=data,columns=cols)
    del data
    df.to_pickle("data/delta/species_disp_{}.pkl".format(f),compression='infer')
    del df

def pool_handler_delta():
  global n_cores
  p = multiprocessing.Pool(n_cores)
  p.map(disp_delta,list_species_count)
  return


species_filenames=[i for i in glob.glob("data/species/*.csv")]
n_molecules=len(species_filenames)

list_species_count=[i for i in range(1,n_molecules+1)]

if __name__ == '__main__':
  pool_handler_delta()




#combine all the delta's file into a single files column wise
print("combining all displacement files into a single file")
def combine(f):
    df_species = pd.read_pickle("data/delta/species_disp_{}.pkl".format(f))
    x=np.array(df_species["xu"].values)
    y=np.array(df_species["yu"].values)
    z=np.array(df_species["zu"].values)
    return x,y,z


ions_filenames=[i for i in glob.glob("data/species/*.csv")]
n_molecules=int(len(ions_filenames)+1)
n_cations=int(len(ions_filenames)*0.5+1)

df = pd.read_pickle("data/delta/species_disp_1.pkl")
delta=(df["delta"].values)


x=[]
y=[]
z=[]
cols=[]
for f in range(1,n_molecules):
  print("Reading: {}".format(f))
  x.append(combine(f)[0])
  y.append(combine(f)[1])
  z.append(combine(f)[2])

x.insert(0,delta)
y.insert(0,delta)
z.insert(0,delta)


#In this section, we are calculating Einstein terms
print("calculating Einstein terms")
x_product=np.zeros((len(df)))
y_product=np.zeros((len(df)))
z_product=np.zeros((len(df)))
for p1 in range(1,len(x)):
    for p2 in range(1,len(x)):
      x_product+=x[p1]*x[p2]
      y_product+=y[p1]*y[p2]
      z_product+=z[p1]*z[p2]

final=x_product+y_product+z_product
cols=["delta","product"]

data=np.vstack((delta,final)).T
pf=pd.DataFrame(data=data,columns=cols)
pf['delta'] = (pf['delta'])*20
del data
sf=pf.groupby(["delta"]).mean()
sf.to_csv("summary/einstein_terms.csv")
del sf

#In this section, we are calculating cross terms
print("calculating cross terms")
x_product=np.zeros((len(df)))
y_product=np.zeros((len(df)))
z_product=np.zeros((len(df)))
for p1 in range(1,len(x)):
    for p2 in range(p1+1,len(x)):
      x_product+=x[p1]*x[p2]
      y_product+=y[p1]*y[p2]
      z_product+=z[p1]*z[p2]

final=x_product+y_product+z_product
cols=["delta","product"]

data=np.vstack((delta,final)).T
pf=pd.DataFrame(data=data,columns=cols)
pf['delta'] = (pf['delta'])*20
del data
sf=pf.groupby(["delta"]).mean()
sf.to_csv("summary/cross_terms.csv")
del sf

#In this section, we are calculating self terms of cations
print("calculating cation self terms")
x_product=np.zeros((len(df)))
y_product=np.zeros((len(df)))
z_product=np.zeros((len(df)))
for p1 in range(1,n_cations):
    x_product+=x[p1]**2
    y_product+=y[p1]**2
    z_product+=z[p1]**2

final=x_product+y_product+z_product
cols=["delta","product"]

data=np.vstack((delta,final)).T
pf=pd.DataFrame(data=data,columns=cols)
pf['delta'] = (pf['delta'])*20
del data
sf=pf.groupby(["delta"]).mean()
sf.to_csv("summary/self_cation.csv")
del sf

#In this section, we are calculating self terms of anions
print("calculating anion self terms")
x_product=np.zeros((len(df)))
y_product=np.zeros((len(df)))
z_product=np.zeros((len(df)))
for p1 in range(n_cations,n_molecules):
    x_product+=x[p1]**2
    y_product+=y[p1]**2
    z_product+=z[p1]**2

final=x_product+y_product+z_product
cols=["delta","product"]

data=np.vstack((delta,final)).T
pf=pd.DataFrame(data=data,columns=cols)
pf['delta'] = (pf['delta'])*20
del data
sf=pf.groupby(["delta"]).mean()
sf.to_csv("summary/self_anion.csv")
del sf

#In this section, we are calculating cross terms of cation-cation
print("calculating cation-cation cross terms")
x_product=np.zeros((len(df)))
y_product=np.zeros((len(df)))
z_product=np.zeros((len(df)))
for p1 in range(1,n_cations):
    for p2 in range(p1+1,n_cations):
      x_product+=x[p1]*x[p2]
      y_product+=y[p1]*y[p2]
      z_product+=z[p1]*z[p2]

final=x_product+y_product+z_product
cols=["delta","product"]

data=np.vstack((delta,final)).T
pf=pd.DataFrame(data=data,columns=cols)
pf['delta'] = (pf['delta'])*20
del data
sf=pf.groupby(["delta"]).mean()
sf.to_csv("summary/cross_cation_cation.csv")
del sf

#In this section, we are calculating cross terms of cation-anion
print("calculating cation-anion cross terms")
x_product=np.zeros((len(df)))
y_product=np.zeros((len(df)))
z_product=np.zeros((len(df)))
for p1 in range(1,n_cations):
    for p2 in range(n_cations,n_molecules):
      x_product+=x[p1]*x[p2]
      y_product+=y[p1]*y[p2]
      z_product+=z[p1]*z[p2]

final=x_product+y_product+z_product
cols=["delta","product"]

data=np.vstack((delta,final)).T
pf=pd.DataFrame(data=data,columns=cols)
pf['delta'] = (pf['delta'])*20
del data
sf=pf.groupby(["delta"]).mean()
sf.to_csv("summary/cross_cation_anion.csv")
del sf

#In this section, we are calculating cross terms of anion-anion
print("calculating anion-anion cross terms")
x_product=np.zeros((len(df)))
y_product=np.zeros((len(df)))
z_product=np.zeros((len(df)))
for p1 in range(n_cations,n_molecules):
    for p2 in range(p1+1,n_molecules):
      x_product+=x[p1]*x[p2]
      y_product+=y[p1]*y[p2]
      z_product+=z[p1]*z[p2]

final=x_product+y_product+z_product
cols=["delta","product"]

data=np.vstack((delta,final)).T
pf=pd.DataFrame(data=data,columns=cols)
pf['delta'] = (pf['delta'])*20
del data
sf=pf.groupby(["delta"]).mean()
sf.to_csv("summary/cross_anion_anion.csv")
del sf

print("deleting all the files")
shutil.rmtree("data")
shutil.rmtree("com_files")
shutil.rmtree("index_files")
