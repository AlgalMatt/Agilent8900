import os
import pandas as pd
import numpy as np
import sympy as sym
from tkinter import filedialog
import tkinter as tk
import csv
import math
import warnings
import statsmodels.api as sm
from scipy.io import loadmat
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.formula.api import ols

####### Functions ################
#function for finding occurrences of substring (or a selection of substrings) 
#within an array of strings
def contains1d(array1, string1, ret_array=True, case_sensitive=False):
    """
    Looks for occurrences of substring(s) within an array of strings, returning
    a boolean array. This is essentially the Pandas series.str.contains method, 
    but slightly adapted.
    
    Parameters
    ----------
    array1 : 1d array (list, numpy array)
        array to search.
    string1 : string, or 1d array of strings (list or numpy array)
        substrings used to search for within array1.
    ret_array : boolean, optional
        If true, then the output will be a 1d array of the same size as array1, 
        providing True/False values for each element that contains/does not 
        contain any of the substrings within string1. 
        If false, then the output will be a matrix of len(array1) by
        len(string1) with each column being a separate boolean array of 
        occurrences of each substring in string1 within array1.
        The default is True.
    case_sensitive : boolean, optional
        If true, the search will be case-sensitive. The default is False.

    Returns
    -------
    retarray : numpy array or matrix of len(array1)
        An array of boolean values where True values indicate the presence of the 
        substring, string1, at the same index of array1. An element-wise 
        string1 in array1.

    """
    
    #vectorize lower cases
    nlower=np.vectorize(str.lower)
    
    retarray=[]
    #if argument string1 is a single string
    if type(string1)==str:
        #lower all cases
        if case_sensitive==False:
            array1=nlower(array1)
            string1=string1.lower()
        for i in array1:
            retarray.append(string1 in i)   
    #if string1 is a list of strings             
    else:
        #lower all cases
        if case_sensitive==False:
            array1=nlower(array1)
            string1=nlower(string1)
        retarray=np.full((len(array1), len(string1)), False)
        #iterate over the list of substrings
        for j, s in enumerate(string1):
            #iterate over the array of strings a check if the iterated 
            #substring is in it            
            for i, a in enumerate(array1):
                retarray[i, j]=s in a
        #if true, return a 1D array, else it returns a len(array1) by 
        #len(string1) matrix of occurrences of each substring within the array
        if ret_array:
            retarray=np.any(retarray, axis=1)           
    return retarray



#function for calling unique values with order retained.
def unique(list1):
    array1=np.array(list1)
    indx=np.unique(array1, return_index=True)[1]
    indx.sort()
    uarray=array1[indx]
    return uarray

















flick=0.003209 #as of 19/7/22
flickbc=0.003
flickbc2023=0.0095
ratioel='Ca48'   
stgvals={'Li7_No gas':21.27, 'B11_No gas': 89.18, 'Na23_No gas':8.41, 
         'Mg24_No gas':5.49, 'Mg25_No gas':5.49, 'Al27_No gas':54.11,
         'Ca43_No gas': 1, 'Ca48_No gas': 1,
         'Mn55_No gas': 58.34, 'Sr88_No gas':1.53, 'Cd111_No gas': 0.143, 
         'Ba138_No gas':3, 'Nd146_No gas':5.611, 'U238_No gas':159.928}
stgunits={'Li7_No gas':'umol/mol', 'B11_No gas': 'umol/mol', 
          'Na23_No gas':'mmol/mol', 'Mg24_No gas':'mmol/mol', 
          'Mg25_No gas':'mmol/mol', 'Al27_No gas':'umol/mol', 
          'Ca43_No gas':'mol/mol', 'Ca48_No gas':'mol/mol', 
         'Mn55_No gas': 'umol/mol', 'Sr88_No gas':'mmol/mol', 
         'Cd111_No gas': 'umol/mol', 'Ba138_No gas':'umol/mol', 
         'Nd146_No gas':'umol/mol', 'U238_No gas':'nmol/mol'}  

stgarray=np.array(list(stgvals.values()))

#############Import data #################

#load in the standard values set 
stndpath=path=(r"C:\Users\mdumo\OneDrive - University of St Andrews"
                    r"\Agilent\Matt\stndvals.csv")
stndvals_df=pd.read_csv(stndpath)  
stndvals_df=stndvals_df.set_index('Element')
stndval_names=stndvals_df.columns[1:]


reppath=(r"C:/Users/mdumo/OneDrive - University of St Andrews/Agilent"
      r"/Matt/")
CPSarchive = pd.read_csv(reppath+'raw_CPS_T.csv')
Narchive = pd.read_csv(reppath+'raw_N_T.csv')
SDarchive = pd.read_csv(reppath+'raw_SD_T.csv')
intTimearchive = pd.read_csv(reppath+'raw_intTime_T.csv')
Repsarchive = pd.read_csv(reppath+'raw_Reps_T.csv')
PAarchive = pd.read_csv(reppath+'raw_PA_T.csv')

#Combine archives

keepcols=['RunName', 'Time', 'Index', 'Runorder', 'Sample', 'Vial', 
          'Li7', 'B11', 'Na23', 'Mg24', 'Mg25', 'Al27', 'Ca48', 
          'Mn55', 'Zn64', 'Sr88', 'Y89', 'Mo92', 'Cd111', 'Ba138', 'La139',
          'Ce140', 'Pr141', 'Nd146', 'Sm147', 'Eu153', 'Gd157', 'Tb159', 
          'Dy164', 'Ho165', 'Er166', 'Tm169', 'Yb174', 'Lu175', 'Tl205', 
          'U238', 
          'K39_39_H2', 'Ca48_48_H2', 'Mn55_55_H2', 'Fe56_56_H2', 'Ni58_58_H2', 
          'Cu63_63_H2', 'Rb85_85_H2', 
          'P31_47_O2', 'S32_48_O2', 'Ca48_48_O2', 'Co59_75_O2' ]

cpsarchmelt=CPSarchive[keepcols].melt(id_vars=keepcols[0:6], 
                                       var_name='Isotope gas', 
                                       value_name='Mean CPS')
Narchmelt=Narchive[keepcols].melt(id_vars=keepcols[0:6], 
                                       var_name='Isotope gas', 
                                       value_name='N')
intTimearchmelt=intTimearchive[keepcols].melt(id_vars=keepcols[0:6], 
                                       var_name='Isotope gas', 
                                       value_name='intTime')
SDarchmelt=SDarchive[keepcols].melt(id_vars=keepcols[0:6], 
                                       var_name='Isotope gas', 
                                       value_name='SD')
PAarchmelt=PAarchive[keepcols].melt(id_vars=keepcols[0:6], 
                                       var_name='Isotope gas', 
                                       value_name='PA')


bigarch=PAarchmelt.copy()
bigarch['N']=Narchmelt['N']
bigarch['intTime']=intTimearchmelt['intTime']
bigarch['Mean CPS']=cpsarchmelt['Mean CPS']
bigarch['SD']=SDarchmelt['SD']


bigarch.insert(7, 'Gas mode', np.repeat(['No gas'], len(bigarch),
                                                  axis=0))
idx=contains1d(bigarch['Isotope gas'], '_H2', case_sensitive=True)
bigarch.loc[idx, 'Gas mode']='H2'
idx=contains1d(bigarch['Isotope gas'], '_O2', case_sensitive=True)
bigarch.loc[idx, 'Gas mode']='O2'

#reformatting the reps archive
newreparch=Repsarchive[['RunName', 'Time', 'Index', 'Sample']]
for el in unique(cpsarchmelt['Isotope gas']):
    ellist=[el+'_1', el+'_2', el+'_3', el+'_4', el+'_5']
    newreparch[el]=list(np.array(Repsarchive[ellist]))
reparchmelt=newreparch.melt(id_vars=['RunName', 'Time', 'Index', 'Sample'], 
                            var_name='Isotope gas', value_name='Rep CPS')


#join reparchmelt and bigarch
bigarch=pd.merge(bigarch, reparchmelt)

bigarch_nonan=bigarch.dropna()

bigarch_nonan['Time']=pd.to_datetime(bigarch_nonan['Time'])
bigarch_nonan['CPC']=bigarch_nonan['Mean CPS']*bigarch_nonan['intTime']

#allrep_df = pd.read_csv('Reps.csv')

runnames=unique(bigarch_nonan['RunName'])




    
#Use symbolic mode to apply data processing
#Define symbols
x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym=sym.symbols(
    'x_sym xb2_sym xb1_sym y_sym yb2_sym yb1_sym')
xs1_sym, xs2_sym, ys1_sym, ys2_sym=sym.symbols(
    'xs1_sym xs2_sym ys1_sym ys2_sym')
Dts_sym, Dts1b_sym, Dts2b_sym, Dtb_sym=sym.symbols(
    'Dts_sym Dts1b_sym Dts2b_sym Dtb_sym')
cov_xy_sym, cov_xs1ys1_sym, cov_xs2ys2_sym, cov_xb1yb1_sym, cov_xb2yb2_sym= \
    sym.symbols('''cov_xy_sym cov_xs1ys1_sym cov_xs2ys2_sym cov_xb1yb1_sym 
    cov_xb2yb2_sym''')
s_x_sym, s_y_sym, s_xb1_sym, s_xb2_sym, s_yb1_sym, s_yb2_sym, s_xs1_sym, \
    s_xs2_sym, s_ys1_sym, s_ys2_sym=sym.symbols(
    '''s_x_sym s_y_sym s_xb1_sym s_xb2_sym s_yb1_sym s_yb2_sym s_xs1_sym 
    s_xs2_sym s_ys1_sym s_ys2_sym''')

#TE/Ca ratio equation
R_sym=(x_sym - Dtb_sym*xb2_sym + xb1_sym*(Dtb_sym - 1))  \
    /(y_sym - Dtb_sym*yb2_sym + yb1_sym*(Dtb_sym - 1))
    
R_f = sym.lambdify((x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym, 
                    Dtb_sym), R_sym)  

#Variance in R
R_var_sym=s_x_sym**2*R_sym.diff(x_sym)**2+s_y_sym**2*R_sym.diff(y_sym)**2\
    +s_xb1_sym**2*R_sym.diff(xb1_sym)**2+s_xb2_sym**2*R_sym.diff(xb2_sym)**2\
    +s_yb1_sym**2*R_sym.diff(yb1_sym)**2+s_yb2_sym**2*R_sym.diff(yb2_sym)**2\
    +2*cov_xy_sym*R_sym.diff(x_sym)*R_sym.diff(y_sym)\
    +2*cov_xb1yb1_sym*R_sym.diff(xb1_sym)*R_sym.diff(yb1_sym)\
    +2*cov_xb2yb2_sym*R_sym.diff(xb2_sym)*R_sym.diff(yb2_sym)

R_var_f = sym.lambdify((x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym, 
                        Dtb_sym, cov_xy_sym, cov_xb1yb1_sym, cov_xb2yb2_sym, 
                        s_x_sym, s_y_sym, s_xb1_sym, s_xb2_sym, s_yb1_sym, 
                        s_yb2_sym), R_var_sym)    


#Bracketed sample equation                    
B_sym=-(x_sym - Dtb_sym*xb2_sym + xb1_sym*(Dtb_sym - 1))/((((Dts_sym - 1)\
    *(xs1_sym - Dts1b_sym*xb2_sym + xb1_sym*(Dts1b_sym - 1)))\
    /(ys1_sym - Dts1b_sym*yb2_sym + yb1_sym*(Dts1b_sym - 1))\
    -(Dts_sym*(xs2_sym - Dts2b_sym*xb2_sym + xb1_sym*(Dts2b_sym - 1)))\
    /(ys2_sym - Dts2b_sym*yb2_sym + yb1_sym*(Dts2b_sym - 1)))\
    *(y_sym - Dtb_sym*yb2_sym + yb1_sym*(Dtb_sym - 1)))

B_f = sym.lambdify((x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym, Dtb_sym, 
              xs1_sym, ys1_sym, Dts1b_sym, xs2_sym, ys2_sym, Dts2b_sym, 
              Dts_sym), B_sym)       
    
#Variance in B    
B_var_sym=s_x_sym**2*B_sym.diff(x_sym)**2+s_y_sym**2*B_sym.diff(y_sym)**2\
    +s_xb1_sym**2*B_sym.diff(xb1_sym)**2+s_xb2_sym**2*B_sym.diff(xb2_sym)**2\
    +s_xs1_sym**2*B_sym.diff(xs1_sym)**2+s_xs2_sym**2*B_sym.diff(xs2_sym)**2\
    +s_yb1_sym**2*B_sym.diff(yb1_sym)**2+s_yb2_sym**2*B_sym.diff(yb2_sym)**2\
    +s_ys1_sym**2*B_sym.diff(ys1_sym)**2+s_ys2_sym**2*B_sym.diff(ys2_sym)**2\
    +2*cov_xy_sym*B_sym.diff(x_sym)*B_sym.diff(y_sym)\
    +2*cov_xb1yb1_sym*B_sym.diff(xb1_sym)*B_sym.diff(yb1_sym)\
    +2*cov_xb2yb2_sym*B_sym.diff(xb2_sym)*B_sym.diff(yb2_sym)\
    +2*cov_xs1ys1_sym*B_sym.diff(xs1_sym)*B_sym.diff(ys1_sym)\
    +2*cov_xs2ys2_sym*B_sym.diff(xs2_sym)*B_sym.diff(ys2_sym)

ec_x_sym=s_x_sym**2*B_sym.diff(x_sym)**2
ec_y_sym=s_y_sym**2*B_sym.diff(y_sym)**2
ec_xb_sym=s_xb1_sym**2*B_sym.diff(xb1_sym)**2+s_xb2_sym**2*B_sym.diff(xb2_sym)**2
ec_yb_sym=s_yb1_sym**2*B_sym.diff(yb1_sym)**2+s_yb2_sym**2*B_sym.diff(yb2_sym)**2
ec_xs_sym=s_xs1_sym**2*B_sym.diff(xs1_sym)**2+s_xs2_sym**2*B_sym.diff(xs2_sym)**2
ec_ys_sym=s_ys1_sym**2*B_sym.diff(ys1_sym)**2+s_ys2_sym**2*B_sym.diff(ys2_sym)**2
ec_covxy_sym=2*cov_xy_sym*B_sym.diff(x_sym)*B_sym.diff(y_sym)
ec_covxbyb_sym=2*cov_xb1yb1_sym*B_sym.diff(xb1_sym)*B_sym.diff(yb1_sym)+2*cov_xb2yb2_sym*B_sym.diff(xb2_sym)*B_sym.diff(yb2_sym)
ec_covxsys_sym=2*cov_xs1ys1_sym*B_sym.diff(xs1_sym)*B_sym.diff(ys1_sym)\
    +2*cov_xs2ys2_sym*B_sym.diff(xs2_sym)*B_sym.diff(ys2_sym)


B_var_f = sym.lambdify((x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym, 
                        Dtb_sym, xs1_sym, ys1_sym, Dts1b_sym, xs2_sym, ys2_sym
                        ,Dts2b_sym, Dts_sym, cov_xy_sym, cov_xb1yb1_sym, 
                        cov_xb2yb2_sym, cov_xs1ys1_sym, cov_xs2ys2_sym, 
                        s_x_sym, s_y_sym, s_xb1_sym, s_xb2_sym,
                        s_yb1_sym, s_yb2_sym, s_xs1_sym, s_ys1_sym, s_xs2_sym, 
                        s_ys2_sym), B_var_sym)


ec_x_f = sym.lambdify((x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym, 
                        Dtb_sym, xs1_sym, ys1_sym, Dts1b_sym, xs2_sym, ys2_sym
                        ,Dts2b_sym, Dts_sym, cov_xy_sym, cov_xb1yb1_sym, 
                        cov_xb2yb2_sym, cov_xs1ys1_sym, cov_xs2ys2_sym, 
                        s_x_sym, s_y_sym, s_xb1_sym, s_xb2_sym,
                        s_yb1_sym, s_yb2_sym, s_xs1_sym, s_ys1_sym, s_xs2_sym, 
                        s_ys2_sym), ec_x_sym)
ec_y_f = sym.lambdify((x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym, 
                        Dtb_sym, xs1_sym, ys1_sym, Dts1b_sym, xs2_sym, ys2_sym
                        ,Dts2b_sym, Dts_sym, cov_xy_sym, cov_xb1yb1_sym, 
                        cov_xb2yb2_sym, cov_xs1ys1_sym, cov_xs2ys2_sym, 
                        s_x_sym, s_y_sym, s_xb1_sym, s_xb2_sym,
                        s_yb1_sym, s_yb2_sym, s_xs1_sym, s_ys1_sym, s_xs2_sym, 
                        s_ys2_sym), ec_y_sym)
ec_xb_f = sym.lambdify((x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym, 
                        Dtb_sym, xs1_sym, ys1_sym, Dts1b_sym, xs2_sym, ys2_sym
                        ,Dts2b_sym, Dts_sym, cov_xy_sym, cov_xb1yb1_sym, 
                        cov_xb2yb2_sym, cov_xs1ys1_sym, cov_xs2ys2_sym, 
                        s_x_sym, s_y_sym, s_xb1_sym, s_xb2_sym,
                        s_yb1_sym, s_yb2_sym, s_xs1_sym, s_ys1_sym, s_xs2_sym, 
                        s_ys2_sym), ec_xb_sym)
ec_yb_f = sym.lambdify((x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym, 
                        Dtb_sym, xs1_sym, ys1_sym, Dts1b_sym, xs2_sym, ys2_sym
                        ,Dts2b_sym, Dts_sym, cov_xy_sym, cov_xb1yb1_sym, 
                        cov_xb2yb2_sym, cov_xs1ys1_sym, cov_xs2ys2_sym, 
                        s_x_sym, s_y_sym, s_xb1_sym, s_xb2_sym,
                        s_yb1_sym, s_yb2_sym, s_xs1_sym, s_ys1_sym, s_xs2_sym, 
                        s_ys2_sym), ec_yb_sym)
ec_xs_f = sym.lambdify((x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym, 
                        Dtb_sym, xs1_sym, ys1_sym, Dts1b_sym, xs2_sym, ys2_sym
                        ,Dts2b_sym, Dts_sym, cov_xy_sym, cov_xb1yb1_sym, 
                        cov_xb2yb2_sym, cov_xs1ys1_sym, cov_xs2ys2_sym, 
                        s_x_sym, s_y_sym, s_xb1_sym, s_xb2_sym,
                        s_yb1_sym, s_yb2_sym, s_xs1_sym, s_ys1_sym, s_xs2_sym, 
                        s_ys2_sym), ec_xs_sym)
ec_ys_f = sym.lambdify((x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym, 
                        Dtb_sym, xs1_sym, ys1_sym, Dts1b_sym, xs2_sym, ys2_sym
                        ,Dts2b_sym, Dts_sym, cov_xy_sym, cov_xb1yb1_sym, 
                        cov_xb2yb2_sym, cov_xs1ys1_sym, cov_xs2ys2_sym, 
                        s_x_sym, s_y_sym, s_xb1_sym, s_xb2_sym,
                        s_yb1_sym, s_yb2_sym, s_xs1_sym, s_ys1_sym, s_xs2_sym, 
                        s_ys2_sym), ec_ys_sym)
ec_covxy_f = sym.lambdify((x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym, 
                        Dtb_sym, xs1_sym, ys1_sym, Dts1b_sym, xs2_sym, ys2_sym
                        ,Dts2b_sym, Dts_sym, cov_xy_sym, cov_xb1yb1_sym, 
                        cov_xb2yb2_sym, cov_xs1ys1_sym, cov_xs2ys2_sym, 
                        s_x_sym, s_y_sym, s_xb1_sym, s_xb2_sym,
                        s_yb1_sym, s_yb2_sym, s_xs1_sym, s_ys1_sym, s_xs2_sym, 
                        s_ys2_sym), ec_covxy_sym)
ec_covxbyb_f = sym.lambdify((x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym, 
                        Dtb_sym, xs1_sym, ys1_sym, Dts1b_sym, xs2_sym, ys2_sym
                        ,Dts2b_sym, Dts_sym, cov_xy_sym, cov_xb1yb1_sym, 
                        cov_xb2yb2_sym, cov_xs1ys1_sym, cov_xs2ys2_sym, 
                        s_x_sym, s_y_sym, s_xb1_sym, s_xb2_sym,
                        s_yb1_sym, s_yb2_sym, s_xs1_sym, s_ys1_sym, s_xs2_sym, 
                        s_ys2_sym), ec_covxbyb_sym)
ec_covxsys_f = sym.lambdify((x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym, 
                        Dtb_sym, xs1_sym, ys1_sym, Dts1b_sym, xs2_sym, ys2_sym
                        ,Dts2b_sym, Dts_sym, cov_xy_sym, cov_xb1yb1_sym, 
                        cov_xb2yb2_sym, cov_xs1ys1_sym, cov_xs2ys2_sym, 
                        s_x_sym, s_y_sym, s_xb1_sym, s_xb2_sym,
                        s_yb1_sym, s_yb2_sym, s_xs1_sym, s_ys1_sym, s_xs2_sym, 
                        s_ys2_sym), ec_covxsys_sym)





                
#load in the standard values set 
stndpath=path=(r"C:\Users\mdumo\OneDrive - University of St Andrews"
                    r"\Agilent\Matt\stndvals.csv")
stndvals_df=pd.read_csv(stndpath)  
stndvals_df=stndvals_df.set_index('Element')
stndval_names=stndvals_df.columns[1:]


Isotopes=unique(bigarch_nonan['Isotope gas'])
bigarch_nonan.insert(7, 'Element', np.repeat(['xx'], len(bigarch_nonan),
                                                  axis=0))
#Create list of elements w/o gas mode and isotope
Elements=[]
for iso in Isotopes:
    isosplit=iso.split('_')[0]
    el="".join([s for s in isosplit if not s.isdigit()])
    Elements.append(el)
    idx=bigarch_nonan['Isotope gas']==iso
    bigarch_nonan.loc[idx, 'Element']=el



   
#make dict for elements and isotopes
isoel_dict=dict(zip(bigarch_nonan['Isotope gas'], bigarch_nonan['Element']))

#Assign measured isotope names to the standard values dataframe
calivals_df=pd.DataFrame()
missing=[]
for i, iso in enumerate(Isotopes):
    #Check if isotope is included in standard spreadsheet, if not then skip
    if all(stndvals_df.index!=isoel_dict[iso]):
        #note the missing isotopes for later
        missing.append(i)
        continue
        
    #Create dataframe with standard values
    row_df=stndvals_df.loc[isoel_dict[iso], :].copy()
    row_df['Isotope']=iso
    row_df=pd.DataFrame(row_df).transpose().reset_index() 
    row_df.rename(columns={'index':'Element'}, inplace=True)
    row_df.set_index("Isotope", inplace=True)  
    calivals_df=pd.concat([calivals_df,row_df])
    
calivals_df=calivals_df.reset_index()
calivals_df=calivals_df.drop_duplicates(subset='Element')
calivals_df=calivals_df.set_index('Element')

big_df=pd.DataFrame()

for run in runnames:
    rundf=bigarch_nonan.loc[bigarch_nonan['RunName']==run].copy()
    
    
    c4=math.gamma(5/2)/math.gamma((5-1)/2)*(
        2/(5-1))**0.5;
    
    sdcps=rundf['SD']
    
    rundf['SD CPS']=sdcps
    
    
    if rundf.empty:
       continue
    isolist=unique(rundf['Isotope gas'])
    numsamples=max(rundf['Runorder'])
   
    runcps=np.reshape(np.array(rundf['Mean CPS']), (-1, numsamples)).T
    runsd=np.reshape(np.array(rundf['SD CPS']), (-1, numsamples)).T
    runreps=np.reshape(np.array(rundf['Rep CPS']), (-1, numsamples)).T
    
    runcpc=np.reshape(np.array(rundf['CPC']), (-1, numsamples)).T
    
    runcps=pd.DataFrame(runcps, columns=isolist)
    runsd=pd.DataFrame(runsd, columns=isolist)
    runreps=pd.DataFrame(runreps, columns=isolist)
    
    
    runinfo=rundf.iloc[0:numsamples, 0:6]
    runinfo=runinfo.reset_index()
    runinfo['Elapse']=[t-runinfo['Time'][0] for t in runinfo['Time']]
    
    if sum(contains1d(runinfo['Sample'], 'stgfrm'))<3:
        continue
    
    if len(runinfo)/sum(contains1d(runinfo['Sample'], 'stgfrm'))>10:
        continue
    
    blkidx=contains1d(runinfo['Sample'], 'blk')
    brktindx=np.array(contains1d(runinfo['Sample'], 'stgfrm')
                      ) & ~np.array(contains1d(runinfo['Sample'], 'stgfrmx'))
    commonstg=runinfo.loc[brktindx, 'Sample'].value_counts().index[0]
    if commonstg.startswith('0.5'):
        brktindx=contains1d(runinfo['Sample'], '0.5stgfrm')
    else:
        omitindx=contains1d(runinfo['Sample'], '0.5stgfrm')
        brktindx=(brktindx) & (~np.array(omitindx))
    
    if ~np.any(blkidx):
        continue
    
    brktdf=runinfo.loc[brktindx]
    blkdf=runinfo.loc[blkidx]
    
    brktorder=brktdf['Runorder']
    blkorder=blkdf['Runorder']
    #Outliers in brkt and blanks
    for el in isolist:
        brktindx2=np.isin(np.array(runinfo['Runorder']), np.array(brktorder))
        blkindx2=np.isin(np.array(runinfo['Runorder']), np.array(blkorder))
        
        brktrat=runcps.loc[brktindx2, el]/runcps.loc[brktindx2, ratioel]
        blkels=runcps.loc[blkindx2, el]
        brktrat_rsd=((runsd.loc[brktindx2, el]/runcps.loc[brktindx2, el])**2
                     +(runsd.loc[brktindx2, 
                                 ratioel]/runcps.loc[brktindx2, ratioel])**2
                     )**0.5
        
        
               
        q75blk, q25blk = np.percentile(blkels, [75 ,25])
        iqr_blk = q75blk - q25blk
        q75brkt, q25brkt = np.percentile(brktrat, [75 ,25])
        iqr_brkt = q75brkt - q25brkt
        outs_blk=(blkels>iqr_blk*2+q75blk) | (blkels<q25blk-iqr_blk*2) 
        
        q75brktrsd, q25brktrsd = np.percentile(brktrat_rsd, [75 ,25])
        iqr_brktrsd = q75brktrsd - q25brktrsd
        
        
        outs_brkt=((brktrat>iqr_brkt*2+q75brkt) | (brktrat<q25brkt-iqr_brkt*2)
                  | (brktrat_rsd<q25brktrsd-iqr_brktrsd*2) | 
                 (brktrat_rsd>q75brktrsd+iqr_brktrsd*2))
        
        
        
        brktorder=np.delete(np.array(brktorder), outs_brkt)
        blkorder=np.delete(np.array(blkorder), outs_blk)
        
    brktindx=np.isin(np.array(runinfo['Runorder']), np.array(brktorder))
    blkindx=np.isin(np.array(runinfo['Runorder']), np.array(blkorder))        
        
    blkrows=np.array(runcps.index[blkindx])
    brktrows=np.array(runcps.index[brktindx])

    
#################Processing ##################
    y_df=pd.DataFrame(np.tile(runcps[ratioel], [len(isolist), 1]).T, 
                         columns=isolist)
    Ca_sd=pd.DataFrame(np.tile(runsd[ratioel], [len(isolist), 1]).T, 
                       columns=isolist)
    repCPS_y=pd.DataFrame(np.tile(runreps[ratioel], [len(isolist), 1]).T, 
                          columns=isolist)
    
    for g in unique(rundf['Gas mode']):
        if g == 'No gas':
            continue
        cols=unique(rundf.loc[rundf['Gas mode']==g, 'Isotope gas'])
        Camode=cols[contains1d(cols, 'Ca48')]
        for c in cols:    
            y_df.loc[:, c]=runcps[Camode]
            Ca_sd.loc[:, c]=runsd[Camode]
            repCPS_y.loc[:, c]=runreps[Camode]
 

    #Generate covariances dataframe    
    cov_df=pd.DataFrame([])
    #cycle through samples with nested loop of isotopes to get covariances
    for i, row in runreps.iterrows():
        cov_array=np.array([])
        for el in isolist:
            #array of numerator isotopes
            rep_x_array=np.array(list(row[el]))
            #array of denominators (Ca)
            rep_y_array=np.array(list(repCPS_y.loc[i, el]))   
            #covariances
            cov_array=np.append(cov_array, 
                             np.cov(np.vstack((rep_x_array, rep_y_array)))[0, 1])
        #make into dataframe   
        cov_row=pd.DataFrame([cov_array], columns=isolist)   
        cov_df=pd.concat([cov_df, cov_row], ignore_index=True)    
    
    
       

    #Initialise ratio and bracketed data frames
    ratio_smpl_df=runinfo.copy()
    brkt_smpl_df=runinfo.copy()
    ratio_smpl_se_df=runinfo.copy()
    brkt_smpl_se_df=runinfo.copy()
    cali_smpl_df=runinfo.copy()
    cali_smpl_se_df=runinfo.copy()
    ec_1_df=runinfo.copy()
    ec_2_df=runinfo.copy()
    ec_3_df=runinfo.copy()
    ec_4_df=runinfo.copy()
    ec_5_df=runinfo.copy()
    ec_6_df=runinfo.copy()
    ec_7_df=runinfo.copy()
    ec_8_df=runinfo.copy()
    ec_9_df=runinfo.copy()
    Neff_bc_df=runinfo.copy()
    NeffRSD_bc_df=runinfo.copy()
    theo_R_rsd_df=runinfo.copy()
    brkt_r_df=runinfo.copy()
    
    ec_alllist=[]
    #Cycle through sample by sample to calculate R and B
    for i, row in runcps.iterrows():
        #If a blank, skip
        if any(i==blkrows):
            continue
        
        
        #Number of blanks and stnds used for the corrections
        #find closest blank(s)
        blkorder=np.vstack((blkrows, np.abs(blkrows-i)))
        blkorder=np.vstack((blkorder, blkrows-i))
        blkorder = blkorder[:, np.argsort(blkorder[1,:], axis=0)]
        #If sample is before or after all blanks, use just one closest blank
        if np.all(blkrows-i>=0) or np.all(blkrows-i<=0):
            blk_r=np.array([blkorder[0, 0],blkorder[0, 0]])        
        #Otherwise use the two closest (braketing blanks)
        #If the first blank is before the sample in the run find one after
        elif blkorder[2, 0]<0:
            #get the next nearest blk that is of opposite sign direction away
            blk2=blkorder[0, np.where(blkorder[2, 1:]>0)[0][0]+1]          
            blk_r=np.array([blkorder[0, 0], blk2])
        #otherwise find a blank before the sample
        else:
            blk2=blkorder[0, np.where(blkorder[2, 1:]<0)[0][0]+1]          
            blk_r=np.array([blk2,blkorder[0, 0]])
               
        #find closest bracketing standard(s)
        brktorder=np.vstack((brktrows, np.abs(brktrows-i)))
        brktorder=np.vstack((brktorder, brktrows-i))
        brktorder = brktorder[:, np.argsort(brktorder[1,:], axis=0)]
        #If sample is before or after all stnds, use just one closest stnd
        if np.all(brktrows-i>=0) or np.all(brktrows-i<=0) or any(brktrows-i==0):
            brkt_r=np.array([brktorder[0, 0],brktorder[0, 0]])
        #Otherwise use the two closest (braketing standards)
        #If the first stnd is before the sample in the run find one after
        elif brktorder[2, 0]<0:
            #get the next nearest stnd that is of opposite sign direction away
            brk2=brktorder[0, np.where(brktorder[2, 1:]>0)[0][0]+1]          
            brkt_r=np.array([brktorder[0, 0], brk2])
        #otherwise find a stnd before the sample
        else:
            brk2=brktorder[0, np.where(brktorder[2, 1:]<0)[0][0]+1]          
            brkt_r=np.array([brk2,brktorder[0, 0]])
              
        #Assign components of processing
        x=np.array(row[isolist]) #TE sample
        y=np.array(y_df.loc[i, isolist]) #Ca sample
        xb1=np.array(runcps.loc[blk_r[0], isolist]) #TE blank 1
        yb1=np.array(y_df.loc[blk_r[0], isolist])#Ca blank 1
        xb2=np.array(runcps.loc[blk_r[1], isolist])#TE blank 2
        yb2=np.array(y_df.loc[blk_r[1], isolist])#Ca blank 2
        
        xs1=np.array(runcps.loc[brkt_r[0], isolist])#TE stnd 1
        ys1=np.array(y_df.loc[brkt_r[0], isolist])#Ca stnd 1
        xs2=np.array(runcps.loc[brkt_r[1], isolist])#TE stnd 2
        ys2=np.array(y_df.loc[brkt_r[1], isolist])#Ca stnd 2
        
                    
        #Assign errors   
        s_x=np.array(runsd.loc[i, isolist])
        s_y=np.array(Ca_sd.loc[i, isolist])
        s_xb1=np.array(runsd.loc[blk_r[0], isolist])
        s_yb1=np.array(Ca_sd.loc[blk_r[0], isolist])
        s_xb2=np.array(runsd.loc[blk_r[1], isolist])
        s_yb2=np.array(Ca_sd.loc[blk_r[1], isolist])
        s_xs1=np.array(runsd.loc[brkt_r[0], isolist])
        s_ys1=np.array(Ca_sd.loc[brkt_r[0], isolist])
        s_xs2=np.array(runsd.loc[brkt_r[1], isolist])
        s_ys2=np.array(Ca_sd.loc[brkt_r[1], isolist])
        
        #Assign covariances        
        cov_xy=np.array(cov_df.loc[i])
        cov_xb1yb1=np.array(cov_df.loc[blk_r[0]])
        cov_xb2yb2=np.array(cov_df.loc[blk_r[1]])
        cov_xs1ys1=np.array(cov_df.loc[brkt_r[0]])
        cov_xs2ys2=np.array(cov_df.loc[brkt_r[1]])      

        #Assign fractional distance between two blanks
        if blk_r[0]==blk_r[1]:
            #If not between two blanks, then distance is 0
            Dtb=0
            Dts1b=0
            Dts2b=0
            xb2=0;
            yb2=0;
            s_yb2=0;
            s_xb2=0;
            s_c_xb2yb2=0;
        else:
            #Sample
            Dtb=(runinfo.loc[i,'Elapse']-runinfo.loc[blk_r[0], 'Elapse'])/(
                runinfo.loc[blk_r[1], 'Elapse']-runinfo.loc[blk_r[0], 'Elapse'])
            #Bracketing standard 1
            Dts1b=(runinfo.loc[brkt_r[0], 'Elapse']
                   -runinfo.loc[blk_r[0], 'Elapse'])/(
                       runinfo.loc[blk_r[1], 'Elapse']
                       -runinfo.loc[blk_r[0], 'Elapse'])
            #Bracketing standard 2
            Dts2b=(runinfo.loc[brkt_r[1], 'Elapse']
                   -runinfo.loc[blk_r[0], 'Elapse'])/(
                       runinfo.loc[blk_r[1], 'Elapse']
                       -runinfo.loc[blk_r[0], 'Elapse'])   
                    
        
        #Assign fractional distance between bracketing standards
        if brkt_r[0]==brkt_r[1]:
            Dts=0
            if blk_r[0]!=blk_r[1]:
                Dts1b=1
                Dts2b=0
            xs2=0
            ys2=0
            s_ys2=0
            s_xs2=0
            s_c_xs2ys2=0
            
        else:
            Dts=(runinfo.loc[i,'Elapse']-runinfo.loc[brkt_r[0], 'Elapse'])/(
                runinfo.loc[brkt_r[1], 'Elapse']-runinfo.loc[brkt_r[0], 'Elapse'])    
              
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            #TE count ratio and variance
            ratio_smpl=R_f(x, xb2, xb1, y, yb2, yb1, Dtb)
            ratio_smpl_var=R_var_f(x, xb2, xb1, y, yb2, yb1, Dtb, cov_xy, 
                                   cov_xb1yb1, cov_xb2yb2, s_x, s_y, s_xb1, s_xb2, 
                                   s_yb1, s_yb2)
            #Bracketed count ratio and variance
            brkt_smpl = B_f(x, xb2, xb1, y, yb2, yb1, Dtb, xs1, ys1, Dts1b, xs2, 
                            ys2, Dts2b, Dts)
            brkt_smpl_var=B_var_f(x, xb2, xb1, y, yb2, yb1, Dtb, xs1, ys1, Dts1b, 
                                  xs2, ys2, Dts2b, Dts, cov_xy, cov_xb1yb1, 
                                  cov_xb2yb2, cov_xs1ys1, cov_xs2ys2, s_x, s_y,
                                  s_xb1, s_xb2, s_yb1, s_yb2, s_xs1, s_ys1, s_xs2, 
                                  s_ys2)
            
            ec_x=ec_x_f(x, xb2, xb1, y, yb2, yb1, Dtb, xs1, ys1, Dts1b, 
                                  xs2, ys2, Dts2b, Dts, cov_xy, cov_xb1yb1, 
                                  cov_xb2yb2, cov_xs1ys1, cov_xs2ys2, s_x, s_y,
                                  s_xb1, s_xb2, s_yb1, s_yb2, s_xs1, s_ys1, s_xs2, 
                                  s_ys2)/brkt_smpl_var
            ec_y=ec_y_f(x, xb2, xb1, y, yb2, yb1, Dtb, xs1, ys1, Dts1b, 
                                  xs2, ys2, Dts2b, Dts, cov_xy, cov_xb1yb1, 
                                  cov_xb2yb2, cov_xs1ys1, cov_xs2ys2, s_x, s_y,
                                  s_xb1, s_xb2, s_yb1, s_yb2, s_xs1, s_ys1, s_xs2, 
                                  s_ys2)/brkt_smpl_var
            ec_xb=ec_xb_f(x, xb2, xb1, y, yb2, yb1, Dtb, xs1, ys1, Dts1b, 
                                  xs2, ys2, Dts2b, Dts, cov_xy, cov_xb1yb1, 
                                  cov_xb2yb2, cov_xs1ys1, cov_xs2ys2, s_x, s_y,
                                  s_xb1, s_xb2, s_yb1, s_yb2, s_xs1, s_ys1, s_xs2, 
                                  s_ys2)/brkt_smpl_var
            ec_yb=ec_yb_f(x, xb2, xb1, y, yb2, yb1, Dtb, xs1, ys1, Dts1b, 
                                  xs2, ys2, Dts2b, Dts, cov_xy, cov_xb1yb1, 
                                  cov_xb2yb2, cov_xs1ys1, cov_xs2ys2, s_x, s_y,
                                  s_xb1, s_xb2, s_yb1, s_yb2, s_xs1, s_ys1, s_xs2, 
                                  s_ys2)/brkt_smpl_var
            ec_xs=ec_xs_f(x, xb2, xb1, y, yb2, yb1, Dtb, xs1, ys1, Dts1b, 
                                  xs2, ys2, Dts2b, Dts, cov_xy, cov_xb1yb1, 
                                  cov_xb2yb2, cov_xs1ys1, cov_xs2ys2, s_x, s_y,
                                  s_xb1, s_xb2, s_yb1, s_yb2, s_xs1, s_ys1, s_xs2, 
                                  s_ys2)/brkt_smpl_var
            ec_ys=ec_ys_f(x, xb2, xb1, y, yb2, yb1, Dtb, xs1, ys1, Dts1b, 
                                  xs2, ys2, Dts2b, Dts, cov_xy, cov_xb1yb1, 
                                  cov_xb2yb2, cov_xs1ys1, cov_xs2ys2, s_x, s_y,
                                  s_xb1, s_xb2, s_yb1, s_yb2, s_xs1, s_ys1, s_xs2, 
                                  s_ys2)/brkt_smpl_var
            ec_covxy=ec_covxy_f(x, xb2, xb1, y, yb2, yb1, Dtb, xs1, ys1, Dts1b, 
                                  xs2, ys2, Dts2b, Dts, cov_xy, cov_xb1yb1, 
                                  cov_xb2yb2, cov_xs1ys1, cov_xs2ys2, s_x, s_y,
                                  s_xb1, s_xb2, s_yb1, s_yb2, s_xs1, s_ys1, s_xs2, 
                                  s_ys2)/brkt_smpl_var
            ec_covxbyb=ec_covxbyb_f(x, xb2, xb1, y, yb2, yb1, Dtb, xs1, ys1, Dts1b, 
                                  xs2, ys2, Dts2b, Dts, cov_xy, cov_xb1yb1, 
                                  cov_xb2yb2, cov_xs1ys1, cov_xs2ys2, s_x, s_y,
                                  s_xb1, s_xb2, s_yb1, s_yb2, s_xs1, s_ys1, s_xs2, 
                                  s_ys2)/brkt_smpl_var
            ec_covxsys=ec_covxsys_f(x, xb2, xb1, y, yb2, yb1, Dtb, xs1, ys1, Dts1b, 
                                  xs2, ys2, Dts2b, Dts, cov_xy, cov_xb1yb1, 
                                  cov_xb2yb2, cov_xs1ys1, cov_xs2ys2, s_x, s_y,
                                  s_xb1, s_xb2, s_yb1, s_yb2, s_xs1, s_ys1, s_xs2, 
                                  s_ys2)/brkt_smpl_var
        
        
        
                 
        Caloc=[k for k, c in enumerate(isolist==ratioel) if c]           
        cpc_x=runcpc[i, :] 
        cpc_y=runcpc[i, Caloc]
        cpc_b=np.mean([runcpc[blk_r[0], :], runcpc[blk_r[1], :]], axis=0)  
        
        
        
        
        
        
            
        ec_list=[ec_x, ec_y, ec_xb, ec_yb, ec_xs, ec_ys, ec_covxy, ec_covxbyb, 
                 ec_covxsys]    
            
        ec_list=list(np.array(ec_list).T)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
        
        #convert variance into standard error
            ratio_smpl_se=ratio_smpl_var**0.5/c4/5**0.5
            brkt_smpl_se=brkt_smpl_var**0.5/c4/5**0.5
            Neff_bc_smpl=((cpc_x-cpc_b)*(cpc_y-cpc_b[Caloc]))/(
                cpc_x+cpc_b+cpc_y+cpc_b[Caloc])
            NeffRSD_bc_smpl=1/Neff_bc_smpl**0.5
            
            Neff_smpl=cpc_x*cpc_y/(cpc_x+cpc_y)
            theorsd_Neff_smpl=(1/Neff_smpl+flick**2)**0.5
            
            if (np.any(contains1d(rundf['Gas mode'], 'H2')) | 
                np.any(contains1d(rundf['Gas mode'], 'O2'))):
                f=flickbc2023
            else:
                f=flickbc
                
            
            theorsd_Neffbc_smpl=(1/Neff_bc_smpl+f**2)**0.5
        
        #put into dataframes
        ratio_smpl_df.loc[i, isolist]=ratio_smpl
        brkt_smpl_df.loc[i, isolist]=brkt_smpl
        ratio_smpl_se_df.loc[i, isolist]=ratio_smpl_se
        brkt_smpl_se_df.loc[i, isolist]=brkt_smpl_se
        Neff_bc_df.loc[i, isolist]=Neff_bc_smpl
        NeffRSD_bc_df.loc[i, isolist]=NeffRSD_bc_smpl
        theo_R_rsd_df.loc[i, isolist]=theorsd_Neffbc_smpl
        
        ec_1_df.loc[i, isolist]=ec_x
        ec_2_df.loc[i, isolist]=ec_y
        ec_3_df.loc[i, isolist]=ec_xb
        ec_4_df.loc[i, isolist]=ec_yb
        ec_5_df.loc[i, isolist]=ec_xs
        ec_6_df.loc[i, isolist]=ec_ys
        ec_7_df.loc[i, isolist]=ec_covxy
        ec_8_df.loc[i, isolist]=ec_covxbyb
        ec_9_df.loc[i, isolist]=ec_covxsys
        
        brkt_r_df.loc[i, ['upper', 'lower']]=brkt_r
        
       
        
        stdarray=[]
        for iso in isolist:
            if isoel_dict[iso]=='Ca':
                stdarray.append(1)
            else:
                stdarray.append(calivals_df.loc[isoel_dict[iso], 'STGFrm'])
        
        stdarray=np.array(stdarray)
        
        
        #calibrate sample to known standard values
        cali_smpl=brkt_smpl*stdarray
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cali_smpl_se=(brkt_smpl_se/brkt_smpl)*cali_smpl
        cali_smpl_df.loc[i, isolist]=cali_smpl
        cali_smpl_se_df.loc[i, isolist]=cali_smpl_se
    
    
    theo_cali_rsd_df=runinfo.copy()
    for i, row in runcps.iterrows():
        if any(i==blkrows):
            continue
        brkt_u=brkt_r_df.loc[i, 'upper']
        brkt_l=brkt_r_df.loc[i, 'lower']
        theo_x=np.array(theo_R_rsd_df.loc[i, isolist])
        theo_s1=np.array(theo_R_rsd_df.loc[brkt_u, isolist])
        theo_s2=np.array(theo_R_rsd_df.loc[brkt_l, isolist])
        theo_s_mean=np.mean([theo_s1, theo_s2], axis=0)
        theo_cali_rsd_df.loc[i, isolist]=(theo_x**2+theo_s_mean**2)**0.5
        
    
    #merge onto bigdf
    idcols=runinfo.columns
    cali_smpl_melt=cali_smpl_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='Cali')
    cali_se_melt=cali_smpl_se_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='Cali se')
    brkt_melt=brkt_smpl_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='Brkted')
    brkt_se_melt=brkt_smpl_se_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='Brkted se')
    ratio_melt=ratio_smpl_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='Ratio')
    ratio_se_melt=ratio_smpl_se_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='Ratio se')
    
    ec1melt=ec_1_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='ec_x')
    ec2melt=ec_2_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='ec_y')
    ec3melt=ec_3_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='ec_xb')
    ec4melt=ec_4_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='ec_yb')
    ec5melt=ec_5_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='ec_xs')
    ec6melt=ec_6_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='ec_ys')
    ec7melt=ec_7_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='ec_covxy')
    ec8melt=ec_8_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='ec_covxbyb')
    ec9melt=ec_9_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='ec_covxsys')
    
   
    Neff_bc_melt=Neff_bc_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='Neff_bc')
    NeffRSD_bc_melt=NeffRSD_bc_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='NeffRSD_bc')
    
    theo_R_rsd_melt=theo_R_rsd_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='theo_R_rsd')
    theo_cali_rsd_melt=theo_cali_rsd_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='theo_cali_rsd')
    
    Caloc=[k for k, c in enumerate(isolist==ratioel) if c]
    cpcCa=runcpc[:, Caloc]
    Neff=runcpc*cpcCa/(runcpc+cpcCa)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        NeffRSD=1/(Neff**0.5)
    
    Neff_df=runinfo.copy()
    NeffRSD_df=runinfo.copy()
    Neff_df[isolist]=Neff
    NeffRSD_df[isolist]=NeffRSD
    
    Neff_melt=Neff_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='Neff')
    NeffRSD_melt=NeffRSD_df.melt(id_vars=idcols, var_name='Isotope gas',
                                     value_name='NeffRSD')
    
    

    rundf=pd.merge(rundf, brkt_melt)
    rundf=pd.merge(rundf, brkt_se_melt)
    rundf=pd.merge(rundf, ratio_melt)
    rundf=pd.merge(rundf, ratio_se_melt)
    rundf=pd.merge(rundf, cali_smpl_melt)
    rundf=pd.merge(rundf, cali_se_melt)
    rundf=pd.merge(rundf, ec1melt)
    rundf=pd.merge(rundf, ec2melt)
    rundf=pd.merge(rundf, ec3melt)
    rundf=pd.merge(rundf, ec4melt)
    rundf=pd.merge(rundf, ec5melt)
    rundf=pd.merge(rundf, ec6melt)
    rundf=pd.merge(rundf, ec7melt)
    rundf=pd.merge(rundf, ec8melt)
    rundf=pd.merge(rundf, ec9melt)
    rundf=pd.merge(rundf, Neff_melt)
    rundf=pd.merge(rundf, NeffRSD_melt)
    rundf=pd.merge(rundf, Neff_bc_melt)
    rundf=pd.merge(rundf, NeffRSD_bc_melt)
    rundf=pd.merge(rundf, theo_R_rsd_melt)
    rundf=pd.merge(rundf, theo_cali_rsd_melt)
    
    big_df=pd.concat([big_df, rundf], axis=0)
    print(run)
    

big_df.to_csv('big_stgfrm_s_df.csv')

