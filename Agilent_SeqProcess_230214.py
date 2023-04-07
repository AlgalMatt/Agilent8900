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



#Checkbox menu options for selecting samples etc.
#Function for creating a checkbox that returns the indexes of the checked boxes
def fancycheckbox(items,  title="", defaults=None, single=False):
    """    
    Creates a pop-up simple checkbox from a list of items. Returns indexes of 
    the checked items.

    Parameters
    ----------
    items : 1d array (list or numpy array)
        list of items to fill the checkbox.
    title : string, optional
        Descriptive title of the checkbox window. The default is "".
    defaults : boolean array, optional
        Boolean array indicating which items to have check boxes ticked by 
        default. The default is None.
    single : boolean, optional
        If true, only one checkbox can be selected. The default is False.

    Returns
    -------
    selected_indexes : numpy.array
        array of indexes of checked items.

    """
    global selected
    #if no defaults used, create a list of False to implement defaults.
    if defaults is None:
        defaults=[False]*len(items)
    # Create the main window
    window = tk.Tk()
    window.title(title)
    #Keep the window at the front of other apps.
    window.lift()
    window.attributes("-topmost", True)
       
    # Create a list of booleans to store the state of each checkbox
    selected = defaults
    
    # Function to update the list of selected items
    def update_selected(var):
        global selected        
        if single:
            for i in range(len(cb_vars)):
                if i != var:
                    cb_vars[i].set(False)
        selected = [cb_vars[i].get() for i in range(len(cb_vars))]
    
    # Create a list to store the checkbox variables
    cb_vars = []
    
    #How many rows should the checkboxes be limited to
    num_rows = 12
    
    #The title of the window
    label=tk.Label(window, text=title, font=("Helvetica", 14))
    label.grid(row=0, column=len(items)//num_rows//2, pady=5)
    
    # Create a 4-column grid of checkboxes
    for i, item in enumerate(items):
        cb_var = tk.BooleanVar()
        cb = tk.Checkbutton(window, text=item, variable=cb_var,
                            command=lambda var=i: update_selected(var), font=("Arial",12),fg="blue", bg="white")
        cb.grid(row=i%num_rows+1, column=i//num_rows, sticky="W")
        cb_vars.append(cb_var)
        if defaults[i]:
            cb.select()

    # Create a "Submit" button
    submit_button = tk.Button(window, text="Submit", command=lambda: window.destroy())
    submit_button.grid(row=num_rows+2, column=len(items)//num_rows//2)
    
    # Run the main loop
    window.mainloop()
    selected_indexes = np.array([i for i, x in enumerate(selected) if x])
    return selected_indexes


#function for calling unique values with order retained.
def unique(list1):
    array1=np.array(list1)
    indx=np.unique(array1, return_index=True)[1]
    indx.sort()
    uarray=array1[indx]
    return uarray



#############Import data #################

#Select directory
root=tk.Tk()
root.withdraw()
root.attributes('-topmost', True)
folder_select=filedialog.askdirectory()
root.destroy()

#reads batch file
folder_list=os.listdir(folder_select)
Batchlogloc=folder_select+'/BatchLog.csv'
batchdf=pd.read_csv(folder_select+'/BatchLog.csv') 
t=pd.to_datetime(batchdf.iloc[:, 1]) #gets time string and converts to datetime
batch_t=batchdf[(batchdf['Acquisition Result']=='Pass') & \
                (batchdf["Sample Type"].str.contains("Tune")==False)]
batch_t["Acq. Date-Time"]=t #put datetime back into df

#Get the runname
fsplit=folder_select.split('/')
runname=fsplit[-1]

#Get a list of the sample names, make them into directories and put them into the main df
values=[]
for i in batch_t["File Name"]:
    sfsplit=i.split('\\')    
    values.append(folder_select+'/'+sfsplit[-1])
   # print(sfsplit[-1])
batch_t['Sample Folder']=values
batch_t=batch_t.reset_index()
#Setup run info table
Run_df=pd.DataFrame(np.repeat(runname, len(batch_t)), columns=['Run Name'])
Run_df=pd.concat([Run_df, batch_t[['Acq. Date-Time', 'Sample Name', 'Vial#']]],
                 axis=1)
Run_df['Elapse']=batch_t['Acq. Date-Time']-batch_t['Acq. Date-Time'][0]

#empty dataframes
repCPS_all_df=pd.DataFrame()
repPA_all_df=pd.DataFrame()
repSD_all_df=pd.DataFrame()



#Start iterating through samples
for folder in batch_t['Sample Folder']:
    #Get list of subfolders containing replicates and gas modes. Exclude quickscan
    csvlist = [s for s in os.listdir(folder) if ".csv" in s and "quickscan" 
               not in s]
    
    #Find out how many repeats and gas modes there are by reading first sample
    testmode=[]
    headlist=[]
    #iterate through all files in sample directory
    for c in csvlist:
        #open each file and take the gas mode by stripping other content
        with open(folder+'/'+c, newline='') as f:
            reader = csv.reader(f)
            row1 = next(reader)[0].rsplit('/')[-1].strip('\n ')
            testmode.append(row1)
    
    #Get the number of repeats and gases by counting occurrences of gas modes
    numrepeats=testmode.count(list(set(testmode))[0])
    numgases=len(set(testmode))
    
    #Extract all data from all samples
    #Create empty lists      
    listCPS=[]
    listPA=[]
    listSD=[]
    #iterate through replicates  
    for i in range(numrepeats):
        #Empty dataframe for each repeat
        allgas_df=pd.DataFrame()
        #iterate through gas modes
        for j in range(numgases):
            fileloc=folder+'/'+csvlist[i+j*numrepeats] #directory of file
            gas_df=pd.read_csv(fileloc, skiprows=list(range(0, 7)), header=0) 
            
            #find the NaN and "Printed..." footer rows and remove
            dropindex=np.where(gas_df.iloc[:, 0].str.contains("Print", na=True))
            gas_df.drop(dropindex[0], inplace=True)
                                             
            gas_df=gas_df.drop(gas_df.tail(1).index) #remove print info
            #get the current gas mode
            gasmodetxt=testmode[i+j*numrepeats].strip(' ') 
            #Make df of current gas mode
            #Combine mass and element to make isotope column
            gas_df["Isotope_gas"]=gas_df["Element"]+\
                gas_df.iloc[:, 0]+"_"+gasmodetxt 
            gas_df["Gas mode"]=gasmodetxt
                        
            #PA column often wrongly named, so need to rename it.
            #first find the column next to CPS            
            idx=np.where(gas_df.columns == 'CPS')[0]+1
            gas_df['PA']=gas_df.iloc[:, idx]
            #Concat all gas modes of this repeat
            allgas_df=pd.concat([allgas_df, gas_df], ignore_index=True)
        #extract cps, PA and stdev info from all gas modes of this repeat
        listCPS.append((np.array(allgas_df['CPS']))) 
        listPA.append((np.array(allgas_df['PA'])))
        listSD.append((np.array(allgas_df['SD'])))
    
    #arrays of isotopes and gas modes used
    isotopes=np.array(allgas_df['Isotope_gas'])    
    Gasmodes=np.array(allgas_df['Gas mode']) 
    
    #Create nested list of CPS, PA and stdev replicates
    repCPS=[list(s) for s in np.vstack(listCPS).T]
    repPA=[list(s) for s in np.vstack(listPA).T]
    repSD=[list(s) for s in np.vstack(listSD).T]
               
    #Form above lists into df
    repCPS_df=pd.DataFrame([repCPS], columns=isotopes)
    repPA_df=pd.DataFrame([repPA], columns=isotopes)
    repSD_df=pd.DataFrame([repSD], columns=isotopes)
        
    #Concatenate all samples
    repCPS_all_df=pd.concat([repCPS_all_df, repCPS_df], ignore_index=True)       
    repPA_all_df=pd.concat([repPA_all_df, repPA_df], ignore_index=True)  
    repSD_all_df=pd.concat([repSD_all_df, repSD_df], ignore_index=True) 
            
# Add in the info    
repCPS_all_df=pd.concat([Run_df, repCPS_all_df], axis=1)            
repPA_all_df=pd.concat([Run_df, repPA_all_df], axis=1)                 
repSD_all_df=pd.concat([Run_df, repSD_all_df], axis=1) 

#Create means and stdevs.
CPSmean=np.array([])
CPSstd=np.array([])
for lab, row in repCPS_all_df.iterrows():    
    bb=np.array([np.array([np.array(b).mean(), np.array(b).std(ddof=1)]) for
                 b in row[isotopes]])
    CPSmean=np.concatenate([CPSmean, bb[:, 0]])
    CPSstd=np.concatenate([CPSstd, bb[:, 1]])
CPSmean=CPSmean.reshape([-1, len(bb)])
CPSstd=CPSstd.reshape([-1, len(bb)])
#Make dataframes
CPSmean_df=pd.DataFrame(CPSmean, columns=isotopes)
CPSstd_df=pd.DataFrame(CPSstd, columns=isotopes)
# Add in the info 
CPSmean_df=pd.concat([Run_df, CPSmean_df], axis=1)            
CPSstd_df=pd.concat([Run_df, CPSstd_df], axis=1) 

#Create average P/A table
PAarray=[]
for lab, row in repPA_all_df[isotopes].iterrows(): #iterate over rows of df
    els=np.array(list(row))
    PAlist=[]
    for x in els: #iterate through elements
        if all(x=='P') | all(x=='A'):
            PAlist.append(str(x[0]))
        else:
            PAlist.append('M')
    PAarray.append(PAlist) #list of lists
    
    
PA_df=pd.DataFrame(np.array(PAarray), columns=isotopes) #make the df
PA_df=pd.concat([Run_df, PA_df], axis=1) #add the info

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
  
B_var_f = sym.lambdify((x_sym, xb2_sym, xb1_sym, y_sym, yb2_sym, yb1_sym, 
                        Dtb_sym, xs1_sym, ys1_sym, Dts1b_sym, xs2_sym, ys2_sym
                        ,Dts2b_sym, Dts_sym, cov_xy_sym, cov_xb1yb1_sym, 
                        cov_xb2yb2_sym, cov_xs1ys1_sym, cov_xs2ys2_sym, 
                        s_x_sym, s_y_sym, s_xb1_sym, s_xb2_sym,
                        s_yb1_sym, s_yb2_sym, s_xs1_sym, s_ys1_sym, s_xs2_sym, 
                        s_ys2_sym), B_var_sym)


######### Options ########

#Assign default blank indexes based off sample names    
blkdefaults=list(batch_t['Sample Name'].str.contains('blk',
                                                     case=False).astype(int))

#Checkbox for selecting blanks
namelist=[str(i+1)+')  '+s for i, s in enumerate(list(batch_t['Sample Name']))]
blkrows=fancycheckbox(namelist, defaults=blkdefaults, title=("Check the blanks"
                                                             " are selected"))

#Assign default bracket standard indexes based off sample names    
brktdefaults=list(batch_t['Sample Name'].str.contains('stgfrm',
                                                     case=False).astype(int))
#set the dummy standard to false
brktdefaults[0]=False
#Checkbox for selecting bracketing standards
brktrows=fancycheckbox(namelist, defaults=brktdefaults, 
                       title=("Check the bracketing standards are selected"))

#iterate through gas modes to select each ratio element and store as dict
ratioels={}
isotopes_bygas={}
isotopes_bygas_element={}
for gas in unique(Gasmodes):
    gasels=isotopes[Gasmodes==gas]
    ratioel_default=contains1d(gasels, 'Ca48')
    ratioel=gasels[fancycheckbox(gasels, defaults=ratioel_default, 
                                 single=True, title=("Select ratio element"
                                                     ))][0]
    #dict of ratio elements
    ratioels[gas]=ratioel
    #dict of measured full isotope names
    isotopes_bygas[gas]=gasels
    #dict of measured isotope names as element only (no mas or gas mode)
    #This is for referencing with the stndvals_df.
    isotopes_bygas_element[gas]=[
        s.split('_')[0].strip('1234567890') for s in gasels] 
    
#Select calibration method
calistyle_list=['Single-point', 'Calibration curve']
calistyle=calistyle_list[fancycheckbox(calistyle_list, defaults=[True, False],
                                       single=True, title=("Select calibration"
                                                           " method"))[0]]

calinames=unique(batch_t['Sample Name'][brktrows])
calirows=[]
if calistyle=='Calibration curve':
    #Select cali standards
    unique_names=unique(batch_t['Sample Name'])
    #default cali standard names
    calidefaults=contains1d(unique_names, ['stgfrm', 'stgcco', 'stglim',
                                        'stgcrl'])
    #remove dummy standard at beginning
    calidefaults[contains1d(unique_names, 'stgfrmx')]=False
    #Choose calibration standards
    calinames=unique_names[fancycheckbox(unique_names, defaults=calidefaults, 
                                         title=("Select names of the "
                                                "calibration standards"))]
    #Find them in the sequence and get the index
    calindx_default=batch_t['Sample Name'].isin(calinames)
    #Removes extra bracketing standards (those that aren't adjacent)
    calirows=[]
    for i, c in enumerate(calindx_default):
        if c and calindx_default[max([i-1, 0])]==False and \
            calindx_default[min([i+1, len(calindx_default)-1])]==False:
            calindx_default[i]=False
            
    #User select the cali stnds from sequence
    calindx=fancycheckbox(namelist, defaults=calindx_default, 
                          title=("Check that the correct calibration standards"
                                 " are selected"))
else:
    calinames=unique(Run_df.loc[brktrows, 'Sample Name'])
        
                
#load in the standard values set 
stndpath=path=(r"C:\Users\mdumo\OneDrive - University of St Andrews"
                    r"\Agilent\Matt\stndvals.csv")
stndvals_df=pd.read_csv(stndpath)  
stndvals_df=stndvals_df.set_index('Element')
stndval_names=stndvals_df.columns[1:]

#make dict for elements and isotopes
isoel_dict=dict(zip(allgas_df['Isotope_gas'], allgas_df['Element']))

#Assign measured isotope names to the standard values dataframe
calivals_df=pd.DataFrame()
missing=[]
for i, iso in enumerate(isotopes):
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
    
  
#Manually associate each bracketing or cali standard to one in stndvals.
stnd_dict={}
for cali in calinames:    
    associate_default=[x in cali for x in stndval_names]
    cal_associate=stndval_names[fancycheckbox(stndval_names, 
                                              defaults=associate_default,
                                              single=True, 
                                              title=("Link {} to the correct " 
                                              "standard name".format(cali)))]   
    stnd_dict[cali]=calivals_df[cal_associate]


#Remove some standards from cali curves
#cycle through elements, list P/A, find majority, reject minority
if calistyle=='Calibration curve':
    PA_df_cali=PA_df.loc[calindx, :]
    for iso in isotopes:             
        PAcounts=PA_df_cali[iso].value_counts()
        PA_minor=PA_df_cali[PA_df_cali[iso]!=PAcounts.index[0]]['Sample Name']
        if PA_minor.shape[0]>0:       
            for st in PA_minor:
                stnd_dict[st].loc[iso]=np.NaN
    
        


#################Processing ##################

#create array of Ca counts using different gases
y_df=Run_df.copy() #mean Ca cps
s_y_df=Run_df.copy() #sd Ca
repCPS_y_df=Run_df.copy() #replicate Ca cps
   
for gas in list(ratioels.keys()):
    #Create matrix of Ca cps values
    rat_mat=np.tile(np.array(CPSmean_df[ratioels[gas]]), 
                    [len(isotopes_bygas[gas]), 1])
    sd_mat=np.tile(np.array(CPSstd_df[ratioels[gas]]), 
                    [len(isotopes_bygas[gas]), 1])
    
    rep_mat=np.tile(np.array(repCPS_all_df[ratioels[gas]]), 
                    [len(isotopes_bygas[gas]), 1])
    #convert to data frame
    rat_df=pd.DataFrame(rat_mat.T, columns=isotopes_bygas[gas])
    y_df=pd.concat([y_df, rat_df], axis=1)
    
    sd_df=pd.DataFrame(sd_mat.T, columns=isotopes_bygas[gas])
    s_y_df=pd.concat([s_y_df, sd_df], axis=1)
    
    rep_df=pd.DataFrame(rep_mat.T, columns=isotopes_bygas[gas])
    repCPS_y_df=pd.concat([repCPS_y_df, rep_df], axis=1)



#Generate covariances dataframe    
cov_df=pd.DataFrame([])
#cycle through samples with nested loop of isotopes to get covariances
for i, row in repCPS_all_df.iterrows():
    cov_array=np.array([])
    for iso in isotopes:
        #array of numerator isotopes
        rep_x_array=np.array(list(row[iso]))
        #array of denominators (Ca)
        rep_y_array=np.array(list(repCPS_y_df.loc[i, iso]))   
        #covariances
        cov_array=np.append(cov_array, 
                         np.cov(np.vstack((rep_x_array, rep_y_array)))[0, 1])
    #make into dataframe   
    cov_row=pd.DataFrame([cov_array], columns=isotopes)   
    cov_df=pd.concat([cov_df, cov_row], ignore_index=True)        

#Initialise ratio and bracketed data frames
ratio_smpl_df=Run_df.copy()
brkt_smpl_df=Run_df.copy()
ratio_smpl_se_df=Run_df.copy()
brkt_smpl_se_df=Run_df.copy()
cali_smpl_df=Run_df.copy()
cali_smpl_se_df=Run_df.copy()

#get the calibration isotopes
cali_isos=isotopes.copy()
cali_isos=np.delete(cali_isos, missing)
   
#Cycle through sample by sample to calculate R and B
for i, row in CPSmean_df.iterrows():
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
    x=np.array(row[isotopes]) #TE sample
    y=np.array(y_df.loc[i, isotopes]) #Ca sample
    xb1=np.array(CPSmean_df.loc[blk_r[0], isotopes]) #TE blank 1
    yb1=np.array(y_df.loc[blk_r[0], isotopes])#Ca blank 1
    xb2=np.array(CPSmean_df.loc[blk_r[1], isotopes])#TE blank 2
    yb2=np.array(y_df.loc[blk_r[1], isotopes])#Ca blank 2
    
    xs1=np.array(CPSmean_df.loc[brkt_r[0], isotopes])#TE stnd 1
    ys1=np.array(y_df.loc[brkt_r[0], isotopes])#Ca stnd 1
    xs2=np.array(CPSmean_df.loc[brkt_r[1], isotopes])#TE stnd 2
    ys2=np.array(y_df.loc[brkt_r[1], isotopes])#Ca stnd 2
    
    #Assign fractional distance between two blanks
    if blk_r[0]==blk_r[1]:
        #If not between two blanks, then distance is 0
        Dtb=0
        Dts1b=0
        Dts2b=0
    else:
        #Sample
        Dtb=(row['Elapse']-Run_df.loc[blk_r[0], 'Elapse'])/(
            Run_df.loc[blk_r[1], 'Elapse']-Run_df.loc[blk_r[0], 'Elapse'])
        #Bracketing standard 1
        Dts1b=(Run_df.loc[brkt_r[0], 'Elapse']
               -Run_df.loc[blk_r[0], 'Elapse'])/(
                   Run_df.loc[blk_r[1], 'Elapse']
                   -Run_df.loc[blk_r[0], 'Elapse'])
        #Bracketing standard 2
        Dts2b=(Run_df.loc[brkt_r[1], 'Elapse']
               -Run_df.loc[blk_r[0], 'Elapse'])/(
                   Run_df.loc[blk_r[1], 'Elapse']
                   -Run_df.loc[blk_r[0], 'Elapse'])   
    
    #Assign fractional distance between bracketing standards
    if brkt_r[0]==brkt_r[1]:
        Dts=0
    else:
        Dtb=(row['Elapse']-Run_df.loc[brkt_r[0], 'Elapse'])/(
            Run_df.loc[brkt_r[1], 'Elapse']-Run_df.loc[brkt_r[0], 'Elapse'])    
    
    #Assign errors   
    s_x=np.array(CPSstd_df.loc[i, isotopes])
    s_y=np.array(s_y_df.loc[i, isotopes])
    s_xb1=np.array(CPSstd_df.loc[blk_r[0], isotopes])
    s_yb1=np.array(s_y_df.loc[blk_r[0], isotopes])
    s_xb2=np.array(CPSstd_df.loc[blk_r[1], isotopes])
    s_yb2=np.array(s_y_df.loc[blk_r[1], isotopes])
    s_xs1=np.array(CPSstd_df.loc[brkt_r[0], isotopes])
    s_ys1=np.array(s_y_df.loc[brkt_r[0], isotopes])
    s_xs2=np.array(CPSstd_df.loc[brkt_r[1], isotopes])
    s_ys2=np.array(s_y_df.loc[brkt_r[1], isotopes])
    
    #Assign covariances        
    cov_xy=np.array(cov_df.loc[i])
    cov_xb1yb1=np.array(cov_df.loc[blk_r[0]])
    cov_xb2yb2=np.array(cov_df.loc[blk_r[1]])
    cov_xs1ys1=np.array(cov_df.loc[brkt_r[0]])
    cov_xs2ys2=np.array(cov_df.loc[brkt_r[1]])
    
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
    
    
     
    #c4 function for adjusting SE for low number of measurements
    c4=math.gamma(numrepeats/2)/math.gamma((numrepeats-1)/2)*(
        2/(numrepeats-1))**0.5;
    #convert variance into standard deviation
    ratio_smpl_se=ratio_smpl_var**0.5/c4/numrepeats**0.5
    brkt_smpl_se=brkt_smpl_var**0.5/c4/numrepeats**0.5
    
    #put into dataframes
    ratio_smpl_df.loc[i, isotopes]=ratio_smpl
    brkt_smpl_df.loc[i, isotopes]=brkt_smpl
    ratio_smpl_se_df.loc[i, isotopes]=ratio_smpl_se
    brkt_smpl_se_df.loc[i, isotopes]=brkt_smpl_se
         
    #calibration (single point)
    if calistyle=='Single-point': 
        #delete the elements missing from the cali standard list
        b2=np.delete(brkt_smpl, missing)
        b2_se=np.delete(brkt_smpl_se, missing)        
        #make standard array based on the bracketing standard
        #Note, this does mean that different bracketing standards can be used
        #throughout the run
        cali_array=np.array(stnd_dict[Run_df.loc[brkt_r[0], 'Sample Name']])
        
        #calibrate sample to known standard values
        cali_smpl=b2*np.squeeze(cali_array)
        cali_smpl_se=(b2_se/b2)*cali_smpl
        cali_smpl_df.loc[i, cali_isos]=cali_smpl
        cali_smpl_se_df.loc[i, cali_isos]=cali_smpl_se
        
    
#calibration (cali curve)   
if calistyle=='Calibration curve':
    stnd_array=np.empty((0, len(calivals_df)))  
    stnd_df=Run_df.loc[calindx]
    for c in calindx:
        s_array=np.array(stnd_dict[Run_df.loc[c, 'Sample Name']]).T
        stnd_array=np.append(stnd_array, s_array, axis=0)
    
    stnd_df[cali_isos]=stnd_array
    
    #perform the cali curve
    for i, iso in enumerate(cali_isos):
       #get the bracketed values of the calibration standards
        brktiso_arr=np.array(brkt_smpl_df.loc[calindx, iso])
        #and their SEs
        brktiso_se_arr=np.array(brkt_smpl_se_df.loc[calindx, iso], 
                                dtype=np.float64)
        #assign the true values to an array
        stndiso_arr=np.array(stnd_df[iso])
        
        #determine weights based on SE
        w=1/brktiso_se_arr**2
        
        #Reshape and add intercept to bracketed values
        X=brktiso_arr
        Xwint=np.empty(shape=(len(X), 2), dtype=np.float64)
        Xwint[:,0]=1
        Xwint[:, 1]=X
        #reshape true values
        Y=stndiso_arr.reshape(-1, 1)
        #Perform weighted least squares
        mdl = sm.WLS(Y, Xwint, weights=w)
        res_wls = mdl.fit()
        
        #get R-squared of WLS
        r_sq = res_wls.rsquared
        #get parameters and their SEs of WLS (0 = intercept, 1 = gradient)
        params=res_wls.params 
        params_se=res_wls.bse 
        
        #Apply regression to samples.
        cali_smpl_df[iso]=brkt_smpl_df[iso]*params[1]+params[0]
        #propagate uncertainty
        cali_smpl_se_df[iso]=(((brkt_smpl_se_df[iso]/brkt_smpl_df[iso])**2 
            + (params_se[1]/params[1])**2)*cali_smpl_df[iso]**2
                              +params_se[0]**2)**0.5
       
        
#Enter name for file
def textinputbox():
    root=tk.Tk()
    #Button function that saves input value and closes the window
    def retrieve_input():
        global inputValue
        inputValue=textBox.get("1.0","end-1c")
        root.destroy()
    #Text input
    textBox=tk.Text(root, height=2, width=10)
    textBox.grid(row=0, column=0)
    #Save button (see function above)
    buttonSave=tk.Button(root, height=1, width=10, text="Save", 
                        command=lambda: retrieve_input())
    buttonSave.grid(row=1, column=0)  
    
    tk.mainloop()   
    return inputValue    
        
savefile=textinputbox()

        
       
        
            
    
    
    





