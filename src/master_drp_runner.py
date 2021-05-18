# -*- coding: utf-8 -*-
"""
Created on Tue May 18 11:29:24 2021

@author: ave41
"""

import sys
import os
# sys.path.insert(1,"C:\\Users\\ave41\\OneDrive - University of Canterbury\\Master's 2021\\ASTR480 Research\\ASTR480 Code\\Data Reduction Pipeline\\ObsData_v4\\ALERT")
path = "C:\\Users\\ave41\\OneDrive - University of Canterbury\\Master's 2021\\ASTR480 Research\\ASTR480 Code\\Data Reduction Pipeline\\DataReductionPipeline\\src"
os.chdir(path)
from drp_funcs import *
path = "C:\\Users\\ave41\\OneDrive - University of Canterbury\\Master's 2021\\ASTR480 Research\\ASTR480 Code\\Data Reduction Pipeline\\ObsData_v4\\ALERT"
os.chdir(path) #from now on, we are in this directory


# pseudocode for DRP:
#     1. read headers of fits files
#     2. make list of fits files whose Target names are the same
#     3. find matching Dark + Flat data
#     4. run DRP as normal for those
#     5. repeat for other targets
    
# ASSUMPTIONS: all files we wish to reduce are in the starting folder (path above)
# and are sorted into ALERT, DARK and FLAT folders

# step one: go to ALERT folder
# step two: read in all fits files
# step three: for each file, read header and put target name, fits filename into lists

# hdul = fits.open(fits_image_filename)

fits_files = []
for entry in os.listdir(path):
    fits_files.append(entry)

targetnames = []
files = []
for file in fits_files:
    hdul = fits.open(file)
    targetnames.append(hdul[0].header['FIELD'])

# refer to drp_funcs.py/def exptime_checker(IMAGElist) for structure idea.
non_duplicated_targetnames = []
for target in targetnames:
    if target not in non_duplicated_targetnames:
    # this condition is to ensure that we don't have a list with 
    # repeating target names. 
        non_duplicated_targetnames.append(target)

for a_target in non_duplicated_targetnames:
    lst_of_files = []
    for file in fits_files:
        this_file = fits.open(file)
        this_file_target = hdul[0].header['FIELD']
        if a_target == this_file_target:
            lst_of_files.append(file)
        else:
            pass
        target_names_dict = {a_target:lst_of_files}  #not giving correct output as expected >:(
      
    
    
    
    
    
    
    
    
    
    
    
    
    