# -*- coding: utf-8 -*-
"""
Created on Tue May 18 11:29:24 2021

@author: ave41
"""
#%%
# initialising timer so we can count elapsed time
from pytictoc import TicToc
t = TicToc() # create TicToc instance
t.tic() # Start timer

# importing packages
import sys
import os

# initialising starting directory
home_path = "C:\\Users\\ave41\\OneDrive - University of Canterbury\\Master's 2021\\" \
            "ASTR480 Research\\ASTR480 Code\\Data Reduction Pipeline\\DataReductionPipeline\\src"
os.chdir(home_path) #from now on, we are in this directory

# importing functions
from drp_funcs import *

# ASSUMPTIONS: all files we wish to reduce are in the starting folder (first path)
# and are sorted into ALERT, DARK and FLAT folders

###############################################################################
#---------------------SECTION ONE: SORTING ALERTS-----------------------------# 
###############################################################################

# changing to ALERT folder
path = "C:\\Users\\ave41\\OneDrive - University of Canterbury\\Master's 2021\\" \
       "ASTR480 Research\\ASTR480 Code\\Data Reduction Pipeline\\ObsData_v4\\ALERT"
os.chdir(path) #from now on, we are in this directory

# making list of all files in ALERT folder
fits_files = []
for entry in os.listdir(path):
    fits_files.append(entry)

# getting target names for all files
targetnames = []
for file in fits_files:
    hdul = fits.open(file)
    targetnames.append(hdul[0].header['FIELD'])

# getting list of non-duplicated target names
# refer to drp_funcs.py/def exptime_checker(IMAGElist) for structure idea.
non_duplicated_targetnames = []
for target in targetnames:
    if target not in non_duplicated_targetnames:
    # this condition is to ensure that we don't have a list with 
    # repeating target names. 
        non_duplicated_targetnames.append(target)

# making dictionary of target name as key and list of corresponding files as value
target_names_dict = dict()
for a_target in non_duplicated_targetnames:
    lst_of_files = []
    for file in fits_files:
        this_file = fits.open(file)
        this_file_target = this_file[0].header['FIELD']
        if a_target == this_file_target:
            lst_of_files.append(file)
        else:
            pass
    target_names_dict.update({a_target: lst_of_files})

# getting exposore times for each target
exptimes = []
for key,value in target_names_dict.items():
    exptimes.append([key,exptime_checker(value)])   

###############################################################################
#---------------------SECTION TWO: DATA REDUCTION-----------------------------# 
###############################################################################
#%%
##---------------------------MAKING MASTER BIASES----------------------------##
# reading in bias files from BIAS folder
BIAS_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/"
                 "Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/"
                 "ObsData_v4/DARK")
# making/checking MBIAS path/folder
MBIAS_path = path_checker(BIAS_path,'Master Biases')

BIAS_imgs = ImageFileCollection(BIAS_path,glob_exclude=['/*-0.fit','/*-99.fit'])

# selecting images
BIAS_files = BIAS_imgs.files_filtered(EXPTIME=0,include_path=True)
BIAS_chips_files = chip_separator(BIAS_files)

# calling mbias function to make master biases for each chip
mbias_maker(BIAS_chips_files,MBIAS_path)

# reading in master bias files from Master Biases folder
MBIAS_imgs = ImageFileCollection(MBIAS_path, keywords='*')
MBIAS_files = MBIAS_imgs.files_filtered(COMBINED=True,
                                        include_path=True)
MBIAS_chips_files = chip_separator(MBIAS_files)

# uncomment the following line if you want image count statistics
# code will take ~6 mins to run
# img_stats(BIAS_files)
#%%
##-----------------------------CALIBRATING DARKS-----------------------------##
# reading in dark files from DARK folder
DARK_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/ObsData_v4/DARK")
# making/checking Calibrated Darks path/folder
DARK_cal_path = path_checker(DARK_path,'Calibrated Darks')
# making/checking MDARK path/folder
MDARK_path = path_checker(DARK_path,'Master Darks')

DARK_imgs = ImageFileCollection(DARK_path,glob_exclude=['/*-0.fit','/*-99.fit'])

# selecting images
DARK_files = DARK_imgs.files_filtered(FIELD='              dark',include_path=True)
DARK_chips_files = chip_separator(DARK_files)

# calling dark_calibrator function to calibrate all the darks 
dark_calibrator(DARK_chips_files,MBIAS_chips_files,DARK_cal_path)

# reading in calibrated dark files from Calibrated Darks folder
# ,glob_exclude=['/*-0.fit','/*-99.fit','/*-1-*.fit']
DARK_cal_imgs = ImageFileCollection(DARK_cal_path)
DARK_cal_files = DARK_cal_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
                                              include_path=True)
DARK_cal_chips_files = chip_separator(DARK_cal_files)

# uncomment the following line if you want image count statistics
# code will take ~6 mins to run
# img_stats(DARK_files)
#%%
##----------------------------MAKING MASTER DARKS----------------------------##
# calling mdark_maker function to make master darks
mdark_maker(DARK_cal_chips_files,MDARK_path)


# reading in master dark files from Master Darks folder
# excluding non-science quicklook images and biases
MDARK_imgs = ImageFileCollection(MDARK_path, keywords='*')
MDARK_files = MDARK_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
                                        include_path=True)
MDARK_ccds = MDARK_imgs.ccds(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>', 
                             combined=True)
MDARK_chips_files = chip_separator(MDARK_files)

#%%
# ##--------------------------------FLATS-----------------------------------##
FLAT_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/ObsData_v3/FLAT")

FLAT_imgs = ImageFileCollection(FLAT_path,glob_exclude=['*[*-0.fit]'])
FLAT_files = FLAT_imgs.files_filtered(FIELD ='              flat',include_path=True)
# making/checking MFLAT path/folder
FLAT_cal_path = path_checker(FLAT_path,'Calibrated Flats')
MFLAT_path = path_checker(FLAT_path,'Master Flats')
# FLAT_exptimes = exptime_checker(FLAT_files)
FLAT_chips_files = chip_separator(FLAT_files)


flat_calibrator(FLAT_chips_files,MDARK_chips_files,FLAT_cal_path)

FLAT_cal_imgs = ImageFileCollection(FLAT_cal_path)
FLAT_cal_files = FLAT_cal_imgs.files_filtered(FIELD   = '              flat' ,
                                              include_path=True)
FLAT_cal_chips_files = chip_separator(FLAT_cal_files)

##----------------------------MAKING MASTER FLATS----------------------------##
mflat_maker(FLAT_cal_chips_files,MFLAT_path)

#%%
###############################################################################
#--------------------SECTION THREE: IMAGE CALIBRATION-------------------------# 
###############################################################################





#%%
#================================ don't touch ================================#

###############################################################################
#-------------------------------END OF CODE-----------------------------------# 
############################################################################### 

t.toc() # Print elapsed time
    
###############################################################################
#-------------------------------END OF CODE-----------------------------------# 
############################################################################### 