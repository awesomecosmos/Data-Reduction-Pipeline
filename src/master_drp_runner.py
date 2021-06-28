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
ALERT_path = "//spcsfs/ave41/astro/ave41/ObsData_v6/ALERT"
os.chdir(ALERT_path) #from now on, we are in this directory

# making list of all files in ALERT folder
# selecting images and excluding non-science images
good_ALERT_files = []
ALERTs_to_include = ['/*-1.fit','/*-2.fit','/*-3.fit','/*-4.fit','/*-5.fit',
              '/*-6.fit','/*-7.fit','/*-8.fit','/*-9.fit','/*-10.fit']
for i in ALERTs_to_include:
    good_ALERT_file = glob.glob(ALERT_path + i)
    good_ALERT_files += good_ALERT_file

# # getting target names for all files
targetnames = []
for file in good_ALERT_files:
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
    for file in good_ALERT_files:
        this_file = fits.open(file)
        this_file_target = this_file[0].header['FIELD']
        if a_target == this_file_target:
            lst_of_files.append(file)
        else:
            pass
    target_names_dict.update({a_target:lst_of_files})

# getting exposore times for each target
exptimes = []
for key,value in target_names_dict.items():
    exptimes.append([key,exptime_checker(value)])   

##--------------------------PREPARATION WORK---------------------------------##
# this is a list of file extensions we wish to keep. 
# It is assumed here we want Chips 1-10 only.
to_include = ['/*-1.fit','/*-2.fit','/*-3.fit','/*-4.fit','/*-5.fit',
              '/*-6.fit','/*-7.fit','/*-8.fit','/*-9.fit','/*-10.fit']

##-------------------------------PATHWORK------------------------------------##
reduced_ALERT_path = path_checker(ALERT_path,'Reduced ALERT')
# reading in bias files from BIAS folder
BIAS_path = Path("//spcsfs/ave41/astro/ave41/ObsData_v6/DARK")
# making/checking MBIAS path/folder
MBIAS_path = path_checker(BIAS_path,'Master Biases')
# reading in dark files from DARK folder
DARK_path = Path("//spcsfs/ave41/astro/ave41/ObsData_v6/DARK")
# making/checking Calibrated Darks path/folder
DARK_cal_path = path_checker(DARK_path,'Calibrated Darks')
# making/checking MDARK path/folder
MDARK_path = path_checker(DARK_path,'Master Darks')

FLAT_path = Path("//spcsfs/ave41/astro/ave41/ObsData_v6/FLAT")

# making/checking Calibrated Flats path/folder
FLAT_cal_path = path_checker(FLAT_path,'Calibrated Flats')
# making/checking MFLAT path/folder
MFLAT_path = path_checker(FLAT_path,'Master Flats')
# making/checking MFLAT_counts path/folder
MFLAT_counts_path = path_checker(MFLAT_path,'Master Flats by Counts')

#%%
###############################################################################
#---------------------SECTION TWO: DATA REDUCTION-----------------------------# 
###############################################################################

#%%
##---------------------------MAKING MASTER BIASES----------------------------##
# selecting images
BIAS_imgs = ImageFileCollection(BIAS_path,glob_exclude=['/*-0.fit','/*-99.fit'])
# make this into a try/else block?
BIAS_files = BIAS_imgs.files_filtered(EXPTIME=1,include_path=True) # EXPTIME is either 0 or 1
BIAS_chips_files = chip_separator(BIAS_files)

BIAS_counts = img_counts(BIAS_files)

# calling mbias function to make master biases for each chip
mbias_files_for_log = mbias_maker(BIAS_chips_files,MBIAS_path)

# reading in master bias files from Master Biases folder
MBIAS_imgs = ImageFileCollection(MBIAS_path, keywords='*')
MBIAS_files = MBIAS_imgs.files_filtered(COMBINED=True,
                                        include_path=True)
MBIAS_chips_files = chip_separator(MBIAS_files)

MBIAS_counts = img_counts(MBIAS_files)

# uncomment the following line if you want image count statistics
# code will take ~6 mins to run
# img_stats(BIAS_files)
#%%
##-----------------------------CALIBRATING DARKS-----------------------------##
# selecting images and excluding non-science images
good_files = []
for i in to_include:
    good_file = glob.glob("//spcsfs/ave41/astro/ave41/ObsData_v6/DARK" 
                          + i)
    good_files += good_file

# selecting images
DARK_imgs = ImageFileCollection(filenames=good_files)
DARK_files = DARK_imgs.files_filtered(FIELD='              dark',include_path=True)

# now we need to get rid of the biases (0s and 1s exposures)
# the >1 condition will take care of it - no darks strictly less than 1s 
# will be calibrated. This is because sometimes the biases can be either 0s or 1s.
good_DARK_files = []
for DARK_file in DARK_files:
    hdu1 = fits.open(DARK_file)
    dark_exptime = hdu1[0].header['EXPTIME']
    if dark_exptime >1:
        good_DARK_files.append(DARK_file)

DARK_chips_files = chip_separator(good_DARK_files)

DARK_counts = img_counts(good_DARK_files)

# calling dark_calibrator function to calibrate all the darks 
dark_calibrator(DARK_chips_files,MBIAS_chips_files,DARK_cal_path)

print('DARK_counts',DARK_counts)
print(" ")


##----------------------------MAKING MASTER DARKS----------------------------##
# reading in calibrated dark files from Calibrated Darks folder
DARK_cal_imgs = ImageFileCollection(DARK_cal_path)
DARK_cal_files = DARK_cal_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
                                              include_path=True)
DARK_cal_chips_files = chip_separator(DARK_cal_files)

DARK_cal_counts = img_counts(DARK_cal_files)

print('DARK_cal_counts',DARK_cal_counts)
print(" ")

# uncomment the following line if you want image count statistics
# code will take ~6 mins to run
# img_stats(DARK_files)
#%%
# calling mdark_maker function to make master darks
mdark_maker(DARK_cal_chips_files,MDARK_path)

# reading in master dark files from Master Darks folder
MDARK_imgs = ImageFileCollection(MDARK_path, keywords='*')
MDARK_files = MDARK_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
                                        include_path=True)
#%%
MDARK_ccds = MDARK_imgs.ccds(COMBINED=True)
MDARK_chips_files = chip_separator(MDARK_files)
# MDARK_counts = img_counts(MDARK_files)

##----------------------------------FLATS------------------------------------##
# selecting images and excluding non-science images
good_files = []
for i in to_include:
    good_file = glob.glob("//spcsfs/ave41/astro/ave41/ObsData_v6/FLAT" 
                          + i)
    good_files += good_file

# selecting images
FLAT_imgs = ImageFileCollection(filenames=good_files)
FLAT_files = FLAT_imgs.files_filtered(FIELD ='              flat',
                                  include_path=True)

# FLAT_counts = img_counts(FLAT_files)

# sorting files appropriately for future use
FLAT_exptimes = exptime_checker(FLAT_files)
FLAT_chips_files = chip_separator(FLAT_files)

# finding closest dark exposure times to flat exposure times
n_combined_dark = len(MDARK_files)
expected_exposure_times = set(FLAT_exptimes) 
actual_exposure_times = set(h['EXPTIME'] for h in MDARK_imgs.headers(combined=True))
combined_darks = {ccd.header['EXPTIME']: ccd for ccd in MDARK_ccds}
#%%
# calling flat_calibrator function to calibrate all the darks 
flat_calibrator(FLAT_chips_files,MDARK_chips_files,FLAT_cal_path,
                actual_exposure_times,combined_darks)
#%%
##----------------------------MAKING MASTER FLATS----------------------------##
# reading in calibrated flat files from Calibrated Flats folder
FLAT_cal_imgs = ImageFileCollection(FLAT_cal_path)
FLAT_cal_files = FLAT_cal_imgs.files_filtered(FIELD   = '              flat' ,
                                              include_path=True)
FLAT_cal_chips_files = chip_separator(FLAT_cal_files)

FLAT_cal_counts = img_counts(FLAT_cal_files)

# calling mflat_maker function to make master darks
mflat_maker(FLAT_cal_chips_files,MFLAT_path)

# reading in master flat files from Master Flats folder
MFLAT_imgs = ImageFileCollection(MFLAT_path, keywords='*')
MFLAT_files = MFLAT_imgs.files_filtered(FIELD   = '              flat',
                                        include_path=True)
MFLAT_ccds = MFLAT_imgs.ccds(FIELD   = '              flat', 
                             combined=True)
MFLAT_chips_files = chip_separator(MFLAT_files)
MFLAT_counts = img_counts(MFLAT_files)
#%%
# categorising master flats by number of counts
hi_counts_flats,ok_counts_flats,lo_counts_flats,last_resort_flats,trash_counts_flats=flats_count_classifier(MFLAT_files)
counts_sep_flats = [hi_counts_flats,ok_counts_flats,lo_counts_flats,last_resort_flats,trash_counts_flats]

#%%
mflat_maker_for_counts(counts_sep_flats,MFLAT_counts_path)

# reading in master flat counts files from Master Flats folder
MFLAT_counts_imgs = ImageFileCollection(MFLAT_counts_path, keywords='*')
MFLAT_counts_files = MFLAT_counts_imgs.files_filtered(FIELD   = '              flat',
                                        include_path=True)
MFLAT_counts_ccds = MFLAT_counts_imgs.ccds(FIELD   = '              flat', 
                             combined=True)
MFLAT_counts_chips_files = chip_separator(MFLAT_counts_files)


#================================ don't touch ================================#
#%%
###############################################################################
#--------------------SECTION THREE: IMAGE CALIBRATION-------------------------# 
###############################################################################

# reducing ALERT data
ALERT_reducer(target_names_dict,reduced_ALERT_path,MDARK_chips_files,
              MFLAT_counts_chips_files,MDARK_imgs,combined_darks)

###############################################################################
#-------------------------------END OF CODE-----------------------------------# 
############################################################################### 

t.toc() # Print elapsed time
    
###############################################################################
#-------------------------------END OF CODE-----------------------------------# 
############################################################################### 