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
ALERT_path = "C:\\Users\\ave41\\OneDrive - University of Canterbury\\Master's 2021\\" \
       "ASTR480 Research\\ASTR480 Code\\Data Reduction Pipeline\\ObsData_v4\\ALERT"
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
BIAS_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/"
                 "Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/"
                 "ObsData_v4/DARK")
# making/checking MBIAS path/folder
MBIAS_path = path_checker(BIAS_path,'Master Biases')
# reading in dark files from DARK folder
DARK_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/"
                 "ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/"
                 "ObsData_v4/DARK")
# making/checking Calibrated Darks path/folder
DARK_cal_path = path_checker(DARK_path,'Calibrated Darks')
# making/checking MDARK path/folder
MDARK_path = path_checker(DARK_path,'Master Darks')

FLAT_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/"
                 "ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/"
                 "ObsData_v4/FLAT")

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
BIAS_files = BIAS_imgs.files_filtered(EXPTIME=0,include_path=True) # EXPTIME is either 0 or 1
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
    good_file = glob.glob("C:\\Users\\ave41\\OneDrive - University of Canterbury\\Master's 2021\\ASTR480 Research\\ASTR480 Code\\Data Reduction Pipeline\\ObsData_v4\\DARK" 
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
# MDARK_ccds = MDARK_imgs.ccds(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>', 
#                              combined=True)
#%%
MDARK_ccds = MDARK_imgs.ccds(COMBINED=True)
MDARK_chips_files = chip_separator(MDARK_files)

MDARK_counts = img_counts(MDARK_files)

# ##--------------------------------FLATS-----------------------------------##
# selecting images and excluding non-science images
good_files = []
for i in to_include:
    good_file = glob.glob("C:\\Users\\ave41\\OneDrive - University of Canterbury\\Master's 2021\\ASTR480 Research\\ASTR480 Code\\Data Reduction Pipeline\\ObsData_v4\\FLAT" 
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
# .headers(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
#                                                                     combined=True)
# SUBBIAS = 'ccd=<CCDData>, master=<CCDData>'
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
def flats_count_classifier(flats_lst):
    hi_counts_flats = ['hi_counts',[]]
    ok_counts_flats = ['ok_counts',[]]
    lo_counts_flats = ['lo_counts',[]]
    last_resort_flats = ['last_resort',[]]
    trash_counts_flats = ['trash_counts',[]]
    # to_inspect_flats = []
    
    for flat_img in flats_lst:
        flat_data = fits.getdata(flat_img)
        # getting counts for flat
        img_mean = np.mean(flat_data)
        
        if img_mean <= 55000 and img_mean > 40000:
            hi_counts_flats[1].append(flat_img)
        elif img_mean <= 40000 and img_mean > 30000:
            ok_counts_flats[1].append(flat_img)
        elif img_mean <= 30000 and img_mean >= 20000:
            lo_counts_flats[1].append(flat_img)
        elif img_mean <= 63000 and img_mean > 55000:
            last_resort_flats[1].append(flat_img)
        elif img_mean < 20000 and img_mean >= 10000:
            last_resort_flats[1].append(flat_img)
        else: 
            trash_counts_flats[1].append(flat_img)
    
    return hi_counts_flats,ok_counts_flats,lo_counts_flats,last_resort_flats,trash_counts_flats
#%%
hi_counts_flats,ok_counts_flats,lo_counts_flats,last_resort_flats,trash_counts_flats=flats_count_classifier(MFLAT_files)

counts_sep_flats = [hi_counts_flats,ok_counts_flats,lo_counts_flats,last_resort_flats,trash_counts_flats]
#%%
def mflat_maker_for_counts(counts_sep_flats,MFLAT_counts_path):
    """
    This function deals witH [ADD STUFF HERE ADND FIX THIS ENTIRE DOCTSTING]
    
    Parameters
    ----------
    cal_flat_chip_sep_files : list
        List of list of filenames of calibrated flats for each chip.
    
    MFLAT_path : WindowsPath object
        Path to directory where Master Flats are to be saved.
    
    Returns
    -------
    Nothing.
    """
    # for each category of counts
    for counts_sep_flats_categories in counts_sep_flats:
        this_category = counts_sep_flats_categories[0]
        counts_sep_flats_category = counts_sep_flats_categories[1]
        
        if len(counts_sep_flats_category) == 0:
            pass
        else:
            # seperating list of flats in this category by chip num
            chip_seperated_files = chip_separator(counts_sep_flats_category)
            # for each array of files for each chip length:
            for chip_files in chip_seperated_files:
                if len(chip_files) == 0:
                    pass
                else:
                    exptimes = []
                    for chip_file in chip_files:
                        # extracting header information for this set of files
                        hdu1 = fits.open(chip_file)
                        exptime = hdu1[0].header['EXPTIME']
                        exptimes.append(exptime)
                        a_hdu1 = fits.open(chip_file)
                        chip_num = a_hdu1[0].header['CHIP']
                
                    # combining all the flats of this set together
                    master_flat = ccdp.combine(chip_files,unit=u.adu,
                                              method='average',
                                              sigma_clip=True, 
                                              sigma_clip_low_thresh=5, 
                                              sigma_clip_high_thresh=5,
                                              sigma_clip_func=np.ma.median, 
                                              sigma_clip_dev_func=mad_std,
                                              mem_limit=350e6)
                        
                    # writing keywords to header
                    master_flat.meta['combined'] = True
                    # master_flat.meta['all_exps'] = exptimes
                    master_flat.meta['CHIP'] = chip_num
                        
                    # writing combined dark as a fits file
                    master_flat.write(MFLAT_counts_path / 
                                          'mflat-{}-chip{}.fit'.format(this_category,chip_num),
                                                                         overwrite=True)

#%%
mflat_maker_for_counts(counts_sep_flats,MFLAT_counts_path)

# reading in master flat counts files from Master Flats folder
MFLAT_counts_imgs = ImageFileCollection(MFLAT_counts_path, keywords='*')
MFLAT_counts_files = MFLAT_counts_imgs.files_filtered(FIELD   = '              flat',
                                        include_path=True)
MFLAT_counts_ccds = MFLAT_counts_imgs.ccds(FIELD   = '              flat', 
                             combined=True)
MFLAT_counts_chips_files = chip_separator(MFLAT_counts_files)

#%%
# print(MBIAS_counts)
# MFLAT_counts_counts = img_counts(MFLAT_counts_files,plots=True)
print('BIAS_counts',BIAS_counts)
print(" ")
print('DARK_counts',DARK_counts)
print(" ")
print('DARK_cal_counts',DARK_cal_counts)
print(" ")
print('MDARK_counts',MDARK_counts)
print(" ")
# print('FLAT_counts',FLAT_counts)
# print(" ")
# print('FLAT_cal_counts',FLAT_cal_counts)
# print(" ")
# print('MFLAT_counts',MFLAT_counts)
# print(" ")
# print('MFLAT_counts_files',MFLAT_counts_counts)

#%%
###############################################################################
#--------------------SECTION THREE: IMAGE CALIBRATION-------------------------# 
###############################################################################

# target_names_dict
# exptimes_lst

# input items for func: ALERT_path,MDARK_chip_sep_files,MFLAT_chip_sep_files
# input items values  : ALERT_path,MDARK_chips_files,MFLAT_chips_files

MDARK_imgs = ImageFileCollection(MDARK_path, keywords='*')
MDARK_ccds = MDARK_imgs.ccds(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>', 
                             combined=True)
combined_darks = {ccd.header['EXPTIME']: ccd for ccd in MDARK_ccds}

MFLAT_imgs = ImageFileCollection(MFLAT_path, keywords='*')
MFLAT_ccds = MFLAT_imgs.ccds(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>', 
                             combined=True)
combined_flats = {ccd.header['EXPTIME']: ccd for ccd in MFLAT_ccds}

#%%
def ALERT_reducer(target_names_dict,reduced_ALERT_path,MDARK_chip_sep_files,MFLAT_chip_sep_files,
                  MDARK_imgs,MFLAT_imgs,combined_darks,combined_flats):
    """
    target_names_dict
    
    reduced_ALERT_path
    
    MDARK_chip_sep_files: MDARK_chips_files
    
    MFLAT_chip_sep_files: MFLAT_chips_files
    
    MDARK_imgs
    
    MFLAT_imgs
    
    combined_darks
    
    combined_flats
    """
    for key, value in target_names_dict.items():
        ALERT_chips_files = chip_separator(value)
        
        for a_index,ALERT_chips in enumerate(ALERT_chips_files):
            # getting master darks  and flats for current chip
            MDARK_chips_file = MDARK_chip_sep_files[a_index]
            MFLAT_chips_file = MFLAT_chip_sep_files[a_index]
            
            ALERT_exptimes = exptime_checker(ALERT_chips)
            
            # finding closest dark exposure times to ALERT exposure times
            d_actual_exposure_times = set(h['EXPTIME'] for h in MDARK_imgs.headers(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
                                                                         combined=True)) 
            for ALERT_file in ALERT_chips:
                al_hdu1 = fits.open(ALERT_file)
                # extracting header data for later saving purposes
                al_file_name = al_hdu1[0].header['RUN'].strip(' ')
                ALERT_exptime = al_hdu1[0].header['EXPTIME']
                al_obs_set = al_hdu1[0].header['SET'].strip(' ')
                al_chip_num = al_hdu1[0].header['CHIP']
                al_filter = al_hdu1[0].header['COLOUR'].strip(' ')
                    
                # making CCDData object for ALERT which we are calibrating
                ALERT_ccd = CCDData.read(ALERT_file,unit=u.adu)
                    
                for mdark in MDARK_chips_file:
                    # finding master dark of matching exp
                    mdark_hdu1 = fits.open(mdark)
                    mdark_ccd = CCDData.read(mdark,unit=u.adu)
                    mdark_exptime = mdark_hdu1[0].header['EXPTIME']
                        
                    # here, we are finding an exact match for the ALERT
                    if mdark_exptime == ALERT_exptime:
                        MDARK_exptime = mdark_exptime
                        MDARK_to_subtract = CCDData.read(mdark,unit=u.adu)
                        print("MDARK_exptime",MDARK_exptime)
                            
                    else:
                        # Find the correct dark exposure
                        MDARK_exptime = find_nearest_dark_exposure(mdark_ccd,d_actual_exposure_times)
                        MDARK_to_subtract = combined_darks[MDARK_exptime]
                        print("MDARK_exptime",MDARK_exptime)
                    
                for mflat in MFLAT_chips_file:
                    # finding master dark of matching exp
                    mflat_hdu1 = fits.open(mflat)
                    mflat_ccd = CCDData.read(mflat,unit=u.adu)
                    mflat_exptime = mflat_hdu1[0].header['EXPTIME']
                        
                    MFLAT_exptime = mflat_exptime
                    MFLAT_to_divide = CCDData(mflat_ccd,unit=u.adu)
                    print("MFLAT_exptime",MFLAT_exptime)
                        
                    # NEED TO CHOOSE A BETTER METHOD FOR CHOOSING THE BEST FLAT
                    # EXP LENGTH!!!!!!!!! BELOW METHOD WILL NOT WORK!
                    # WILL LIKELY NEED TO DO SOMETHING WITH COUNTS!
            
                reduced_ALERT = ccdp.ccd_process(ALERT_ccd,
                                             dark_frame=MDARK_to_subtract,
                                             master_flat=MFLAT_to_divide,
                                             data_exposure=ALERT_exptime*u.second,
                                             dark_exposure=MDARK_exptime*u.second)
                # writing keywords to header
                reduced_ALERT.meta['reduced'] = True
                reduced_ALERT.meta['EXPTIME'] = ALERT_exptime
                reduced_ALERT.meta['CHIP'] = al_chip_num
                reduced_ALERT.meta['RUN'] = al_file_name
                reduced_ALERT.meta['SET'] = al_obs_set
                reduced_ALERT.meta['COLOUR'] = al_filter
                
                # Save the result
                reduced_ALERT.write(reduced_ALERT_path / 
                              "reduced-{}-{}-{}-{}-{}-{}.fit".format(key.strip(' '),
                                                                     al_file_name,
                                                                     ALERT_exptime,
                                                                     al_filter,
                                                                     al_obs_set,
                                                                     al_chip_num),
                                                                     overwrite=True)
        

#%%
ALERT_reducer(target_names_dict,reduced_ALERT_path,MDARK_chips_files,MFLAT_chips_files,
              MDARK_imgs,MFLAT_imgs,combined_darks,combined_flats)

####
# median_count = [np.median(data) for data in FLAT_imgs.data()]
# mean_count = [np.mean(data) for data in FLAT_imgs.data()]
# #%%

# # make this so that it generates 1 plot per chip, so a total of 10 plots
# plt.plot(median_count[20:30], label='median',color="darkviolet")
# plt.plot(mean_count[20:30], label='mean',color="palevioletred")
# plt.xlabel('Image number')
# plt.ylabel('Count (ADU)')
# plt.title('Pixel value in calibrated flat frames')
# plt.legend()
# print(median_count)

#%%
#================================ don't touch ================================#

###############################################################################
#-------------------------------END OF CODE-----------------------------------# 
############################################################################### 

t.toc() # Print elapsed time
    
###############################################################################
#-------------------------------END OF CODE-----------------------------------# 
############################################################################### 