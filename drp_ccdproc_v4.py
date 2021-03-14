# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 13:45:48 2021

@author: ave41
"""

# drp_ccdproc_v2.py serves as a testing sandbox for drp_ccdproc_v1.py.

###############################################################################
#-------------------SECTION ONE: IMPORTING PACKAGES---------------------------# 
###############################################################################

# basic Python packages
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# path-type packages
import os
import glob
from pathlib import Path

# warnings
import warnings
warnings.filterwarnings('ignore')

# Astropy packages
from astropy.io import fits
from astropy import units as u
from astropy.stats import mad_std
from astropy.nddata import CCDData
from astropy.utils.data import get_pkg_data_filename

# ccdproc packages
import ccdproc as ccdp
from ccdproc import Combiner
from ccdproc import ImageFileCollection
from ccdproc.utils.sample_directory import sample_directory_with_files

# user-defined packages
from convenience_functions import show_image

###############################################################################
#-------------------SECTION TWO: HELPER FUNCTIONS-----------------------------# 
###############################################################################

def find_nearest_dark_exposure(image, dark_exposure_times, tolerance=0.5):
    """
    Find the nearest exposure time of a dark frame to the exposure time of the image,
    raising an error if the difference in exposure time is more than tolerance.
    
    Source of this function: 
    https://mwcraig.github.io/ccd-as-book/05-03-Calibrating-the-flats.html
    
    Parameters
    ----------
    
    image : astropy.nddata.CCDData
        Image for which a matching dark is needed.
    
    dark_exposure_times : list
        Exposure times for which there are darks.
    
    tolerance : float or ``None``, optional
        Maximum difference, in seconds, between the image and the closest dark. Set
        to ``None`` to skip the tolerance test.
    
    Returns
    -------
    
    float
        Closest dark exposure time to the image.
    """

    dark_exposures = np.array(list(dark_exposure_times))
    idx = np.argmin(np.abs(dark_exposures - image.header['exptime']))
    closest_dark_exposure = dark_exposures[idx]

    if (tolerance is not None and 
        np.abs(image.header['exptime'] - closest_dark_exposure) > tolerance):
        
        raise RuntimeError('Closest dark exposure time is {} for flat of exposure '
                           'time {}.'.format(closest_dark_exposure, a_flat.header['exptime']))
        
    
    return closest_dark_exposure

#-----------------------------------------------------------------------------
def chip_num_extractor(img):
    """ This function extracts the chip number from the filename.
    
    Input parameter(s): 
    * img - file path.
      dtype: string
      
    Output parameter(s):
    * number of chip
      dtype: int 
    """
    return abs(int(img[-6:-4]))


def chip_separator(IMAGElist, filetype):
    """ This function creates an array of arrays containing the filenames for 
    each chip.
    
    Input parameter(s): 
    * IMAGElist - list of filenames.
      dtype: list
    * filetype - type of image, i.e. ALERT, DARKS or FLATS.
      dtype: string
      
    Output parameter(s):
    * chip_names_lst - array of arrays containing filenames for each chip.
      dtype: list 
    """
    chip_names_lst = []    
    for i in range(1,11):
        chip_names = []
        for j in IMAGElist:
            if chip_num_extractor(j) == i:
                chip_names.append(j)
        chip_names_lst.append(chip_names)
    return chip_names_lst


def path_checker(origin_path,folder_name):
    """
    This function checks whether a subdirectory exists. If it doesn't, it will 
    be created.
    
    Parameters
    ----------
    origin_path : WindowsPath
        Path of origin, i.e. the directory in which the subdirectory is to be 
        created.
    
    folder_name : str
        Name of subidrectory.
    
    Returns
    -------
    new_path: WindowsPath
        Path of subdirectory.
    """
    new_dir = os.path.join(origin_path,folder_name)
    if os.path.isdir(new_dir) == False:
        os.mkdir(new_dir)
    else:
        pass
    new_path = Path(new_dir)
    return new_path

###############################################################################
#-------------------SECTION THREE: DATA REDUCTION-----------------------------# 
###############################################################################
#%%
##---------------------------MAKING MASTER BIASES----------------------------##

# reading in bias files from BIAS folder
BIAS_path = Path('Obs Data/DARK')
BIAS_imgs = ImageFileCollection(BIAS_path,glob_exclude=['/*-0.fit','/*-99.fit'])

# making/checking MBIAS path/folder
MBIAS_path = path_checker(BIAS_path,'Master Biases')

# selecting images
BIAS_files = BIAS_imgs.files_filtered(EXPTIME=1,include_path=True)
BIAS_chips_files = chip_separator(BIAS_files,'BIAS')

for index,BIAS_chips in enumerate(BIAS_chips_files):
    # getting some numbers for display/saving purposes later
    chip_num = chip_num_extractor(BIAS_chips_files[index][0])
    num_of_biases = len(BIAS_chips_files[0])
    
    # converting each BIAS file to a fits array
    BIAS_fits = [fits.getdata(BIAS_file) for BIAS_file in BIAS_chips]
    # converting each BIAS array to a CCDData object
    BIAS_ccd = [CCDData(BIAS_fit,unit=u.adu) for BIAS_fit in BIAS_fits]

    # combining all the biases together
    master_bias = ccdp.combine(BIAS_ccd,unit=u.adu,
                              method='average',
                              sigma_clip=True, 
                              sigma_clip_low_thresh=5, 
                              sigma_clip_high_thresh=5,
                              sigma_clip_func=np.ma.median, 
                              sigma_clip_dev_func=mad_std,
                              mem_limit=350e6)
    # setting combined bias header
    master_bias.meta['combined'] = True
    # writing combined bias as a fits file
    master_bias.write(MBIAS_path / 'mbias_chip{}.fit'.format(chip_num),overwrite=True)

    # plotting single bias compared to combined bias
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    show_image(BIAS_ccd[0], cmap='gray', ax=ax1, fig=fig, percl=90)
    ax1.set_title('Single calibrated bias for Chip {}'.format(chip_num))
    show_image(master_bias.data, cmap='gray', ax=ax2, fig=fig, percl=90)
    ax2.set_title('{} bias images combined for Chip {}'.format(num_of_biases,chip_num))

#%%
##-----------------------------CALIBRATING DARKS-----------------------------##

# reading in dark files from DARK folder
DARK_path = Path('Obs Data/DARK')
DARK_imgs = ImageFileCollection(DARK_path,glob_exclude=['/*-0.fit','/*-99.fit'])

# making/checking MBIAS path/folder
DARK_cal_path = path_checker(DARK_path,'Calibrated Darks')

DARK_files = DARK_imgs.files_filtered(FIELD='              dark',include_path=True)
DARK_chips_files = chip_separator(DARK_files, 'DARKS')

for DARK_chips in DARK_chips_files:
    for DARK_file in DARK_chips:
        # extracting header data for display/saving purposes later
        hdu1 = fits.open(DARK_file)
        file_name = hdu1[0].header['RUN'].strip(' ')
        exptime = hdu1[0].header['EXPTIME']
        obs_set = hdu1[0].header['SET'].strip(' ')
        chip_num = hdu1[0].header['CHIP']
    
        # converting each DARKS file to a fits array
        DARK_fits = fits.getdata(DARK_file) 
        # converting each DARK array to a CCDData object
        DARK_ccd = CCDData(DARK_fits,unit=u.adu)

        # Subtract bias from each Dark
        ccd = ccdp.subtract_bias(DARK_ccd, master_bias)
        # setting combined dark header
        #ccd.meta['calibrated'] = True
        ccd.meta['RUN'] = file_name
        ccd.meta['EXPTIME'] = exptime
        ccd.meta['SET'] = obs_set
        ccd.meta['CHIP'] = chip_num
        # Save the result
        ccd.write(DARK_cal_path 
                  / "calibrated_dark-{}-{}-{}-{}.fit".format(file_name,exptime,
                                                            obs_set,chip_num),
                                                            overwrite=True)

#%%
##----------------------------MAKING MASTER DARKS----------------------------##
# reading in calibrated dark files from Calibrated Darks folder
DARK_cal_imgs = ImageFileCollection(DARK_cal_path,glob_exclude=['/*-0.fit','/*-99.fit'])
DARK_cal_files = DARK_cal_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',include_path=True)
DARK_cal_chips_files = chip_separator(DARK_cal_files,'DARKS')

# making/checking MDARK path/folder
MDARK_path = path_checker(DARK_path,'Master Darks')

for DARK_chips in DARK_cal_chips_files:
    # making/re-setting folder to store CCD darks for each chip
    chip_n_darks_lst = []
    
    # for each file for each chip number:
    for DARK_file in DARK_chips:
        # extracting header data for display/saving purposes later
        hdu1 = fits.open(DARK_file)
        file_name = hdu1[0].header['RUN'].strip(' ')
        exptime = hdu1[0].header['EXPTIME']
        obs_set = hdu1[0].header['SET'].strip(' ')
        chip_num = hdu1[0].header['CHIP']
        
        # converting each DARKS file to a fits array
        DARK_cal_fits = fits.getdata(DARK_file) 
        # converting each DARK array to a CCDData object
        DARK_cal_ccd = CCDData(DARK_cal_fits,unit=u.adu)
        # appending each DARK CCD file to the list (for each chip)
        chip_n_darks_lst.append(DARK_cal_ccd)

    # combining all the darks together
    master_dark = ccdp.combine(chip_n_darks_lst,unit=u.adu,
                              method='average',
                              sigma_clip=True, 
                              sigma_clip_low_thresh=5, 
                              sigma_clip_high_thresh=5,
                              sigma_clip_func=np.ma.median, 
                              sigma_clip_dev_func=mad_std,
                              mem_limit=350e6)
    # setting combined dark header
    master_dark.meta['combined'] = True
    # writing combined dark as a fits file
    master_dark.write(MDARK_path / 'mdark-{}-chip{}.fit'.format(exptime,
                                                                chip_num),
                                                          overwrite=True)
    # writing keywords to header
    ccd.meta['RUN'] = file_name
    ccd.meta['EXPTIME'] = exptime
    ccd.meta['SET'] = obs_set
    ccd.meta['CHIP'] = chip_num

    # plotting single dark compared to combined dark
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        #plotting single dark for current chip
    show_image(chip_n_darks_lst[0], cmap='gray', ax=ax1, fig=fig, percl=90)
    ax1.set_title('Single calibrated dark for Chip {}'.format(chip_num))
        # plotting combined dark for current chip
    show_image(master_dark.data, cmap='gray', ax=ax2, fig=fig, percl=90)
    ax2.set_title('{}s Master Dark for Chip {}'.format(exptime,chip_num))



#%%
# ##--------------------------------FLATS-----------------------------------##
# # (8/3/21) currently this section combines all the chips
# # will need to sort this out
# ##------------------------------##
# # FLATS_path = Path('Obs Data/FLATS')
# # FLATS_imgs = ImageFileCollection(FLATS_path, keywords='*')

# # # (8/3/21) will need to add an earlier section which stores the exptime of 
# # # the science image to put in here, instead of 300
# # FLATS_files = FLATS_imgs.files_filtered(FIELD ='              flat',include_path=True)

# # # converting each DARKS file to a fits array
# # FLATS_fits = [fits.getdata(FLATS_file) for FLATS_file in FLATS_files]
# # # converting each BIAS array to a CCDData object
# # FLATS_ccd = [CCDData(FLATS_fit,unit=u.adu) for FLATS_fit in FLATS_fits]


# pseudocode for first section of FLATS:
# flats_exp_times = []
# for FLAT in FLATS:
    # exptime = exptime(from header)
    # flats_exp_times.append(exptime)

# change dir to DARKS
# use calibrated darks path
# for DARK in DARKS:
    # exptime_darks = exptime(from header)
    # for flat_exp in flats_exp_times:
        # if flat_exp == exptime_darks:
            # save DARK to new folder in FLATS folder called 'Cal Darks for Flats'
        # else:
            # pass
#https://mwcraig.github.io/ccd-as-book/05-03-Calibrating-the-flats.html        
        



#%%

# ifc_reduced = ccdp.ImageFileCollection(DARKS_cal_path)
# # combined_dark_files = ifc_reduced.files_filtered(EXPTIME=300,combined=True)
# # combined_dark_files = ifc_reduced.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',combined=True)
# combined_dark_files = ifc_reduced.files_filtered(SIMPLE=True,combined=True)
# flat_image_type = '              flat'

# #%%

# n_combined_dark = len(combined_dark_files)
# expected_exposure_times = set([5, 10, 20, 30])

# if n_combined_dark < 10:
#     raise RuntimeError('One or more combined dark is missing. Please re-run the dark notebook.')
# elif n_combined_dark > 10:
#     raise RuntimeError('There are more combined dark frames than expected.')
    
# actual_exposure_times = set(h['exptime'] for h in ifc_reduced.headers(imagetyp='dark', combined=True))

# if (expected_exposure_times - actual_exposure_times):
#     raise RuntimeError('Encountered unexpected exposure time in combined darks. '
#                        'The unexpected times are {}'.format(actual_exposure_times - expected_exposure_times))





