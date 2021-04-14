# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 13:45:48 2021

@author: ave41
"""

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
import astropy
from astropy.io import fits
import astropy.units as u
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

# def find_nearest_dark_exposure(image, dark_exposure_times, tolerance=0.5):
#     """
#     Find the nearest exposure time of a dark frame to the exposure time of the image,
#     raising an error if the difference in exposure time is more than tolerance.
    
#     Source of this function: 
#     https://mwcraig.github.io/ccd-as-book/05-03-Calibrating-the-flats.html
    
#     Parameters
#     ----------
    
#     image : astropy.nddata.CCDData
#         Image for which a matching dark is needed.
    
#     dark_exposure_times : list
#         Exposure times for which there are darks.
    
#     tolerance : float or ``None``, optional
#         Maximum difference, in seconds, between the image and the closest dark. Set
#         to ``None`` to skip the tolerance test.
    
#     Returns
#     -------
    
#     float
#         Closest dark exposure time to the image.
#     """

#     dark_exposures = np.array(list(dark_exposure_times))
#     idx = np.argmin(np.abs(dark_exposures - image.header['EXPTIME']))
#     closest_dark_exposure = dark_exposures[idx]

#     if (tolerance is not None and 
#         np.abs(image.header['EXPTIME'] - closest_dark_exposure) > tolerance):
        
#         raise RuntimeError('Closest dark exposure time is {} for flat of exposure '
#                            'time {}.'.format(closest_dark_exposure, a_flat.header['EXPTIME']))
        
    
#     return closest_dark_exposure

#-----------------------------------------------------------------------------
def chip_num_extractor(img):
    """ This function extracts the chip number from the filename.
    
    Parameters
    ----------  
    img : string
        file path.
      
    Returns
    -------
    chip_num : int
        number of chip 
    """
    hdu1 = fits.open(img)
    chip_num = hdu1[0].header['CHIP']
    return chip_num
    #return abs(int(img[-6:-4]))


def chip_separator(IMAGElist):
    """ This function creates an array of arrays containing the filenames for 
    each chip.
    
    Parameters
    ---------- 
    IMAGElist : list
        list of filenames.
      
    Returns
    -------
    chip_names_lst : list 
        array of arrays containing filenames for each chip. 
    """
    chip_names_lst = []    
    for chip in range(1,11):
        # range(1,11) is range of the 10 chips
        chip_names = []
        for image in IMAGElist:
            if chip_num_extractor(image) == chip:
                chip_names.append(image)
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

def exptime_checker(IMAGElist):
    """
    This function creates a list of all exposure times from a list of FITS files.
    
    Parameters
    ----------
    IMAGElist : list
        List of FITS filenames from which the exposure times are to be extracted.
        FITS files must have header keyword 'EXPTIME'.
    
    Returns
    -------
    exptimes: list
        List of exposure times for the IMAGElist.
    """
    exptimes = []
    for image in IMAGElist:
        hdu1 = fits.open(image)
        exptime = hdu1[0].header['EXPTIME']
        if exptime not in exptimes:
            # this condition is to ensure that we don't have a list with 
            # repeating exposure times.
            exptimes.append(exptime)
        else:
            pass
    return exptimes

def exptime_separator(IMAGElist):
    """
    This function creates a array of arrays of FITS files, where each sub-array
    contains the filenames for an exposure length from the function 
    exptime_checker(IMAGElist).
    
    Parameters
    ----------
    IMAGElist : list
        List of FITS filenames from which the exposure times are to be extracted.
        FITS files must have header keyword 'EXPTIME'.
    
    Returns
    -------
    exptimes_lst: list
        List of lists of filenames for each exposure length.
    """
    exptimes_lst = []
    # getting the exposure lengths from the image list
    exptimes = exptime_checker(IMAGElist)
    for exptime in range(len(exptimes)):
        exptime_files = []
        for image in IMAGElist:
             hdu1 = fits.open(image)
             exptime_of_file = hdu1[0].header['EXPTIME']
             # checking if each file matches the exposure length
             if exptime_of_file == exptimes[exptime]:
                 exptime_files.append(image)
        exptimes_lst.append(exptime_files)    
    return exptimes_lst
#%%
###############################################################################
#----------------SECTION THREE: DATA REDUCTION FUNCTIONS----------------------#
###############################################################################

def mbias_maker(bias_chip_sep_files,MBIAS_path):
    """
    This function deals with making a master bias for each chip. It creates
    FITS master bias files for each chip, and can also show image comparisons
    of a single bias vs a master bias, if that block of code is uncommented.
    
    Parameters
    ----------
    bias_chip_sep_files : list
        List of list of filenames of biases for each chip.
    
    MBIAS_path : WindowsPath object
        Path to directory where Master Biases are to be saved.
    
    Returns
    -------
    Nothing.
    """
    for index,BIAS_chips in enumerate(bias_chip_sep_files):
        # getting some numbers for display/saving purposes later
        chip_num = chip_num_extractor(bias_chip_sep_files[index][0])
        num_of_biases = len(bias_chip_sep_files[0])
    
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
        master_bias.meta['CHIP'] = chip_num
        
        # writing combined bias as a fits file
        master_bias.write(MBIAS_path / 'mbias_chip{}.fit'.format(chip_num),overwrite=True)

        # # plotting single bias compared to combined bias
        # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        # show_image(BIAS_ccd[0], cmap='gray', ax=ax1, fig=fig, percl=90)
        # ax1.set_title('Single calibrated bias for Chip {}'.format(chip_num))
        # show_image(master_bias.data, cmap='gray', ax=ax2, fig=fig, percl=90)
        # ax2.set_title('{} bias images combined for Chip {}'.format(num_of_biases,chip_num))


def dark_calibrator(dark_chip_sep_files,MBIAS_chip_sep_files,DARK_cal_path):
    """
    This function deals with calibrating the darks for each chip. It calibrates
    each dark (which has been chip-seprated) for each exposure using the master
    biases, and can also show image comparisons of a single dark vs a calibrated 
    dark, if that block of code is uncommented.
    
    Parameters
    ----------
    dark_chip_sep_files : list
        List of list of filenames of darks for each chip.
        
    MBIAS_chip_sep_files : list
        List list of master biases for each chip.
    
    DARK_cal_path : WindowsPath object
        Path to directory where Calibrated Darks are to be saved.
    
    Returns
    -------
    Nothing.
    """
    for d_index,DARK_chips in enumerate(dark_chip_sep_files):
        # getting master biases for current chip
        mbias_for_chip = MBIAS_chip_sep_files[d_index][0]
        # converting each DARKS file to a fits array
        mbias_fits = fits.getdata(mbias_for_chip) 
        # converting each DARK array to a CCDData object
        mbias_ccd = CCDData(mbias_fits,unit=u.adu)
        
        for DARK_file in DARK_chips:
            # extracting header data for display/saving purposes later
            hdu1 = fits.open(DARK_file)
            d_file_name = hdu1[0].header['RUN'].strip(' ')
            dark_exptime = hdu1[0].header['EXPTIME']
            d_obs_set = hdu1[0].header['SET'].strip(' ')
            d_chip_num = hdu1[0].header['CHIP']
            
            img_name = '{}-{}-{}-{}.fit'.format(d_file_name,dark_exptime,
                                                d_obs_set,d_chip_num)
    
            # converting each DARKS file to a fits array
            DARK_fits = fits.getdata(DARK_file) 
            # converting each DARK array to a CCDData object
            DARK_ccd = CCDData(DARK_fits,unit=u.adu)

            # Subtract bias from each Dark
            ccd = ccdp.subtract_bias(DARK_ccd,mbias_ccd)
            # setting combined dark header
            ccd.meta['calibrated'] = 'True'
            ccd.meta['RUN'] = d_file_name
            ccd.meta['EXPTIME'] = dark_exptime
            ccd.meta['SET'] = d_obs_set
            ccd.meta['CHIP'] = d_chip_num
            # Save the result
            ccd.write(DARK_cal_path 
                  / "calibrated_dark-{}".format(img_name),overwrite=True)
            
            
def mdark_maker(dark_chip_sep_files,MDARK_path):
    """
    This function deals with making a master dark for each chip. It creates
    FITS master dark files for each chip for each exposure time, and can also 
    show image comparisons of a single dark vs a master dark, if that block of 
    code is uncommented.
    
    Parameters
    ----------
    dark_chip_sep_files : list
        List of list of filenames of calibrated darks for each chip.
    
    MDARK_path : WindowsPath object
        Path to directory where Master Darks are to be saved.
    
    Returns
    -------
    Nothing.
    """
    for DARK_chips in dark_chip_sep_files:
        # seperating list of files by their exposure lengths
        exptime_seperated_files = exptime_separator(DARK_chips)
    
        # for each array of files for each exposure length:
        for exptime_seperated_exps in exptime_seperated_files:
            # extracting header information for this set of files
            hdu1 = fits.open(exptime_seperated_exps[0])
            exptime = hdu1[0].header['EXPTIME']
            chip_num = hdu1[0].header['CHIP']
        
            # getting CCD image for plotting purposes
            dark_fits = fits.getdata(exptime_seperated_exps[0]) 
            dark_ccd = CCDData(dark_fits,unit=u.adu) 

            # combining all the darks of this set together
            master_dark = ccdp.combine(exptime_seperated_exps,unit=u.adu,
                                       method='average',
                                       sigma_clip=True, 
                                       sigma_clip_low_thresh=5, 
                                       sigma_clip_high_thresh=5,
                                       sigma_clip_func=np.ma.median, 
                                       sigma_clip_dev_func=mad_std,
                                       mem_limit=350e6)
        
            # writing keywords to header
            master_dark.meta['combined'] = True
            master_dark.meta['EXPTIME'] = exptime
            master_dark.meta['CHIP'] = chip_num
        
            # writing combined dark as a fits file
            master_dark.write(MDARK_path / 'mdark-{}-chip{}.fit'.format(exptime,
                                                                  chip_num),
                                                             overwrite=True)

            # # plotting single dark compared to combined dark
            # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
            # #plotting single dark for current chip
            # show_image(dark_ccd, cmap='gray', ax=ax1, fig=fig, percl=90)
            # ax1.set_title('Single calibrated dark for Chip {}'.format(chip_num))
            # # plotting combined dark for current chip
            # show_image(master_dark.data, cmap='gray', ax=ax2, fig=fig, percl=90)
            # ax2.set_title('{}s Master Dark for Chip {}'.format(exptime,chip_num))

