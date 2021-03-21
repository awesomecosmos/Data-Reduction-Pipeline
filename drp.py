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
    idx = np.argmin(np.abs(dark_exposures - image.header['EXPTIME']))
    closest_dark_exposure = dark_exposures[idx]

    if (tolerance is not None and 
        np.abs(image.header['EXPTIME'] - closest_dark_exposure) > tolerance):
        
        raise RuntimeError('Closest dark exposure time is {} for flat of exposure '
                           'time {}.'.format(closest_dark_exposure, a_flat.header['EXPTIME']))
        
    
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
    hdu1 = fits.open(img)
    chip_num = hdu1[0].header['CHIP']
    return chip_num
    #return abs(int(img[-6:-4]))


def chip_separator(IMAGElist):
    """ This function creates an array of arrays containing the filenames for 
    each chip.
    
    Input parameter(s): 
    * IMAGElist - list of filenames.
      dtype: list
      
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
    for i in range(len(exptimes)):
        exptime_files = []
        for j in IMAGElist:
             hdu1 = fits.open(j)
             exptime_of_file = hdu1[0].header['EXPTIME']
             # checking if each file matches the exposure length
             if exptime_of_file == exptimes[i]:
                 exptime_files.append(j)
        exptimes_lst.append(exptime_files)    
    return exptimes_lst

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
BIAS_chips_files = chip_separator(BIAS_files)

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

    # # plotting single bias compared to combined bias
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    # show_image(BIAS_ccd[0], cmap='gray', ax=ax1, fig=fig, percl=90)
    # ax1.set_title('Single calibrated bias for Chip {}'.format(chip_num))
    # show_image(master_bias.data, cmap='gray', ax=ax2, fig=fig, percl=90)
    # ax2.set_title('{} bias images combined for Chip {}'.format(num_of_biases,chip_num))

#%%
##-----------------------------CALIBRATING DARKS-----------------------------##

# reading in dark files from DARK folder
DARK_path = Path('Obs Data/DARK')
DARK_imgs = ImageFileCollection(DARK_path,glob_exclude=['/*-0.fit','/*-99.fit'])

# making/checking MBIAS path/folder
DARK_cal_path = path_checker(DARK_path,'Calibrated Darks')

DARK_files = DARK_imgs.files_filtered(FIELD='              dark',include_path=True)
DARK_chips_files = chip_separator(DARK_files)

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
        ccd.meta['calibrated'] = 'True'
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
# excluding non-science quicklook images and biases
DARK_cal_imgs = ImageFileCollection(DARK_cal_path,
                                    glob_exclude=['/*-0.fit','/*-99.fit','/*-1-*.fit'])
DARK_cal_files = DARK_cal_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
                                              include_path=True)
DARK_cal_chips_files = chip_separator(DARK_cal_files)

# making/checking MDARK path/folder
MDARK_path = path_checker(DARK_path,'Master Darks')

for DARK_chips in DARK_cal_chips_files:
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


#%%
# ##--------------------------------FLATS-----------------------------------##

FLAT_path = Path('Obs Data/FLAT')
# FLAT_path = Path('Obs Data_old/FLATS')
to_include = ['/*-1.fit','/*-2.fit','/*-3.fit','/*-4.fit','/*-5.fit',
              '/*-6.fit','/*-7.fit','/*-8.fit','/*-9.fit','/*-10.fit']
good_files = []
for i in to_include:
    good_file = glob.glob('Obs Data/FLAT' + i)
    # good_file = glob.glob('Obs Data_old/FLATS' + i)
    good_files += good_file

FLAT_imgs = ImageFileCollection(filenames=good_files)

# FLAT_files = FLAT_imgs.files_filtered(FIELD ='        flat_round',
#                                  include_path=True)
FLAT_files = FLAT_imgs.files_filtered(FIELD ='              flat',
                                  include_path=True)
# making/checking MFLAT path/folder
FLAT_cal_path = path_checker(FLAT_path,'Calibrated Flats')
MFLAT_path = path_checker(FLAT_path,'Master Flats')
FLAT_exptimes = exptime_checker(FLAT_files)
FLAT_chips_files = chip_separator(FLAT_files)

# FLAT_ccds = FLAT_imgs.ccds(FIELD ='        flat_round',            
#                            ccd_kwargs={'unit': 'adu'},                            
#                            return_fname=True)

FLAT_ccds = FLAT_imgs.ccds(FIELD ='              flat',            
                            ccd_kwargs={'unit': 'adu'},                            
                            return_fname=True)

#%%
# reading in master dark files from Master Darks folder
# excluding non-science quicklook images and biases
MDARK_imgs = ImageFileCollection(MDARK_path, keywords='*')
MDARK_files = MDARK_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
                                        include_path=True)
MDARK_ccds = MDARK_imgs.ccds(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>', 
                             combined=True)
MDARK_chips_files = chip_separator(MDARK_files)

n_combined_dark = len(MDARK_files)
expected_exposure_times = set(FLAT_exptimes)
actual_exposure_times = set(h['EXPTIME'] for h in MDARK_imgs.headers(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>', combined=True))

if (expected_exposure_times - actual_exposure_times):
    raise RuntimeError('Encountered unexpected exposure time in combined darks. '
                       'The unexpected times are {}'.format(actual_exposure_times - expected_exposure_times))

combined_darks = {ccd.header['EXPTIME']: ccd for ccd in MDARK_ccds}

#%%
for f_index,FLAT_chips in enumerate(FLAT_chips_files):
    # getting master darks for current chip
    MDARK_chips_file = MDARK_chips_files[f_index]
    
    for FLAT_file in FLAT_chips:
        hdu1 = fits.open(FLAT_file)
        f_file_name = hdu1[0].header['RUN'].strip(' ')
        flat_exptime = hdu1[0].header['EXPTIME']
        f_obs_set = hdu1[0].header['SET'].strip(' ')
        f_chip_num = hdu1[0].header['CHIP']
        
        FLAT_ccd = CCDData(FLAT_file,unit=u.adu)
        
        for mdark in MDARK_chips_file:
            mdark_hdu1 = fits.open(mdark)
            mdark_exptime = mdark_hdu1[0].header['EXPTIME']
            if mdark_exptime == flat_exptime:
                MDARK_exptime = mdark_exptime
                MDARK_to_subtract = CCDData(mdark,unit=u.adu)
        
        # closest_dark = find_nearest_dark_exposure(FLAT_ccd, actual_exposure_times)

        # MDARK_exptime_u = u.Quantity(MDARK_exptime,u.second)
        # flat_exptime_u = u.Quantity(flat_exptime,u.second)
        MDARK_exptime_u = MDARK_exptime*u.second
        flat_exptime_u = flat_exptime*u.second
        # MDARK_to_subtract.header['EXPTIME'] = MDARK_exptime_u
        
        # header_exptime = ccdp.Keyword('EXPTIME', unit=u.second)
        
        calibrated_flat = ccdp.subtract_dark(ccd=FLAT_ccd, 
                              master=MDARK_to_subtract,
                              dark_exposure=MDARK_exptime_u,
                              data_exposure=flat_exptime_u,
                              # exposure_time='EXPTIME',
                              exposure_unit=u.second,
                              scale=False)
        
        # calibrated_flat = ccdp.flat_correct(MDARK_to_subtract,FLAT_ccd)
                              

#%%

for flat_ccd, file_name in FLAT_ccds:   
    
    # extracting header data for display/saving purposes later
    hdu1 = fits.open(flat_ccd)
    file_name = hdu1[0].header['RUN'].strip(' ')
    exptime = hdu1[0].header['EXPTIME']
    obs_set = hdu1[0].header['SET'].strip(' ')
    chip_num = hdu1[0].header['CHIP']
    
    # Find the correct dark exposure
    closest_dark = find_nearest_dark_exposure(flat_ccd, actual_exposure_times)
    
    # Subtract the dark current 
    calibrated_flat = ccdp.subtract_dark(flat_ccd, 
                             combined_darks[closest_dark],
                             exposure_time='exptime', 
                             exposure_unit=u.second)

    # Save the result; there are some duplicate file names so pre-pend "flat"
    calibrated_flat.write(FLAT_cal_path / 
              "calibrated_flat-{}-{}-{}-{}.fit".format(file_name,exptime,
                                                       obs_set,chip_num),
                                                       overwrite=True)


#%%

for FLAT_chips in FLAT_chips_files:
    # seperating list of files by their exposure lengths
    exptime_seperated_files = exptime_separator(FLAT_chips)
    
    # for each array of files for each exposure length:
    for exptime_seperated_exps in exptime_seperated_files:
        # extracting header information for this set of files
        hdu1 = fits.open(exptime_seperated_exps[0])
        exptime = hdu1[0].header['EXPTIME']
        chip_num = hdu1[0].header['CHIP']
        
        # getting CCD image for plotting purposes
        flat_fits = fits.getdata(exptime_seperated_exps[0]) 
        flat_ccd = CCDData(flat_fits,unit=u.adu) 

        # combining all the darks of this set together
        master_flat = ccdp.subtract_dark(a_flat_reduced, 
                                         combined_darks[closest_dark], 
                                         exposure_time='exptime', 
                                         exposure_unit=u.second, 
                                         scale=False)
        
        # writing keywords to header
        master_flat.meta['combined'] = True
        master_flat.meta['EXPTIME'] = exptime
        master_flat.meta['CHIP'] = chip_num
        
        # writing combined dark as a fits file
        master_flat.write(MFLAT_path / 'mflat-{}-chip{}.fit'.format(exptime,
                                                                  chip_num),
                                                              overwrite=True)




# for FLAT_chips in FLAT_chips_files:
#     # seperating list of files by their exposure lengths
#     exptime_seperated_files = exptime_separator(FLAT_chips)
    
#     # for each array of files for each exposure length:
#     for exptime_seperated_exps in exptime_seperated_files:
#         # extracting header information for this set of files
#         hdu1 = fits.open(exptime_seperated_exps[0])
#         exptime = hdu1[0].header['EXPTIME']
#         chip_num = hdu1[0].header['CHIP']
        
#         # getting CCD image for plotting purposes
#         flat_fits = fits.getdata(exptime_seperated_exps[0]) 
#         flat_ccd = CCDData(flat_fits,unit=u.adu) 

#         # combining all the darks of this set together
#         master_flat = ccdp.combine(exptime_seperated_exps,unit=u.adu,
#                               method='average',
#                               sigma_clip=True, 
#                               sigma_clip_low_thresh=5, 
#                               sigma_clip_high_thresh=5,
#                               sigma_clip_func=np.ma.median, 
#                               sigma_clip_dev_func=mad_std,
#                               mem_limit=350e6)
        
#         # writing keywords to header
#         master_flat.meta['combined'] = True
#         master_flat.meta['EXPTIME'] = exptime
#         master_flat.meta['CHIP'] = chip_num
        
#         # writing combined dark as a fits file
#         master_flat.write(MFLAT_path / 'mflat-{}-chip{}.fit'.format(exptime,
#                                                                   chip_num),
#                                                              overwrite=True)

#         # plotting single dark compared to combined dark
#         fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
#         #plotting single dark for current chip
#         show_image(flat_ccd, cmap='gray', ax=ax1, fig=fig, percl=90)
#         ax1.set_title('Single calibrated flat for Chip {}'.format(chip_num))
#         # plotting combined dark for current chip
#         show_image(master_flat.data, cmap='gray', ax=ax2, fig=fig, percl=90)
#         ax2.set_title('{}s Master Flat for Chip {}'.format(exptime,chip_num))



































