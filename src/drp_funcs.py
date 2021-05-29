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
from astropy.visualization import hist
from astropy.utils.data import get_pkg_data_filename

# ccdproc packages
import ccdproc as ccdp
from ccdproc import Combiner
from ccdproc import ImageFileCollection
from ccdproc.utils.sample_directory import sample_directory_with_files

# Seaborn packages 
import seaborn as sns

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
    
    # a_flat = CCDData.read(image[0], unit='adu')

    if (tolerance is not None and 
        np.abs(image.header['EXPTIME'] - closest_dark_exposure) > tolerance):
        
        raise RuntimeError('Closest dark exposure time is {} for flat of exposure '
                           'time {}.'.format(closest_dark_exposure, image.header['EXPTIME']))
    
    return closest_dark_exposure

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
def img_stats(img_list):
    """
    This function provides image count statistics for FITS files in an image 
    list. It uses a violin plot to display the distribution of the counts.
    
     Parameters
    ----------
    img_list : list
        List of filenames of FITS files.
    
    Returns
    -------
    Violin plots, saved.
    
    """
    for img in img_list:
        image_data = fits.getdata(img)
        
        # extracting data from header for display purposes
        hdu1 = fits.open(img)
        file_name = hdu1[0].header['RUN'].strip(' ')
        exptime = hdu1[0].header['EXPTIME']
        obs_set = hdu1[0].header['SET'].strip(' ')
        chip_num = hdu1[0].header['CHIP']
        
        img_name = '{}-{}-{}-{}.fit'.format(file_name,exptime,obs_set,chip_num)
        
        # plotting figure
        plt.figure()
        sns.set_theme(style="whitegrid")
        ax = sns.violinplot(x=image_data.flatten(),color="mediumorchid")
        ax.set_title('Distribution of counts of {}'.format(img_name))
        ax.set_xlabel('Counts')
        plt.savefig("violin_{}.jpg".format(img_name),dpi=900)
        plt.show()
        
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
        

def mbias_maker2(bias_chip_sep_files,MBIAS_path):
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
    bias_list_to_return = []
    for BIAS_chips in bias_chip_sep_files:
        # seperating list of files by their exposure lengths
        exptime_seperated_files = exptime_separator(BIAS_chips)
        bias_list_to_return.append(exptime_seperated_files)
        # for each array of files for each exposure length:
        for exptime_seperated_exps in exptime_seperated_files:
            # extracting header information for this set of files
            hdu1 = fits.open(exptime_seperated_exps[0])
            exptime = hdu1[0].header['EXPTIME']
            chip_num = hdu1[0].header['CHIP']
            

            # combining all the biases together
            master_bias = ccdp.combine(exptime_seperated_exps,unit=u.adu,
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
            master_bias.meta['EXPTIME'] = exptime
        
            # writing combined bias as a fits file
            master_bias.write(MBIAS_path / 'mbias-{}-chip{}.fit'.format(exptime,chip_num),
                                                                          overwrite=True)
        

        # # getting CCD image for plotting purposes
        # bias_fits = fits.getdata(exptime_seperated_exps[0]) 
        # bias_ccd = CCDData(bias_fits,unit=u.adu) 

        # # plotting single bias compared to combined bias
        # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        # show_image(BIAS_ccd[0], cmap='gray', ax=ax1, fig=fig, percl=90)
        # ax1.set_title('Single calibrated bias for Chip {}'.format(chip_num))
        # show_image(master_bias.data, cmap='gray', ax=ax2, fig=fig, percl=90)
        # ax2.set_title('{} bias images combined for Chip {}'.format(num_of_biases,chip_num))
        
    return bias_list_to_return #to add to calibration log


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
        mbias_for_chip = MBIAS_chip_sep_files[d_index]#[0]
        # converting each DARKS file to a fits array
        mbias_fits = fits.getdata(mbias_for_chip[0]) 
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
            cal_dark = ccdp.subtract_bias(DARK_ccd,mbias_ccd)
            # setting combined dark header
            cal_dark.meta['calibrated'] = 'True'
            cal_dark.meta['RUN'] = d_file_name
            cal_dark.meta['EXPTIME'] = dark_exptime
            cal_dark.meta['SET'] = d_obs_set
            cal_dark.meta['CHIP'] = d_chip_num
            # Save the result
            cal_dark.write(DARK_cal_path 
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

#%%
def flat_calibrator(flat_chip_sep_files,MDARK_chip_sep_files,FLAT_cal_path,
                    actual_exposure_times,combined_darks): 
    """
    This function deals with calibrating the flats for each chip. It calibrates
    each flat (which has been chip-separated) for each exposure using the master
    darks, and can also show image comparisons of a single flat vs a calibrated 
    flat, if that block of code is uncommented.
    
    Parameters
    ----------
    flat_chip_sep_files : list
        List of list of filenames of flats for each chip.
        
    MDARK_chip_sep_files : list
        List of list of master darks for each chip.
    
    FLAT_cal_path : WindowsPath object
        Path to directory where Calibrated Flats are to be saved.
    
    actual_exposure_times : set
        Set of exposure times found.
    
    combined_dark : dict
        Dictionary of darks for each exposure time.
    
    Returns
    -------
    Nothing.
    """
    for f_index,FLAT_chips in enumerate(flat_chip_sep_files):
            # getting master darks for current chip
            MDARK_chips_file = MDARK_chip_sep_files[f_index]
            
            for FLAT_file in FLAT_chips:
                hdu1 = fits.open(FLAT_file)
                # extracting header data for later saving purposes
                f_file_name = hdu1[0].header['RUN'].strip(' ')
                flat_exptime = hdu1[0].header['EXPTIME']
                f_obs_set = hdu1[0].header['SET'].strip(' ')
                f_chip_num = hdu1[0].header['CHIP']
                
                # making CCDData object for flat which we are calibrating
                FLAT_ccd = CCDData.read(FLAT_file,unit=u.adu)
                
                for mdark in MDARK_chips_file:
                    # finding master dark of matching exp
                    mdark_hdu1 = fits.open(mdark)
                    mdark_ccd = CCDData.read(mdark,unit=u.adu)
                    mdark_exptime = mdark_hdu1[0].header['EXPTIME']
                    
                    # here, we are finding an exact match for the flat
                    if mdark_exptime == flat_exptime:
                        MDARK_exptime = mdark_exptime
                        MDARK_to_subtract = CCDData.read(mdark,unit=u.adu)
                        
                    else:
                        # Find the correct dark exposure
                        MDARK_exptime = find_nearest_dark_exposure(mdark_ccd, actual_exposure_times)
                        MDARK_to_subtract = combined_darks[MDARK_exptime]
    
                MDARK_exptime_u = MDARK_exptime*u.second #produces a Quantity object
                flat_exptime_u = flat_exptime*u.second   #produces a Quantity object
                
                cal_flat = ccdp.subtract_dark(ccd=FLAT_ccd, #CCD array of flat
                                                     master=MDARK_to_subtract, #CCD array of master dark
                                                     dark_exposure=MDARK_exptime_u,
                                                     data_exposure=flat_exptime_u,
                                                     scale=False)
                # Save the result
                cal_flat.write(FLAT_cal_path / 
                      "calibrated_flat-{}-{}-{}-{}.fit".format(f_file_name,flat_exptime,
                                                                f_obs_set,f_chip_num),
                                                                overwrite=True) 

def mflat_maker(cal_flat_chip_sep_files,MFLAT_path):
    """
    This function deals with making a master flat for each chip. It creates
    FITS master flat files for each chip for each exposure time, and can also 
    show image comparisons of a single flat vs a master flat, if that block of 
    code is uncommented.
    
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

    for FLAT_chips in cal_flat_chip_sep_files:
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
            master_flat = ccdp.combine(exptime_seperated_exps,unit=u.adu,
                                  method='average',
                                  sigma_clip=True, 
                                  sigma_clip_low_thresh=5, 
                                  sigma_clip_high_thresh=5,
                                  sigma_clip_func=np.ma.median, 
                                  sigma_clip_dev_func=mad_std,
                                  mem_limit=350e6)
            
            # writing keywords to header
            master_flat.meta['combined'] = True
            master_flat.meta['EXPTIME'] = exptime
            master_flat.meta['CHIP'] = chip_num
            
            # writing combined dark as a fits file
            master_flat.write(MFLAT_path / 'mflat-{}-chip{}.fit'.format(exptime,
                                                                      chip_num),
                                                                  overwrite=True)
            
            
            
def img_counts(img_list,plots=False):
    avg_counts = []
    # avg_counts_dict=dict()
    
    for img in img_list:
        image_data = fits.getdata(img)
        
        # extracting data from header for display purposes
        hdu1 = fits.open(img)
        # file_name = hdu1[0].header['RUN'].strip(' ')
        # exptime = hdu1[0].header['EXPTIME']
        # obs_set = hdu1[0].header['SET'].strip(' ')
        chip_num = hdu1[0].header['CHIP']
        img_name = 'chip{}.fit'.format(chip_num)

        # getting statistical data
        img_min = np.min(image_data)
        img_max = np.max(image_data)
        img_mean = np.mean(image_data)
        # img_std = np.std(image_data)
        
        avg_counts.append('Chip'+str(chip_num)+':'+str(img_mean.round(2)))
        
        # setting number of bins 
        NBINS = 100
        
        if plots==True:
            plt.figure()
            plt.hist(image_data.flatten(),bins=NBINS,label='counts')
            plt.axvline(x=img_min,linestyle='--',label='min {}'.format(img_min),alpha=0.5)
            plt.axvline(x=img_max,linestyle='--',label='max {}'.format(img_max),alpha=0.5)
            plt.axvline(x=img_mean,linestyle='-',linewidth=0.5,color='b',label='mean {:.2f}'.format(img_mean),alpha=1)
            plt.legend()
            plt.grid()
            plt.xlabel('Count level in image')
            plt.ylabel('Number of pixels with that count')
            plt.title('Histogram of counts of {}'.format(img_name))
            plt.savefig("hist_{}.jpg".format(img_name))
            plt.show() 
        
    return avg_counts