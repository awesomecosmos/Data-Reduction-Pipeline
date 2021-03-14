# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 10:39:13 2021

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


def plot_grid(datacube,imagenames):
    """ function to plot an image cube
        for this code a "cube" means a stack of image data arrays
    """
    no_A = len(datacube) ## number of individual images in the cube
    xplots = int(np.around(np.sqrt(no_A))) ## number of image grid columns
    yplots = xplots + 1 ## number of image grid rows--sometimes there are one fewer, but that's okay

#     print no_A, xplots, yplots ## this line is for troubleshooting
    
    gs = gridspec.GridSpec(yplots, xplots) ## define the image grid
    plt.figure(figsize=(15,15)) ## set the figure size
    for i in range(no_A): 
        ## plots each individual image within the image grid: 
        B = datacube[i]
        plt.subplot(gs[i])
        plt.imshow(np.log10(B), origin='lower', cmap='gray')
        plt.title(imagenames[i])
        
#-----------------------------------------------------------------------------# 

def path_stripper(IMAGElist):
    """ Function to strip preceding pathname of file, such that we only have the filename.
    
    Input parameter(s):
    * IMAGElist - list of IMAGE paths. IMAGE refers to either: DARKS, FLATS or ALERTS.
      dtype: list
    
    Output parameter(s):
    * filenames_lst - list of stripped filenames, with no preceding path.
      dtype: list 
    """
    filenames_lst = []
    for file in IMAGElist:
        filenames_lst.append(os.path.basename(file))
    return filenames_lst


def nonscienceimg_filtering(path): #quicklook
    """ This function filters out the non-science images in a typical MOA dataset.
    These non-science images are 'quick-look' images and end in '*-0.fit'and '*-99.fit'.
    We do not want to include these in our science image dataset.
    The science images are images for each MOA chip and end in '*(1-10).fit'.
    We only want the science images.
    
    Input parameter(s): 
    * path - looks something like this: '.../DARKS'. 
      dtype: string
 
    In lieu of 'DARKS': 'FLATS' and 'ALERT' can also be used.
    This corresponds to MOA file/folder-naming protocol.
    """
    all_data_files = glob.glob(path + '/*.fit')
    science_img_files = np.setdiff1d(all_data_files, glob.glob((path + '/*-0.fit')))
    science_img_files = np.setdiff1d(science_img_files, glob.glob((path + '/*-99.fit')))
    return science_img_files


def imagelist_creator(path):
    """ This function creates an image list for each filetype.
    
    Input parameter(s): 
    * path - looks something like this: '.../DARKS'. 
      dtype: string
    """
    IMAGElist = nonscienceimg_filtering(path).tolist()
    return IMAGElist


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

   
def medianing(raw_image_data,chips):
    """ This function medians the images for each chip and makes a master for each.
    
    Input parameter(s): 
    * raw_image_data - raw image data arrays.
      dtype: dict
    * chips - array of array of chips.
      dtype: list
      
    Output parameter(s):
    * master_img_lst - list of master images.
      dtype: int 
    """
    to_median = []
    master_img_lst = []
    for chip in chips:
        for img in chip:
            img_array = raw_image_data[img]
            to_median.append(img_array)
        master_img = np.median(to_median, axis=0)
        master_img_lst.append(master_img)
    return master_img_lst


#%%

###############################################################################
#-------------------SECTION THREE: DATA REDUCTION-----------------------------# 
###############################################################################

##------------BIASES------------##
# (8/3/21) currently this section combines all the chips
# will need to sort this out
##------------------------------##

# reading in bias files from BIAS folder
# ensure there are no quicklooks (-0 or -99) otherwise will not work
BIAS_path = Path('Obs Data/BIAS')
BIAS_imgs = ImageFileCollection(BIAS_path, keywords='*')
BIAS_files = BIAS_imgs.files_filtered(EXPTIME=1,include_path=True)

# converting each BIAS file to a fits array
BIAS_fits = [fits.getdata(BIAS_file) for BIAS_file in BIAS_files]
# converting each BIAS array to a CCDData object
BIAS_ccd = [CCDData(BIAS_fit,unit=u.adu) for BIAS_fit in BIAS_fits]

# combining all the biases together
combined_bias = ccdp.combine(BIAS_ccd,unit=u.adu,
                              method='average',
                              sigma_clip=True, 
                              sigma_clip_low_thresh=5, 
                              sigma_clip_high_thresh=5,
                              sigma_clip_func=np.ma.median, 
                              sigma_clip_dev_func=mad_std,
                              mem_limit=350e6)
# setting combined bias header
combined_bias.meta['combined'] = True
# writing combined bias as a fits file
combined_bias.write(BIAS_path / 'combined_bias.fit',overwrite=True)

# plotting single bias compared to combined bias
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
show_image(BIAS_ccd[0], cmap='gray', ax=ax1, fig=fig, percl=90)
ax1.set_title('Single calibrated bias')
show_image(combined_bias.data, cmap='gray', ax=ax2, fig=fig, percl=90)
ax2.set_title('{} bias images combined'.format(len(BIAS_files)))

#%%

##------------DARKS------------##
# (8/3/21) currently this section combines all the chips
# will need to sort this out
##------------------------------##

# reading in dark files from DARKS folder
# ensure there are no quicklooks (-0 or -99) otherwise will not work
DARKS_path = Path('Obs Data/DARKS')
DARKS_imgs = ImageFileCollection(DARKS_path, keywords='*')

# (8/3/21) will need to add an earlier section which stores the exptime of 
# the science image to put in here, instead of 300
DARKS_files = DARKS_imgs.files_filtered(EXPTIME=300,include_path=True)

# converting each DARKS file to a fits array
DARKS_fits = [fits.getdata(DARKS_file) for DARKS_file in DARKS_files]
# converting each BIAS array to a CCDData object
DARKS_ccd = [CCDData(DARKS_fit,unit=u.adu) for DARKS_fit in DARKS_fits]

# need to add a block here wich checks for existence of Calib Darks folder
# if it doens't exist, then it makes the dir
# else it just goes ahead and does the path thing
calibrated_darks_dir = os.path.join(DARKS_path,'Calibrated Darks')
os.mkdir(calibrated_darks_dir) 
DARKS_cal_path = Path(calibrated_darks_dir)

for i in range(len(DARKS_ccd)):
    # Subtract bias
    ccd = ccdp.subtract_bias(DARKS_ccd[i], combined_bias)
    # Save the result
    filename = DARKS_files[i][-18:-4]
    ccd.write(DARKS_cal_path / "calibrated_dark_{}.fit".format(filename),overwrite=True)

#%%
# reading in calibrated dark files from DARKS folder
# ensure there are no quicklooks (-0 or -99) otherwise will not work
DARKS_cal_imgs = ImageFileCollection(DARKS_cal_path, keywords='*')

# (8/3/21) will need to add an earlier section which stores the exptime of 
# the science image to put in here, instead of 300
DARKS_cal_files = DARKS_cal_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',include_path=True)

# converting each DARKS file to a fits array
DARKS_cal_fits = [fits.getdata(DARKS_cal_file) for DARKS_cal_file in DARKS_cal_files]
# converting each BIAS array to a CCDData object
DARKS_cal_ccd = [CCDData(DARKS_cal_fit,unit=u.adu) for DARKS_cal_fit in DARKS_cal_fits]

# combining all the biases together
combined_dark = ccdp.combine(DARKS_cal_ccd,unit=u.adu,
                              method='average',
                              sigma_clip=True, 
                              sigma_clip_low_thresh=5, 
                              sigma_clip_high_thresh=5,
                              sigma_clip_func=np.ma.median, 
                              sigma_clip_dev_func=mad_std,
                              mem_limit=350e6)
# setting combined bias header
combined_dark.meta['combined'] = True
# writing combined bias as a fits file
combined_dark.write(DARKS_cal_path / 'combined_dark.fit',overwrite=True)

# plotting single bias compared to combined bias
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
show_image(DARKS_cal_ccd[0], cmap='gray', ax=ax1, fig=fig, percl=90)
ax1.set_title('Single calibrated dark')
show_image(combined_dark.data, cmap='gray', ax=ax2, fig=fig, percl=90)
ax2.set_title('{} dark images combined'.format(len(DARKS_cal_files)))


#%%
##------------FLATS------------##
# (8/3/21) currently this section combines all the chips
# will need to sort this out
##------------------------------##
# FLATS_path = Path('Obs Data/FLATS')
# FLATS_imgs = ImageFileCollection(FLATS_path, keywords='*')

# # (8/3/21) will need to add an earlier section which stores the exptime of 
# # the science image to put in here, instead of 300
# FLATS_files = FLATS_imgs.files_filtered(FIELD ='              flat',include_path=True)

# # converting each DARKS file to a fits array
# FLATS_fits = [fits.getdata(FLATS_file) for FLATS_file in FLATS_files]
# # converting each BIAS array to a CCDData object
# FLATS_ccd = [CCDData(FLATS_fit,unit=u.adu) for FLATS_fit in FLATS_fits]












