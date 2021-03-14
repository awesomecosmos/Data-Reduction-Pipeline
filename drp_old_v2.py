# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 10:53:53 2020

@author: Aayushi Verma
"""

###############################################################################
#-------------------SECTION ONE: IMPORTING PACKAGES---------------------------# 
###############################################################################

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import glob
import os

from astropy.io import fits
#from scipy.ndimage import interpolation as interp

#from skimage.feature.register_translation import (register_translation, _upsampled_dft)

import warnings
warnings.filterwarnings('ignore')

from pathlib import Path

from astropy.nddata import CCDData
from ccdproc import ImageFileCollection
import ccdproc as ccdp

 
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


###############################################################################
#-------------------SECTION THREE: DATA REDUCTION-----------------------------# 
###############################################################################
    
## creating a list of files for each type of image
DARKSlist = imagelist_creator('Obs Data_old/DARKS')
FLATSlist = imagelist_creator('Obs Data_old/FLATS')
ALERTlist = imagelist_creator('Obs Data_old/ALERT')

## combining the lists together for a masterlist of all the images
all_images_list = DARKSlist + FLATSlist + ALERTlist

## creating a dictionary for raw image arrays
raw_image_data = {}
for image_name in all_images_list: 
    raw_image_data[image_name] = fits.getdata(image_name)

#----------DARKS----------#
    #write a func tht checks the counts for each dark for ech chip according to that table
    
    
    
## making a master dark for each chip
DARKS_chips = chip_separator(DARKSlist,'DARKS')
mdark_lst = medianing(raw_image_data,DARKS_chips)

for mdark in mdark_lst:
    plt.figure(figsize=(15,15)) 
    plt.imshow(np.log10(mdark), origin='lower', cmap='gray');
    plt.title('Master DARK')

#-----------------
#%%
# filneames of flats and science frames that have not yet been bias-subtracted: 
dark_subtracted_list_in = ALERTlist + FLATSlist
#print debiaslist_in ## for troubleshooting

## filenames for the corresponding bias-subtracted images:
dark_subtracted_list_out = ['dark-subtracted_' + img for img in dark_subtracted_list_in]



## subtract the master bias from each of the raw science & flat frames: 

dark_subtracted_data_out = {} ## dictionary for the debiased images

for i in range(len(dark_subtracted_list_in)):  
    dark_subtracted_data_out[dark_subtracted_list_out[i]] = raw_image_data[dark_subtracted_list_in[i]] - master_bias

## python note: we're iterating over an integer, through lists, and the lists were defined in the same order. 
## we wouldn't want to iterate through dictionaries this way because dictionaries are unordered. 





































