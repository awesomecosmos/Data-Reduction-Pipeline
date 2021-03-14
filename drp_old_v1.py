# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 14:23:56 2021

@author: Aayushi Verma
"""

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import glob
import os
from astropy.io import fits
from scipy.ndimage import interpolation as interp

from skimage.feature.register_translation import (register_translation, _upsampled_dft)
#import skimage.feature.register_translation 
#import skimage.feature
#import skimage.feature.register_translation.register_translation, skimage.feature.register_translation._upsampled_dft

## This turns off warnings: not a great way to code
## But when we show the images, sometimes we're taking the logarithm of zero and it doesn't like that
## Which would matter if we were doing math, but we're just taking a look at images, so we can ignore it. 
import warnings
warnings.filterwarnings('ignore')

## function to plot an image cube
## for this code a "cube" means a stack of image data arrays

def plot_grid(datacube,imagenames):
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
        
# def path_stripper(IMAGElist):
#     """ Function to strip preceding pathname of file, such that we only have the filename.
    
#     Input parameter(s):
#     * IMAGElist - list of IMAGE paths. IMAGE refers to either: DARKS, FLATS or ALERTS.
#       dtype: list
    
#     Output parameter(s):
#     * filenames_lst - list of stripped filenames, with no preceding path.
#       dtype: list 
#     """
#     filenames_lst = []
#     for file in IMAGElist:
#         filenames_lst.append(os.path.basename(file))
#     return filenames_lst

def nonscienceimg_filtering(path):
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

def datacube_creator(IMAGElist,filetype,raw_image_data):
    """ This function creates a 'datacube' for each image in IMAGElist.
    
    Input parameter(s): 
    * IMAGElist - list of IMAGE paths. IMAGE refers to either: DARKS, FLATS or ALERTS.
      dtype: list
    * filetype - corresponds to type of file. Could be either DARKS, FLATS or ALERT.
      dtype: string
    * raw_image_data - dictionary of path:data array.
      dtype: dictionary
    """
    IMAGEcube = np.stack([raw_image_data[IMAGEframe] for IMAGEframe in IMAGElist],axis=0)
    return IMAGEcube

DARKSlist = imagelist_creator('Obs Data/DARKS')
FLATSlist = imagelist_creator('Obs Data/FLATS')
ALERTlist = imagelist_creator('Obs Data/ALERT')

## now put all the lists together for a masterlist of all the images:
all_images_list = DARKSlist + FLATSlist + ALERTlist
# print(all_images_list)

raw_image_data = {}
for image_name in all_images_list: 
    raw_image_data[image_name] = fits.getdata(image_name)

#for image in all_images_list: print(raw_image_data[image].shape) ##for troubleshooting
## (check to make sure they're all the same size)


DARKScube = datacube_creator(DARKSlist,'DARKS',raw_image_data)
FLATScube = datacube_creator(FLATSlist,'FLATS',raw_image_data)
ALERTcube = datacube_creator(ALERTlist,'ALERT',raw_image_data)

plot_grid(DARKScube,DARKSlist)

master_dark = np.median(DARKScube, axis=0) ## to median combine them instead

plt.figure(figsize=(15,15)) 
plt.imshow(np.log10(master_dark), origin='lower', cmap='gray');
plt.title('Master Dark')



# filneames of flats and science frames that have not yet been bias-subtracted: 
dark_subtracted_list_in = ALERTlist + FLATSlist
#print debiaslist_in ## for troubleshooting

## filenames for the corresponding bias-subtracted images:
dark_subtracted_list_out = ['dark_subtracted_' + img for img in dark_subtracted_list_in]

## subtract the master bias from each of the raw science & flat frames: 
#debias_data_out = {} ## dictionary for the debiased images

dark_subtracted_data_out = {} ## dictionary for the debiased images

for i in range(len(dark_subtracted_list_in)):  
    dark_subtracted_data_out[dark_subtracted_list_out[i]] = raw_image_data[dark_subtracted_list_in[i]] - master_dark

## python note: we're iterating over an integer, through lists, and the lists were defined in the same order. 
## we wouldn't want to iterate through dictionaries this way because dictionaries are unordered. 

## create an array of debiased images
dark_subtracted_cube = np.stack([dark_subtracted_data_out[image] for image in dark_subtracted_list_out],axis=0)

## show the images: 
plot_grid(dark_subtracted_cube,dark_subtracted_list_out)

## this results in "debiased" science & flat frames
## for one of the science frames, the difference is shown here--top is the raw image, 
## bottom is the debiased image

im = 3  ## the image we're looking at. M42 is 0 through 4, the flats are 5 through 8

## these values are the maximum and minumum pixel values mapped to different grays 
## where min gets shown as black and max as white
## all the plots shown are in log base 10, so the actual pixel values aren't really 1.5-3 or so
## but those are the values after we take the logarithm in this case. 
graymin = 1
graymax = 4

plt.figure(1);
plt.figure(figsize=(15,15)) ;
plt.imshow(np.log10(raw_image_data[dark_subtracted_list_in[im]]), origin='lower', cmap='gray', vmin=graymin, vmax=graymax);
plt.title(dark_subtracted_list_in[im])

plt.figure(2);
plt.figure(figsize=(15,15)) ;
plt.imshow(np.log10(dark_subtracted_data_out[dark_subtracted_list_out[im]]), origin='lower', cmap='gray', vmin=graymin, vmax=graymax);
plt.title(dark_subtracted_list_out[im])

## first we need a list of JUST the debiased flat images to work with: 
dark_subtracted_flat_list = ['dark_subtracted_' + image for image in FLATSlist] 

## create an array of debiased v-flat images 
flatcube = np.stack([dark_subtracted_data_out[flat_frame] for flat_frame in dark_subtracted_flat_list],axis=0)

## average the images in the stack
master_flat = np.average(flatcube, axis=0)


## here's the master flat: 
our_vmin = 4000
our_vmax = 8000
plt.figure(figsize=(15,15))
#plt.imshow((master_flat), origin='lower', cmap='gray', vmin=our_vmin, vmax=our_vmax)
plt.imshow((master_flat), origin='lower', cmap='gray')
plt.title('Master Flat')

## it's not shown in the same log scale as everything else so that we can see the flatfield variations better. 

print('master flat median: ' + str(np.median(master_flat)) + " counts" )
print('master flat mean: ' + str(np.mean(master_flat)) + " counts" )
print('master flat max value: ' + str(np.max(master_flat)) + " counts" )
print('master flat min value: ' + str(np.min(master_flat)) + " counts" )


normalized_master_flat = master_flat/np.mean(master_flat)


print('normalized master flat median: ' + str(np.median(normalized_master_flat)) )
print('normalized master flat mean: ' + str(np.mean(normalized_master_flat)) )
print('normalized master flat max value: ' + str(np.max(normalized_master_flat)) )
print('normalized master flat min value: ' + str(np.min(normalized_master_flat)) )

## normalized master flat: 
plt.figure(figsize=(15,15))
plt.imshow((normalized_master_flat), origin='lower', cmap='gray', vmin=.95, vmax=1.1)
plt.title('Normalized Master Flat')

## we'll start with a list of the debiased M42 images: 
debias_m42_list = ['debiased_' + im for im in ALERTlist]
print(len(debias_m42_list))
print(debias_m42_list[0]) ## this line is for troubleshooting

## and we'll make a corresponding list to name the flattened images: 
flat_debias_m42_list = ['flattened_' + im for im in debias_m42_list]
print(flat_debias_m42_list[0]) ## this line is for troubleshooting
print(len(flat_debias_m42_list))
## create an empty dictionary to populate with the completely corrected science frames: 
flat_debias_data_out = {} 

## and populate the dictionary with each corrected image
## where the dictionary keys = the images in flat_debias_m42_list
## we're iterating over an integer here again because the lists match up
for i in range(len(debias_m42_list)): 
    flat_debias_data_out[flat_debias_m42_list[i]] = dark_subtracted_data_out[debias_m42_list[i]] / normalized_master_flat
    
## create an array of corrected M42 images
m42cube = np.stack([flat_debias_data_out[science_frame] for science_frame in flat_debias_m42_list],axis=0)
# print m42cube ## this line is for troubleshooting
# m42cube.shape ## this line is for troubleshooting

## show the images: 
plot_grid(m42cube,ALERTlist)


## so if we compare the initial raw science frame (1) to a debiased (2) and a flatfield & bias-corrected frame (3) 
## we can see the progression in lessening noise

## grayscale mapping min/max values: 
graymin = 1.25
graymax = 3

im = 3  ## the image we're looking at. there are 5 images of M42, indexed 0 through 4

## (1) raw image
im1 = ALERTlist[im]
data_im1 = raw_image_data[im1]

## (2) debiased image
im1_d = debias_m42_list[im]
data_im1_d = dark_subtracted_data_out[im1_d]

## (3) flattened and debiased image
im1_d_f = flat_debias_m42_list[im]
data_im1_d_f = flat_debias_data_out[im1_d_f]

plt.figure(1)
plt.figure(figsize=(15,15))
plt.imshow(np.log10(data_im1), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
plt.title('(1) Raw Science Image: ' + im1)

plt.figure(2)
plt.figure(figsize=(15,15))
plt.imshow(np.log10(data_im1_d), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
plt.title('(2) Debiased Image: ' + im1_d)

plt.figure(3)
plt.figure(figsize=(15,15))
plt.imshow(np.log10(data_im1_d_f), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
plt.title('(3) Flatfielded & Debiased Image: ' + im1_d_f)


## array of images + average combine: 
m42_cube = np.stack(flat_debias_data_out.values(),axis=0)
m42_stacked = np.average(m42_cube, axis=0)

## plotting: 
plt.figure(1);
plt.figure(figsize=(15,15));
plt.title('V-Band M42: Stacked But Not Aligned');
plt.imshow(np.log10(m42_stacked), origin='lower', cmap='gray', vmin=1.5, vmax=3);

## choose an image to define as zero shift:
zero_shift_image = flat_debias_m42_list[3]

## find all shifts for other images: 
imshifts = {} # dictionary to hold the x and y shift pairs for each image
for image in flat_debias_m42_list: 
    ## register_translation is a function that calculates shifts by comparing 2-D arrays
    result, error, diffphase = register_translation(
        flat_debias_data_out[zero_shift_image], 
        flat_debias_data_out[image], 1000)
    imshifts[image] = result
    
# print imshifts ## for troubleshooting

## new list for shifted image names: 
shifted_m42_list = ['shifted_' + im for im in flat_debias_m42_list]

## new dictionary for shifted image data: 
shifted_m42_data = {}
for i in range(len(shifted_m42_list)):
    ## interp.shift is the function doing the heavy lifting here,
    ## it's reinterpolating each array into the new, shifted one
    shifted_m42_data[shifted_m42_list[i]] = interp.shift(
        flat_debias_data_out[flat_debias_m42_list[i]], 
        imshifts[flat_debias_m42_list[i]])
    
## array of aligned arrays: 
m42cube  = np.stack(shifted_m42_data.values(),axis=0)

## average combined final image: 
m42_stacked = np.average(m42cube, axis=0)
# m42_cube.shape

## show the final image array as an image: 
plt.figure(1)
plt.figure(figsize=(15,15));
plt.title('V-Band M42: Aligned and Stacked');
plt.imshow(np.log10(m42_stacked), origin='lower', cmap='gray', vmin=1.5, vmax=3)

## 1st image is raw science frame
## 2nd image has been debiased (science - bias)
## 3rd image has been flatfield-corrected ((science - master_bias)/(normalized_master_flat))
## 4th image has been aligned and stacked (sum(shifted((science - master_bias)/(normalized_master_flat)))

## grayscale mapping min/max values: 
graymin = 1.45
graymax = 3.25

im = 1  ## the image we're looking at. there are 5 images of M42, indexed 0 through 4

## (1) raw image
im1 = ALERTlist[im]
data_im1 = raw_image_data[im1]

## (2) debiased image
im1_d = debias_m42_list[im]
data_im1_d = dark_subtracted_data_out[im1_d]

## (3) flattened and debiased image
im1_d_f = flat_debias_m42_list[im]
data_im1_d_f = flat_debias_data_out[im1_d_f]

plt.figure(1)
plt.figure(figsize=(15,15))
plt.imshow(np.log10(data_im1), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
plt.title('(1) Raw Science Image: ' + im1)

plt.figure(2)
plt.figure(figsize=(15,15))
plt.imshow(np.log10(data_im1_d), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
plt.title('(2) Debiased Image: ' + im1_d)

plt.figure(3)
plt.figure(figsize=(15,15))
plt.imshow(np.log10(data_im1_d_f), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
plt.title('(3) Flatfielded & Debiased Image: ' + im1_d_f)

plt.figure(4)
plt.figure(figsize=(15,15))
plt.imshow(np.log10(m42_stacked), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
plt.title('(4) Final Aligned & Stacked Image')


























