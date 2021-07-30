# -*- coding: utf-8 -*-

###############################################################################
#-------------------SECTION ZERO: IMPORTING PACKAGES--------------------------#
###############################################################################

# initialising timer so we can count elapsed time
from pytictoc import TicToc
t = TicToc() # create TicToc instance
t.tic() # Start timer

# importing packages
import sys
import os

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
from astropy.io.fits import Header
import astropy.units as u
from astropy.stats import mad_std
from astropy.nddata import CCDData
from astropy.visualization import hist
from astropy.utils.data import get_pkg_data_filename

# Astroquery packages 
from astroquery.astrometry_net import AstrometryNet

# Misc packages
import astroalign as aa

# initialising starting directory
code_home_path = "C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/01 Data Reduction Pipeline/DataReductionPipeline/src"
os.chdir(code_home_path) #from now on, we are in this directory

# importing functions
from drp_funcs import *

#%%
###############################################################################
#---------------------SECTION ONE: INITIALISATION-----------------------------#
###############################################################################

ALERT_path = "//spcsfs/ave41/astro/ave41/ObsData_v6/ALERT"
os.chdir(ALERT_path) #from now on, we are in this directory
reduced_ALERT_path = path_checker(ALERT_path,'Reduced ALERT')
# WCS_cal_path = path_checker(reduced_ALERT_path,'WCS Calibrated')
WCS_cal_path = path_checker(code_home_path,'WCS Calibrated')

#//spcsfs/ave41/astro/ave41/ObsData_v6/ALERT/Reduced ALERT/reduced-2020_UB5-A4228-60-R-a-1.fit

# reading in reduced ALERT files from Reduced ALERTS folder
reduced_ALERT_imgs = ImageFileCollection(reduced_ALERT_path, keywords='*')
reduced_ALERT_files = reduced_ALERT_imgs.files_filtered(REDUCED = True,
                                        include_path=True)

# initialising starting directory
os.chdir(code_home_path) #from now on, we are in this directory

###############################################################################
#--------------------------SECTION TWO: ASTROMETRY----------------------------#
###############################################################################
#%%

ast = AstrometryNet()
ast.api_key = "kbhqokfxlzyezitf"

single_test_img = "//spcsfs//ave41//astro//ave41//ASP_TestData_v1//reduced-C2021_A6-A4213-60-R-a-3.fit"
# # single_test_img = "C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/01 Data Reduction Pipeline/DataReductionPipeline/reduced-C2021_A6-A4213-60-R-a-3.fit"
# single_test_img_hdul = fits.open(single_test_img)
# # single_test_img_hdul = fits.open(single_test_img,mode="update")
# single_test_img_hdr1 = single_test_img_hdul[0].header

# run_filename = single_test_img_hdul[0].header['RUN'].strip(' ')
# exptime = single_test_img_hdul[0].header['EXPTIME']
# obs_set = single_test_img_hdul[0].header['SET'].strip(' ')
# chip_num = single_test_img_hdul[0].header['CHIP']
# filter_colour = single_test_img_hdul[0].header['COLOUR'].strip(' ')

# filename_to_write = "wcs_cal-{}-{}-{}-{}-{}.fits".format(run_filename,exptime,
#                                                          filter_colour,obs_set,
#                                                          chip_num)

# #%%
# test_wcs_header = ast.solve_from_image(single_test_img,solve_timeout=1000,
#                                        force_image_upload=False)
# #%%
# single_test_img_hdul.close(single_test_img)

# #%%
# a_single_test_img_hdul = fits.open(single_test_img,mode="update")
# # single_test_img_hdr1 = single_test_img_hdul[0].header

# fits.writeto(single_test_img,data=None,header=test_wcs_header,overwrite=True)

# test_wcs_hdr_items = test_wcs_header.items() #produces a Generator object
# it = iter(test_wcs_hdr_items) #produces a Generator object
# while True:
#      try:
#          my_items = next(it) #produces a Tuple
#          hdr_key = my_items[0] 
#          hdr_val = my_items[1]
#          single_test_img_hdr1.set(hdr_key,hdr_val)
#      except StopIteration:
#          break

#single_test_img_hdul.close()
#%%
# single_test_img_hdul.writeto(single_test_img,overwrite=True)
#%%

# test_keys, test_values = test_wcs_header.items()
# for i in range(len(test_wcs_header)):
#     hdr1.set(test_keys[i],test_values[i])
#%%
# hdul.append(fits.ImageHDU(data=None, header=test_wcs_header, name='wcs_hdr'))
# hdul.writeto(WCS_cal_path/filename_to_write,overwrite=True)

with fits.open(single_test_img, "append") as single_test_img_hdul:
    single_test_img_hdr1 = single_test_img_hdul[0].header
    single_test_img_data = single_test_img_hdul[0].data

    run_filename = single_test_img_hdul[0].header['RUN'].strip(' ')
    exptime = single_test_img_hdul[0].header['EXPTIME']
    obs_set = single_test_img_hdul[0].header['SET'].strip(' ')
    chip_num = single_test_img_hdul[0].header['CHIP']
    filter_colour = single_test_img_hdul[0].header['COLOUR'].strip(' ')
    
    filename_to_write = "wcs_cal-{}-{}-{}-{}-{}.fits".format(run_filename,exptime,
                                                             filter_colour,obs_set,
                                                             chip_num)
    
    #%%
    test_wcs_header = ast.solve_from_image(single_test_img,solve_timeout=1000,
                                           force_image_upload=False)
    #%%
    single_test_img_new = "//spcsfs//ave41//astro//ave41//ASP_TestData_v1//reduced-C2021_A6-A4213-60-R-a-3-NEW.fit"
    test_wcs_hdr_items = test_wcs_header.items() #produces a Generator object
    it = iter(test_wcs_hdr_items) #produces a Generator object
    while True:
          try:
              my_items = next(it) #produces a Tuple
              hdr_key = my_items[0] 
              hdr_val = my_items[1]
              single_test_img_hdr1.set(hdr_key,hdr_val)
          except StopIteration:
              break
    
    fits.writeto(single_test_img_new, data=single_test_img_data, header=single_test_img_hdr1, 
                 output_verify="fix", overwrite=True)

#%%

# #%%

for reduced_ALERT_file in reduced_ALERT_files:
    try_again = True
    submission_id = None
    
    
    
    hdul = fits.open(reduced_ALERT_file)
    hdr1 = hdul[0].header  
    
    run_filename = hdul[0].header['RUN'].strip(' ')
    exptime = hdul[0].header['EXPTIME']
    obs_set = hdul[0].header['SET'].strip(' ')
    chip_num = hdul[0].header['CHIP']
    filter_colour = hdul[0].header['COLOUR'].strip(' ')
    
    filename_to_write = "wcs_cal-{}-{}-{}-{}-{}.fits".format(run_filename,exptime,
                                                              filter_colour,obs_set,
                                                              chip_num)
    
    while try_again:
        try:
            if not submission_id:
                wcs_header = ast.solve_from_image(reduced_ALERT_file,
                                                  submission_id=submission_id,
                                                  solve_timeout=10000,
                                                  force_image_upload=False)
            else:
                wcs_header = ast.monitor_submission(submission_id,
                                                    solve_timeout=100000)
        except TimeoutError as e:
            submission_id = e.args[1]
        else:
            # got a result, so terminate
            try_again = False



#%%
for reduced_ALERT_file in reduced_ALERT_files:
    hdul = fits.open(reduced_ALERT_file)
    hdr1 = hdul[0].header
    
    run_filename = hdul[0].header['RUN'].strip(' ')
    exptime = hdul[0].header['EXPTIME']
    obs_set = hdul[0].header['SET'].strip(' ')
    chip_num = hdul[0].header['CHIP']
    filter_colour = hdul[0].header['COLOUR'].strip(' ')
    
    filename_to_write = "wcs_cal-{}-{}-{}-{}-{}.fits".format(run_filename,exptime,
                                                              filter_colour,obs_set,
                                                              chip_num)
    
    wcs_header = ast.solve_from_image(reduced_ALERT_file,solve_timeout=10000,force_image_upload=False)
    
    hdul.append(fits.ImageHDU(data=None, header=wcs_header, name='wcs_hdr'))
    hdul.writeto(WCS_cal_path/filename_to_write,clobber=True)
    
# # #================================ don't touch ================================#

# # ###############################################################################
# # #-------------------------------END OF CODE-----------------------------------#
# # ###############################################################################

# # t.toc() # Print elapsed time

# # ###############################################################################
# # #-------------------------------END OF CODE-----------------------------------#
# # ###############################################################################