# -*- coding: utf-8 -*-

# The aim of this code is to use astrometry.net to solve from a list of images.

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
from astropy.wcs import WCS
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
from asp_funcs import *

#%%
###############################################################################
#---------------------SECTION ONE: INITIALISATION-----------------------------#
###############################################################################

# reading in reduced ALERT files from Reduced ALERTS folder
ALERT_path = Path("//spcsfs/ave41/astro/ave41/ObsData_v6/ALERT")

"//spcsfs/ave41/astro/ave41/ObsData-2021-02-16/ALERT"


reduced_ALERT_path = path_checker(ALERT_path,'Reduced ALERT')
reduced_ALERT_imgs = ImageFileCollection(reduced_ALERT_path, keywords='*')
reduced_ALERT_files = reduced_ALERT_imgs.files_filtered(REDUCED = True,
                                        include_path=True)

# defining paths/reading in Test Data files
ASP_test_data_path = Path("//spcsfs/ave41/astro/ave41/ASP_TestData_v1/")
WCS_test_cal_path = path_checker(ASP_test_data_path,'WCS Calibrated')

# THIS IMAGE SOLVES SUCCESSFULLY
single_test_img = "//spcsfs//ave41//astro//ave41//ASP_TestData_v1//reduced-C2021_A6-A4213-60-R-a-3.fit"

# THE NEXT THREE IMAGES DO NOT SOLVE SUCCESSFULLY
another_test_img = "//spcsfs//ave41//astro//ave41//ASP_TestData_v1//reduced-C2021_A7-A4266-300-R-a-3.fit"
another_test_img2 = "//spcsfs//ave41//astro//ave41//ASP_TestData_v1//reduced-C2021_A7-A4265-300-R-a-3.fit"
another_test_img3 = "//spcsfs//ave41//astro//ave41//ASP_TestData_v1//reduced-C2021_A7-A4266-300-R-a-3.fit"

# testdata = [str(ASP_test_data_path) +"/" + n for n in os.listdir(ASP_test_data_path) if (n.endswith('fit') and n.__contains__('-60-')) ]
testdata = [str(ASP_test_data_path) +"/" + n for n in os.listdir(ASP_test_data_path) if (n.endswith('fit'))]





# actual_data = [str(reduced_ALERT_path) +"/" + n for n in os.listdir(reduced_ALERT_path) if (n.endswith('fit') and n.__contains__('-300-') and n.__contains__('-3'))]







# new_testdata = [single_test_img,single_test_img]

# test_lst = []
# test_lst.append(another_test_img2)
# test_lst.append(another_test_img)
# test_lst2 = reduced_ALERT_files[:3]

###############################################################################
#--------------------------SECTION TWO: ASTROMETRY----------------------------#
###############################################################################
#%%
# creating an AstrometryNet instance
ast = AstrometryNet()

# insert your API key here
# for more information, see:
# https://nova.astrometry.net/api_help
ast.api_key = "kbhqokfxlzyezitf" 

#%%
#--------------------------------------
# Testing for 1 single image
#
test_wcs_header = ast.solve_from_image(single_test_img,solve_timeout=1000,
                                   force_image_upload=True)
#%%
# wcs_writer is a user-defined function to write the WCS header to the original image.
test_ccd_obj = wcs_writer(test_wcs_header, single_test_img, WCS_test_cal_path)

#%%
#--------------------------------------
# Loop to read in a file from a list of files, and solve that file.

for reduced_ALERT_file in testdata:
    
    try_again = True
    submission_id = None
    
    while try_again:
        try:
            if not submission_id:
                print("ACCEPTED for {}!".format(reduced_ALERT_file))
                wcs_header = ast.solve_from_image(reduced_ALERT_file,
                                              solve_timeout=1000,
                                              submission_id=submission_id,
                                              force_image_upload=True)
                                              # ra_key="RA      ",
                                              # dec_key="DEC     ")
    
            else:
                print("FAIL for {}!".format(reduced_ALERT_file))
                wcs_header = ast.monitor_submission(submission_id,
                                                    solve_timeout=10000)
        except TimeoutError as e:
            print("FAIL for {}!".format(reduced_ALERT_file))
            submission_id = e.args[1]
        else:
            # got a result, so terminate
            try_again = False
    
    if wcs_header:
        print("SUCCESS for {}!".format(reduced_ALERT_file))
        ccd_obj = wcs_writer(wcs_header, reduced_ALERT_file, WCS_test_cal_path)
    else:
        # Code to execute when solve fails
        print("FAIL")

#%%    
#================================ don't touch ================================#

###############################################################################
#-------------------------------END OF CODE-----------------------------------#
###############################################################################

t.toc() # Print elapsed time

###############################################################################
#-------------------------------END OF CODE-----------------------------------#
###############################################################################