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
date_of_calibration = datetime.today().strftime('%Y-%m-%d')
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++CHANGES++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cal_log_date = "ObsData-2021-02-16"

# reading in reduced ALERT files from Reduced ALERTS folder
reduced_ALERT_path = Path("//spcsfs/ave41/astro/ave41/ObsData-2021-02-12/ALERT/Reduced ALERT")

# filtering data files to only choose 300-s exposures for Chip 3
data_to_calibrate = [str(reduced_ALERT_path) +"/" + n 
                     for n in os.listdir(reduced_ALERT_path) if (n.endswith('fit') 
                     and n.__contains__('-300-') and n.__contains__('-3'))]

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++CHANGES++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

test_cal_log_date = "ASP_TestData_v1"

# defining paths/reading in Test Data files
ASP_test_data_path = Path("//spcsfs/ave41/astro/ave41/ASP_TestData_v1/")
WCS_test_cal_path = path_checker(ASP_test_data_path,'WCS Calibrated')
outputs_path = path_checker(WCS_test_cal_path,'Outputs')

testdata = [str(ASP_test_data_path) +"/" + n for n in os.listdir(ASP_test_data_path) if (n.endswith('fit'))]

# creating Calibration Log
log_filename = "asp_log-{}.txt".format(test_cal_log_date)
log_path = Path(str(outputs_path) + "\\" + log_filename)
calibration_log = open(log_path,"w")

calibration_log.write("Calibration Log for Astrometrical Calibration for {}".format(test_cal_log_date)+"\n")
calibration_log.write("Date of calibration: {}".format(date_of_calibration)+"\n")
calibration_log.write("Reduced ALERTS to calibrate astrometrically:"+"\n")
calibration_log.write(str(testdata)+"\n")

spam = t.tocvalue()
calibration_log.write(str(spam)+"\n")
calibration_log.close()

#------------------------------------------------------------------------------
# # defining paths to store outputs of ASP 
# WCS_cal_path = path_checker(reduced_ALERT_path,'WCS Calibrated')
# outputs_path = path_checker(WCS_cal_path,'Outputs')

# # creating Calibration Log
# log_filename = "asp_log-{}.txt".format(cal_log_date)
# log_path = Path(str(outputs_path) + "\\" + log_filename)
# calibration_log = open(log_path,"w")

# calibration_log.write("Calibration Log for Astrometrical Calibration for {}".format(cal_log_date)+"\n")
# calibration_log.write("Reduced ALERTS to calibrate astrometrically:"+"\n")
# calibration_log.write(str(data_to_calibrate)+"\n")

# spam = t.tocvalue()
# calibration_log.write(str(spam)+"\n")
# calibration_log.close()

###############################################################################
#--------------------------SECTION TWO: ASTROMETRY----------------------------#
###############################################################################

# creating an AstrometryNet instance
ast = AstrometryNet()

# insert your API key here
# for more information, see:
# https://nova.astrometry.net/api_help
ast.api_key = "kbhqokfxlzyezitf" 

#%%
#--------------------------------------
# Loop to read in a file from a list of files, and solve that file.

asp(testdata,log_path,WCS_test_cal_path,t,ast)

#%%    
#================================ don't touch ================================#

###############################################################################
#-------------------------------END OF CODE-----------------------------------#
###############################################################################

t.toc() # Print elapsed time

###############################################################################
#-------------------------------END OF CODE-----------------------------------#
###############################################################################