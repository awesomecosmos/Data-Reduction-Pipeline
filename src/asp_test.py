# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 18:53:06 2021

@author: ave41
"""

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

test_cal_log_filename = "ASP_TestData_v1"

#%%
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++CHANGES++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#------------------------------- TEST DATA-------------------------------------

# defining paths/reading in Test Data files
ASP_test_data_path = Path("//spcsfs/ave41/astro/ave41/ASP_TestData_v1")
WCS_test_cal_path = path_checker(ASP_test_data_path,'WCS Calibrated')
outputs_path = path_checker(WCS_test_cal_path,'Outputs')

testdata = [str(ASP_test_data_path) +"/" + n for n in os.listdir(ASP_test_data_path) if (n.endswith('fit'))]

# creating Calibration Log
log_filename = "asp_log-{}.txt".format(test_cal_log_filename)
log_path = Path(str(outputs_path) + "\\" + log_filename)
calibration_log = open(log_path,"w")
calibration_log.write("Calibration Log for Astrometrical Calibration for {}".format(test_cal_log_filename)+"\n")
calibration_log.write("Date of calibration: {}".format(date_of_calibration)+"\n")
calibration_log.write("Reduced ALERTS to calibrate astrometrically:"+"\n")
calibration_log.write(str(testdata)+"\n")
spam = t.tocvalue()
calibration_log.write(str(spam)+"\n")
calibration_log.close()
#%%
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
asp(data_to_calibrate,log_path,WCS_cal_path,t,ast)

#%%    
#================================ don't touch ================================#

###############################################################################
#-------------------------------END OF CODE-----------------------------------#
###############################################################################

t.toc() # Print elapsed time
# writing Calibration Log
calibration_log = open(log_path,"a")
calibration_log.write("----- end of astrometrical calibration -----"+"\n")
spam = t.tocvalue()
calibration_log.write(str(spam)+"\n")
calibration_log.close()

###############################################################################
#-------------------------------END OF CODE-----------------------------------#
###############################################################################