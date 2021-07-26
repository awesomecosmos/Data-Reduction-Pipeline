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

# initialising starting directory
code_home_path = "C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/01 Data Reduction Pipeline/DataReductionPipeline/src"
os.chdir(code_home_path) #from now on, we are in this directory

# importing functions
from drp_funcs import *

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

#%%

###############################################################################
#---------------------SECTION ONE: INITIALISATION-----------------------------#
###############################################################################

ALERT_path = "//spcsfs/ave41/astro/ave41/ObsData_v6/ALERT"
os.chdir(ALERT_path) #from now on, we are in this directory
reduced_ALERT_path = path_checker(ALERT_path,'Reduced ALERT')

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
hdul = fits.open(single_test_img)
hdr1 = hdul[0].header
# extract stuff from header to write in filename later
my_ra_key = str(hdul[0].header['RA-NOW  '])
my_dec_key = str(hdul[0].header['DEC-NOW '])
# my_ra_key = hdul[0].header['RA-NOW ']
# my_dec_key = hdul[0].header['DEC-NOW']
radius = 1 #degree
# my_ra_dec_units = (u.arcsecond,u.arcsecond)
my_ra_dec_units = ('second', 'arcsecond')
#%%
test_wcs_header = ast.solve_from_image(single_test_img,solve_timeout=1000,force_image_upload=False)
                                       # parity=2,
                                       # ra_key='RA-NOW  ',dec_key='DEC-NOW ')
                                       # ra_dec_units=my_ra_dec_units)

#%%
copied_img = "reduced-C2021_A6-A4213-60-R-a-3 - Copy.fit"
# og_header = hdul.copy()
# og_header.tofile(fileobj='new_file.fits', overwrite=False)
test_wcs_header.tofile(fileobj=copied_img, overwrite=True)

#%%

new_hdu = fits.ImageHDU(test_wcs_header)

og_header = hdul[0].header
primary_hdu = fits.PrimaryHDU(header=og_header)

new_hdul = fits.HDUList([primary_hdu, test_wcs_header])

#%%

new_hdul = fits.HDUList()
new_hdul.append(fits.PrimaryHDU(header=hdr1))
new_hdul.append(fits.ImageHDU(data=None, header=test_wcs_header, name='wcs_hdr'))

new_hdul.writeto('new_test.fits', clobber=True)

#%%

#hdul.append(fits.PrimaryHDU(header=hdr1))
hdul.append(fits.ImageHDU(data=None, header=test_wcs_header, name='wcs_hdr'))

hdul.writeto('new_test.fits', clobber=True)

#%%
from astropy.coordinates import SkyCoord
my_ra_key = hdul[0].header['RA-NOW ']
my_dec_key = hdul[0].header['DEC-NOW']
my_ra_dec_units = (u.arcsecond,u.arcsecond)

yeet = SkyCoord(my_ra_key, my_dec_key, unit=my_ra_dec_units)
print(yeet)

#fwhm=3,

#%%
wcs_headers_lst = []
for reduced_ALERT_file in reduced_ALERT_files:
    wcs_header = ast.solve_from_image(reduced_ALERT_file,solve_timeout=10000,force_image_upload=False)
    #wcs_headers_lst.append(wcs_header)


#================================ don't touch ================================#

###############################################################################
#-------------------------------END OF CODE-----------------------------------#
###############################################################################

t.toc() # Print elapsed time

###############################################################################
#-------------------------------END OF CODE-----------------------------------#
###############################################################################