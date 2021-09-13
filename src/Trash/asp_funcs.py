# -*- coding: utf-8 -*-

# This file contains all necessary packages and functions required for the 
# Astrometric Calibration of the images. This file is required to run the file 
# asp.py.
# Import it as: 
# from asp_funcs import *

###############################################################################
#-------------------SECTION ONE: IMPORTING PACKAGES---------------------------#
###############################################################################

# warnings
import warnings
warnings.filterwarnings('ignore')

# Astropy packages
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.nddata import CCDData

#%%
###############################################################################
#--------------------SECTION ONE: HELPER FUNCTIONS----------------------------#
###############################################################################

def wcs_writer(wcs_header, image, WCS_cal_path):
    """
    This function writes the WCS header object information to the original
    FITS file.
    
    Parameters
    ----------
    wcs_header : io.fits.header.Header
        An Astropy Header object returned after solving using astrometry.net.
    
    image : str
        Path of FITS image file to which WCS information is to be appended.
    
    WCS_cal_path : WindowsPath object
        Path to directory where astrometrically-calibrated image is to be saved.
    
    Returns
    -------
    image_ccd_with_wcs : astropy.nddata.ccddata.CCDData
        CCDData object of astrometrically-calibrated image.
    """
    with fits.open(image, "append") as img_hdul:
        img_hdr1 = img_hdul[0].header
        img_ccd = CCDData.read(image,unit=u.adu)
    
        run_filename = img_hdr1['RUN'].strip(' ')
        exptime = img_hdr1['EXPTIME']
        filter_colour = img_hdr1['COLOUR'].strip(' ')
        obs_set = img_hdr1['SET'].strip(' ')
        chip_num = img_hdr1['CHIP']
        
        filename_to_write = "wcs_cal-{}-{}-{}-{}-{}.fit".format(run_filename,exptime,
                                                                filter_colour,obs_set,
                                                                chip_num)
        
        img_ccd_object = CCDData(img_ccd, wcs=WCS(wcs_header))
        img_ccd_object.write(WCS_cal_path/filename_to_write, overwrite=True)
    
    return img_ccd_object

###############################################################################
#------------------SECTION TWO: ASTROMETRY FUNCTIONS-------------------------#
###############################################################################


# ###############################################################################
# #-------------------------------END OF CODE-----------------------------------#
# ###############################################################################