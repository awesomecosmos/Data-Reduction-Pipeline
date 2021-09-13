# -*- coding: utf-8 -*-

# This file contains all necessary packages and functions required for the 
# Astrometric Calibration of the images. This file is required to run the file 
# asp.py.
# Import it as: 
# from asp_funcs import *

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
# import astropy
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

def asp(lst_of_files,cal_log_path,WCS_cal_path,t,ast):
    """
    This function deals with solving for the astrometry for images.
    
    Parameters
    ----------
    lst_of_files : list
        List of filenames for images.
    
    cal_log_path : WindowsPath object
        Path to directory where calibration log is to be saved.
    
    WCS_cal_path : WindowsPath object
        Path to directory where astrometrically-calibrated image is to be saved.
    
    t : TicToc object
        Tictoc instance, for calibration log purposes.
        
    ast : AstrometryNet object
        AstrometryNet instance, for solving purposes.
    
    Returns
    -------
    Nothing.
    """
    for reduced_ALERT_file in lst_of_files:
        
        try_again = True
        submission_id = None
        
        while try_again:
            try:
                if not submission_id:
                    output_msg = "ACCEPTED for {}!".format(reduced_ALERT_file)
                    print(output_msg)
                    
                    # writing Calibration Log
                    calibration_log = open(cal_log_path,"a")
                    calibration_log.write(output_msg+"\n")
                    spam = t.tocvalue()
                    calibration_log.write(str(spam)+"\n")
                    calibration_log.close()
                    
                    # solving the image
                    wcs_header = ast.solve_from_image(reduced_ALERT_file,
                                                      solve_timeout=1000,
                                                      submission_id=submission_id,
                                                      force_image_upload=True)
                                                      # ra_key="RA      ",
                                                      # dec_key="DEC     ")
        
                else:
                    output_msg = "FAIL for {}!".format(reduced_ALERT_file)
                    print(output_msg)
                    
                    # writing Calibration Log
                    calibration_log = open(cal_log_path,"a")
                    calibration_log.write(output_msg+"\n")
                    spam = t.tocvalue()
                    calibration_log.write(str(spam)+"\n")
                    calibration_log.close()
                    
                    # monitors submission 
                    wcs_header = ast.monitor_submission(submission_id,
                                                        solve_timeout=10000)
            except TimeoutError as e:
                output_msg = "FAIL for {}!".format(reduced_ALERT_file)
                print(output_msg)
                
                # writing Calibration Log
                calibration_log = open(cal_log_path,"a")
                calibration_log.write(output_msg+"\n")
                spam = t.tocvalue()
                calibration_log.write(str(spam)+"\n")
                calibration_log.close()
                
                # sets submission ID
                submission_id = e.args[1]
            else:
                # got a result, so terminate
                try_again = False
        
        if wcs_header:
            output_msg = "SUCCESS for {}!".format(reduced_ALERT_file)
            print(output_msg)
            
            # writing Calibration Log
            calibration_log = open(cal_log_path,"a")
            calibration_log.write(output_msg+"\n")
            spam = t.tocvalue()
            calibration_log.write(str(spam)+"\n")
            calibration_log.close()
            
            # writing the new header to the original file
            ccd_obj = wcs_writer(wcs_header, reduced_ALERT_file, WCS_cal_path)
        else:
            # Code to execute when solve fails
            output_msg = "FAIL for {}! Something has gone wrong.".format(reduced_ALERT_file)
            print(output_msg)
            
            # writing Calibration Log
            calibration_log = open(cal_log_path,"a")
            calibration_log.write(output_msg+"\n")
            spam = t.tocvalue()
            calibration_log.write(str(spam)+"\n")
            calibration_log.close()

# ###############################################################################
# #-------------------------------END OF CODE-----------------------------------#
# ###############################################################################