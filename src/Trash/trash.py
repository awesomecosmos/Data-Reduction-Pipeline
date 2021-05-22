# -*- coding: utf-8 -*-
"""
Created on Thu May 20 13:51:48 2021

@author: ave41
"""
# sys.path.insert(1,"C:\\Users\\ave41\\OneDrive - University of Canterbury\\Master's 2021\\ASTR480 Research\\ASTR480 Code\\Data Reduction Pipeline\\ObsData_v4\\ALERT")


###############################################################################
#-------------------SECTION ONE: IMPORTING PACKAGES---------------------------# 
###############################################################################

# importing drp functions from /src/ folder
import sys
sys.path.insert(1,"C:\\Users\\ave41\\OneDrive - University of Canterbury\\Master's 2021\\ASTR480 Research\\ASTR480 Code\\Data Reduction Pipeline\\DataReductionPipeline\\src")
from drp_funcs import *

# # basic Python packages
# import numpy as np
# import matplotlib
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec



# # Astropy packages
# from astropy.io import fits
# # from astropy import stats
# from astropy.visualization import hist

# # basic Python packages
# import numpy as np
# import matplotlib
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec

# # path-type packages
# import os
# import glob
# from pathlib import Path

# # warnings
# import warnings
# warnings.filterwarnings('ignore')

# # Astropy packages
# from astropy.io import fits
# import astropy.units as u
# from astropy.stats import mad_std
# from astropy.nddata import CCDData
# from astropy.utils.data import get_pkg_data_filename

# # ccdproc packages
# import ccdproc as ccdp
# from ccdproc import Combiner
# from ccdproc import ImageFileCollection
# from ccdproc.utils.sample_directory import sample_directory_with_files

# import seaborn as sns

# # user-defined packages
# from convenience_functions import show_image

#----------------------------------------------

# def img_stats(img_list):
#     for img in img_list:
#         image_data = fits.getdata(img)
        
#         # extracting data from header for display purposes
#         hdu1 = fits.open(img)
#         file_name = hdu1[0].header['RUN'].strip(' ')
#         exptime = hdu1[0].header['EXPTIME']
#         obs_set = hdu1[0].header['SET'].strip(' ')
#         chip_num = hdu1[0].header['CHIP']
        
#         img_name = '{}-{}-{}-{}.fit'.format(file_name,exptime,obs_set,chip_num)

#         # # getting statistical data
#         # img_min = np.min(image_data)
#         # img_max = np.max(image_data)
#         # img_mean = np.mean(image_data)
#         # # img_std = np.std(image_data)
        
#         # # setting number of bins 
#         # NBINS = 100
        
#         # plotting figure
#         plt.figure()
#         sns.set_theme(style="whitegrid")
#         ax = sns.violinplot(x=image_data.flatten(),color="mediumorchid")
#         ax.set_title('Distribution of counts of {}'.format(img_name))
#         ax.set_xlabel('Counts')
#         plt.savefig("violin_{}.jpg".format(img_name),dpi=900)
#         plt.show()
        
        # plt.figure()
        # plt.hist(image_data.flatten(),bins=NBINS,label='counts')
        # plt.axvline(x=img_min,linestyle='--',label='min {}'.format(img_min),alpha=0.5)
        # plt.axvline(x=img_max,linestyle='--',label='max {}'.format(img_max),alpha=0.5)
        # plt.axvline(x=img_mean,linestyle='-',linewidth=0.5,color='b',label='mean {:.2f}'.format(img_mean),alpha=1)
        # plt.legend()
        # plt.grid()
        # plt.xlabel('Count level in image')
        # plt.ylabel('Number of pixels with that count')
        # plt.title('Histogram of counts of {}'.format(img_name))
        # plt.savefig("hist_{}.jpg".format(img_name))
        # plt.show()

# reading in bias files from BIAS folder
# BIAS_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/ObsData_v3/DARK")
# BIAS_imgs = ImageFileCollection(BIAS_path,glob_exclude=['/*-0.fit','/*-99.fit'])

# # selecting images
# BIAS_files = BIAS_imgs.files_filtered(EXPTIME=1,include_path=True)

# img_stats(BIAS_files)



# exptimes = []
# for value in target_names_dict.values():
#     exptimes.append(exptime_checker(value))



#===============F L A T S=====================#

# FLAT_path = Path('Obs Data/FLAT')
# # FLAT_path = Path('Obs Data_old/FLATS')
# to_include = ['/*-1.fit','/*-2.fit','/*-3.fit','/*-4.fit','/*-5.fit',
#               '/*-6.fit','/*-7.fit','/*-8.fit','/*-9.fit','/*-10.fit']
# good_files = []
# for i in to_include:
#     good_file = glob.glob('Obs Data/FLAT' + i)
#     good_files += good_file

# FLAT_imgs = ImageFileCollection(filenames=good_files)

# FLAT_files = FLAT_imgs.files_filtered(FIELD ='              flat',
#                                   include_path=True)
# # making/checking MFLAT path/folder
# FLAT_cal_path = path_checker(FLAT_path,'Calibrated Flats')
# MFLAT_path = path_checker(FLAT_path,'Master Flats')
# FLAT_exptimes = exptime_checker(FLAT_files)
# FLAT_chips_files = chip_separator(FLAT_files)

# FLAT_ccds = FLAT_imgs.ccds(FIELD ='              flat',            
#                             ccd_kwargs={'unit': 'adu'},                            
#                             return_fname=True)


# FLAT_imgs = ImageFileCollection(FLAT_path,glob_exclude=['*-0.fit','*-99.fit'])

# to_include = ['/*-1.fit','/*-2.fit','/*-3.fit','/*-4.fit','/*-5.fit',
#               '/*-6.fit','/*-7.fit','/*-8.fit','/*-9.fit','/*-10.fit']
# to_exclude = ['/*-0.fit','/*-99.fit']
# good_files = []
# for i in to_include:
#     good_file = glob.glob('Obs Data/FLAT' + i)
#     good_files += good_file


# n_combined_dark = len(MDARK_files)
# expected_exposure_times = set(FLAT_exptimes)
# actual_exposure_times = set(h['EXPTIME'] for h in MDARK_imgs.headers(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>', combined=True))

# if (expected_exposure_times - actual_exposure_times):
#     raise RuntimeError('Encountered unexpected exposure time in combined darks. '
#                        'The unexpected times are {}'.format(actual_exposure_times - expected_exposure_times))

# combined_darks = {ccd.header['EXPTIME']: ccd for ccd in MDARK_ccds}


# making/checking MDARK path/folder
# MFLAT_path = path_checker(FLAT_path,'Master Flats')

# calling mdark_maker function to make master darks
# mdark_maker(FLAT_cal_chips_files,MFLAT_path)




            
            
            # # Save the result
            # calibrated_flat.write(FLAT_cal_path / 
            #       "calibrated_flat-{}-{}-{}-{}.fit".format(f_file_name,flat_exptime,
            #                                                 f_obs_set,f_chip_num),
            #                                                 overwrite=True)
        

# #---------------DRAFT CODE---------------------------#
# for f_index,FLAT_chips in enumerate(FLAT_chips_files):
#     # getting master darks for current chip
#     MDARK_chips_file = MDARK_chips_files[f_index]
    
#     for FLAT_file in FLAT_chips:
#         hdu1 = fits.open(FLAT_file)
#         f_file_name = hdu1[0].header['RUN'].strip(' ')
#         flat_exptime = hdu1[0].header['EXPTIME']
#         f_obs_set = hdu1[0].header['SET'].strip(' ')
#         f_chip_num = hdu1[0].header['CHIP']
        
#         flat_fits = fits.getdata(FLAT_file)
#         FLAT_ccd = CCDData(flat_fits,unit=u.adu)
#         fits.info(FLAT_ccd)
        
#         for mdark in MDARK_chips_file:
#             mdark_hdu1 = fits.open(mdark)
#             mdark_exptime = mdark_hdu1[0].header['EXPTIME']
#             if mdark_exptime == flat_exptime:
#                 MDARK_exptime = mdark_exptime
#                 MDARK_to_subtract = CCDData.read(mdark,unit=u.adu)
                
#         # print(FLAT_ccd.unit())
#         # print(MDARK_to_subtract.unit())
        
#         MDARK_exptime_u = MDARK_exptime*u.second #produces a Quantity object
#         flat_exptime_u = flat_exptime*u.second   #produces a Quantity object
        
#         calibrated_flat = ccdp.subtract_dark(ccd=FLAT_ccd, #CCD array of flat
#                                              master=MDARK_to_subtract, #CCD array of master dark
#                                              dark_exposure=MDARK_exptime_u,
#                                              data_exposure=flat_exptime_u,
#                                              # exposure_unit=u.second,
#                                              scale=False)
        
        
#         # Save the result
#         calibrated_flat.write(FLAT_cal_path / 
#               "calibrated_flat-{}-{}-{}-{}.fit".format(f_file_name,flat_exptime,
#                                                         f_obs_set,f_chip_num),
#                                                         overwrite=True)                          

#%%

# for flat_ccd, file_name in FLAT_ccds:   
    
#     # extracting header data for display/saving purposes later
#     hdu1 = fits.open(flat_ccd)
#     file_name = hdu1[0].header['RUN'].strip(' ')
#     exptime = hdu1[0].header['EXPTIME']
#     obs_set = hdu1[0].header['SET'].strip(' ')
#     chip_num = hdu1[0].header['CHIP']
    
#     # Find the correct dark exposure
#     closest_dark = find_nearest_dark_exposure(flat_ccd, actual_exposure_times)
    
#     # Subtract the dark current 
#     calibrated_flat = ccdp.subtract_dark(flat_ccd, 
#                              combined_darks[closest_dark],
#                              exposure_time='exptime', 
#                              exposure_unit=u.second)

#     # Save the result; there are some duplicate file names so pre-pend "flat"
#     calibrated_flat.write(FLAT_cal_path / 
#               "calibrated_flat-{}-{}-{}-{}.fit".format(file_name,exptime,
#                                                         obs_set,chip_num),
#                                                         overwrite=True)


#%%

# for FLAT_chips in FLAT_chips_files:
#     # seperating list of files by their exposure lengths
#     exptime_seperated_files = exptime_separator(FLAT_chips)
    
#     # for each array of files for each exposure length:
#     for exptime_seperated_exps in exptime_seperated_files:
#         # extracting header information for this set of files
#         hdu1 = fits.open(exptime_seperated_exps[0])
#         exptime = hdu1[0].header['EXPTIME']
#         chip_num = hdu1[0].header['CHIP']
        
#         # getting CCD image for plotting purposes
#         flat_fits = fits.getdata(exptime_seperated_exps[0]) 
#         flat_ccd = CCDData(flat_fits,unit=u.adu) 

#         # combining all the darks of this set together
#         master_flat = ccdp.subtract_dark(flat_ccd, 
#                                          combined_darks[closest_dark], 
#                                          exposure_time='exptime', 
#                                          exposure_unit=u.second, 
#                                          scale=False)
        
#         # writing keywords to header
#         master_flat.meta['combined'] = True
#         master_flat.meta['EXPTIME'] = exptime
#         master_flat.meta['CHIP'] = chip_num
        
#         # writing combined dark as a fits file
#         master_flat.write(MFLAT_path / 'mflat-{}-chip{}.fit'.format(exptime,
#                                                                   chip_num),
#                                                               overwrite=True)






#         # plotting single dark compared to combined dark
#         fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
#         #plotting single dark for current chip
#         show_image(flat_ccd, cmap='gray', ax=ax1, fig=fig, percl=90)
#         ax1.set_title('Single calibrated flat for Chip {}'.format(chip_num))
#         # plotting combined dark for current chip
#         show_image(master_flat.data, cmap='gray', ax=ax2, fig=fig, percl=90)
#         ax2.set_title('{}s Master Flat for Chip {}'.format(exptime,chip_num))



# the structure of the code is such that we need to calibrate flats of 
# all exposure lengths (in this case: 5s, 10s, 20s and 30s) for each chip of
# the camera (there are 10 chips).
# FLAT_chips_files is a 10x10 array of 10 flats of different exposures
# per chip.
# this block of code aims to calibrate all flats by subtracting the
# approporiate master dark (matching the exposure length and chip number of 
# the flat).






































