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

def img_stats(img_list):
    for img in img_list:
        image_data = fits.getdata(img)
        
        # extracting data from header for display purposes
        hdu1 = fits.open(img)
        file_name = hdu1[0].header['RUN'].strip(' ')
        exptime = hdu1[0].header['EXPTIME']
        obs_set = hdu1[0].header['SET'].strip(' ')
        chip_num = hdu1[0].header['CHIP']
        
        img_name = '{}-{}-{}-{}.fit'.format(file_name,exptime,obs_set,chip_num)

        # # getting statistical data
        # img_min = np.min(image_data)
        # img_max = np.max(image_data)
        # img_mean = np.mean(image_data)
        # # img_std = np.std(image_data)
        
        # # setting number of bins 
        # NBINS = 100
        
        # plotting figure
        plt.figure()
        sns.set_theme(style="whitegrid")
        ax = sns.violinplot(x=image_data.flatten(),color="mediumorchid")
        ax.set_title('Distribution of counts of {}'.format(img_name))
        ax.set_xlabel('Counts')
        plt.savefig("violin_{}.jpg".format(img_name),dpi=900)
        plt.show()
        
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
BIAS_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/ObsData_v3/DARK")
BIAS_imgs = ImageFileCollection(BIAS_path,glob_exclude=['/*-0.fit','/*-99.fit'])

# selecting images
BIAS_files = BIAS_imgs.files_filtered(EXPTIME=1,include_path=True)

img_stats(BIAS_files)



# exptimes = []
# for value in target_names_dict.values():
#     exptimes.append(exptime_checker(value))










