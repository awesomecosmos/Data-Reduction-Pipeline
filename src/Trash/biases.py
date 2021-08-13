# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 14:09:10 2021

@author: ave41
"""

# # reading in bias files from BIAS folder
# BIAS_path = Path("//spcsfs/ave41/astro/ave41/ObsData-2021-02-17/DARK")
# BIAS_path = Path("//spcsfs/ave41/astro/ave41/ObsData_v6/DARK")
# BIAS_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/01 Data Reduction Pipeline/ObsData_v3/DARK")
# 
# # making/checking MBIAS path/folder
# MBIAS_path = path_checker(BIAS_path,'Master Biases')

##---------------------------MAKING MASTER BIASES----------------------------##
# # selecting images
# BIAS_imgs = ImageFileCollection(BIAS_path,glob_exclude=['/*-0.fit','/*-99.fit'])
# # make this into a try/else block?
# BIAS_files = BIAS_imgs.files_filtered(EXPTIME=1,include_path=True) # EXPTIME is either 0 or 1

# # separating the biases by chip number
# BIAS_chips_files = chip_separator(BIAS_files)

# # getting the number of counts for each bias
# BIAS_counts = img_counts(BIAS_files)

# # calling mbias function to make master biases for each chip
# mbias_maker(BIAS_chips_files,MBIAS_path)

# # reading in master bias files from Master Biases folder
# MBIAS_imgs = ImageFileCollection(MBIAS_path, keywords='*')
# MBIAS_files = MBIAS_imgs.files_filtered(COMBINED=True,
#                                         include_path=True)

# # separating the master biases by chip number
# MBIAS_chips_files = chip_separator(MBIAS_files)

# # getting the number of counts for each master bias
# MBIAS_counts = img_counts(MBIAS_files)

# spam = t.tocvalue()

# # writing calibration info to calibration log
# calibration_log = open(log_path,"a")
# calibration_log.write(str(spam)+"\n")
# calibration_log.write("Bias Files"+"\n")
# for BIAS_chips_file in BIAS_chips_files:
#     calibration_log.write(str(BIAS_chips_file)+"\n")
# calibration_log.write("Bias Counts"+"\n")
# calibration_log.write(str(BIAS_counts)+"\n")
# calibration_log.write("Master Bias Files"+"\n")
# for MBIAS_chips_file in MBIAS_chips_files:
#     calibration_log.write(str(MBIAS_chips_file)+"\n")
# calibration_log.write("Master Bias Counts"+"\n")
# calibration_log.write(str(MBIAS_counts)+"\n")
# calibration_log.close()

# uncomment the following line if you want image count statistics
# code will take ~6 mins to run
# img_stats(MBIAS_files)

# # reading in calibrated dark files from Calibrated Darks folder
# DARK_cal_imgs = ImageFileCollection(DARK_cal_path)
# DARK_cal_files = DARK_cal_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
#                                               include_path=True)

# # separating the calibrated darks by chip number
# DARK_cal_chips_files = chip_separator(DARK_cal_files)

# # getting the number of counts for each calibrated dark
# DARK_cal_counts = img_counts(DARK_cal_files)

# # calling mdark_maker function to make master darks
# mdark_maker(DARK_cal_chips_files,MDARK_path)

# calling mdark_maker function to make master darks


# MDARK_files = MDARK_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
#                                         include_path=True)


# MDARK_ccds = MDARK_imgs.ccds(BUNIT   = 'adu     ')



# calibration_log.write("Calibrated Dark Files"+"\n")
# for DARK_cal_chips_file in DARK_cal_chips_files:
#     calibration_log.write(str(DARK_cal_chips_file)+"\n")
# calibration_log.write("Calibrated Dark Counts"+"\n")
# calibration_log.write(str(DARK_cal_counts)+"\n")