# -*- coding: utf-8 -*-

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
import astropy
from astropy.io import fits
import astropy.units as u
from astropy.stats import mad_std
from astropy.nddata import CCDData
from astropy.visualization import hist
from astropy.utils.data import get_pkg_data_filename

# ccdproc packages
import ccdproc as ccdp
from ccdproc import Combiner
from ccdproc import ImageFileCollection
from ccdproc.utils.sample_directory import sample_directory_with_files

# Seaborn packages
import seaborn as sns

# user-defined packages
from convenience_functions import show_image

###############################################################################
#-------------------SECTION TWO: HELPER FUNCTIONS-----------------------------#
###############################################################################

def find_nearest_dark_exposure(image, dark_exposure_times, tolerance=60):
    """
    Find the nearest exposure time of a dark frame to the exposure time of the image,
    raising an error if the difference in exposure time is more than tolerance.

    Source of this function:
    https://mwcraig.github.io/ccd-as-book/05-03-Calibrating-the-flats.html

    Parameters
    ----------

    image : astropy.nddata.CCDData
        Image for which a matching dark is needed.

    dark_exposure_times : list
        Exposure times for which there are darks.

    tolerance : float or ``None``, optional
        Maximum difference, in seconds, between the image and the closest dark. Set
        to ``None`` to skip the tolerance test.

    Returns
    -------

    float
        Closest dark exposure time to the image.
    """

    dark_exposures = np.array(list(dark_exposure_times))
    idx = np.argmin(np.abs(dark_exposures - image.header['EXPTIME']))
    closest_dark_exposure = dark_exposures[idx]

    if (tolerance is not None and
        np.abs(image.header['EXPTIME'] - closest_dark_exposure) > tolerance):

        raise RuntimeError('Closest dark exposure time is {} for flat of exposure '
                           'time {}.'.format(closest_dark_exposure, image.header['EXPTIME']))

    return closest_dark_exposure

#-----------------------------------------------------------------------------
def chip_num_extractor(img):
    """ This function extracts the chip number from the filename.

    Parameters
    ----------
    img : string
        file path.

    Returns
    -------
    chip_num : int
        number of chip
    """
    hdu1 = fits.open(img)
    chip_num = hdu1[0].header['CHIP']
    return chip_num
    #return abs(int(img[-6:-4]))


def chip_separator(IMAGElist):
    """ This function creates an array of arrays containing the filenames for
    each chip.

    Parameters
    ----------
    IMAGElist : list
        list of filenames.

    Returns
    -------
    chip_names_lst : list
        array of arrays containing filenames for each chip.
    """
    chip_names_lst = []
    for chip in range(1,11):
        # range(1,11) is range of the 10 chips
        chip_names = []
        for image in IMAGElist:
            if chip_num_extractor(image) == chip:
                chip_names.append(image)
        chip_names_lst.append(chip_names)
    return chip_names_lst


def path_checker(origin_path,folder_name):
    """
    This function checks whether a subdirectory exists. If it doesn't, it will
    be created.

    Parameters
    ----------
    origin_path : WindowsPath
        Path of origin, i.e. the directory in which the subdirectory is to be
        created.

    folder_name : str
        Name of subidrectory.

    Returns
    -------
    new_path: WindowsPath
        Path of subdirectory.
    """
    new_dir = os.path.join(origin_path,folder_name)
    if os.path.isdir(new_dir) == False:
        os.mkdir(new_dir)
    else:
        pass
    new_path = Path(new_dir)
    return new_path

def exptime_checker(IMAGElist):
    """
    This function creates a list of all exposure times from a list of FITS files.

    Parameters
    ----------
    IMAGElist : list
        List of FITS filenames from which the exposure times are to be extracted.
        FITS files must have header keyword 'EXPTIME'.

    Returns
    -------
    exptimes: list
        List of exposure times for the IMAGElist.
    """
    exptimes = []
    for image in IMAGElist:
        hdu1 = fits.open(image)
        exptime = hdu1[0].header['EXPTIME']
        if exptime not in exptimes:
            # this condition is to ensure that we don't have a list with
            # repeating exposure times.
            exptimes.append(exptime)
        else:
            pass
    return exptimes

def exptime_separator(IMAGElist):
    """
    This function creates a array of arrays of FITS files, where each sub-array
    contains the filenames for an exposure length from the function
    exptime_checker(IMAGElist).

    Parameters
    ----------
    IMAGElist : list
        List of FITS filenames from which the exposure times are to be extracted.
        FITS files must have header keyword 'EXPTIME'.

    Returns
    -------
    exptimes_lst: list
        List of lists of filenames for each exposure length.
    """
    exptimes_lst = []
    # getting the exposure lengths from the image list
    exptimes = exptime_checker(IMAGElist)
    for exptime in range(len(exptimes)):
        exptime_files = []
        for image in IMAGElist:
             hdu1 = fits.open(image)
             exptime_of_file = hdu1[0].header['EXPTIME']
             # checking if each file matches the exposure length
             if exptime_of_file == exptimes[exptime]:
                 exptime_files.append(image)
        exptimes_lst.append(exptime_files)
    return exptimes_lst

#%%
###############################################################################
#----------------SECTION THREE: STATISTICAL FUNCTIONS-------------------------#
###############################################################################

def img_stats(img_list,plots_path):
    """
    This function provides image count statistics for FITS files in an image
    list. It uses a violin plot to display the distribution of the counts.

     Parameters
    ----------
    img_list : list
        List of filenames of FITS files.
    
    plots_path : WindowsPath object
        Path to directory where Plots are to be saved.

    Returns
    -------
    Violin plots, saved.

    """
    for img in img_list:
        image_data = fits.getdata(img)

        # extracting data from header for display purposes
        hdu1 = fits.open(img)
        file_name = hdu1[0].header['RUN'].strip(' ')
        exptime = hdu1[0].header['EXPTIME']
        obs_set = hdu1[0].header['SET'].strip(' ')
        chip_num = hdu1[0].header['CHIP']

        img_name = '{}-{}-{}-{}.fit'.format(file_name,exptime,obs_set,chip_num)

        # plotting figure
        plt.figure()
        sns.set_theme(style="whitegrid")
        ax = sns.violinplot(x=image_data.flatten(),color="mediumorchid")
        ax.set_title('Distribution of counts of {}'.format(img_name))
        ax.set_xlabel('Counts')
        plt.savefig(plots_path/"violin_{}.jpg".format(img_name),dpi=900)
        plt.show()

def img_counts(img_list,plots_path=None,plots=False):
    """
    This function finds the counts of the images.

    Parameters
    ----------
    img_list : list
        List of images.
    
    plots_path : WindowsPath object
        Path to directory where Plots are to be saved.

    plots : bool, optional
        Boolean (True/False) of plot display option.

    Returns
    -------
    avg_counts_for_display : list
        List of number of average counts.

    """
    avg_counts_for_display = []
    tot_avg_counts = []
    img_exptimes = exptime_separator(img_list)

    for img_exptime in img_exptimes:
        img_exp_chips = chip_separator(img_exptime)
        for img_exp_chip in img_exp_chips:
            avg_counts = []
            for img in img_exp_chip:
                image_data = fits.getdata(img)
                # extracting data from header for display purposes
                hdu1 = fits.open(img)
                exptime = hdu1[0].header['EXPTIME']
                chip_num = hdu1[0].header['CHIP']
                img_name = 'chip{}.fit'.format(chip_num)

                # getting statistical data
                img_min = np.min(image_data)
                img_max = np.max(image_data)
                img_mean = np.mean(image_data)

                avg_counts.append(img_mean)

                # setting number of bins
                NBINS = 100

                if plots==True:
                    plt.figure()
                    plt.hist(image_data.flatten(),bins=NBINS,label='counts')
                    plt.axvline(x=img_min,linestyle='--',label='min {}'.format(img_min),alpha=0.5)
                    plt.axvline(x=img_max,linestyle='--',label='max {}'.format(img_max),alpha=0.5)
                    plt.axvline(x=img_mean,linestyle='-',linewidth=0.5,color='b',label='mean {:.2f}'.format(img_mean),alpha=1)
                    plt.legend()
                    plt.grid()
                    plt.xlabel('Count level in image')
                    plt.ylabel('Number of pixels with that count')
                    plt.title('Histogram of counts of {}'.format(img_name))
                    plt.savefig(plots_path/"hist_{}.jpg".format(img_name))
                    plt.show()

                avg_counts_to_return = np.mean(avg_counts)
                chip_num_to_return = str(chip_num)
                avg_counts_for_display.append('Chip'+str(chip_num_to_return)+':'+str(avg_counts_to_return.round(2)))
    return avg_counts_for_display

#%%
###############################################################################
#----------------SECTION FOUR: DATA REDUCTION FUNCTIONS-----------------------#
###############################################################################

def mbias_maker(bias_chip_sep_files,MBIAS_path,plots_path=None,plots=False):
    """
    This function deals with making a master bias for each chip. It creates
    FITS master bias files for each chip, and can also show image comparisons
    of a single bias vs a master bias, if that block of code is uncommented.

    Parameters
    ----------
    bias_chip_sep_files : list
        List of list of filenames of biases for each chip.

    MBIAS_path : WindowsPath object
        Path to directory where Master Biases are to be saved.
    
    plots_path : WindowsPath object
        Path to directory where Plots are to be saved.
    
    plots : bool, optional
        Boolean (True/False) of plot display option.

    Returns
    -------
    Nothing.
    """
    bias_list_to_return = []
    for BIAS_chips in bias_chip_sep_files:
        # seperating list of files by their exposure lengths
        exptime_seperated_files = exptime_separator(BIAS_chips)
        bias_list_to_return.append(exptime_seperated_files)
        # for each array of files for each exposure length:
        for exptime_seperated_exps in exptime_seperated_files:
            # extracting header information for this set of files
            hdu1 = fits.open(exptime_seperated_exps[0])
            exptime = hdu1[0].header['EXPTIME']
            chip_num = hdu1[0].header['CHIP']
            BIAS_ccd = CCDData(exptime_seperated_exps[0],unit=u.adu)

            # combining all the biases together
            master_bias = ccdp.combine(exptime_seperated_exps,unit=u.adu,
                                       method='average',
                                       sigma_clip=True,
                                       sigma_clip_low_thresh=5,
                                       sigma_clip_high_thresh=5,
                                       sigma_clip_func=np.ma.median,
                                       sigma_clip_dev_func=mad_std,
                                       mem_limit=350e6)

            # setting combined bias header
            master_bias.meta['combined'] = True
            master_bias.meta['CHIP'] = chip_num
            master_bias.meta['EXPTIME'] = exptime

            # writing combined bias as a fits file
            master_bias.write(MBIAS_path / 'mbias-{}-chip{}.fit'.format(exptime,chip_num),
                                                                          overwrite=True)
        if plots==True:
            # getting CCD image for plotting purposes
            bias_fits = fits.getdata(exptime_seperated_exps[0])
            bias_ccd = CCDData(bias_fits,unit=u.adu)

            # plotting single bias compared to combined bias
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
            show_image(BIAS_ccd, cmap='gray', ax=ax1, fig=fig, percl=90)
            ax1.set_title('Single calibrated bias for Chip {}'.format(chip_num))
            show_image(master_bias.data, cmap='gray', ax=ax2, fig=fig, percl=90)
            ax2.set_title('Master Bias for Chip {}'.format(chip_num))
            plt.savefig(plots_path/"bias_vs_mbias-{}-chip{}.jpg".format(exptime,chip_num))

    return bias_list_to_return #to add to calibration log


def dark_calibrator(dark_chip_sep_files,MBIAS_chip_sep_files,DARK_cal_path):
    """
    This function deals with calibrating the darks for each chip. It calibrates
    each dark (which has been chip-seprated) for each exposure using the master
    biases, and can also show image comparisons of a single dark vs a calibrated
    dark, if that block of code is uncommented.

    Parameters
    ----------
    dark_chip_sep_files : list
        List of list of filenames of darks for each chip.

    MBIAS_chip_sep_files : list
        List list of master biases for each chip.

    DARK_cal_path : WindowsPath object
        Path to directory where Calibrated Darks are to be saved.

    Returns
    -------
    Nothing.
    """
    for d_index,DARK_chips in enumerate(dark_chip_sep_files):
        # getting master biases for current chip
        mbias_for_chip = MBIAS_chip_sep_files[d_index]#[0]
        # converting each DARKS file to a fits array
        mbias_fits = fits.getdata(mbias_for_chip[0])
        # converting each DARK array to a CCDData object
        mbias_ccd = CCDData(mbias_fits,unit=u.adu)

        for DARK_file in DARK_chips:
            # extracting header data for display/saving purposes later
            hdu1 = fits.open(DARK_file)
            d_file_name = hdu1[0].header['RUN'].strip(' ')
            dark_exptime = hdu1[0].header['EXPTIME']
            d_obs_set = hdu1[0].header['SET'].strip(' ')
            d_chip_num = hdu1[0].header['CHIP']

            img_name = '{}-{}-{}-{}.fit'.format(d_file_name,dark_exptime,
                                                d_obs_set,d_chip_num)

            # converting each DARKS file to a fits array
            DARK_fits = fits.getdata(DARK_file)
            # converting each DARK array to a CCDData object
            DARK_ccd = CCDData(DARK_fits,unit=u.adu)

            # Subtract bias from each Dark
            cal_dark = ccdp.subtract_bias(DARK_ccd,mbias_ccd)
            # setting combined dark header
            cal_dark.meta['calibrated'] = 'True'
            cal_dark.meta['RUN'] = d_file_name
            cal_dark.meta['EXPTIME'] = dark_exptime
            cal_dark.meta['SET'] = d_obs_set
            cal_dark.meta['CHIP'] = d_chip_num
            # Save the result
            cal_dark.write(DARK_cal_path
                  / "calibrated_dark-{}".format(img_name),overwrite=True)


def mdark_maker(dark_chip_sep_files,MDARK_path,plots_path=None,plots=False):
    """
    This function deals with making a master dark for each chip. It creates
    FITS master dark files for each chip for each exposure time, and can also
    show image comparisons of a single dark vs a master dark, if that block of
    code is uncommented.

    Parameters
    ----------
    dark_chip_sep_files : list
        List of list of filenames of calibrated darks for each chip.

    MDARK_path : WindowsPath object
        Path to directory where Master Darks are to be saved.
    
    plots_path : WindowsPath object
        Path to directory where Plots are to be saved.
    
    plots : bool, optional
        Boolean (True/False) of plot display option.

    Returns
    -------
    Nothing.
    """
    for DARK_chips in dark_chip_sep_files:
        # seperating list of files by their exposure lengths
        exptime_seperated_files = exptime_separator(DARK_chips)

        # for each array of files for each exposure length:
        for exptime_seperated_exps in exptime_seperated_files:
            # extracting header information for this set of files
            hdu1 = fits.open(exptime_seperated_exps[0])
            exptime = hdu1[0].header['EXPTIME']
            chip_num = hdu1[0].header['CHIP']

            # combining all the darks of this set together
            master_dark = ccdp.combine(exptime_seperated_exps,unit=u.adu,
                                       method='average',
                                       sigma_clip=True,
                                       sigma_clip_low_thresh=5,
                                       sigma_clip_high_thresh=5,
                                       sigma_clip_func=np.ma.median,
                                       sigma_clip_dev_func=mad_std,
                                       mem_limit=350e6)

            # writing keywords to header
            master_dark.meta['combined'] = True
            master_dark.meta['EXPTIME'] = exptime
            master_dark.meta['CHIP'] = chip_num

            # writing combined dark as a fits file
            master_dark.write(MDARK_path / 'mdark-{}-chip{}.fit'.format(exptime,
                                                                  chip_num),
                                                             overwrite=True)
            if plots==True:
                # getting CCD image for plotting purposes
                dark_fits = fits.getdata(exptime_seperated_exps[0])
                dark_ccd = CCDData(dark_fits,unit=u.adu)
                # plotting single dark compared to combined dark
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
                #plotting single dark for current chip
                show_image(dark_ccd, cmap='gray', ax=ax1, fig=fig, percl=90)
                ax1.set_title('Single calibrated dark for Chip {}'.format(chip_num))
                # plotting combined dark for current chip
                show_image(master_dark.data, cmap='gray', ax=ax2, fig=fig, percl=90)
                ax2.set_title('{}s Master Dark for Chip {}'.format(exptime,chip_num))
                plt.savefig(plots_path/"cal-dark_vs_mdark-{}-chip{}.jpg".format(exptime,chip_num))

#%%

def flats_selector(flats_txt_path,FLAT_path_str,science_files,include=True):
    """
    This function selects good or bad flats and filters them from the flats in
    the directory.

    Parameters
    ----------
    flats_txt_path : str
        Name of text file containing filtering criteria for flats.
        Example: 'flats.txt'

    include : bool
        True if need to include the filtering criteria. 
        False if need to exclude the filtering criteria.
        Default = True.
    
    FLAT_path_str : str
        String of path to directory where Flats are stored.
    
    science_files : list
        List of science files to use.

    Returns
    -------
    good_flat_files : list
        List of filenames to ues for further flats processing.
    """
    with open(flats_txt_path) as f:
        list_of_selected_flats = f.read().splitlines() 
    
    good_flat_files = []
    
    if include is True:
        for a_selected_flat in list_of_selected_flats:
            good_flat_file = glob.glob(FLAT_path_str + a_selected_flat)
            good_flat_files += good_flat_file
    
    else:
        bad_flat_files = []
        for a_selected_flat in list_of_selected_flats:
            bad_flat_file = glob.glob(FLAT_path_str + a_selected_flat)
            bad_flat_files += bad_flat_file

        for science_flat in science_files:
            if science_flat not in bad_flat_files:
                good_flat_files.append(science_flat)
    
    return good_flat_files


def flat_calibrator(flat_chip_sep_files,MDARK_chip_sep_files,FLAT_cal_path,
                    actual_exposure_times,combined_darks):
    """
    This function deals with calibrating the flats for each chip. It calibrates
    each flat (which has been chip-separated) for each exposure using the master
    darks, and can also show image comparisons of a single flat vs a calibrated
    flat, if that block of code is uncommented.

    Parameters
    ----------
    flat_chip_sep_files : list
        List of list of filenames of flats for each chip.

    MDARK_chip_sep_files : list
        List of list of master darks for each chip.

    FLAT_cal_path : WindowsPath object
        Path to directory where Calibrated Flats are to be saved.

    actual_exposure_times : set
        Set of exposure times found.

    combined_dark : dict
        Dictionary of darks for each exposure time.

    Returns
    -------
    Nothing.
    """
    for f_index,FLAT_chips in enumerate(flat_chip_sep_files):
    # for FLAT_chips in flat_chip_sep_files:
            # getting master darks for current chip
            MDARK_chips_file = MDARK_chip_sep_files[f_index]
            # print(MDARK_chips_file)
            for FLAT_file in FLAT_chips:
                # print(FLAT_file)
                hdu1 = fits.open(FLAT_file)
                # extracting header data for later saving purposes
                f_file_name = hdu1[0].header['RUN'].strip(' ')
                flat_exptime = hdu1[0].header['EXPTIME']
                f_obs_set = hdu1[0].header['SET'].strip(' ')
                f_chip_num = hdu1[0].header['CHIP']

                # making CCDData object for flat which we are calibrating
                FLAT_ccd = CCDData.read(FLAT_file,unit='adu')

                for mdark in MDARK_chips_file:
                    # finding master dark of matching exp
                    mdark_hdu1 = fits.open(mdark)
                    mdark_ccd = CCDData.read(mdark,unit=u.adu)
                    mdark_exptime = mdark_hdu1[0].header['EXPTIME']

                    # here, we are finding an exact match for the flat
                    if mdark_exptime == flat_exptime:
                        MDARK_exptime = mdark_exptime
                        MDARK_to_subtract = CCDData.read(mdark,unit=u.adu)

                    else:
                        # Find the correct dark exposure
                        MDARK_exptime = find_nearest_dark_exposure(mdark_ccd,actual_exposure_times)
                        MDARK_to_subtract = combined_darks[MDARK_exptime]

                MDARK_exptime_u = MDARK_exptime*u.second #produces a Quantity object
                flat_exptime_u = flat_exptime*u.second   #produces a Quantity object

                cal_flat = ccdp.subtract_dark(ccd=FLAT_ccd, #CCD array of flat
                                              master=MDARK_to_subtract, #CCD array of master dark
                                              dark_exposure=MDARK_exptime_u,
                                              data_exposure=flat_exptime_u,
                                              scale=False)
                # Save the result
                cal_flat.write(FLAT_cal_path /
                      "calibrated_flat-{}-{}-{}-{}.fit".format(f_file_name,flat_exptime,
                                                                f_obs_set,f_chip_num),
                                                                overwrite=True)

def mflat_maker(cal_flat_chip_sep_files,MFLAT_path):
    """
    This function deals with making a master flat for each chip. It creates
    FITS master flat files for each chip for each exposure time, and can also
    show image comparisons of a single flat vs a master flat, if that block of
    code is uncommented.

    Parameters
    ----------
    cal_flat_chip_sep_files : list
        List of list of filenames of calibrated flats for each chip.

    MFLAT_path : WindowsPath object
        Path to directory where Master Flats are to be saved.

    Returns
    -------
    Nothing.
    """

    for FLAT_chips in cal_flat_chip_sep_files:
        # seperating list of files by their exposure lengths
        exptime_seperated_files = exptime_separator(FLAT_chips)

        # for each array of files for each exposure length:
        for exptime_seperated_exps in exptime_seperated_files:
                # extracting header information for this set of files
                hdu1 = fits.open(exptime_seperated_exps[0])
                exptime = hdu1[0].header['EXPTIME']
                chip_num = hdu1[0].header['CHIP']

                # getting CCD image for plotting purposes
                flat_fits = fits.getdata(exptime_seperated_exps[0])
                flat_ccd = CCDData(flat_fits,unit=u.adu)

                # combining all the flats of this set together
                master_flat = ccdp.combine(exptime_seperated_exps,unit=u.adu,
                                      method='average',
                                      sigma_clip=True,
                                      sigma_clip_low_thresh=5,
                                      sigma_clip_high_thresh=5,
                                      sigma_clip_func=np.ma.median,
                                      sigma_clip_dev_func=mad_std,
                                      mem_limit=350e6)

                # writing keywords to header
                master_flat.meta['combined'] = True
                master_flat.meta['EXPTIME'] = exptime
                master_flat.meta['CHIP'] = chip_num

                # writing combined dark as a fits file
                master_flat.write(MFLAT_path /
                                  'mflat-{}-chip{}.fit'.format(exptime,chip_num),
                                                                 overwrite=True)


def flats_count_classifier(flats_lst):
    """
    This function categorises each master flat according to its counts.

    Parameters (1)
    --------------
    flats_lst : list
        List of master flats.

    Returns (5)
    -----------
    hi_counts_flats : list
        List of flats with counts between 40,000 and 55,000.

    ok_counts_flats : list
        List of flats with counts between 30,000 and 40,000.

    lo_counts_flats : list
        List of flats with counts between 20,000 and 30,000.

    last_resort_flats : list
        List of flats with counts between 55,000 and 63,000, or 10,000 and 20,000.

    trash_counts_flats : list
        List of flats with counts above 63,000 or below 10,000.

    """
    hi_counts_flats = ['hi_counts',[]]
    ok_counts_flats = ['ok_counts',[]]
    lo_counts_flats = ['lo_counts',[]]
    last_resort_flats = ['last_resort',[]]
    trash_counts_flats = ['trash_counts',[]]

    for flat_img in flats_lst:
        flat_data = fits.getdata(flat_img)
        # getting counts for flat
        img_mean = np.mean(flat_data)

        if img_mean <= 55000 and img_mean > 40000:
            hi_counts_flats[1].append(flat_img)
        elif img_mean <= 40000 and img_mean > 30000:
            ok_counts_flats[1].append(flat_img)
        elif img_mean <= 30000 and img_mean >= 20000:
            lo_counts_flats[1].append(flat_img)
        elif img_mean <= 63000 and img_mean > 55000:
            last_resort_flats[1].append(flat_img)
        elif img_mean < 20000 and img_mean >= 10000:
            last_resort_flats[1].append(flat_img)
        else:
            trash_counts_flats[1].append(flat_img)

    return hi_counts_flats,ok_counts_flats,lo_counts_flats,last_resort_flats,trash_counts_flats


def mflat_maker_for_counts(counts_sep_flats,MFLAT_counts_path):
    """
    This function makes master flats for each chip for each category of counts.

    Parameters
    ----------
    counts_sep_flats : list
        List of lists of catgories and their associated master flat filenames.

    MFLAT_counts_path : WindowsPath object
        Path to directory where Master Flats by Counts are to be saved.

    Returns
    -------
    Nothing.

    """
    # for each category of counts
    for counts_sep_flats_categories in counts_sep_flats:
        this_category = counts_sep_flats_categories[0] #extracting category as a str
        counts_sep_flats_category = counts_sep_flats_categories[1]

        if len(counts_sep_flats_category) == 0:
            pass
        else:
            # seperating list of flats in this category by chip num
            chip_seperated_files = chip_separator(counts_sep_flats_category)
            # for each array of files for each chip length:
            for chip_files in chip_seperated_files:
                if len(chip_files) == 0:
                    pass
                else:
                    exptimes = []
                    for chip_file in chip_files:
                        # extracting header information for this set of files
                        hdu1 = fits.open(chip_file)
                        chip_num = hdu1[0].header['CHIP']

                    # combining all the flats of this set together
                    master_flat = ccdp.combine(chip_files,unit=u.adu,
                                              method='average',
                                              sigma_clip=True,
                                              sigma_clip_low_thresh=5,
                                              sigma_clip_high_thresh=5,
                                              sigma_clip_func=np.ma.median,
                                              sigma_clip_dev_func=mad_std,
                                              mem_limit=350e6)

                    # writing keywords to header
                    master_flat.meta['combined'] = True
                    master_flat.meta['CATEGORY'] = this_category
                    master_flat.meta['CHIP'] = chip_num

                    # writing combined dark as a fits file
                    master_flat.write(MFLAT_counts_path /
                                          'mflat-{}-chip{}.fit'.format(this_category,chip_num),
                                                                         overwrite=True)


def ALERT_reducer(target_names_dict,reduced_ALERT_path,MDARK_chip_sep_files,
                  MFLAT_counts_chips_files,MDARK_imgs,combined_darks,plots_path,plots=False):
    """
    This function reduces the ALERT data using the appropriate darks and flats.

    Parameters (8)
    ----------
    target_names_dict : dict
        Dictionary of target names and their filenames.

    reduced_ALERT_path : WindowsPath object
        Path to directory where Reduced ALERTs are to be saved.

    MDARK_chip_sep_files : list
        List of list of master darks for each chip.

    MFLAT_counts_chips_files : list
        List of list of master flats for each chip.

    MDARK_imgs : ImageFileCollection (ccdproc.image_collection.ImageFileCollection)
        Image collection of master darks.

    combined_darks : dict
        Dictionary of combined darks.
    
    plots_path : WindowsPath object
        Path to directory where Plots are to be saved.
    
    plots : bool, optional
        Boolean (True/False) of plot display option.

    Returns
    -------
    Nothing.

    """
    for key, value in target_names_dict.items():
        ALERT_chips_files = chip_separator(value)

        for a_index,ALERT_chips in enumerate(ALERT_chips_files):
            # getting master darks  and flats for current chip
            MDARK_chips_file = MDARK_chip_sep_files[a_index]
            MFLAT_chips_file = MFLAT_counts_chips_files[a_index]

            ALERT_exptimes = exptime_checker(ALERT_chips)

            # finding closest dark exposure times to ALERT exposure times
            d_actual_exposure_times = set(h['EXPTIME'] for h in MDARK_imgs.headers(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
                                                                         combined=True))
            for ALERT_file in ALERT_chips:
                al_hdu1 = fits.open(ALERT_file)
                # extracting header data for later saving purposes
                al_file_name = al_hdu1[0].header['RUN'].strip(' ')
                ALERT_exptime = al_hdu1[0].header['EXPTIME']
                al_obs_set = al_hdu1[0].header['SET'].strip(' ')
                al_chip_num = al_hdu1[0].header['CHIP']
                al_filter = al_hdu1[0].header['COLOUR'].strip(' ')

                # making CCDData object for ALERT which we are calibrating
                ALERT_ccd = CCDData.read(ALERT_file,unit=u.adu)

                for mdark in MDARK_chips_file:
                    # finding master dark of matching exp
                    mdark_hdu1 = fits.open(mdark)
                    mdark_ccd = CCDData.read(mdark,unit=u.adu)
                    mdark_exptime = mdark_hdu1[0].header['EXPTIME']

                    # here, we are finding an exact match for the ALERT
                    if mdark_exptime == ALERT_exptime:
                        MDARK_exptime = mdark_exptime
                        MDARK_to_subtract = CCDData.read(mdark,unit=u.adu)
                        print("MDARK_exptime",MDARK_exptime)

                    else:
                        # Find the correct dark exposure
                        MDARK_exptime = find_nearest_dark_exposure(mdark_ccd,d_actual_exposure_times)
                        MDARK_to_subtract = combined_darks[MDARK_exptime]
                        print("MDARK_exptime",MDARK_exptime)
                try:
                    # for mflat_lst in MFLAT_chips_file:
                    for mflat in MFLAT_chips_file:
                        mflat_hdu1 = fits.open(mflat)
                        mflat_ccd = CCDData.read(mflat,unit=u.adu)
                        mflat_exptime = mflat_hdu1[0].header['EXPTIME']
                        mflat_category = mflat_hdu1[0].header['CATEGORY']
    
                        if mflat_category == 'hi_counts':
                            MFLAT_exptime = mflat_exptime
                            MFLAT_to_divide = CCDData(mflat_ccd,unit=u.adu)
                        elif mflat_category == 'ok_counts':
                            MFLAT_exptime = mflat_exptime
                            MFLAT_to_divide = CCDData(mflat_ccd,unit=u.adu)
                        elif mflat_category == 'lo_counts':
                            MFLAT_exptime = mflat_exptime
                            MFLAT_to_divide = CCDData(mflat_ccd,unit=u.adu)
                        elif mflat_category == 'last_resort':
                            MFLAT_exptime = mflat_exptime
                            MFLAT_to_divide = CCDData(mflat_ccd,unit=u.adu)
                        else: #trash_counts
                            MFLAT_exptime = mflat_exptime
                            MFLAT_to_divide = CCDData(mflat_ccd,unit=u.adu)
                            
                    
                    reduced_ALERT = ccdp.ccd_process(ALERT_ccd,
                                                     dark_frame=MDARK_to_subtract,
                                                     master_flat=MFLAT_to_divide,
                                                     data_exposure=ALERT_exptime*u.second,
                                                     dark_exposure=MDARK_exptime*u.second)
                except UnboundLocalError:
                    # for mflat_lst in MFLAT_chips_file:
                    for mflat in MFLAT_chips_file:
                        mflat_hdu1 = fits.open(mflat)
                        mflat_ccd = CCDData.read(mflat,unit=u.adu)
                        mflat_exptime = mflat_hdu1[0].header['EXPTIME']
                        mflat_category = mflat_hdu1[0].header['CATEGORY']
    
                        if mflat_category == 'hi_counts':
                            MFLAT_exptime = mflat_exptime
                            MFLAT_to_divide = CCDData(mflat_ccd,unit=u.adu)
                        elif mflat_category == 'ok_counts':
                            MFLAT_exptime = mflat_exptime
                            MFLAT_to_divide = CCDData(mflat_ccd,unit=u.adu)
                        elif mflat_category == 'lo_counts':
                            MFLAT_exptime = mflat_exptime
                            MFLAT_to_divide = CCDData(mflat_ccd,unit=u.adu)
                        elif mflat_category == 'last_resort':
                            MFLAT_exptime = mflat_exptime
                            MFLAT_to_divide = CCDData(mflat_ccd,unit=u.adu)
                        else: #trash_counts
                            pass
                            
                    reduced_ALERT = ccdp.ccd_process(ALERT_ccd,
                                                     dark_frame=MDARK_to_subtract,
                                                     master_flat=MFLAT_to_divide,
                                                     data_exposure=ALERT_exptime*u.second,
                                                     dark_exposure=MDARK_exptime*u.second)
                    
                
                filename_to_write = "reduced-{}-{}-{}-{}-{}-{}.fit".format(key.strip(' '),
                                                                     al_file_name,
                                                                     ALERT_exptime,
                                                                     al_filter,
                                                                     al_obs_set,
                                                                     al_chip_num)
                # writing keywords to header
                reduced_ALERT.meta['reduced'] = True
                reduced_ALERT.meta['EXPTIME'] = ALERT_exptime
                reduced_ALERT.meta['CHIP'] = al_chip_num
                reduced_ALERT.meta['RUN'] = al_file_name
                reduced_ALERT.meta['SET'] = al_obs_set
                reduced_ALERT.meta['COLOUR'] = al_filter
                reduced_ALERT.meta['FILENAME'] = filename_to_write

                # Save the result
                reduced_ALERT.write(reduced_ALERT_path/filename_to_write,overwrite=True)
                if plots==True:
                    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
                    show_image(ALERT_ccd, cmap='gray', ax=ax1, fig=fig, percl=90)
                    ax1.set_title('Raw ALERT Image for Chip {} for {}'.format(al_chip_num,key))
                    show_image(reduced_ALERT.data, cmap='gray', ax=ax2, fig=fig, percl=90)
                    ax2.set_title('Reduced ALERT Image for Chip {} for {}'.format(al_chip_num,key))
                    plt.savefig(plots_path/"raw_vs_reduced_ALERT-{}-{}-{}-{}-{}-{}.jpg".format(key.strip(' '),
                                                                                          al_file_name,
                                                                                          ALERT_exptime,
                                                                                          al_filter,
                                                                                          al_obs_set,
                                                                                          al_chip_num))

###############################################################################
#-------------------------------END OF CODE-----------------------------------#
###############################################################################