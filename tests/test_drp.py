# -*- coding: utf-8 -*-
#! /usr/bin/env python
"""
Created on Tue Apr  6 12:15:37 2021

@author: ave41
"""
# importing drp functions from /src/ folder
import sys
sys.path.insert(1,"C:\\Users\\ave41\\OneDrive - University of Canterbury\\Master's 2021\\ASTR480 Research\\ASTR480 Code\\Data Reduction Pipeline\\DataReductionPipeline\\src")
from drp_funcs import *
import unittest

###############################################################################
#---------------------SECTION ONE: READING IN DATA----------------------------# 
###############################################################################

# reduced_ALERT_path = path_checker(ALERT_path,'Reduced ALERT')
to_include = ['/*-1.fit','/*-2.fit','/*-3.fit','/*-4.fit','/*-5.fit',
              '/*-6.fit','/*-7.fit','/*-8.fit','/*-9.fit','/*-10.fit']
# reading in bias files from BIAS folder
BIAS_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/"
                 "Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/"
                 "ObsData_v4/DARK")
BIAS_imgs = ImageFileCollection(BIAS_path,glob_exclude=['/*-0.fit','/*-99.fit'])
BIAS_files = BIAS_imgs.files_filtered(EXPTIME=0,include_path=True) # EXPTIME is either 0 or 1
# making/checking MBIAS path/folder
MBIAS_path = path_checker(BIAS_path,'Master Biases')
MBIAS_imgs = ImageFileCollection(MBIAS_path, keywords='*')
MBIAS_files = MBIAS_imgs.files_filtered(COMBINED=True,
                                        include_path=True)
# reading in dark files from DARK folder
DARK_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/"
                 "ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/"
                 "ObsData_v4/DARK")
good_files = []
for i in to_include:
    good_file = glob.glob("C:\\Users\\ave41\\OneDrive - University of Canterbury\\Master's 2021\\ASTR480 Research\\ASTR480 Code\\Data Reduction Pipeline\\ObsData_v4\\DARK" 
                          + i)
    good_files += good_file

# selecting images
DARK_imgs = ImageFileCollection(filenames=good_files)
DARK_files = DARK_imgs.files_filtered(FIELD='              dark',include_path=True)

# now we need to get rid of the biases (0s and 1s exposures)
# the >1 condition will take care of it - no darks strictly less than 1s 
# will be calibrated. This is because sometimes the biases can be either 0s or 1s.
good_DARK_files = []
for DARK_file in DARK_files:
    hdu1 = fits.open(DARK_file)
    dark_exptime = hdu1[0].header['EXPTIME']
    if dark_exptime >1:
        good_DARK_files.append(DARK_file)
# making/checking Calibrated Darks path/folder
DARK_cal_path = path_checker(DARK_path,'Calibrated Darks')
DARK_cal_imgs = ImageFileCollection(DARK_cal_path)
DARK_cal_files = DARK_cal_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
                                              include_path=True)
# making/checking MDARK path/folder
MDARK_path = path_checker(DARK_path,'Master Darks')
MDARK_imgs = ImageFileCollection(MDARK_path, keywords='*')
MDARK_files = MDARK_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
                                        include_path=True)
FLAT_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/"
                 "ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/"
                 "ObsData_v3/FLAT")
good_files = []
for i in to_include:
    good_file = glob.glob("C:\\Users\\ave41\\OneDrive - University of Canterbury\\Master's 2021\\ASTR480 Research\\ASTR480 Code\\Data Reduction Pipeline\\ObsData_v3\\FLAT" 
                          + i)
    good_files += good_file

# selecting images
FLAT_imgs = ImageFileCollection(filenames=good_files)
FLAT_files = FLAT_imgs.files_filtered(FIELD ='              flat',
                                  include_path=True)
# making/checking Calibrated Flats path/folder
FLAT_cal_path = path_checker(FLAT_path,'Calibrated Flats')
FLAT_cal_imgs = ImageFileCollection(FLAT_cal_path)
FLAT_cal_files = FLAT_cal_imgs.files_filtered(FIELD   = '              flat' ,
                                              include_path=True)
# making/checking MFLAT path/folder
MFLAT_path = path_checker(FLAT_path,'Master Flats')
MFLAT_imgs = ImageFileCollection(MFLAT_path, keywords='*')
MFLAT_files = MFLAT_imgs.files_filtered(FIELD   = '              flat',
                                        include_path=True)


test_fits_file = "C://Users//ave41//OneDrive - University of Canterbury//Master's 2021//ASTR480 Research//ASTR480 Code//Data Reduction Pipeline//DataReductionPipeline//tests//D21350-60-a-5.fit"
BIAS_chips_files = chip_separator(BIAS_files)
DARK_cal_chips_files = chip_separator(DARK_cal_files)

###############################################################################
#-------------------------SECTION TWO: UNIT TESTS-----------------------------# 
###############################################################################
    
class Test(unittest.TestCase):  
 #------------------------------------------------------------------------------ 
    def test_keyword_checker(self):
        """
        Tests if image headers contains certain keywords.
        This works in levels, with an image from each level being checked.
        Unprocessed data > Calibrated data > Master data > Reduced ALERTs
        """
        expected_keywords = ['EXPTIME','CHIP','RUN','SET']
        lst_of_levels_to_check = [BIAS_files[0],MBIAS_files[0],DARK_files[0],DARK_cal_files[0],
                                  MDARK_files[0],FLAT_files[0],FLAT_cal_files[0],MFLAT_files[0]]
        for level in lst_of_levels_to_check:
            for keyword in expected_keywords:
                try:
                    hdu1 = fits.open(level) 
                    hdr = hdu1[0].header
                    assert keyword in hdr
                except:
                    print("The keyword {} is not contained within the header of {}. Please inspect the header.".format(keyword,level))
#------------------------------------------------------------------------------
    
    def test_chipNumExtractor1(self):
        """
        Tests if the original function provides the expected output on a 
        test fits file.
        """
        try:
            actual = chip_num_extractor(test_fits_file)
            expected = 5
            self.assertEqual(actual,expected)
        except:
            print("There's something wrong with chip_num_extractor. Not providing expected output.")
    
    def test_chipNumExtractor2(self):
        """
        Tests if the image header has the keyword 'CHIP'.
        """
        try:
            hdu1 = fits.open(test_fits_file)
            actual = hdu1[0].header['CHIP']
            expected = 5
            self.assertEqual(actual,expected)
        except:
            print("FITS header is missing keyword 'CHIP'. Needed for function to work.")
        
    def test_chipNumExtractor3(self):
        """
        Tests if the image input is a string.
        """
        try:
            actual = type(test_fits_file)
            expected = str
            self.assertEqual(actual,expected)
        except:
            print("The input is not a string. Please enter the image path as a string.")
        
    def test_chipNumExtractor4(self):
        """
        Tests if the image is a CCD-type object when opened as a FITS file.
        """
        try:
            open_good_fits = fits.getdata(test_fits_file) 
            ccd = CCDData(open_good_fits,unit=u.adu)
            actual = type(ccd)
            expected = astropy.nddata.ccddata.CCDData
            self.assertEqual(actual,expected)
        except:
            print("The FITS image has not been able to convert to a astropy.nddata.ccddata.CCDData object.")
    
    def test_chipNumExtractor5(self):
        """
        Tests if the image header contains the keyword 'CHIP'.
        """
        try:
            hdu1 = fits.open(test_fits_file) 
            hdr = hdu1[0].header
            assert 'CHIP' in hdr
        except:
            print("The keyword 'CHIP' is not contained within the header. Please inspect the header.")
#------------------------------------------------------------------------------  
    
    def test_chipSeparator1(self):
        """
        Tests if the input is a list.
        """
        try:
            actual = type(BIAS_chips_files)
            expected = list
            self.assertEqual(actual,expected)
        except:
            print("The input is not a list.")
            
    def test_chipSeparator2(self):
        """
        Tests if the content in the input list is a string.
        """
        try:
            func = chip_separator(BIAS_files)
            actual = type(func[0][0])
            expected = str
            self.assertEqual(actual,expected)
        except:
            print("The information inside the input list is not of type 'string'.")
#------------------------------------------------------------------------------  

    def test_path_checker2(self):
        """
        Tests if the input folder_name is a string.
        """
        folder_name = 'Master Biases'
        try:
            actual = type(folder_name)
            expected = str
            self.assertEqual(actual,expected)
        except:
            print("The input folder_name is not a string. It should be a string.")
#------------------------------------------------------------------------------  
    
    def test_exptime_checker1(self):
        """
        Tests if IMAGElist is a list.
        """
        IMAGElist = DARK_cal_chips_files[0]
        try:
            actual = type(IMAGElist)
            expected = list
            self.assertEqual(actual,expected)
        except:
            print("The input IMAGElist is not a list. It should be a list.")
            
    def test_exptime_checker2(self):
        """
        Tests if IMAGElist[i] is a string.
        """
        IMAGElist_component = DARK_cal_chips_files[0][0]
        try:
            actual = type(IMAGElist_component)
            expected = str
            self.assertEqual(actual,expected)
        except:
            print("The input IMAGElist[0] is not a str. It should be a str.")
            
    def test_exptime_checker3(self):
        """
        Tests if EXPTIME is an int.
        """
        IMAGElist_component = DARK_cal_chips_files[0][0]
        hdu1 = fits.open(IMAGElist_component)
        exptime = hdu1[0].header['EXPTIME']
        try:
            actual = type(exptime)
            expected = int
            self.assertEqual(actual,expected)
        except:
            print("The EXPTIME in the header is not an int. It should be an int.")
#------------------------------------------------------------------------------ 
    
    def test_exptime_separator1(self):
        """
        Tests if the input IMAGElist is a list.
        """
        IMAGElist = DARK_cal_chips_files[0]
        try:
            actual = type(IMAGElist)
            expected = list
            self.assertEqual(actual,expected)
        except:
            print("The input IMAGElist is not a list. It should be a list.")
    
    def test_exptime_separator2(self):
        """
        Tests if IMAGElist[i] is a string.
        """
        IMAGElist_component = DARK_cal_chips_files[0][0]
        try:
            actual = type(IMAGElist_component)
            expected = str
            self.assertEqual(actual,expected)
        except:
            print("The input IMAGElist[0] is not a str. It should be a str.")
    
    def test_exptime_separator3(self):
        """
        Tests if EXPTIME is an int.
        """

        IMAGElist_component = DARK_cal_chips_files[0][0]
        hdu1 = fits.open(IMAGElist_component)
        exptime = hdu1[0].header['EXPTIME']
        try:
            actual = type(exptime)
            expected = int
            self.assertEqual(actual,expected)
        except:
            print("The EXPTIME in the header is not an int. It should be an int.")
 #------------------------------------------------------------------------------ 
    
    def test_img_stats1(self):
        """
        Tests if the input IMAGElist is a list.
        """
        try:
            actual = type(DARK_cal_files)
            expected = list
            self.assertEqual(actual,expected)
        except:
            print("The input IMAGElist is not a list. It should be a list.")   
#------------------------------------------------------------------------------ 

    def test_mbias_maker1(self):
        """
        Tests if bias_chip_sep_files is a list.
        """
        try:
            actual = type(BIAS_chips_files)
            expected = list
            self.assertEqual(actual,expected)
        except:
            print("The input bias_chip_sep_files is not a list. It should be a list.") 
    
    def test_mbias_maker2(self):
        """
        Tests if bias_chip_sep_files[0] is a list.
        """
        try:
            actual = type(BIAS_chips_files[0])
            expected = list
            self.assertEqual(actual,expected)
        except:
            print("The bias_chip_sep_files[0] is not a list. It should be a list.") 
    
    def test_mbias_maker3(self):
        """
        Tests if bias_chip_sep_files[0][0] is a str.
        """

        try:
            actual = type(BIAS_chips_files[0][0])
            expected = str
            self.assertEqual(actual,expected)
        except:
            print("The bias_chip_sep_files[0][0] is not a str. It should be a str.") 

#------------------------------------------------------------------------------ 

    # def test_dark_calibrator1(self):
    #     """
    #     Tests if dark_chip_sep_files is a list.
    #     """
    
    
    
    
    
    
#================================ don't touch ================================#

###############################################################################
#-------------------------------END OF CODE-----------------------------------# 
############################################################################### 

if __name__ == '__main__':
    unittest.main()
    
###############################################################################
#-------------------------------END OF CODE-----------------------------------# 
############################################################################### 

    
# tests to write:
# path_checker(origin_path,folder_name)
    # 1. origin_path is Pathobject [done - but not working]
    # 2. folder_name is string [done]
    # 3. os is imported?
    # 4. the os in which we are running this code is Windows?
# exptime_checker(IMAGElist)
    # 1. IMAGElist is list [done]
    # 2. IMAGElist[i] is str [done]
    # 3. IMAGElist has dimensions 10 [not needed]
    

# already written:
# chip_num_extractor(img)
    # 1. img is string [done]
    # 2. img header has keyword CHIP [done]
    # 3/ img is np.array [done]
# chip_separator(IMAGElist)
    # 1. IMAGElist is list [done]
    # 2. IMAGElist[i] is string [done]
    
    
# not needed:
    # def test_exptime_checker3(self):
    #     """
    #     Tests if len(IMAGElist) == 10, one for each chip.
    #     """
    #     DARK_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/ObsData_v3/DARK")
    #     DARK_cal_path = path_checker(DARK_path,'Calibrated Darks')
    #     DARK_cal_imgs = ImageFileCollection(DARK_cal_path)
    #     DARK_cal_files = DARK_cal_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
    #                                           include_path=True)
    #     DARK_cal_chips_files = chip_separator(DARK_cal_files)
    #     IMAGElist = DARK_cal_chips_files[0]
    #     try:
    #         actual = len(IMAGElist)
    #         expected = 10
    #         self.assertEqual(actual,expected)
    #     except:
    #         print("The length of IMAGElist is not 10. It should be 10, one for each chip.")   
    
    
    
## COME BACK TO THIS UNIT TEST LATER:
    
        # def test_path_checker1(self):
    #     """
    #     Tests if the input Path is a PathObject.
    #     """
    #     BIAS_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/ObsData_v3/DARK")
    #     try:
    #         # actual = str(type(BIAS_path))
    #         actual = str(type(BIAS_path))
    #         expected = 'pathlib.WindowsPath'
    #         self.assertEqual(actual,expected)
    #     except:
    #         print("There's something wrong with path_checker. Not providing expected output.")
    
    # def test_path_checker1(self):
    #     """
    #     Tests if the input Path is a PathObject.
    #     """
    #     BIAS_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/ObsData_v3/DARK")
    #     try:
    #         # actual = str(type(BIAS_path))
    #         actual = BIAS_path
    #         expected = 'pathlib.WindowsPath'
    #         self.isinstance(actual,pathlib.WindowsPath)
        # except:
        #     print("There's something wrong with path_checker. Not providing expected output.")
    
    
    
    