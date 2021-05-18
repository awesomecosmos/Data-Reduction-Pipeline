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
    
class Test(unittest.TestCase):  
#------------------------------------------------------------------------------
    def test_chipNumExtractor1(self):
        """
        Tests if the original function provides the expected output on a 
        test fits file.
        """
        good_fits = "C://Users//ave41//OneDrive - University of Canterbury//Master's 2021//ASTR480 Research//ASTR480 Code//Data Reduction Pipeline//DataReductionPipeline//tests//D21350-60-a-5.fit"
        try:
            actual = chip_num_extractor(good_fits)
            expected = 5
            self.assertEqual(actual,expected)
        except:
            print("There's something wrong with chip_num_extractor. Not providing expected output.")
    
    def test_chipNumExtractor2(self):
        """
        Tests if the image header has the keyword 'CHIP'.
        """
        good_fits = "C://Users//ave41//OneDrive - University of Canterbury//Master's 2021//ASTR480 Research//ASTR480 Code//Data Reduction Pipeline//DataReductionPipeline//tests//D21350-60-a-5.fit"
        try:
            hdu1 = fits.open(good_fits)
            actual = hdu1[0].header['CHIP']
            expected = 5
            self.assertEqual(actual,expected)
        except:
            print("FITS header is missing keyword 'CHIP'. Needed for function to work.")
        
    def test_chipNumExtractor3(self):
        """
        Tests if the image input is a string.
        """
        good_fits = "C://Users//ave41//OneDrive - University of Canterbury//Master's 2021//ASTR480 Research//ASTR480 Code//Data Reduction Pipeline//DataReductionPipeline//tests//D21350-60-a-5.fit"
        try:
            actual = type(good_fits)
            expected = str
            self.assertEqual(actual,expected)
        except:
            print("The input is not a string. Please enter the image path as a string.")
        
    def test_chipNumExtractor4(self):
        """
        Tests if the image is a CCD-type object when opened as a FITS file.
        """
        good_fits = "C://Users//ave41//OneDrive - University of Canterbury//Master's 2021//ASTR480 Research//ASTR480 Code//Data Reduction Pipeline//DataReductionPipeline//tests//D21350-60-a-5.fit"
        try:
            open_good_fits = fits.getdata(good_fits) 
            ccd = CCDData(open_good_fits,unit=u.adu)
            actual = type(ccd)
            expected = astropy.nddata.ccddata.CCDData
            self.assertEqual(actual,expected)
        except:
            print("The FITS image has not been able to convert to a astropy.nddata.ccddata.CCDData object.")
#------------------------------------------------------------------------------  
    
    def test_chipSeparator1(self):
        """
        Tests if the input is a list.
        """
        BIAS_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/ObsData_v3/DARK")
        BIAS_imgs = ImageFileCollection(BIAS_path,glob_exclude=['/*-0.fit','/*-99.fit'])
        BIAS_files = BIAS_imgs.files_filtered(EXPTIME=1,include_path=True)
        try:
            test_lst = BIAS_files
            func = chip_separator(test_lst)
            actual = type(func)
            expected = list
            self.assertEqual(actual,expected)
        except:
            print("The input is not a list.")
            
    def test_chipSeparator2(self):
        """
        Tests if the content in the input list is a string.
        """
        BIAS_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/ObsData_v3/DARK")
        BIAS_imgs = ImageFileCollection(BIAS_path,glob_exclude=['/*-0.fit','/*-99.fit'])
        BIAS_files = BIAS_imgs.files_filtered(EXPTIME=1,include_path=True)
        try:
            test_lst = BIAS_files
            func = chip_separator(test_lst)
            actual = type(func[0][0])
            expected = str
            self.assertEqual(actual,expected)
        except:
            print("The information inside the input list is not of type 'string'.")
#------------------------------------------------------------------------------  

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
    
    def test_path_checker1(self):
        """
        Tests if the input Path is a PathObject.
        """
        BIAS_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/ObsData_v3/DARK")
        try:
            # actual = str(type(BIAS_path))
            actual = BIAS_path
            expected = 'pathlib.WindowsPath'
            self.isinstance(actual,pathlib.WindowsPath)
        except:
            print("There's something wrong with path_checker. Not providing expected output.")
            
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
        DARK_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/ObsData_v3/DARK")
        DARK_cal_path = path_checker(DARK_path,'Calibrated Darks')
        DARK_cal_imgs = ImageFileCollection(DARK_cal_path)
        DARK_cal_files = DARK_cal_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
                                              include_path=True)
        DARK_cal_chips_files = chip_separator(DARK_cal_files)
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
        DARK_path = Path("C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/Data Reduction Pipeline/ObsData_v3/DARK")
        DARK_cal_path = path_checker(DARK_path,'Calibrated Darks')
        DARK_cal_imgs = ImageFileCollection(DARK_cal_path)
        DARK_cal_files = DARK_cal_imgs.files_filtered(SUBBIAS = 'ccd=<CCDData>, master=<CCDData>',
                                              include_path=True)
        DARK_cal_chips_files = chip_separator(DARK_cal_files)
        IMAGElist_component = DARK_cal_chips_files[0][0]
        try:
            actual = type(IMAGElist_component)
            expected = str
            self.assertEqual(actual,expected)
        except:
            print("The input IMAGElist[0] is not a str. It should be a str.")
    



#=============================================================================#
if __name__ == '__main__':
    unittest.main()
#=============================================================================#

    
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
    
    
    
    
    
    
    
    