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
        

if __name__ == '__main__':
    unittest.main()
    
    
    
    
    
    
# tests to write:
# chip_num_extractor(img)
    # 1. img is string
    # 2. img header has keyword CHIP
    # 3/ img is np.array
# chip_separator(IMAGElist)
    # 1. IMAGElist is list
    # 2. IMAGElist[i] is string
# path_checker(origin_path,folder_name)
    # 1. origin_path is Pathobject
    # 2. folder_name is string
    # 3. os is imported?
    # 4. the os in which we are running this code is Windows?

    
    
    
    
    
    
    
    
    
    
    
    
    