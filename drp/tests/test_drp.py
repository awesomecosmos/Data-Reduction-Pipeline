# -*- coding: utf-8 -*-
#! /usr/bin/env python
"""
Created on Tue Apr  6 12:15:37 2021

@author: ave41
"""

import sys
sys.path.insert(1,"C:\\Users\\ave41\\OneDrive - University of Canterbury\\Master's 2021\\ASTR480 Research\\ASTR480 Code\\Data Reduction Pipeline\\DataReductionPipeline\\drp\\src")
from drp import chip_num_extractor

import unittest
from astropy.io import fits
    
class Test(unittest.TestCase):
    def test_chipNumExtractor1_success(self):
        good_fits = "C://Users//ave41//OneDrive - University of Canterbury//Master's 2021//ASTR480 Research//ASTR480 Code//Data Reduction Pipeline//DataReductionPipeline//test//D21350-60-a-5.fit"
        actual = chip_num_extractor(good_fits)
        expected = 5
        self.assertEqual(actual,expected)
    
    # def test_chipNumExtractor1_exception(self):
    #     bad_fits = 
        

if __name__ == '__main__':
    unittest.main()
    
    
    
    
    
    
# tests to write:
# chip_num_extractor(img)
    # 1. img is string
    # 2. img header has keyword CHIP
# chip_separator(IMAGElist)
    # 1. IMAGElist is list
    # 2. IMAGElist[i] is string
# path_checker(origin_path,folder_name)
    # 1. origin_path is Pathobject
    # 2. folder_name is string
    # 3. os is imported?
    # 4. the os in which we are running this code is Windows?

    
    
    
    
    
    
    
    
    
    
    
    
    