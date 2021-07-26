# DRP > src
This folder contains the executable files for running the program.

Main (Important) Running Files
------------------------------

- ```master_drp_runner.py``` is the file required to do the data reduction. In this file, you may change the file directories as required. Shell scripting for this pipeline is a future project, for the time being, please manually run this file.

- ```drp_funcs.py``` is the file containing the required functions to run the program. This is automatically read into ```master_drp_runner.py```. Do not modify unless necessary.

- ```asp.py``` is the file for doing astrometric calibration. Currently it functions as an independent piece of code, however will eventually be merged into ```master_drp_runner.py``` and ```drp_funcs.py```. "ASP" stands for Astrometry Solving Pipeline.

Miscellaneous Files
------------------------------

- ```convenience_functions.py``` is a file containing useful functions for some functions in the program. Do not modify.

- ```calibration_log.txt``` is a text file which logs the data reduction output from ```master_drp_runner.py```. This file needs to be present (can be empty) in the repository in order for ```master_drp_runner.py``` to run, i.e. ```master_drp_runner.py``` does not create this file, it assumes its existence in the folder.

- ```flats.txt``` is a text file which contains filtering keywords to select certain flat files. This file needs to be present in the repository if its associated function, ```flats_selector``` is to be run.

- ```Trash``` is a folder containing files that are no longer required. Ignore this folder.
