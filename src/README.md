# DRP > src
This folder contains the executable files for running the program.

```master_drp_runner.py``` is the file required to do the data reduction. Here, you may change the file directories as required.

```drp_funcs.py``` is the file containing the required functions to run the program. This is automatically read into ```master_drp_runner.py```. Do not modify unless necessary.

```convenience_functions.py``` is a file containing useful functions for some functions in the program. Do not modify.

```calibration_log.txt``` is a text file which logs the data reduction output from ```master_drp_runner.py```. This file needs to be present in the repository in order for ```master_drp_runner.py``` to run.

```flats.txt``` is a text file which contains filtering keywords to select certain flat files. This file needs to be present in the repository if its associated function, ```flats_selector``` is to be run.

```Trash``` is a folder containing files that are no longer required. Ignore this folder.
