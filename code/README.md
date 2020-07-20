# Folders

 - archive : Where old modules and notebooks are stored.

 - temp : where temporary files are stored. Usually, results outputted by scripts are put there.

 - cpp : where the source code for the c++ modules are located.


# Modules

 - corrfunc.py : The main module where the correlation function and the helper functions are located.

 - test_corrfunc.py : Where unit tests on the corrfunc module are located.

 - processing.py : This file contains functions that process the dustmap data. As soon as the consistency tests are done most of the function are going into the archive.

 - plot_corr_script.py : This module is used to make the correlation plots with the correlation function in corrfunc.py, and it dumps the results in the temp folder. This script usually runs in the background, as for high resolution it can take time to produce results and so notebooks are not so much an option.

# Notebooks

 - Two-point-correlation-function.ipynb : This file is used to run interactive tests for either debugging, or to deepen one's understanding of the correlation function. 
