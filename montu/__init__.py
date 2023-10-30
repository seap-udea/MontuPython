###############################################################
# Import montu modules
###############################################################
from montu.version import *
from montu.util import *
from montu.time import *
from montu.stars import *
from montu.observer import *
from montu.sebau import *

###############################################################
# Aliases
###############################################################
D2H = Util.dec2hex
VPRINT = Util.vprint
PRINTDF = Util.print_df
TABLEDF = Util.table_df

###############################################################
# External modules
###############################################################
import warnings
import numpy as np

###############################################################
# Constants
###############################################################
# Numerical Constants
RAD = 180/np.pi
DEG = 1/RAD
MARCSEC = 1e-3/3600 # milliarcsec in degrees 


###############################################################
# Initialization
###############################################################
# Avoid warnings
warnings.filterwarnings("ignore")

# Load basic SPICE kernels
Util.load_kernels()

# Showing version 
print(f"Running MontuPython version {version}")
