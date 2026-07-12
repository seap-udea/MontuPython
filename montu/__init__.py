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
D2H = Util.dec2sex
S2D = Util.sex2dec
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

# Load planetary data
ALLPLANETS = Util.load_planets()

# Showing version 
def welcome_version():
    print(f"MontuPython version {version}. 𓋹 𓍘 𓋴 𓎛 𓂡 𓁘 (ii-ti m Htp, HkAx Hn'-k)")

def welcome_translate():
    print(f"𓋹 𓍘 𓋴 (ii-ti m Htp: Ii-ti em hotep = Welcome in peace) 𓎛 𓂡 𓁘 (HkAx Hn'-k: Heka hen-ek = May Heka (magic) be with you)")

welcome_version()
