###############################################################
# Import montu modules
###############################################################
from montu.version import *
from montu.util import *
from montu.physics import *
from montu.time import *
from montu.observer import *
from montu.horizon import *
from montu.heka import *
from montu.sebau import *
from montu.stars import *
from montu.maps import *
from montu.phenomena import *

###############################################################
# Aliases
###############################################################
D2S = Util.dec2sex
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

# Load planetary data
ALLPLANETS = load_planets()

# Showing version 
def welcome_version():
    print(f"MontuPython version {version}. 𓇍𓇋𓇋𓏏𓅓𓊵 𓎛𓎡𓄿𓀭𓎛𓈖𓂝𓎡 (ii-ti m Htp, HkAx Hn'-k)")

def welcome_translate():
    print(f"𓇍𓇋𓇋𓏏𓅓𓊵 (ii-ti m Htp: Ii-ti em hotep = Welcome in peace) 𓎛𓎡𓄿𓀭𓎛𓈖𓂝𓎡 (HkAx Hn'-k: Heka hen-ek = May Heka (magic) be with you)")

welcome_version()
