
###############################################################
# Import montu modules
###############################################################
from montu.version import *
from montu.util import *
from montu.time import *
from montu.stars import *

###############################################################
# Aliases
###############################################################
D2H = Util.dec2hex
VPRINT = Util.vprint
PRINTDF = Util.print_df
TABLEDF = Util.table_df

###############################################################
# Import other modules
###############################################################
import warnings

###############################################################
# Initialization
###############################################################
# Avoid warnings
warnings.filterwarnings("ignore")

# Load basic SPICE kernels
Util.load_kernels()
