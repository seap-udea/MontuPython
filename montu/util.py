###############################################################
# Montu interdependencies
###############################################################
import montu

###############################################################
# Required packages
###############################################################
import inspect
import os
import requests
import tqdm

import numpy as np
import spiceypy as spy
from tabulate import tabulate

###############################################################
# Module constants
###############################################################
BASIC_KERNELS = {
    'naif0012.tls':'',
    'frame.tk':'',
    'pck00011.tpc':'',
    'earth_assoc_itrf93.tf':''
}
PRECISION_KERNELS = {
    'latest_leapseconds.tls':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls',
    'de441_part-1.bsp':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp',
    'de441_part-2.bsp':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp',
}
KERNELS_LOADED = dict()

from tqdm import tqdm
def GENERATOR():
    """This routine is intended to create a while True loop for tqdm 
    counter
    """
    while True:yield
WHILE_TRUE = lambda:tqdm(GENERATOR())
PROGRESS = lambda iterable:tqdm(iterable)

###############################################################
# Montu Python Util Class
###############################################################
class Util(object):

    def vprint(verbose,*args):
        """Print messages in verbose mode
        """
        if verbose:
            print(*args)

    def arange(start, stop, step=1, endpoint=True):
        """Same as np.arange but including endpoint

        Source: https://stackoverflow.com/a/68551927
        """
        arr = np.arange(start, stop, step)
        if endpoint and arr[-1]+step==stop:
            arr = np.concatenate([arr,[stop]])
        return arr
    
    def print_df(df):
        """Print DataFrame.
        
        Parameters:
            df: Pandas DataFrame:
                DataFrame to print.
        """
        from IPython.display import display,HTML
        display(HTML(df.to_html()))

    def table_df(df,format='github'):
        """Present a DataFrame in a tabular form

        format: string, default = 'github':
            Format of the table.

            Other formats: “plain”,“simple”,“github”,“grid”,“fancy_grid”,“pipe”,
            “orgtbl”,“jira”,“presto”,“pretty”,“psql”,“rst”,“mediawiki”,“moinmoin”,
            “youtrack”,“html”,“latex”,“latex_raw”,“latex_booktabs”,“textile”, 
        """
        print(tabulate(df,headers='keys',tablefmt=format))

    def dt2cal(dt,bce=False):
        """Convert array of datetime64 to a calendar array of year, month, day, hour,
        minute, seconds, microsecond with these quantites indexed on the last axis.

        Parameters
            dt : datetime64 array (...)
                numpy.ndarray of datetimes of arbitrary shape

        Returns
            cal : uint32 array (..., 7)
            
                calendar array with last axis representing year, month, day, hour,
                minute, second, microsecond

        Adapted from: https://stackoverflow.com/a/56260054
        """
        out = np.empty(dt.shape + (7,), dtype="u4")
        Y, M, D, h, m, s = [dt.astype(f"M8[{x}]") for x in "YMDhms"]
        out[..., 0] = Y + 1970 # Gregorian Year
        
        out[..., 1] = (M - Y) + 1 # month
        out[..., 2] = (D - M) + 1 # dat
        out[..., 3] = (dt - D).astype("m8[h]") # hour
        out[..., 4] = (dt - h).astype("m8[m]") # minute
        out[..., 5] = (dt - m).astype("m8[s]") # second
        out[..., 6] = (dt - s).astype("m8[us]") # microsecond

        #out = np.array([float(o) for o in out])
        out = [int(o) for o in out]
        out = [-1]+out if bce else [1]+out
        return out

    def load_kernels(kernels=BASIC_KERNELS,dir='montmp/',verbose=False):
        
        # Check if dir exists
        if not os.path.exists(dir):
            os.system(f"mkdir -p {dir}")

        # Load kernel
        for kernel,item in kernels.items():

            #
            if kernel in KERNELS_LOADED.keys():
                Util.vprint(verbose,f"Kernel {kernel} already loaded, skipping")
                continue

            # Local kernel
            if len(item) == 0:
                if verbose:print(f"Loading local kernel {kernel}")
                kernel_file = Util._data_path(kernel)
                if os.path.isfile(kernel_file):
                    spy.furnsh(kernel_file)
                else:
                    raise AssertionError(f"Kernel file '{kernel}' not found in data directory")
                KERNELS_LOADED[kernel] = True
                continue

            # Remote kernel
            kernel_path = dir+"/"+kernel
            
            if not os.path.exists(kernel_path):
                # Download kernel if it is not yet downloaded
                if verbose:print(f"Downloading '{kernel}'...")
                Util._wget(item,kernel_path)
            
            # Once downloaded furnish kernel
            if verbose:print(f"Loading kernel {kernel}")
            spy.furnsh(kernel_path)
            KERNELS_LOADED[kernel] = True

    def dec2hex(dec,string=True):

        dec = float(dec)
        sgn = np.sign(dec)
        dec = abs(dec)
        h = int(dec)
        mf = 60*(dec - int(dec))
        m = int(mf)
        s = 60*(mf - m)
        if string:
            ret = f"{int(sgn*h):02d}:{int(m):02d}:{s:.3f}"
        else:
            ret = sgn*h,m,s
        return ret

    def string_difference(string1, string2):
        """Calculate the difference between two strings
        """
        A = set(string1.split()) 
        B = set(string2.split()) 
        str_diff = A.symmetric_difference(B)
        isEmpty = (len(str_diff) == 0)
        return str_diff
    
    def haversine_distance(lat1, lon1, lat2, lon2):
        """Compute angular distance between two points
        """
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        return c
    
    def montu_mark(ax):
        """Add a water mark to a 2d or 3d plot.
        
        Parameters:
        
            ax: Class axes: 
                Axe where the watermark will be placed.
        """
        #Get the height of axe
        axh=ax.get_window_extent().transformed(ax.get_figure().dpi_scale_trans.inverted()).height
        fig_factor=axh/4
        
        #Options of the water mark
        args=dict(
            rotation=270,ha='left',va='top',
            transform=ax.transAxes,color='pink',fontsize=6*fig_factor,zorder=100
        )
        
        #Text of the water mark
        mark=f"MontuPython {montu.version}"
        
        #Choose the according to the fact it is a 2d or 3d plot
        try:
            ax.add_collection3d
            plt_text=ax.text2D
        except:
            plt_text=ax.text
            
        text=plt_text(1,1,mark,**args);
        return text

    def get_methods(my_class):
        """Get a list of the methods for class my_class
        """
        return sorted([member[0] for member in inspect.getmembers(my_class) if '__' not in member[0]])

    def _data_path(filename,check=False):
        """Get the full path of the `datafile` which is one of the datafiles provided with the package.
        
        Parameters:
            filename: Name of the data file, string.
            
        Return:
            Full path to package datafile in the python environment.
            
        """
        file_path = os.path.join(os.path.dirname(__file__),'data',filename)
        if check and (not os.path.isfile(file_path)):
            raise ValueError(f"File '{filename}' does not exist in data directory")
        return file_path

    def _wget(url, filename, verbose=False):
        """Get a file from a url and store it as filename

        Source: ChatGPT
        """
        response = requests.get(url, stream=True)
        total_size = int(response.headers.get('content-length', 0))
        Util.vprint(verbose,f"Downloading {filename} from {url} [size = {total_size}]")

        # Initialize the progress bar
        progress_bar = tqdm.tqdm(total=total_size, unit='B', unit_scale=True)

        with open(filename, 'wb') as file:
            for data in response.iter_content(chunk_size=1024):
                file.write(data)
                progress_bar.update(len(data))

        progress_bar.close()

    def _linear_map(mapped,observed):
        a = (observed[1]-observed[0])/(mapped[1]-mapped[0])
        b = observed[0] - a*mapped[0]
        map = lambda x:a*x+b
        return map

class Dictobj(object):
    """Convert a dictionary to an object

    Examples:
        ob = Dictobj(a=2,b=3)
        print(ob.a,ob.b)
        ob = Dictobj(dict=dict(a=2,b=3))
        print(ob.a,ob.b)
        ob = Dictobj(dict={'a':2,'b':3})
        print(ob.a,ob.b)
    """

    def __init__(self, **kwargs):
        if 'dict' in kwargs.keys():
            kwargs.update(kwargs['dict'])
        for key, value in kwargs.items():
            if key == 'dict':continue
            setattr(self, key, value)

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return self.__str__()
