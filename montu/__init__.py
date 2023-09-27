import numpy as np
import spiceypy as spy
import os
import requests
from tqdm import tqdm 

# Version
from montu.version import *

# Global variables
KERNELS = {
    'latest_leapseconds.tls':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls',
    'de441_part-1.bsp':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp',
    'pck00011.tpc':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc',
    'de440.bsp':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp',
}

class Montu(object):
    # Load kernels
    def load_kernels(dir='/tmp/'):
        # Check if dir exists
        if not os.path.exists(dir):
            os.system(f"mkdir -p {dir}")    
        for kernel,item in KERNELS.items():
            kernel_path = dir+"/"+kernel
            if not os.path.exists(kernel_path):
                print(f"Downloading '{kernel}'...")
                Montu._wget(item,kernel_path)
            else:
                print(f"Loading kernel {kernel}")
                spy.furnsh(kernel_path)

    # Util routines
    def _wget(url, file_name):
        """
        Source: ChatGPT
        """
        response = requests.get(url, stream=True)
        total_size = int(response.headers.get('content-length', 0))

        # Initialize the progress bar
        progress_bar = tqdm(total=total_size, unit='B', unit_scale=True)

        with open(file_name, 'wb') as file:
            for data in response.iter_content(chunk_size=1024):
                file.write(data)
                progress_bar.update(len(data))

        progress_bar.close()

    def _data_path(filename):
        """
        Get the full path of the `datafile` which is one of the datafiles provided with the package.
        
        Parameters:
            filename: Name of the data file, string.
            
        Return:
            Full path to package datafile in the python environment.
            
        """
        return os.path.join(os.path.dirname(__file__),'data',filename);