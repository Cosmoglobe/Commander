import sys
import os
import numpy as np
import healpy as hp        
import h5py
import math
import time
import pathlib
# Getting full path to Mathew's library as object
commander_tools_path = pathlib.Path(__file__).absolute().parents[2].joinpath('python','commander_tools').resolve()
# Appending the path to PYTHONPATH
sys.path.append(str(commander_tools_path))
# Importing necessary modules from Mathew's library 
from tod_tools import commander_tod as tod


'''
 The `add_field(self, fieldName, data, compression=None)` is what you need to create a new field inside the file.
'''
def main():
    pass


def make_od():
    pass

if __name__ == '__main__':
    start_time = time.time()
    main()
    end_time = time.time()
    run_time = end_time - start_time
    print(f"Script run time: {run_time}s")



