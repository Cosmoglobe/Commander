#---------------------------------------------
import sys
import os
import numpy as np
import healpy as hp        
import h5py
import math
import time
import argparse
import re
from pathlib import Path
# Getting full path to Mathew's library as an object
commander_tools_path = Path(__file__).absolute().parents[2].joinpath('python','commander_tools').resolve()
# Appending the path to PYTHONPATH
sys.path.append(str(commander_tools_path))
# Importing necessary modules from Mathew's library 
from tod_tools import commander_tod as tod
#---------------------------------------------


'''
The `commander_tod(self, outPath, name, version=None, dicts=None, overwrite=False)` needs to be instantiated.
The `init_file(self, freq, od, mode='r')` is what you need to create a file
The `add_field(self, fieldName, data, compression=None)` is what you need to create a new field inside the file.
'''
def main():
    print("Script has started")
    parser = argparse.ArgumentParser(description='Adjust Level3 QUIET data for Commander3')
    # 
    make_od()

def make_od():
    # 
    level3_dir = Path('/mn/stornext/d16/cmbco/bp/maksym/quiet/data/Q/ces/patch_gc') 
    output_dir = Path('/mn/stornext/d16/cmbco/bp/maksym/quiet/data/Q/ces/patch_gc/output')
    if not Path.is_dir(output_dir):
        Path.mkdir(output_dir)
    version = '0.0.1'
    # Initialising tod object
    comm_tod = tod.commander_tod(output_dir, 'QUIET', version, dicts=None, overwrite=False)
    # Creating new file
    comm_tod.init_file(1, 1, mode='w')
    # Getting file names inside specified directory and removing the path component
    level3_data_files = sorted(level3_dir.rglob('*.hdf')) 
    level3_data_files = [data_file.name for data_file in level3_data_files]
    # Retrieving CES values from the file names 
    compiled_pattern = re.compile('[\d]')
    level3_ces = [int("".join(compiled_pattern.findall(data_file))) for data_file in level3_data_files] 
    #

if __name__ == '__main__':
    start_time = time.time()
    main()
    end_time = time.time()
    run_time = end_time - start_time
    print(f"Script run time: {run_time}s")



