#---------------------------------------------
import numpy as np
import healpy as hp
import h5py
import time
import math
import re
from tqdm import tqdm
from astropy.io import fits
from joblib import cpu_count, Parallel, delayed, parallel_backend
from pathlib import Path
# Getting full path to Mathew's library as an object
commander_tools_path = Path(__file__).absolute().parents[2].joinpath('python','commander_tools').resolve()
# Appending the path to `PYTHONPATH`, so no need to 
# modify it externally (in your `.bashrc` etc.)
import sys
sys.path.append(str(commander_tools_path))
# Importing necessary modules from Mathew's library 
from tod_tools import commander_instrument as comm_inst
#---------------------------------------------
"""
This script creates QUIET instrument file and 
places it inside `output` folder. It uses com-
mon Commander3 structure to produce a file:

<file name>/
           |-<freq 1>/
                     |- bandpass
                     |- bandpassx
           |-<freq 2>
           | ...
           |-<freq N>
           |-<det 1>/
                    |- bandpass
                    |- bandpassx
                    |- beam/
                           |- B
                           |- E
                           |- T
                    |- beamlmax                 Dataset {1}
                    |- beammmax                 Dataset {1}
                    |- centFreq                 Dataset {1}
                    |- elip                     Dataset {1}
                    |- fwhm                     Dataset {1}
                    |- mbeam_eff                Dataset {1}
                    |- psi_ell                  Dataset {1}
                    |- sl                       Group
                    |- sllmax                   Dataset {1}
                    |- slmmax                   Dataset {1}
           |-<det 2>
           |  ...
           |-<det M>
           |- common

where
<freq 1> = Q
<freq 2> = W
"""
#---------------------------------------------
"""
We need beam, bandpass and sielobes (for each diode?). 

The instrument file should be one file with all the information below (LFI example)

top level -- 30, 44, 70 <= frequencies for each of these there is a bandpass
030                      Group
044/                     Group
    bandpass                 Dataset {319}
    bandpassx                Dataset {319}
070                      Group
18M                      Group
18S                      Group
19M                      Group
19S                      Group
20M                      Group
20S                      Group
21M                      Group
21S                      Group
22M                      Group
22S                      Group
23M                      Group
23S                      Group
24M                      Group
24S                      Group
25M                      Group
25S                      Group
26M                      Group
26S                      Group
27M                      Group
27S                      Group
28M/                     Group
    bandpass                 Dataset {310} <= response 
    bandpassx                Dataset {310} <= independent variable (frequency), array of all the freqs bandpass was measured at
    beam/                    Group <= a_lm representation of the beams
        B                        Dataset {9006001}
        E                        Dataset {9006001}
        T                        Dataset {9006001}
    beamlmax                 Dataset {1}
    beammmax                 Dataset {1}
    centFreq                 Dataset {1}
    elip                     Dataset {1}
    fwhm                     Dataset {1}
    mbeam_eff                Dataset {1}
    psi_ell                  Dataset {1}
    sl                       Group
    sllmax                   Dataset {1}
    slmmax                   Dataset {1}
28S                      Group
common/                  Group
    version

Q
W
<det1>
<det2>
...
common/
       version


Need to use `commander_instrument.py` with following fields:
    __init__(self, path, fileName, version, mode)
    add_field(self, fieldName, data)
"""
def main():
    """
    Main method of the script
    """
    version = np.string_('0.0.2')


if __name__ == '__main__':
    filepath = "/mn/stornext/u3/hke/quiet_data/auxilliary"
    data = fits.open(f'{filepath}/quiet_qband_temp_beam_140910.fits')
    print(data[0].header)
    print(data[1].header)



