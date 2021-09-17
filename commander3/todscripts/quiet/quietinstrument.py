#---------------------------------------------
import numpy as np
from scipy.integrate import trapz
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
Beam files etc. for QUIET:
https://lambda.gsfc.nasa.gov/product/quiet/quiet_band_get.cfm
Similar stuff for WMAP for comparison:
https://lambda.gsfc.nasa.gov/product/map/dr5/bandpass_get.cfm

QUIET product description (on Lambda):
This product contains ancillary data for QUIET's 1st+2nd season galactic fields data release. The first set of ancillary data are the beam files. One is labeled 'instrumental' and the other is labelel 'effective'. The instrumental beam corresponds to the instantaneous beam profile of the instrument and is useful for calibration, observations of planets, etc. The effective beam corresponds to the beam after map-making. This is different from the instrumental beam due to pointing uncertainties. The telescope mount is slightly flexible, so it produces an additional smoothing term when averaging over pointing uncertainties and must be accounted for in the final map and any derived products such as the power spectrum estimation.

While the Q-band bandpass profile corresponds to a properly co-added average over detectors, the W-band profile is simply a representative profile for a typical detector. It therefore has an additional uncertainty.

From this presentation we get the FWHM:
https://www-group.slac.stanford.edu/ais/publicDocs/presentation123.pdf
FWHM ~ 28arcmin @ Q-band
FWHM ~ 13arcmin @ W-band
http://bccp.berkeley.edu/beach_program/presentations12/Buder.pdf
FWHM ~ 27arcmin @ Q-band
FWHM ~ 12arcmin @ W-band


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
"""
QUIET just used bl, so there will no be blm

- sllmax = slmmax <= for now at least
- beam/
    |- B <= some data from Lambda
    |- E <= some stuff from Lambda
    |- T <= some data from Lambda
  The T, E and B will be the same (technically incorrect but will do for now).
- beamlamx = 0
- beammmax = 0 #np.zeros(len(beamlmax))
- bandpass and bandpassx are the same for all detectors.
"""
"""
RIMO -- reduced instrumental model
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
    output_dir = Path('/mn/stornext/d16/cmbco/bp/maksym/quiet/data/Q/ces/patch_gc/output')
    if not Path.is_dir(output_dir):
        Path.mkdir(output_dir)
    version = "0.0.1" #np.string_('0.0.1')
    freqs = ['Q']
    #---------------------------------------------
    # Retrieving data
    #---------------------------------------------
    # Read-in the Bandpass data
    bandpass_file = "quiet_qband_bandpass_v1.txt"
    #with open("quiet_qband_bandpass_v1.txt") as bandpass_file:
    #    bandpass_data = bandpass_file.readlines()
    nu, response = np.loadtxt(bandpass_file).T
    #for n, r in nu, response:
    #print(f"{nu}, {response}")
    #---------------------------------------------
    instrument_filename = "QUIET_instrument_v" + str(version) + ".h5"
    # Instantiating Commander instrument Class
    cinst = comm_inst.commander_instrument(output_dir, instrument_filename, version, mode="w")

    # TODO: add frequency loop here
    # Bandpass and Bandpassx for each frequency
    cinst.add_bandpass(freqs[0], nu, response)
    # TODO: add the detector loop here 
    det = "01"
    # Bandpass and Bandpassx for each detector 
    cinst.add_bandpass(det, nu, response)
    # B_l
    beam = fits.open("quiet_qband_beam_instrumental_v1.fits")
    beam_temperature = beam[1].data['temperature']
    #print(beam[1].header)
    #print(beam[1].data['temperature'])

    cinst.add_field(det+"/beam/E", data=beam_temperature)
    cinst.add_field(det+"/beam/B", data=beam_temperature)
    cinst.add_field(det+"/beam/T", data=beam_temperature)
    #
    beammmax = 0
    cinst.add_field(det+"/beammmax", data=beammmax)
    beamlmax = len(beam_temperature) - 1 
    cinst.add_field(det+"/beamlmax", data=beamlmax)
    # The central frequency -- we integrate (or do numerical intergral solving)
    # via all frequencies
    cent_freq = trapz(nu*response, nu)/trapz(response, nu)
    cinst.add_field(det+"/centFreq", data=cent_freq)
    # FWHM -- look above
    # ...
    # Eliptisity -- assume beams are circular
    cinst.add_field(det+"/elip", data=0)
    # What fraction of the beam will be included in the map?
    cinst.add_field(det+"/mbeam_eff", data=1)
    # Rotation angle relative to each detector <= Planck specific, so set to 0
    cinst.add_field(det+"/psi_ell", data=0)
    # Sidelobes
    sl = np.zeros_like(beam_temperature)
    cinst.add_field(det+"/sl/E", data=sl)
    cinst.add_field(det+"/sl/B", data=sl)
    cinst.add_field(det+"/sl/T", data=sl)
    sllmax = beamlmax
    cinst.add_field(det+"/sllmax", data=sllmax)
    slmmax = beammmax
    cinst.add_field(det+"/slmmax", data=slmmax)
    # common/version
    cinst.add_field("common/version", data=np.string_(version))




def make_instrument():
    """
    Method for processing instrumental file
    """
    pass
    

if __name__ == '__main__':
    start_time = time.time()
    print("Script has started!")
    main()
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Script run time: {total_time:.2f} s")
    #filepath = "/mn/stornext/u3/hke/quiet_data/auxilliary"
    #data = fits.open(f'{filepath}/quiet_qband_temp_beam_140910.fits')
    #print(data[0].header)
    #print(data[1].header)



