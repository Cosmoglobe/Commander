import glob
import fitsio
import numpy as np
import pickle
from tqdm import tqdm
from itertools import chain
from time import time
import astropy.units as u
from astropy.time import Time
"""
This script is essentially an adaption of the script
    - Commander/commander3/python/commander_tools/tod_tools/akari_native_tod_reader.py
written by Eirik Gjerløw for AKARI. 
"""


class IRASTODReader:
    """
    An IRASTODReader instance will, upon request, fetch a continuous sequence of data from the
    IRAS FITS TOD files, joined together across files. The data to be fetched is defined by the
    user through a formatting function.

    Arguments:

    iras_npz_dir (str): The directory in which the .npz tod files reside (should not be
        compressed)
    band (str): The band for which to load the TODs. Should be '012', '025', '060' or '100'.
    fits2output_formatter (function): A function defining the output format of the TOD. The function
        signature should be (f, f_start, f_end, band), where f is an open fitsio.FITS instance
        pointing to one of the TOD files, f_start is the starting index of the current file, f_end
        is the ending index of the current file, and band is the same type of string as above. The
        function should then fetch the desired data corresponding to the range [f_start:f_end] (in
        Python slicing format), as well as any other data desired, in a dictionary. The dictionary
        should not contain sub-dictionaries. Each entry should either be a single scalar or a numpy
        array of length (f_end - f_start).
    load_idx_file_mapping (bool): The TODReader creates an internal mapping between files and the
        global TOD index represented by that file. This option allows the user to load this mapping
        from an existing .pkl file rather than to compute them. The directory from which to load the
        file is given by the 'mapping_dir' argument.
    save_idx_file_mapping(bool): Whether to save the internal mapping to file (not necessary to use
        together with 'load_idx_file_mapping'). The 'mapping_dir' argument is where the file will be
        saved.
    mapping_dir (str): The directory in which a potential mapping file will be saved.
    """
    def __init__(
            self, 
            iras_npz_dir, # /mn/stornext/d5/data/vetleav/IRAS/tod_data/IPAC_level1
            band, 
            band_dets,
            testing = False,
            ):
        
        
        # self.BAND_DET_NUMBERS = {
        #     '012': [23, 24, 25, 26, 27, 28, 29, 30, 47, 48, 49, 50, 51, 52, 53, 54],    #  12 micron / Band 1
        #     ### All detectors in 25 and 60
        #     # '025': [16, 17, 18, 19, 20, 21, 22, 39, 40, 41, 42, 43, 44, 45, 46],        #  25 micron / Band 2
        #     # '060': [8, 9, 10, 11, 12, 13, 14, 15, 31, 32, 33, 34, 35, 36, 37, 38],      #  60 micron / Band 3
        #     ### Dead detectors in 25 and 60 omitted, 17&20 and 36, respectively
        #     '025': [16, 18, 19, 21, 22, 39, 40, 41, 42, 43, 44, 45, 46],        #  25 micron / Band 2
        #     '060': [8, 9, 10, 11, 12, 13, 14, 15, 31, 32, 33, 34, 35, 37, 38],      #  60 micron / Band 3
        #     '100': [1, 2, 3, 4, 5, 6, 7, 55, 56, 57, 58, 59, 60, 61, 62],               # 100 micron / Band 4
        # }
        # self.band_det_numbers = {}
        # for band, det_name_list in band_dets.items():
        #     self.band_det_numbers[band] = [int(det_name.split("-")[-1]) for det_name in det_name_list]

        self.det_name2num = {}
        for band, det_name_list in band_dets.items():
            for det_name in det_name_list:
                self.det_name2num[det_name] = int(det_name.split("-")[-1])# for det_name in det_name_list]

        self.det_num2name = {
            v: k for k, v in self.det_name2num.items()
        }

        # for k,v in self.det_name2num.items():
            # self.det_num2name[v] = k
        # print(self.det_name2num)
        # print()
        # print(self.det_num2name)

        # exit()
        # for k,v in self.band_det_numbers.items():
            # print(f"{k:2}, {v}")
        # print(self.band_det_numbers)
        # exit()

        valid_bands = ['012','025','060','100']
        if type(band) in (int, float): 
            # Allow for band to be specified as integer or float 
            band = f"{band:03.0f}"
        if band not in valid_bands:
            raise ValueError(
                f"Unrecognized band '{band}'. Should be one of {valid_bands}.")
        
        self.band_dets  = band_dets
        self.band       = band
        self.data_dir   = iras_npz_dir
        self.testing    = testing

        ### Should probably move to other script
        self.detector_areas = self.get_detector_area_dict()
        self.detector_freqs = self.get_detector_freq_dict()
        self.flux_to_MJy_sr = self.get_flux_to_MJy_sr_conversion_factor_dict()

        self.filelist   = self.get_det_filelist()


    def test_data(self):
        # file = glob.glob(f"{self.data_dir}/sop543/obs013/det30_continuous.npz")#[0]
        file = np.load(f"{self.data_dir}/sop534/obs013/det30_continuous.npz")
        ### KEYS: t, ra, dec, flux, n_merged, dedup_seconds, sop, obs, det
        tt = file["t"]
        self.t0 = Time('1981-01-01', scale='utc')
        self.epoch_1981      = Time("1981-01-01T00:00:00", scale="utc") # IRAS's t=0
        print('test data exiting')
        exit()

    def get_filelist(self):
        """
        MAYBE OUTDATED
        all_files        = glob.glob(f"{self.data_dir}/*/*/*.npz")
        det_numbers = self.band_det_numbers[self.band]
        if self.testing: 
            # Reduce total number of files to improve speed.  
            all_files   = all_files[::50]
            det_numbers = det_numbers[::2]

        # file_identifiers = set(map(str, [f"det{i:02}_continuous.npz" for i in det_numbers]))
        file_identifiers = set([f"det{i:02}_continuous.npz" for i in det_numbers])

        self.filelist = [
            f for f in all_files
            if any(i in f for i in file_identifiers)
        ]
        if not self.testing:
            # Make sure that number of files found is correct.
            NUM_FILES_PER_BAND = {
                # Numbers below were calculated with the bash script 'self.data_dir/check_nfiles.sh'
                # by using 'find' while manually looping over detectors.
                '012': 92592, # 16 detectors * 5787 files/det
                '025': 75231, # 13 detectors * 5787 files/det (2/15 detectors dead)
                '060': 86805, # 15 detectors * 5787 files/det (1/16 detectors dead)
                '100': 86805, # 15 detectors * 5787 files/det
            }
            err_msg = f"Error! Number of files found ({len(self.filelist)})" 
            err_msg += f" not matching expected values ({NUM_FILES_PER_BAND[self.band]})"
            assert NUM_FILES_PER_BAND[self.band] == len(self.filelist), err_msg
        
        return sorted(self.filelist)
        """
        raise NotImplementedError


    def get_file_index_range(self, filename):
        # OBSOLETE???
        file_idx = self.filelist.index(filename)
        print('get_file_index_range exiting')

        # exit()
        raise NotImplementedError

    def get_det_filelist(self):
        # Return dict of npz files
        #  key: detector number [int] 
        #  val: sorted list of all npz files for said detector.
        # This dictionary is used by the adapter to load the correct files.
        filelist_dir = {}
        # for det_num in tqdm(self.band_det_numbers[self.band]):
        for det_name in tqdm(self.band_dets[self.band]):
            # filelist_dir[f"{det_num:03d}"] = sorted(glob.glob(f"{self.data_dir}/*/*/det{det_num}_continuous.npz"))
            det_num =  self.det_name2num[det_name]
            filelist_dir[det_name] = sorted(glob.glob(f"{self.data_dir}/*/*/det{det_num}_continuous.npz"))
        return filelist_dir


    def get_detector_area_dict(self):
        # Return dictionary containing area of each detector: 
        #   key: detector number (int)
        #   val: Area in steradian. (Astropy quantity.) 
        det_file = "/mn/stornext/d5/data/vetleav/IRAS/detector_area.npy"
        det_area = np.load(det_file, allow_pickle=True).item()

        # print(det_area)
        for det_num, det_name in self.det_num2name.items():
        # for det_num in det_area.keys():
            # det_area[self.det_num2name[det_num]] = det_area.pop(det_num)
            det_area[det_name] = det_area.pop(det_num)

        return det_area
    
    def get_detector_freq_dict(self):
        # Return dictionary containing the frequency of all detectors.
        # Keys of self.band_det_numbers are the obs. wavelength of the band in micron. 
        detector_freq_dict = {}
        # for band, det_list in self.band_det_numbers.items():
        for band, det_list in self.band_dets.items():
            band_wl = float(band) * u.micron
            for det_name in det_list:
                detector_freq_dict[det_name] = band_wl.to(u.Hz, equivalencies=u.spectral())
        return detector_freq_dict
    
    def get_flux_to_MJy_sr_conversion_factor_dict(self):
        # In IPAC level 1 data, TODs are given in units of W/m^2
        # This method creates the 'conv_factor_dict' dictionary
        # consisting of multiplicative factors to go from W/m^2 to MJy/sr 
        # TOD[det] (MJy/sr) = TOD[det] (W/m^2) * conv_factor_dict[det]
        conv_factor_dict = {}
        for det_name in self.detector_freqs.keys():
            det_freq = self.detector_freqs[det_name]
            det_area = self.detector_areas[det_name]
            conv_factor_dict[det_name] = (1 * u.W / u.m**2 / det_freq / det_area).to(u.MJy / u.sr).value
        return conv_factor_dict



if __name__ == '__main__':
    TOD_DIR = "/mn/stornext/d5/data/vetleav/IRAS/tod_data/IPAC_level1"
    import sys
    if len(sys.argv) == 2:
        band = f"{sys.argv[1].zfill(3)}"
    else:
        band = "100"

    TESTING = False
    # TESTING = True

    ITR = IRASTODReader(
        iras_npz_dir=TOD_DIR,
        band=band,
        testing=TESTING,
        )
        