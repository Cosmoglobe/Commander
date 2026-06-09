import glob
import numpy as np
import pickle
from tqdm import tqdm
from itertools import chain
from time import time
import astropy.units as u
from astropy.time import Time
from pathlib import Path
"""
This script is essentially an adaption of the script
    - Commander/commander3/python/commander_tools/tod_tools/akari_native_tod_reader.py
written by Eirik Gjerløw for AKARI. 
"""


# Valid band names for iras 
BANDS = [
    '012', #  12 micron / Band 1
    '025', #  25 micron / Band 2
    '060', #  60 micron / Band 3
    '100'  # 100 micron / Band 4
]

# Detectors in each band. 
BAND_DET_NUMBERS = {
    '012': [23, 24, 25, 26, 27, 28, 29, 30, 47, 48, 49, 50, 51, 52, 53, 54],    #  12 micron / Band 1
    '025': [16, 17, 18, 19, 20, 21, 22, 39, 40, 41, 42, 43, 44, 45, 46],        #  25 micron / Band 2
    '060': [8, 9, 10, 11, 12, 13, 14, 15, 31, 32, 33, 34, 35, 36, 37, 38],      #  60 micron / Band 3
    '100': [1, 2, 3, 4, 5, 6, 7, 55, 56, 57, 58, 59, 60, 61, 62],               # 100 micron / Band 4
}

# Kept for reference. Probably useful later.
DEAD_DETECTOR_LIST = [17, 20, 36]
DEGRADED_PERFORMANCE_DETECTORS = [25, 26, 28, 42] # See expl.suppl. Chapter IV.A.2

DEAD_DETECTORS = {}

# Compute dictionary with:
#   - keys: BANDS
#   - vals: list(str), with "IRAS_xxx-yy", where xxx=band and yy=det
# Detectors in "DEAD_DETECTOR_LIST" are omitted from BAND_DETS  
BAND_DETS = {}
for band in BANDS:
    DEAD_DETECTORS[band] = []
    BAND_DETS[band] = []
    for det_num in BAND_DET_NUMBERS[band]:
        if det_num in DEAD_DETECTOR_LIST:
            DEAD_DETECTORS[band].append(f'IRAS_{band}-{det_num:02d}')
        # elif det_num in DEGRADED_PERFORMANCE_DETECTORS:
            # DEAD_DETECTORS[band].append(f'IRAS_{band}-{det_num:02d}')
        else:
            BAND_DETS[band].append(f'IRAS_{band}-{det_num:02d}')

def get_chunk_min_length(*args):
    """
    Module-level function
    Used for running IRASTODReader._store_chunk_size_dictionary in parallel
    """
    ii, filelist_values = args
    chunk_lengths = np.empty(len(filelist_values), dtype=int)
    for jj, flist in enumerate(filelist_values):
        chunk_lengths[jj] = len(
            np.load(flist[ii], allow_pickle=True)["t"]
        )
    return ii + 1, np.min(chunk_lengths)


class IRASTODReader:
    """
    An IRASTODReader instance will, upon request, fetch a continuous sequence of data from the
    IRAS FITS TOD files, joined together across files. The data to be fetched is defined by the
    user through a formatting function.

    Arguments:

    iras_npz_dir (str): The directory in which the .npz tod files reside (should not be
        compressed)
    band (str): The band for which to load the TODs. Should be on of BANDS, i.e. '012', '025', '060' or '100'.
    band_dets:  Dictionary with list of detectors for each band. 
    load_filelist: If True (default), loads precomputed filelist and ntod for said band.  
    """
    def __init__(
            self, 
            iras_npz_dir, # /mn/stornext/d5/data/vetleav/IRAS/tod_data/IPAC_level1
            band,
            band_dets       = BAND_DETS,
            load_filelist   = True
            ):
        
        if band not in BANDS:
            raise ValueError(
                f"Unrecognized band '{band}'. Should be one of {BANDS}.")
        if type(band) is not str:
            raise TypeError
        
        self.data_dir   = iras_npz_dir # Location of TODs
        self.band       = band
        self.band_dets  = band_dets
        
        # Loc. of dicts containing filelists and and chunk sizes for each band.
        self.aux_dir    = Path("/mn/stornext/d5/data/vetleav/IRAS/tod_data/hdf_commander/aux_files") 

        if load_filelist:
            self.filelist   = self._load_det_filelist_dictionary_band(band)
            self.chunksizes = self._load_chunk_size_dictionary(band)
        else:
            self.filelist   = self.get_det_filelist()
            # self.chunksizes = self._load_chunk_size_dictionary(band)
            # self.chunksizes = self._store_chunk_size_dictionary(band)

    def _get_det_eff_solid_angles(self):
        # Get effective solid angle of IRAS detectors
        # From Table IV.A.1 in IRAS expl. suppl.
        table_file = "/mn/stornext/d5/data/vetleav/IRAS/aux_files/solid_angles_detectors.txt"
        data = np.loadtxt(
            fname=table_file,
            # delimiter=" ",
            # dtype=(int, int, float),
            usecols=(0,1,2),
            unpack=True
        )
        band_numbers = data[0].astype(int) 
        det_numbers  = data[1].astype(int) 
        solid_angles = data[2].astype(float) 

        det_area_dict = {}

        # band_number = np.a
        # for i in range(len(band_number)):
        for band_num, det_num, sol_ang in zip(band_numbers, det_numbers, solid_angles): 
            band = str(band_num).zfill(3)
            # print(band_number[i], det_number[i], solid_angle[i])
            band_det_label = f'IRAS_{band}-{det_num:02d}'
            if band_det_label not in self.band_dets[band]:
                raise RuntimeError("Wth??")
            det_area_dict[band_det_label] = sol_ang * 1e-7 * u.sr
        return det_area_dict

    def store_flux2MJy_sr_conversion_factors(self):
        # Store dictionary containing each detector's multiplicative factor 
        # to convert tod from W/m^2 to MJy/sr.
        # Returns dict with
        #   key: detector name [str]
        #   val: conversion factor [float] 

        # Conversion is based on the effective solid angle of each detector,
        # according to Table IV.A.1 in IRAS expl. suppl. 
        outfile = f"{self.aux_dir}/flux2MJy_conversion_factors.pkl"
        det_area_dict = self._get_det_eff_solid_angles()
        conv_factor_dict = {}
        for band, det_list in self.band_dets.items():
            band_wl = float(band) * u.micron
            band_freq = band_wl.to(u.Hz, equivalencies=u.spectral())
            for det_name in det_list:
                # det_num = int(det_name.split("-")[-1])
                det_area = det_area_dict[det_name]
                conv_factor_dict[det_name] = (1 * u.W / u.m**2 / band_freq / det_area).to(u.MJy / u.sr).value

        with open(outfile, "wb") as pkl_file:
            pickle.dump(conv_factor_dict, pkl_file)
        return 

    def load_flux2MJy_sr_conversion_factors(self):
        # Load dictionary containing each detector's multiplicative factor 
        # to convert tod from W/m^2 to MJy/sr.
        # Returns dict with
        #   key: detector name [str]
        #   val: conversion factor [float] 

        # Conversion is based on the effective solid angle of each detector,
        # according to Table IV.A.1 in IRAS expl. suppl. 

        conv_factor_dict_fname = f"{self.aux_dir}/flux2MJy_conversion_factors.pkl"
        # if not Path(conv_factor_dict_fname).exists():

        with open(conv_factor_dict_fname, "rb") as pkl_file:
            conv_factor_dict = pickle.load(pkl_file)
        return conv_factor_dict


    def get_det_filelist(self):
        """
        Return dict of npz files
         key: detector number [int] 
         val: sorted list of all npz files for said detector.
        This dictionary is used by the adapter to load the correct files.
        """
        filelist_dict = {}
        for det_name in tqdm(self.band_dets[self.band]):
            det_num =  int(det_name.split("-")[-1])
            filelist_dict[det_name] = sorted(glob.glob(f"{self.data_dir}/*/*/det{det_num:02}_continuous.npz"))
            if len(filelist_dict[det_name]) != 5787:
                raise RuntimeError(
                    f"""
                    ERROR! Found {len(filelist_dict[det_name])} files for detector {det_name}.
                    Only valid number is 5787. Aborting"""
                )
        return filelist_dict

    def _store_det_filelist_dictionary(self):
        """
        Store all files for a given band as a dictionary.
        - keys: det_name, e.g. "IRAS_060-XX" for detectors in band 60.
        - vals: list[filepaths], where filepaths is a sorted list of paths to
                all 5878 filenames for a given detector.  
        The stored dictionaries are loaded by the "_load_det_filelist_dictionary_band" method. 
        """
        outdir = self.aux_dir
        for band, det_name_list in tqdm(self.band_dets.items()):
            outfname_band = f"{outdir}/filelist_{band}.pkl"
            if Path(outfname_band).exists:
                # Prevent overwriting
                print(f"{outfname_band} already exists. Skipping")
                continue
            filelist_band = {}
            for det_name in tqdm(det_name_list):
                det_num =  int(det_name.split("-")[-1])
                # Get sorted list of all files for a given detector. 
                filelist_band[det_name] = sorted(glob.glob(f"{self.data_dir}/*/*/det{det_num:02}_continuous.npz"))

                if len(filelist_band[det_name]) != 5787:
                    # All detectors in all bands should have 5787 associated npz files. 
                    raise RuntimeError(
                        f"""
                        ERROR! Found {len(filelist_band[det_name])} files for detector {det_name}.
                        Only valid number is 5787. Aborting"""
                    )
            # Store dictionary as pkl file. 
            with open(outfname_band, "wb") as pkl_file:
                pickle.dump(filelist_band, pkl_file)
        return

    def _load_det_filelist_dictionary_band(self, band):
        """
        Return dictionray containing all npz files for a given band
        computed by the "_store_det_filelist_dictionary_band" method. 
        - keys: det_name, e.g. "IRAS_060-XX" for detectors in band 60.
        - vals: list[filepaths], where filepaths is a sorted list of paths to
                all 5878 filenames for a given detector.  
        
        """
        filelist_fname = f"{self.aux_dir}/filelist_{band}.pkl"
        with open(filelist_fname, "rb") as pkl_file:
            filelist = pickle.load(pkl_file)
        return filelist

    def _store_chunk_size_dictionary(self, band, num_workers=None):
        """
        In bands 060 and 100, some chunks (files) have a 
        varying number of tod samples between detectors, causing
        errors with huffman compression when creating h5 files 
        due to varying file array sizes. 
        Update: A few samples in bands 012 and 025 had this issue, causing
        errors in chunks 201-300 for these bands. 

        Here, we go through all chunks in all detectors, and store 
        the *lowest* number of tod samples for a given  chunk in each band. 
        Tods are then loaded as data[:min_ntod].
        For bands 012 and 025 we set min_ntod=-1,
        as all chunks have equal ntod. 
        Method is brute force and inefficient, but only needs to be run once.  

        Store min_ntod dictionary for each band as pickle file, with
        - keys: chunk_idx, integer in [1,2,...,5786,5787]
        - vals: lowest ntod for the chunk across detectors in band.   
        Output: chunk_sizes_{band}.pkl
        """

        if band not in BANDS:
            raise ValueError
        outfile = f"{self.aux_dir}/chunk_sizes_{band:03}.pkl"
        if Path(outfile).exists():
            # Prevent overwriting
            raise FileExistsError(outfile)
        
        filelist = self._load_det_filelist_dictionary_band(band)

        print(f"Computing min chunk sizes for band {band}")
        N_chunks            = len(list(filelist.values())[0]) # 5787
        
        if num_workers is None:
            # Slooow 
            chunksizes_dict    = {}
            N_dets          = len(filelist.keys())
            chunk_lengths   = np.zeros(N_dets, dtype=int)
            for ii in tqdm(range(N_chunks)):
                # Iterate over all chunks
                chunk_idx = ii + 1
                for jj, flist in enumerate(filelist.values()):
                    # Get number of tods in all detectors for said chunk. 
                    chunk_lengths[jj] = len(np.load(flist[ii], allow_pickle=True)["t"])
                # Find lowest ntod for chunk.
                chunksizes_dict[chunk_idx] = np.min(chunk_lengths)
        else:
            from concurrent.futures import ProcessPoolExecutor, as_completed
            filelist_values = list(filelist.values())
            # OLD METHOD. New method adds progress bar
            # with ProcessPoolExecutor(max_workers=num_workers) as executor:
            #     results = executor.map(
            #         get_chunk_min_length,
            #         [(ii, filelist_values) for ii in range(N_chunks)]
            #     )
            chunksizes_dict = {}

            with ProcessPoolExecutor(max_workers=num_workers) as executor:
                futures = {
                    executor.submit(
                        get_chunk_min_length,
                        ii,
                        filelist_values,
                    ): ii
                    for ii in range(N_chunks)
                }

                for future in tqdm(as_completed(futures), total=N_chunks):
                    chunk_idx, min_length = future.result()
                    chunksizes_dict[chunk_idx] = min_length

        print(f"Done! Storing to {outfile}")
        with open(outfile, "wb") as pkl_file:
            pickle.dump(chunksizes_dict, pkl_file)

    def _load_chunk_size_dictionary(self, band):
        """
        Loads number of tods to use in every chunk for a given band. 
        See docstring in "_store_chunk_size_dictionary" for details. 
        """
        outdir = "/mn/stornext/d5/data/vetleav/IRAS/tod_data"
        outfile = f"{self.aux_dir}/chunk_sizes_{band:03}.pkl"
        with open(outfile, "rb") as pkl_file:
            chunksizes_dict = pickle.load(pkl_file)
        return chunksizes_dict




if __name__ == '__main__':
    TOD_DIR = "/mn/stornext/d5/data/vetleav/IRAS/tod_data/IPAC_level1"

    ITR = IRASTODReader(
        iras_npz_dir    = TOD_DIR,
        band            = band,
        load_filelist   = True
        )
    
    # ITR._get_det_eff_solid_angles()
    # ITR.store_flux2MJy_sr_conversion_factors()
    # ITR.test()
    import os
    print(f"{os.cpu_count()=}")
    # ITR._store_chunk_size_dictionary(band="012")
    # print()
    # ITR._store_chunk_size_dictionary(band="025", num_workers=64)
