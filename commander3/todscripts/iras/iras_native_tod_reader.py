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



BANDS = [
    '012',
    '025',
    '060',
    '100'
]

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

BAND_DETS = {}
for band in BANDS:
    DEAD_DETECTORS[band] = []
    BAND_DETS[band] = []
    for det_num in BAND_DET_NUMBERS[band]:
        if det_num in DEAD_DETECTOR_LIST:
            DEAD_DETECTORS[band].append(f'IRAS_{band}-{det_num:02d}')
        else:
            BAND_DETS[band].append(f'IRAS_{band}-{det_num:02d}')

class IRASTODReader:
    """
    An IRASTODReader instance will, upon request, fetch a continuous sequence of data from the
    IRAS FITS TOD files, joined together across files. The data to be fetched is defined by the
    user through a formatting function.

    Arguments:

    iras_npz_dir (str): The directory in which the .npz tod files reside (should not be
        compressed)
    band (str): The band for which to load the TODs. Should be '012', '025', '060' or '100'.
    """
    def __init__(
            self, 
            iras_npz_dir, # /mn/stornext/d5/data/vetleav/IRAS/tod_data/IPAC_level1
            band,
            band_dets = BAND_DETS,
            load_filelist = False
            ):
        
        if band not in BANDS:
            raise ValueError(
                f"Unrecognized band '{band}'. Should be one of {BANDS}.")
        
        self.data_dir   = iras_npz_dir
        self.band       = band
        self.band_dets  = band_dets


        if load_filelist:
            self.filelist   = self._load_det_filelist_dictionary_band(band)
            self.chunksizes = self._load_chunk_size_dictionary(band)
        else:
            self.filelist   = self.get_det_filelist()
            self.chunksizes = self._load_chunk_size_dictionary(band)


    def get_det_filelist(self):
        # Return dict of npz files
        #  key: detector number [int] 
        #  val: sorted list of all npz files for said detector.
        # This dictionary is used by the adapter to load the correct files.
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

        outdir = "/mn/stornext/d5/data/vetleav/IRAS/tod_data"
        outfname_all = f"{outdir}/filelist_all.pkl"
        filelist_dict = {}
        for band, det_name_list in tqdm(self.band_dets.items()):
            outfname_band = f"{outdir}/filelist_{band}.pkl"
            if Path(outfname_band).exists:
                print(f"{outfname_band} already exists. Skipping")
                continue
            filelist_band = {}
            for det_name in tqdm(det_name_list):
                det_num =  int(det_name.split("-")[-1])
                filelist_band[det_name] = sorted(glob.glob(f"{self.data_dir}/*/*/det{det_num:02}_continuous.npz"))
                if len(filelist_band[det_name]) != 5787:
                    raise RuntimeError(
                        f"""
                        ERROR! Found {len(filelist_dict[det_name])} files for detector {det_name}.
                        Only valid number is 5787. Aborting"""
                    )
            with open(f"{outdir}/filelist_{band}.pkl", "wb") as pkl_file:
                pickle.dump(filelist_band, pkl_file)
            filelist_dict[band] = filelist_band
        if Path(outfname_all).exists:
            print(f"{outfname_band} already exists. Returning")
            return

        with open(outfname_all, "wb") as pkl_file:
            pickle.dump(filelist_dict, pkl_file)

        return filelist_dict

    def _load_det_filelist_dictionary_band(self, band):
        outdir = "/mn/stornext/d5/data/vetleav/IRAS/tod_data"
        filelist_fname = f"{outdir}/filelist_{band}.pkl"
        with open(filelist_fname, "rb") as pkl_file:
            filelist = pickle.load(pkl_file)
        return filelist


    def _load_det_filelist_dictionary_full(self):
        ### OBSOLETE
        outdir = "/mn/stornext/d5/data/vetleav/IRAS/tod_data"
        filelist_fname = f"{outdir}/filelist_all.pkl"
        if not Path(filelist_fname).exists():
            raise FileNotFoundError
        
        with open(f"{outdir}/filelist_all.pkl", "rb") as pkl_file:
            self.filelist_all = pickle.load(pkl_file)
        return self.filelist_all


    def _store_chunk_size_dictionary(self, band):
        if band not in BANDS:
            raise ValueError
        outdir = "/mn/stornext/d5/data/vetleav/IRAS/tod_data"
        filelist_fname = f"{outdir}/filelist_{band:03}.pkl"
        outfile = f"{outdir}/chunk_sizes_{band:03}.pkl"
        if Path(outfile).exists():
            raise FileExistsError
        
        with open(filelist_fname, "rb") as pkl_file:
            filelist = pickle.load(pkl_file)

        print(f"Computing min chunk sizes for band {band}")
        
        N_chunks            = len(list(filelist.values())[0])
        chunksizes_dict    = {}
        if band == "012" or band == "025":
            chunksizes_dict    = {chunk_idx: -1 for chunk_idx in range(1, N_chunks + 1)}
            print(f"Done! Storing to {outfile}")
            with open(outfile, "wb") as pkl_file:
                pickle.dump(chunksizes_dict, pkl_file)
            return 
        
        N_dets          = len(filelist.keys())
        chunk_lengths   = np.zeros(N_dets, dtype=int)
        for ii in tqdm(range(N_chunks)):
        # for ii in range(N_chunks):
            chunk_idx = ii + 1
            for jj, flist in enumerate(filelist.values()):
                chunk_lengths[jj] = len(np.load(flist[ii], allow_pickle=True)["t"])

            chunksizes_dict[chunk_idx] = np.min(chunk_lengths)
        
        print(f"Done! Storing to {outfile}")
        with open(outfile, "wb") as pkl_file:
            pickle.dump(chunksizes_dict, pkl_file)

    def _load_chunk_size_dictionary(self, band):
        print("Loading")
        outdir = "/mn/stornext/d5/data/vetleav/IRAS/tod_data"
        outfile = f"{outdir}/chunk_sizes_{band:03}.pkl"
        with open(outfile, "rb") as pkl_file:
            chunksizes_dict = pickle.load(pkl_file)
        return chunksizes_dict




if __name__ == '__main__':
    TOD_DIR = "/mn/stornext/d5/data/vetleav/IRAS/tod_data/IPAC_level1"
    import sys
    if len(sys.argv) == 2:
        band = f"{sys.argv[1].zfill(3)}"
    else:
        band = "100"

    # TESTING = False
    TESTING = True

    ITR = IRASTODReader(
        iras_npz_dir=TOD_DIR,
        band=band,
        load_filelist=True
        )
    print("STORING")
    # ITR._store_chunk_size_dictionary("100")
    # ITR._store_chunk_size_dictionary("012")
    # ITR._store_chunk_size_dictionary("025")


    ch = ITR._load_chunk_size_dictionary("025")
    print("res:")
    keys = list(ch.keys())
    for k in keys[::100]:
    # for k,v in ch.items()[::100]:
        print(f"{k:4}, {ch[k]}")