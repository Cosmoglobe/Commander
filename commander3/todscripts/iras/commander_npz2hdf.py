from iras_commander_data_adapter import IRASCommanderDataAdapter
import iras_native_tod_reader

from cosmoglobe.tod_tools.commander_hdf_writer import CommanderHDFWriter
from time import time
from pathlib import Path
import os 

# Main script - run to create the IRAS HDF files
IRAS_NPZ_DIR    = "/mn/stornext/d5/data/vetleav/IRAS/tod_data/IPAC_level1"
OUTPATH         = '/mn/stornext/d5/data/vetleav/IRAS/tod_data/hdf_commander'
NSIDE           = 512
NUM_PROCESSES   = 64
OVERWRITE       = True

BANDS       = iras_native_tod_reader.BANDS
BAND_DETS   = iras_native_tod_reader.BAND_DETS
NUM_FILES_PER_SEGMENT = 100


def set_bands_and_band_dets(testing=False, bands=None, verbose=True):
    ### Set band and detectors to be used. 
    # Wrote this method for simpler debugging/testing 
    # If testting = False, uses all detectors in bands
    # If testing = True, uses bands/dets given explicitly below   
    band_dets = iras_native_tod_reader.BAND_DETS
    if not testing:
        if bands is None:
            bands = iras_native_tod_reader.BANDS
        else:
            band_dets_copy = band_dets.copy()
            band_dets = {band: band_dets_copy[band] for band in bands}
    else:
        # bands = ['025']
        # bands = ['060', '100']
        bands = ['060']
        # bands = ['100']
        # print("="*50, "\n","TEST RUN! Reducing the number of detectors.", end="\n"+"="*50+"\n")
        print("="*50)
        print("TEST RUN! Reducing the number of detectors.")
        print("="*50)

        bdets = band_dets.copy()
        for band, dets in bdets.items():
            if band in bands:
                # band_dets[band] = dets[:7]
                band_dets[band] = dets[3:5]
                # band_dets[band] = [dets[6]]#, dets[6]]#, dets[3]]
                # band_dets[band] = [dets[3]]#, dets[6]]#, dets[3]]
                # band_dets[band] = [dets[4], dets[6]]#, dets[3]]
                # band_dets[band] = [dets[3], dets[6]]#, dets[3]]
                # band_dets[band] = dets
                print(f"{band=}:")# {len(dets)} det -->> {len(band_dets[band])} detectors")
                det_nums = [int(det.split("-")[-1]) for det in band_dets[band]]
                print(f"  Ndets: {len(dets)} -->> {len(band_dets[band])}.  det = {det_nums}")
            else:
                band_dets.pop(band)
        

    if verbose:
        print()
        print("#"*50)
        print("Running hdf writer with the following parameters:")
        print(f"    {NSIDE                  = }")
        print(f"    {NUM_FILES_PER_SEGMENT  = }")
        print(f"    {NUM_PROCESSES          = }")
        print(f"    {bands                  = }")
        print(f"    band_dets              = ", end="\n      ")
        print(*[vv for vv in band_dets.items()], sep="\n      ")
        print("#"*50)
        print("\n"*2)
    return bands, band_dets


def run_single_band_write(version):
    # Store files for this version in separate subdir.  
    outpath = f"{OUTPATH}/v{version}"
    # Create outdir. 
    # Prevent unintentional overwriting.
    Path(outpath).mkdir(exist_ok=False) 

    N_dets = [len(BAND_DETS[band]) for band in BANDS]
    print(f"Storing h5 files for band {BANDS} with {N_dets} detectors")
    t0 = time()
    iras_comm_data_adapter = IRASCommanderDataAdapter(
        data_dir                = IRAS_NPZ_DIR, 
        nside                   = NSIDE, 
        version                 = version,
        num_files_per_segment   = NUM_FILES_PER_SEGMENT,
        bands                   = BANDS,
        band_dets               = BAND_DETS
        )
    

    # Write hdf files.
    comm_todwriter = CommanderHDFWriter(iras_comm_data_adapter)
    comm_todwriter.write_hdf_files(
        hdf_output_dir  = outpath,
        overwrite       = True, 
        num_processes   = NUM_PROCESSES,
        bands           = BANDS,
        verbosity       = 2
        )
    t1 = time()
    print(f"DONE! dur: {t1-t0:.3f} sec")

    

def run_multiple_band_write(version, outpath=None, num_processes=NUM_PROCESSES):
    # Store files for this version in separate subdir.  
    if outpath is None:
        outpath = f"{OUTPATH}/v{version}"
    # Create outdir. 
    # Prevent unintentional overwriting.
    exist_ok = True if "test" in outpath else False
    Path(outpath).mkdir(exist_ok=exist_ok) 

    N_dets = [len(BAND_DETS[band]) for band in BANDS]
    print(f"Storing h5 files for band {BANDS} with {N_dets} detectors")

    t0 = time()
    iras_comm_data_adapter = IRASCommanderDataAdapter(
        data_dir                = IRAS_NPZ_DIR, 
        nside                   = NSIDE, 
        version                 = version,
        num_files_per_segment   = NUM_FILES_PER_SEGMENT,
        bands                   = BANDS,
        band_dets               = BAND_DETS
        )

    # Write hdf files.
    comm_todwriter = CommanderHDFWriter(iras_comm_data_adapter)
    comm_todwriter.write_hdf_files(
        hdf_output_dir  = outpath, 
        overwrite       = True, 
        num_processes   = num_processes,
        bands           = BANDS,
        verbosity       = 2
        )
    t1 = time()
    print(f"DONE! dur: {t1-t0:.3f} sec")

def run_multiple_band_write_external(
        outpath,
        version, 
        nside,
        overwrite,
        num_processes           = NUM_PROCESSES,
        num_files_per_segment   = NUM_FILES_PER_SEGMENT,
        bands                   = BANDS,
        band_dets               = BAND_DETS,
        ):

    if Path(outpath).exists() and not overwrite:
        # Prevent overwrite when called from external script
        print(f"Output dir for version{version} already exists and overwrite=False. Aborting")
        return
    # Double protection in case of maximum stupidity.
    Path(outpath).mkdir(parents=True, exist_ok=False) 

    N_dets = [len(BAND_DETS[band]) for band in BANDS]
    print(f"Storing h5 files for band {BANDS} with {N_dets} detectors")

    t0 = time()
    iras_comm_data_adapter = IRASCommanderDataAdapter(
        data_dir                = IRAS_NPZ_DIR, 
        nside                   = nside, 
        version                 = version,
        num_files_per_segment   = num_files_per_segment,
        bands                   = bands,
        band_dets               = band_dets
        )

    # Write hdf files.
    comm_todwriter = CommanderHDFWriter(iras_comm_data_adapter)
    comm_todwriter.write_hdf_files(
        hdf_output_dir  = outpath, 
        overwrite       = overwrite, 
        num_processes   = num_processes,
        bands           = bands,
        verbosity       = 2
        )
    t1 = time()
    print(f"DONE! dur: {t1-t0:.3f} sec")


def main():
    global BANDS
    global BAND_DETS
    global NUM_PROCESSES
    global OVERWRITE
    global NSIDE
    """
    Set value of various parameters here 
    before calling the actual functions.
    This setup was adapted/written for simplified 
    testing/debugging purposes 
    """
    NUM_PROCESSES   = 256
    NSIDE           = 512
    OVERWRITE       = True
    # BANDS, BAND_DETS = set_bands_and_band_dets(testing=True, bands=["060"])
    # BANDS, BAND_DETS = set_bands_and_band_dets(testing=False, bands=["012"])
    # BANDS, BAND_DETS = set_bands_and_band_dets(testing=False, bands=["025"])
    # BANDS, BAND_DETS = set_bands_and_band_dets(testing=False, bands=["060"])
    # BANDS, BAND_DETS = set_bands_and_band_dets(testing=False, bands=["100"])
    # version = None # Integer. Which version of the hdf files we're writing. 
    # run_single_band_write(version=version)

    BANDS, BAND_DETS = set_bands_and_band_dets(bands=["025"], testing=True)
    version = 1
    outpath = f"{OUTPATH}/test_v{version}"

    run_multiple_band_write(version=version, outpath=outpath, num_processes=NUM_PROCESSES)
    # run_multiple_band_write(version=version, outpath=outpath, num_processes=None)



if __name__ == '__main__':
    main()
