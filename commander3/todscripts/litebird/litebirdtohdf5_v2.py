#---------------------------------------------
# Native Python libs 
#---------------------------------------------
import os 
import re
import sys 
import time
import pathlib 
import glob 
import multiprocessing as mp
#---------------------------------------------
# Third-party libraries
#---------------------------------------------
import h5py
import numpy as np 
import healpy as hp
import litebird_sim as lbs
import joblib
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.coordinates import BarycentricMeanEcliptic
# to see progress bar
from tqdm import tqdm
#
from collections import deque
#---------------------------------------------
# Local Development Libs
#---------------------------------------------
# TODO: In the future version you should use `cosmoglobe` and so eliminating
# the need to do this workaround in the first place 
# 
# Getting full path to Mathew's library as an object
commander_tools_path = pathlib.Path(__file__).absolute().parents[2].joinpath(
        'python','commander_tools').resolve()
# Appending the path to `PYTHONPATH`, so no need to 
# modify it externally (in your `.bashrc` etc.)
sys.path.append(str(commander_tools_path))
# Importing necessary modules from Mathew's library 
from tod_tools import commander_tod as comm_tod
#from testlib import LitebirdTodReader3
from tod_tools.litebird_imo import LitebirdImo
#from tod_tools.litebird import litebird
#from tod_tools.litebird_native_tod_reader import LitebirdTodReader#, LitebirdTodReader2
import litebird_sim as lbs

"""
This script converts LiteBIRD sims produced by Giuseppe at el into Commander3
suitable format 
"""


def get_data(nside, k, dfile, scale_loss):

    #print(f"Test statement {k}")
    with h5py.File(dfile, 'r') as readin_file:
        # 
        psi       = np.array(readin_file.get("psi"))
        # working with tods 
        tod_cmb   = np.array(readin_file.get("tod_cmb"))
        tod_dip   = np.array(readin_file.get("tod_dip")) 
        tod_fg    = np.array(readin_file.get("tod_fg")) 
        tod_wn    = np.array(readin_file.get("tod_wn")) 
        
        #Scale noise to number of detectors
        nr_det_out= 4 #Nr of detectors in output
        nr_det    = tod_cmb.shape[0]

        # Scale white noise level to reduced nr of detectors
        scale_nrdet = np.sqrt(nr_det_out/nr_det)
        # Scale whote noise level due to data loss (not done in LB TOD sims)
        tod_wn    = tod_wn * scale_nrdet * scale_loss
        
        #tod_wn_1f_30mHz = np.array(readin_file.get("tod_wn_1f_30mHz"))
        tod_cadd  = tod_cmb + tod_dip + tod_fg + tod_wn
        # TODO: Add calculation of pixels <= pointings
        pointings = np.array(readin_file.get("pointings"))
        theta     = pointings[:,:,0]
        phi       = pointings[:,:,1]
        pixels    = hp.ang2pix(nside, theta, phi)


    # Keeping track of current index
    #return {k: tod_cadd}
    return (tod_cadd, pixels, psi)

def append_remnant(arr, remnant):
    """
    Method appends the remnant array to the list from the left 
    """

    arr = deque(arr)
    arr.appendleft(remnant)
    arr = list(arr)
    
    return arr


def main():
    """
    Main method of the script 
    """
    # ----------------------------------
    # User-Specified Parameters 
    # ----------------------------------

    nprocs = 48 #64 #joblib.cpu_count() - 28 #16
    version = np.string_('0.0.1')
    # The size of one scan in a new file (to be created/output) 
    scan_size = 2**16 #~1 hr
    # The amount of scans to have in a given (to be created/output) file
    scan_num  = 26 
    # Number of psi bins used in Huffman compression. Form 2^n
    npsi = 1 # BeyondPlanck: 4096
    # Start time for observation
    start_time = '2030-04-01T00:00:00'

    # ----------------------------------
    # Retrieving Instrument Data
    # ----------------------------------

    # To access the database etc.
    freqs = ['L1-060']

    imo_db_interface = lbs.Imo()
    imo_version = 'v1.3'
    instrument = 'LFT'
    imo_db_datapath = f"/releases/{imo_version}/satellite"
    # Getting Data From LB Database
    channel_info = imo_db_interface.query(
            f"{imo_db_datapath}/{instrument}/{freqs[0]}/channel_info"
            )
    metadata = channel_info.metadata
    #print(metadata)
    instrument_info = imo_db_interface.query(
        f"/releases/{imo_version}/satellite/LFT/instrument_info"
        )
    #print(instrument_info.metadata)
    # Number of Detectos for a given Channel
    ndets      = metadata["number_of_detectors"]
    # Detector Lables for a given Instrument 
    det_labels = metadata["detector_names"]
    # Sampling rate in Hz
    fsamp      = metadata["sampling_rate_hz"]
    # Knee frequency in MHz 
    #fknee      = metadata["fknee_mhz"]
    # Alpha 
    #alpha      = metadata["alpha"]
    # FWHM in arcmin 
    fwhm       = metadata["fwhm_arcmin"]

    nside = 512 
    
    # Get parameters for data loss (not taken into account in LB TOD sims)
    scan_params = imo_db_interface.query(
                "/releases/v1.3/satellite/scanning_parameters"
                )
    metadata = scan_params.metadata
    #print("scanning parameters")
    #print(scan_params.metadata)
    margin    = metadata["margin"]
    det_yield = metadata["detector_yield"]     # detector yield
    cos_ray   = metadata["cosmic_ray_loss"]    # Cosmic ray loss
    duty_cycl = metadata["observation_duty_cycle"]    # Observation duty cycle
    # Scale white noise level due to data loss (not done in LB TOD sims)
    scale_loss = 1/np.sqrt(det_yield*margin*cos_ray*duty_cycl)

    # ----------------------------------
    # Retrieving Simulations' Data
    # ----------------------------------

    lbdata_dir = pathlib.Path(
            "/mn/stornext/d22/cmbco/litebird/e2e_ns512/sim0000/detectors_LFT_L1-060_T+B/tods"
            )
    output_dir = pathlib.Path(
            "/mn/stornext/u3/ragnaur/data/tut/Commander3_LB_TOD/TODS"
            )
    if not pathlib.Path.is_dir(output_dir):
        pathlib.Path.mkdir(output_dir)

    # Getting file names inside specified directory and removing the path component
    lbdata_files = sorted(lbdata_dir.rglob('*_rank????.hdf5')) 
    #lbdata_files = [data_file.name for data_file in lbdata_files]

    # Getting all the data sizes so later will split everything into chunks
    with h5py.File(lbdata_files[0], 'r') as readin_file:
        tod_cmb  = np.array(readin_file.get("tod_cmb"))
        # The remaining TOD chunk to stitch to a new chunk  
        remnant_tod = np.zeros_like(tod_cmb)
        remnant_pix = np.zeros_like(tod_cmb)
        remnant_psi = np.zeros_like(tod_cmb)

    # Getting the length of the single tod to calculate minimnum number of
    # files to open 
    todlen = len(tod_cmb[0])

    nfiles = len(lbdata_files)
    
    # Freeing up the memory by deleting unnecessary object(s)
    del tod_cmb#, tod_pix, tod_psi

    # The amount of time (in sec) the given scan will have 
    # (this value is unused and it just for us to see)
    scan_time = scan_size / fsamp #19.0

    # Number of scans of length scan_size (out of the entire length todlen) to
    # obtain from a given simulation file 
    N_scans_from_file = todlen / scan_size 
    # If this number is less then scan_num, then we open additional file
    # and to do that we need to know a number of files to open 
    #if N_scans_from_file < scan_num:
    # Minimum number of files to open to get a single output file (with a
    # specified number of scans and length of a single scan). Getting int
    # number and adding one to ensure that we do not go less then 1 file. 
    # [**NOTE**]: Each processor will get this value to work with in parallel. 
    nfiles_to_open = scan_num // N_scans_from_file + 1
    #elif N_scans_from_file >= scan_num:
    #    nfiles_to_open = 1
        # 
    print("#------------------")
    print(f"Total number of simulated files to work with: {nfiles}")
    print(f"The length of a single TOD: {todlen}")
    print(f"The scan size in a new file: {scan_size} ({scan_time:.2f} s or {scan_time/3600:.2f} hrs)")
    print(f"Total Number of scans in a new file: {scan_num}")
    print(f"Number of scans to obtain from simulated input file: {N_scans_from_file:.2f}")
    print(f"The number ratio of new scan to old scan: {scan_num / N_scans_from_file:.2f}")
    print(f"Number of files to open (per CPU process): {nfiles_to_open:.2f}")
    print("#------------------")


    # Ensuring we get to open integer number of files
    nfiles_to_open = int(nfiles_to_open)
    # Splitting it the workload (number of files) into equal batches.
    # So each processor of nprocs will work with nfiles_to_open in a 
    # given batch. The more processors involved, the bigger one batch 
    # and the smaller the end loop (the faster the calculations).
    batches = [ int(nfiles_to_open * nprocs)
            for i in range( nfiles // int(nfiles_to_open * nprocs) ) ]

    # Appending the remnant batch to the list of batches
    # [**NOTE**]: To avoid situation when the split was done equally and so wew
    # have 0 in the end, we do the if statement, i.e.
    # Batches are: [144, 144, 144,144, 144, 144, 144, 144, 144, 0]
    if nfiles > np.sum(batches):
        batches.append(nfiles - np.sum(batches))

    file_ranges = [0]
    for i in range(len(batches)):
        if i == len(batches)-1:
            file_ranges.append(nfiles)
        else:
            file_ranges.append((i+1)*batches[i])

    print(f"Batches are: {batches}")
    print(f"File ranges to process: {file_ranges}")
    print("#------------------")

    # ----------------------------------
    # Starting Calculations/File Creations
    # ----------------------------------

    manager = mp.Manager()
    dicts = {freqs[0]:manager.dict()}#, 44:manager.dict(), 70:manager.dict()}
    ctod = comm_tod.commander_tod(output_dir, 'LB', version, dicts=dicts, overwrite=True)

    # The Operational Day in Full Analogy with Planck
    #od = 1 
    ods = [0]
    # To trace the global scan in all files, i.e. each scan in each file will
    # have unique identifier
    #global_scan_id = 0 

    #remnant_scans = 0
    # Huffmann compression
    huffman = ['huffman', {'dictNum':1}]

    print(f"The remnant tods shape is: {remnant_tod.shape} and the values are:\n{remnant_tod}")

    for i in range(1, len(file_ranges)):
    #for i in tqdm(range(1, 2)):#len(batches)):
        print(f"Working with i = {i}:  {file_ranges[i-1]} -- {file_ranges[i]}")

        # This will pass e.g. 192 files (1 batch) if scan_size = 2**16,
        # scan_num = 20, nprocs = 64. In such a way each core will get 3 files
        # to work with. The resulting array will be of size (192, 48, 462334) 
        superTOD = joblib.Parallel(n_jobs=nprocs, backend="multiprocessing", verbose=2)(joblib.delayed(get_data)
                (nside, k, dfile, scale_loss) 
                for k, dfile in enumerate(lbdata_files[file_ranges[i-1]:file_ranges[i]])) 
                #(nside, k, dfile) for k, dfile in enumerate(lbdata_files))

        superTOD = list(map(list, zip(*superTOD)))
        superPsi = superTOD[2]
        superPix = superTOD[1]
        superTOD = superTOD[0]

        # Adding the remnant array to the super array from the left
        # (supposed to be very fast)
        if i != 1:         
            superTOD = append_remnant(superTOD, remnant_tod)
            superPix = append_remnant(superPix, remnant_pix)
            superPsi = append_remnant(superPsi, remnant_psi)
            # Freeing memory up 
            del remnant_tod, remnant_psi, remnant_pix

        ## Stitching subarrays together
        superTOD = np.concatenate(superTOD, axis=1)
        superPix = np.concatenate(superPix, axis=1)
        superPsi = np.concatenate(superPsi, axis=1)
        # Number of scans from the combined TOD
        nscansTOD = len(superTOD[0]) // scan_size + 1 
        # Number of files (in current iteration, i) to open and write data into  
        ods_shift = nscansTOD // scan_num
        ods.append(ods[i-1] + ods_shift) 
        print(ods)
        
        remnant_tod = []
        remnant_pix = []
        remnant_psi = []

        #print(len(superTOD[0]))
        print(f"For i = {i}: scansTOD = {nscansTOD} => new files = {ods_shift}")
        # Method creates the ods number of files 
        results = joblib.Parallel(n_jobs=nprocs, backend="multiprocessing", verbose=2)(joblib.delayed(make_ods)
                (ctod, imo_db_interface, imo_db_datapath, instrument, freqs[0], 
                    nside, fsamp, ndets, det_labels, scan_size, scan_num, ods[i-1], od, 
                    superTOD, superPix, superPsi, huffman, start_time, npsi) 
                for od in range(ods_shift)) 
                #for od in range(ods[i-1], ods[i], 1)) 

        # Getting the remainder to stitch to the left 
        # 33 * 20  
        remnant_scan_id = ods_shift * scan_num #+ scan_num
        #print(remnant_scan_id)
        #print(len(superTOD[0][remnant_scan_id*scan_size:]))

        for det_idx, det_label in enumerate(det_labels):
            remnant_tod.append(superTOD[det_idx][remnant_scan_id*scan_size:])
            remnant_pix.append(superPix[det_idx][remnant_scan_id*scan_size:])
            remnant_psi.append(superPsi[det_idx][remnant_scan_id*scan_size:])

        # Freeing memory up
        del superPsi, superPix, superTOD
        
        remnant_tod = np.array(remnant_tod)
        remnant_pix = np.array(remnant_pix)
        remnant_psi = np.array(remnant_psi)



    # Writing into last file whatever was left from the main loop (given that
    # it is of a power of two) 
    nscansRem = len(remnant_tod[0]) // scan_size
    print(ods[-1])
    if nscansRem > 0:
        #nscansTOD = len(remnant_tod[0]) // scan_size + 1 
        remnant_scans = nscansTOD - ods_shift * scan_num # <= do I need this one?
        print(f"For i = {i}: scansRem = {nscansRem} => new files = {1}")
        # Here, instead of scan_num, we should use whatever amount of scans we can
        # include into a file
        results = make_ods(ctod, imo_db_interface, imo_db_datapath, instrument, 
                freqs[0], nside, fsamp, ndets, det_labels, scan_size, scan_num, #nscansRem+scan_num, 
                           ods[-1], 0, remnant_tod, remnant_pix, remnant_psi, huffman, start_time, npsi)
    ctod.make_filelists()


def make_ods(ctod, imo_db_interface, imo_db_datapath, instrument, freq, nside, fsamp, 
        ndets, det_labels, scan_size, scan_num, ods, od, superTOD, superPix, superPsi, 
             huffman, start_time, npsi): 

    # Initialising new file 
    ctod.init_file(freq, ods + od + 1, mode='w')

    for local_scan_id in range(scan_num): 
        
        scan_id = od * scan_num + local_scan_id
        #print(f"od = {od}: {global_scan_id}:{global_scan_id + 1}")
        global_scan_id = (ods + od) * scan_num + local_scan_id 
        
        #Add fields common for each scan
        #Size of scan
        ctod.add_field(f"{global_scan_id}".zfill(6) + "/common/ntod", scan_size)
        
        #Time is given by file number being read in; one file is one day
        #The time at the start of this chunk. Space is given for 3 different units if desired.
        time_now = global_scan_id * scan_size/86400/fsamp + Time(start_time, format="isot")
        #print("time: ", time_now, "scan ID: ", global_scan_id, "mjd: ", time_now.mjd)
        ctod.add_field(f"{global_scan_id}".zfill(6) + "/common/time", [time_now.mjd,0,0])
        
        #Satelite position. TODO: Needs to be calculated
        orbit = lbs.SpacecraftOrbit(time_now)
        #spacecraft position and velocity
        pos_vel = lbs.spacecraft_pos_and_vel(orbit, start_time=time_now,
             time_span_s=scan_size/fsamp, delta_time_s=86400.0)
        #print(pos_vel)
        #print(pos_vel.positions_km,pos_vel.velocities_km_s)
        #print(pos_vel.positions_km[0])
        #print(SkyCoord(x=pos_vel.positions_km[0,0], y=pos_vel.positions_km[0,1], z=pos_vel.positions_km[0,2],
        #               frame=BarycentricMeanEcliptic().name, representation_type="cartesian")
        #      .transform_to("galactic").data.to_cartesian().get_xyz().value)
        ctod.add_field(f"{global_scan_id}".zfill(6) + "/common/satpos", [0,0,0])
        
        #vsun: Satellite velocity: The x,y,z velocity of the satellite relative to the sun. 
        #TODO: needs to be calculated 
        #print(pos_vel.velocities_km_s)
        #print(SkyCoord(x=pos_vel.positions_km[1,0], y=pos_vel.positions_km[1,1], z=pos_vel.positions_km[1,2],
        #               v_x=pos_vel.velocities_km_s[0,0], v_y=pos_vel.velocities_km_s[0,1], 
        #               v_z=pos_vel.velocities_km_s[0,2], frame=BarycentricMeanEcliptic().name, 
        #               differential_type='cartesian', representation_type="cartesian")
        #      .transform_to("galactic").data.to_cartesian().get_xyz().value)
        ctod.add_field(f"{global_scan_id}".zfill(6) + "/common/vsun", [0,0,0])

        for det_idx, det_label in enumerate(det_labels):
            # Ensuring the name gets unique indentifier (scan id)
            #prefix = f"{global_scan_id}".zfill(6) + "/" + det_label
            prefix = f"{global_scan_id}".zfill(6) + "/" + det_label

            arr = superTOD[det_idx][scan_id*scan_size:(scan_id+1)*scan_size]

            # TODO: add this logic since otherwise the latest scans will
            # otherwrite existing ones. 
            # If array is empty (mainly viable for the end scans) or its length is less than
            # acceptable chunk, we break the loop
            if len(arr) != scan_size or np.all(arr==0): # <= if array is empty or all are zeros (the scan is non-existent)
                break
            ctod.add_field(prefix + "/tod", arr)

            arr = superPix[det_idx][scan_id*scan_size:(scan_id+1)*scan_size]
            ctod.add_field(prefix + "/pix", arr, huffman)

            arr = superPsi[det_idx][scan_id*scan_size:(scan_id+1)*scan_size]
            ctod.add_field(prefix + "/psi", arr, huffman)
            # Getting Scalars, namely:
            # gain, sigma0, fknee, alpha
            detector_info = imo_db_interface.query(
                    f"{imo_db_datapath}/{instrument}/{freq}/{det_label}/detector_info"
                    )
            metadata = detector_info.metadata
            fknee    = metadata["fknee_mhz"]
            # Alpha 
            alpha    = metadata["alpha"]
            # Making up some values
            gain     = 1.0 
            sigma0   = 1.0
            scalars  = np.array([gain, sigma0, fknee, alpha])
            ctod.add_field(prefix + '/scalars', scalars)

        if len(arr) != scan_size or np.all(arr==0): # <= if array is empty or all are zeros (the scan is non-existent)
            break

        # Finilising Data for Each Scan
        ctod.finalize_chunk(f'{global_scan_id}'.zfill(6))

    # Things common for each scan in a given file
    prefix = 'common'
    #print(det_labels)
    # Detector Labels
    ctod.add_field(prefix + '/det',      np.string_(', '.join(det_labels)))
    #ctod.add_field(prefix + '/det',      np.string_(det_labels))
    #ctod.add_field(prefix + '/datatype', np.string_(datatype))
    ctod.add_field(prefix + '/npsi',     npsi)
    ctod.add_field(prefix + '/nside',    nside)
    # Polarization angles
    if instrument in ('LFT', 'HFT'):
        polang = [0 + np.pi/4 * int(it.split('_')[3][0] == 'U') + np.pi/2 *
                int(it.split('_')[3][1] == 'B') for it in det_labels]
    else:
        polang = [int(it.split('_')[3][:1])*np.pi/180 + np.pi/2 *
                int(it.split('_')[3][2] == 'B') for it in det_labels]
    #print(polang)
    ctod.add_field(prefix + '/polang',   polang)
    # Main Beam
    mbang  =  np.array([0.0] * ndets)
    ctod.add_field(prefix + '/mbang',    mbang)
    # Sampling Rate
    ctod.add_field(prefix + '/fsamp',    fsamp)

    # Closing the file
    ctod.finalize_file()

    return 0 
    



    """
    End data format for a single file

    scan 1/det 1
           ...
          /det N/pix 
                /psi 
                /flag 
                /outP
                /scalars 
                /tod 
          /common/fsamp 
                 /nside
                 /det 
                 /polang 
                 /mbang 
                 /npsi
    ...
    scan M/det 1 
           ...
           /det N 
    common/det
          /nside 
          /npsi 
          /fsamp
          /polang <= optional (?)
          /mbang  <= optional (?)

    where N = 48, and scan is determined by us 
    """


if __name__ == '__main__':
    start_time = time.time()
    print("Script has started!")
    main()
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Script run time: {total_time:.2f} s")




