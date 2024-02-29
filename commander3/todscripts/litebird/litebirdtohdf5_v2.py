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
from astropy.time import TimeDelta
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import BarycentricMeanEcliptic
import warnings 
from astropy.utils.exceptions import AstropyWarning
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
from tod_tools.lfi import lfi
#from testlib import LitebirdTodReader3
from tod_tools.litebird_imo import LitebirdImo
#from tod_tools.litebird import litebird
#from tod_tools.litebird_native_tod_reader import LitebirdTodReader#, LitebirdTodReader2
import litebird_sim as lbs

"""
This script converts LiteBIRD sims produced by Giuseppe at el into Commander3
suitable format 
Need to start server for litebird instrument and activate litebird IMO database to run the code

The script needs user specification of 
- LiteBIRD TOD simulation files directory
- List of detectors to save the data from from each frequency band
- location to output new Commander3 formated simulations

"""

def get_detectors(det_file, detector_list):
    # Takes inn a list of all detectors in the LB TOD simulations and
    # a file where list of detectors that will be used
    # Returns array of detector number (detectors) and new list of detector labels
    detectors = []
    det_labels = []
    file = open(det_file, 'r')
    while True:
        line = file.readline()
        if not line:
            break
        det_labels.append(str.strip(line))
        i = 0
        for l in detector_list: 
            if str.strip(line) == str.strip(l):
                detectors.append(i)
            i += 1
    file.close()
    print(det_labels, detectors)
    return (detectors, det_labels)


def get_data(nside, k, dfile, scale_loss, dets_nr, det_labels):
    """
    In: 
    nside      : nside of map
    k          : ? Not used
    dfile      : date file containing TODS to be read in
                 In baryonic mean ecliptic coordinates
    scale_loss : Data loss due to cosmic ray etc... 
                 Not taken into account by LB sim team
    dets_nr    : list of detector numbers of which to keep data from
                 Used to reduce the data volume to fewer detectors

    Out:
    tod_cadd   : Coadded TOD containing time stream for for instance 
                 cmb, white noise and foregrounds
                 In galactic coordinates
    pixels     : List of pixels corresponding to each tod sample in tod_coadd 
    psi_gal    : Polarization angle (theta, phi) for each sample in tod_coadd
    """
    #Reading from data file:
    dt = np.dtype('f8')
    with h5py.File(dfile, 'r') as readin_file:
        #Polarization angle in ecliptic coordinates
        psi_all   = np.mod(np.array(readin_file.get("psi")), 2 * np.pi) # Project onto the circle
        #Pointings in ecliptic
        pointings_ecl_all = np.array(readin_file.get("pointings"))
        #TODS:
        tod_cmb   = np.array(readin_file.get("tod_cmb"), dtype=dt)
        tod_dip   = np.array(readin_file.get("tod_dip"), dtype=dt)
        tod_fg    = np.array(readin_file.get("tod_fg"), dtype=dt)
        #tod_wn    = np.array(readin_file.get("tod_wn"), dtype=dt)
        tod_wn_1f_30mHz = np.array(readin_file.get("tod_wn_1f_30mHz"), dtype=dt)
 
    # Number of samples in one TOD
    tod_len      = psi_all.shape[1]
    #Scale noise to number of detectors
    nr_det_out   = len(dets_nr) #Nr of detectors in output
    nr_det       = tod_cmb.shape[0]
    # Scale white noise level to reduced nr of detectors
    scale_nrdet  = np.sqrt(nr_det_out/nr_det)
    # Scale white noise level due to data loss (not done in LB TOD sims)
    #tod_wn       = tod_wn * scale_nrdet * scale_loss
    # Scalw wn and 1/f noise
    tod_wn_1f_30mHz = tod_wn_1f_30mHz * scale_nrdet * scale_loss

    # Add TODS
    #tod_cadd_all = tod_cmb + tod_wn + tod_fg #+ tod_dip
    tod_cadd_all = tod_cmb + tod_wn_1f_30mHz + tod_fg + tod_dip
    
    #sigma_0 = np.diff(tod_wn).std() / 2**0.5  # Using Eqn 20 of BP06
        
    #Reduce number of detectors: 
    psi_ecl      = np.zeros(shape=(nr_det_out,tod_len))
    tod_cadd     = np.zeros(shape=(nr_det_out,tod_len))
    pointings_ecl= np.zeros(shape=(nr_det_out,tod_len,2))
    # Keep only detectors from input list
    for d in range(0,nr_det_out):
        psi_ecl[d,:]       = psi_all[int(dets_nr[d]), :]
        tod_cadd[d,:]      = tod_cadd_all[int(dets_nr[d]), :]
        pointings_ecl[d,:,:] = pointings_ecl_all[int(dets_nr[d]),:,:]
    #change from barycentricmeanecliptic --> galactic
    psi_gal   = np.zeros_like(psi_ecl)
    pointings_gal = np.zeros_like(pointings_ecl)
    for d in range(0,nr_det_out):
        pointings_gal[d,:,:], psi_gal[d,:] = lbs.coordinates.rotate_coordinates_e2g(pointings_ecl[d,:,:], pol_angle_ecl=psi_ecl[d,:])
        #if det_labels[d].split('_')[3][-1] == 'B':
        #psi_gal =  2*np.pi - psi_gal
    theta     = pointings_gal[:,:,0]
    phi       = pointings_gal[:,:,1]
    pixels    = hp.ang2pix(nside, theta, phi)#, lonlat=True)
    return (tod_cadd, pixels, psi_gal)


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

    nprocs = 64 #joblib.cpu_count() - 28 #16
    version = np.string_('0.0.1')
    # The size of one scan in a new file (to be created/output) 
    scan_size = 2**16 #~1 hr
    # The amount of scans to have in a given (to be created/output) file
    scan_num  = 25 # ~1 day
    # Number of psi bins used in Huffman compression. Form 2^n
    npsi = 4096 # BeyondPlanck: 4096
    # Start time for observation
    start_time = '2030-04-01T00:00:00'
#    nside = 512 # Defined further down - dependent on FWHM
    output_dir = pathlib.Path(
#         "/mn/stornext/u3/ragnaur/data/tut/Commander3_LB_TOD/TODS/"
#         "/mn/stornext/u3/eirikgje/data/litebird_tods/"
         "/mn/stornext/d16/cmbco/litebird/TOD_analysis/TODS_eirik_cmb_fg_wn_ncorr30_dipol_v3/")
         
    lbdata_dir = pathlib.Path(
         "/mn/stornext/d22/cmbco/litebird/e2e_ns512/sim0000/"
         )
    det_dir   = pathlib.Path(
        "/mn/stornext/u3/ragnaur/data/tut/Commander3_LB_TOD/data_LB/"
    )

    # ----------------------------------
    # Initializing common to all files
    # ----------------------------------
    if not pathlib.Path.is_dir(output_dir):
        pathlib.Path.mkdir(output_dir)

    # Huffman compression
    compArray = [lfi.psiDigitize, lfi.huffman]


    # ----------------------------------
    # Retrieving Instrument Data
    # ----------------------------------    
    imo_db_interface = lbs.Imo(flatfile_location='/mn/stornext/u3/eirikgje/src/litebird_imo/IMO/')
    imo_version = 'v1.3'

    # Get parameters for data loss (not taken into account in LB TOD sims)
    scan_params = imo_db_interface.query(
        f"/releases/{imo_version}/satellite/scanning_parameters"
        )
    metadata = scan_params.metadata
    margin    = metadata["margin"]
    det_yield = metadata["detector_yield"]     # detector yield
    cos_ray   = metadata["cosmic_ray_loss"]    # Cosmic ray loss
    duty_cycl = metadata["observation_duty_cycle"]    # Observation duty cycle
    # Scale white noise level due to data loss (not done in LB TOD sims)
    scale_loss = 1/np.sqrt(det_yield*margin*cos_ray*duty_cycl)

    instrument = ['LFT', 'MFT', 'HFT']
#    instrument = ['MFT', 'HFT']
    imo_db_datapath = f"/releases/{imo_version}/satellite"
    for inst in instrument:
        print("Working with detectors on ", inst)

        # Getting Data From LB Database
        instrument_info = imo_db_interface.query(
            f"{imo_db_datapath}/{inst}/instrument_info"
            )
        freqs = instrument_info.metadata['channel_names']
        # or specify detectors manually:
        #freqs = ['M2-119', 'M1-140', 'M2-166', 'M1-195']
        
        for freq in freqs:
            print("Working with frequency band ", freq)
            channel_info = imo_db_interface.query(
                f"{imo_db_datapath}/{inst}/{freq}/channel_info"
                )
            metadata = channel_info.metadata
            # Number of Detectos for a given Channel
            #ndets      = metadata["number_of_detectors"]
            # Detector Lables for a given Instrument 
            #det_labels = metadata["detector_names"]
            #Read detector labels from file from e2e sim team
            det_labels = []
            all_dets   = "detectors_" + inst + "_" + freq + "_T+B.txt"
            file = open(det_dir / all_dets, 'r')
            while True:
                line = file.readline()
                if not line:
                    break
                det_labels.append(str.strip(line))
            file.close()
            #Get detector numbers for detectors used in current setup
            det_name   = "detectors_" + inst + "_" + freq + "_T+B_4.txt"
            det_file   = det_dir / det_name
            det_nr_used, det_labels= get_detectors(det_file, det_labels)
            dets_nr    = len(det_nr_used)
            # Sampling rate in Hz
            fsamp      = metadata["sampling_rate_hz"]
            # FWHM in arcmin 
            fwhm       = metadata["fwhm_arcmin"]
            if fwhm > 30:
                nside = 512
            else:
                nside = 1024
            
            # ----------------------------------
            # Retrieving Simulations' Data
            # ----------------------------------
            folder = "detectors_" + inst + "_" + freq + "_T+B/tods"
            lbfreq_dir = lbdata_dir / folder
            output_freqname = 'LB_' + '_'.join(freq.split('-')[::-1])
            det_dir_out = output_dir / output_freqname
            
            if not pathlib.Path.is_dir(det_dir_out):
                pathlib.Path.mkdir(det_dir_out)

            # Getting file names inside specified directory and removing the path component
            lbdata_files = sorted(lbfreq_dir.rglob('*_rank????.hdf5')) 
            
            # Getting all the data sizes so later will split everything into chunks
            with h5py.File(lbdata_files[0], 'r') as readin_file:
                tod_wn  = np.array(readin_file.get("tod_wn"))
            
            # The remaining TOD chunk to stitch to a new chunk  
            remnant_tod = np.zeros(shape=(dets_nr,tod_wn.shape[1]))
            remnant_pix = np.zeros(shape=(dets_nr,tod_wn.shape[1]))
            remnant_psi = np.zeros(shape=(dets_nr,tod_wn.shape[1]))

            # Getting the length of the single tod to calculate minimnum number of files
            todlen = len(tod_wn[0])
            nfiles = len(lbdata_files)

            # Calculating sigma0
            scale_nrdet  = np.sqrt(len(tod_wn[:,0])/dets_nr)
            print("Total det nr:", len(tod_wn[:,0]), "det nr used", dets_nr)
            # Scale white noise level due to data loss (not done in LB TOD sims)
            tod_wn       = tod_wn * scale_nrdet * scale_loss
            sigma0       = np.diff(tod_wn).std() / 2**0.5  # Using Eqn 20 of BP06
    
            # Freeing up the memory by deleting unnecessary object(s)
            del tod_wn

            # The amount of time (in sec) the given scan will have 
            # (this value is unused and it just for us to see)
            scan_time = scan_size / fsamp #19.0

            # Number of scans of length scan_size (out of the entire length todlen) to
            # obtain from a given simulation file 
            N_scans_from_file = todlen / scan_size 
            # Minimum number of files to open to get a single output file (with a
            # specified number of scans and length of a single scan). Getting int
            # number and adding one to ensure that we do not go less then 1 file. 
            # [**NOTE**]: Each processor will get this value to work with in parallel. 
            nfiles_to_open = scan_num // N_scans_from_file + 1

            # Ensuring we get to open integer number of files
            nfiles_to_open = int(nfiles_to_open)
            # Splitting it the workload (number of files) into equal batches.
            # So each processor of nprocs will work with nfiles_to_open in a 
            # given batch. The more processors involved, the bigger one batch 
            # and the smaller the end loop (the faster the calculations).
            batches = [ int(nfiles_to_open * nprocs)
                        for i in range( nfiles // int(nfiles_to_open * nprocs) ) ]

            # Appending the remnant batch to the list of batches
            # [**NOTE**]: To avoid situation when the split was done equally and so we
            # have 0 in the end, we do the if statement, i.e.
            if nfiles > np.sum(batches):
                batches.append(nfiles - np.sum(batches))

            file_ranges = [0]
            for i in range(len(batches)):
                if i == len(batches)-1:
                    file_ranges.append(nfiles)
                else:
                    file_ranges.append((i+1)*batches[i])
            print("#------------------")
            print(f"Total number of simulated files to work with: {nfiles}")
            print(f"The length of a single TOD: {todlen}")
            print(f"The scan size in a new file: {scan_size} ({scan_time:.2f} s or {scan_time/3600:.2f} hrs)")
            print(f"Total Number of scans in a new file: {scan_num}")
            print(f"Number of scans to obtain from simulated input file: {N_scans_from_file:.2f}")
            print(f"The number ratio of new scan to old scan: {scan_num / N_scans_from_file:.2f}")
            print(f"Number of files to open (per CPU process): {nfiles_to_open:.2f}")
            print(f"Number of detectors to save: {dets_nr}")
            print(f"Batches are: {batches}")
            print(f"File ranges to process: {file_ranges}")
            print("#------------------")

            # ----------------------------------
            # Starting Calculations/File Creations
            # ----------------------------------

            manager = mp.Manager()
#            dicts = {freq:manager.dict()}
            dicts = {output_freqname:manager.dict()}
#            ctod = comm_tod.commander_tod(det_dir_out, 'LB', version, dicts=dicts, overwrite=True)
            ctod = comm_tod.commander_tod(det_dir_out, None, version, dicts=dicts, overwrite=True)

            # The Operational Day in Full Analogy with Planck
            ods = [0]
            # To trace the global scan in all files, i.e. each scan in each file will
            # have unique identifier
            
            print(f"The remnant tods shape is: {remnant_tod.shape} and the values are:\n{remnant_tod}")
            for i in range(1, len(file_ranges)):
                print(f"Working with i = {i}:  {file_ranges[i-1]} -- {file_ranges[i]}")

                # This will pass e.g. 192 files (1 batch) if scan_size = 2**16,
                # scan_num = 20, nprocs = 64. In such a way each core will get 3 files
                # to work with. The resulting array will be of size (192, 48, 462334) 
                superTOD = joblib.Parallel(n_jobs=nprocs, backend="multiprocessing", verbose=2)(joblib.delayed(get_data)
                                                                                                (nside, k, dfile, scale_loss, det_nr_used, det_labels) 
                                for k, dfile in enumerate(lbdata_files[file_ranges[i-1]:file_ranges[i]])) 

                superTOD = list(map(list, zip(*superTOD)))
                superPsi = superTOD[2]
                superPix = superTOD[1]
                superTOD = superTOD[0]

                # Adding the remnant array to the super array from the left
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
                print(f"For i = {i}: scansTOD = {nscansTOD} => new files = {ods_shift}")
                # Method creates the ods number of files 
                results = joblib.Parallel(n_jobs=nprocs, backend="multiprocessing", verbose=2)(joblib.delayed(make_ods)
                                 (ctod, imo_db_interface, imo_db_datapath, inst, freq, 
                                 nside, fsamp, det_labels, scan_size, scan_num, ods[i-1], od, 
                                 superTOD, superPix, superPsi, compArray, start_time, npsi, sigma0) 
                                 for od in range(ods_shift)) 
                                
                # Getting the remainder to stitch to the left 
                remnant_scan_id = ods_shift * scan_num 
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
                remnant_scans = nscansTOD - ods_shift * scan_num # <= do I need this one?
                print(f"For i = {i}: scansRem = {nscansRem} => new files = {1}")
                # Here, instead of scan_num, we should use whatever amount of scans we can
                # include into a file
                results = make_ods(ctod, imo_db_interface, imo_db_datapath, inst, 
                                   freq, nside, fsamp, det_labels, scan_size, scan_num,
                                   ods[-1], 0, remnant_tod, remnant_pix, remnant_psi, compArray, 
                                   start_time, npsi, sigma0)
            ctod.make_filelists()
        del remnant_tod, remnant_pix, remnant_psi


def make_ods(ctod, imo_db_interface, imo_db_datapath, instrument, freq, nside, fsamp, 
        det_labels, scan_size, scan_num, ods, od, superTOD, superPix, superPsi, 
             compArray, start_time, npsi, sigma0): 
    """
    In: 
    ctod
    imo_db_interface
    imo_db_datapath
    instrument
    freq
    nside
    fsamp
    det_labels
    scan_size
    scan_num
    ods
    od
    superTOD
    superPix
    superPsi
    compArray
    start_time
    npsi
    sigma0
    
    Out: 
    """
    
    output_freqname = 'LB_' + '_'.join(freq.split('-')[::-1])
    # Initialising new file 
    ctod.init_file(output_freqname, ods + od + 1, mode='w')
    ndets = len(det_labels)
    
    for local_scan_id in range(scan_num): 
        scan_id = od * scan_num + local_scan_id
        #print(f"od = {od}: {global_scan_id}:{global_scan_id + 1}")
        global_scan_id = (ods + od) * scan_num + local_scan_id +1
        #print("global_scan_id", global_scan_id)
        #Add fields common for each scan
        ctod.add_field(f"{global_scan_id}".zfill(6) + "/common/ntod", scan_size)
        
        #Time is given by file number being read in; one file is one day
        #The time at the start of this chunk. Space is given for 3 different units if desired.
        #with warnings.catch_warnings(): #Ignore warnings of dubious year
        #warnings.simplefilter('ignore', AstropyWarning)
        warnings.filterwarnings('ignore', category=UserWarning, append=True)
        time_passed = TimeDelta((global_scan_id-1) * scan_size/fsamp, format='sec')
        time_now = time_passed + Time(start_time, format="isot")
        #print("time: ", time_now, "scan ID: ", global_scan_id, "mjd: ", time_now.mjd)
        ctod.add_field(f"{global_scan_id}".zfill(6) + "/common/time", [time_now.mjd,0,0])
        
        #Satelite position. TODO: Needs to be calculated
        orbit = lbs.SpacecraftOrbit(time_now)
        #spacecraft position and velocity
        pos_vel = lbs.spacecraft_pos_and_vel(orbit, start_time=time_now,
             time_span_s=scan_size/fsamp, delta_time_s=86400.0)

        # linear velocity of the spacecraft in the Barycentric Ecliptic reference frame (in km/s)
        vel_ecl = pos_vel.velocities_km_s 
        # position of the spacecraft in the Barycentric Ecliptic reference frame (in kilometers)
        pos_ecl = pos_vel.positions_km 
        pos_vel_gal = (SkyCoord(x=pos_ecl[1,0]*u.km, y=pos_ecl[1,1]*u.km, z=pos_ecl[1,2]*u.km,
                       v_x=vel_ecl[0,0]*u.km/u.s, v_y=vel_ecl[0,1]*u.km/u.s, v_z=vel_ecl[0,2]*u.km/u.s, 
                        frame=BarycentricMeanEcliptic().name,
                       differential_type='cartesian', representation_type="cartesian")
                       .transform_to("galactic")) 
        pos_gal = pos_vel_gal.data.to_cartesian().get_xyz().value*1000
        vel_gal = [pos_vel_gal.velocity.d_x.value*1000, pos_vel_gal.velocity.d_y.value*1000, pos_vel_gal.velocity.d_z.value*1000]

        ctod.add_field(f"{global_scan_id}".zfill(6) + "/common/satpos", pos_gal)
        ctod.add_field(f"{global_scan_id}".zfill(6) + "/common/vsun", vel_gal)

        for det_idx, det_label in enumerate(det_labels):
            # Ensuring the name gets unique indentifier (scan id)
            #prefix = f"{global_scan_id}".zfill(6) + "/" + det_label
            prefix = f"{global_scan_id}".zfill(6) + "/" + det_label

            tod_arr   = superTOD[det_idx][scan_id*scan_size:(scan_id+1)*scan_size]

            # If array is empty (mainly viable for the end scans) or its length is less than
            # acceptable chunk, we break the loop
            if len(tod_arr) != scan_size or np.all(tod_arr==0): 
                break
            ctod.add_field(prefix + "/tod", tod_arr)

            pix_arr   = superPix[det_idx][scan_id*scan_size:(scan_id+1)*scan_size]
            ctod.add_field(prefix + "/pix", pix_arr, [lfi.huffman])

            compArray = [lfi.psiDigitize, lfi.huffman]
            psi_arr   = superPsi[det_idx][scan_id*scan_size:(scan_id+1)*scan_size]
            ctod.add_field(prefix + "/psi", psi_arr, compArray)
            
            flag_arr = np.zeros(shape=pix_arr.shape)
            ctod.add_field(prefix + '/flag', flag_arr, compArray) # [lfi.huffman])

            # Getting Scalars: gain, sigma0, fknee, alpha
            detector_info = imo_db_interface.query(
                    f"{imo_db_datapath}/{instrument}/{freq}/{det_label}/detector_info"
                    )
            metadata  = detector_info.metadata
            # fknee in Hz
            fknee     = metadata["fknee_mhz"] /1000 #from mHz -> Hz
            # Alpha 
            alpha     = metadata["alpha"]*-1
            # Ideal gain
            gain      = 1.0 
            scalars   = np.array([gain, sigma0, fknee, alpha])
            ctod.add_field(prefix + '/scalars', scalars)
            
            # average outer product - estimate of sattellite spinn axis
            # used for load balancing. Set to zero for now
            # outAng   = lfi.ring_outer_product(Theta, Phi)
            outAng    = np.array([0,0])
            ctod.add_field(prefix + '/outP', outAng)

        if len(tod_arr) != scan_size or np.all(tod_arr==0):
            break
        # Finilising Data for Each Scan
        ctod.finalize_chunk(f'{global_scan_id}'.zfill(6))

    # Things common for each scan in a given file
    prefix = 'common'
    #print(det_labels)
    # Detector Labels - det_labels can be given as .txt 
    ctod.add_field(prefix + '/det',      np.string_(', '.join(det_labels)))
    #ctod.add_field(prefix + '/det',      np.string_(det_labels))
    #ctod.add_field(prefix + '/datatype', np.string_(datatype))
    ctod.add_field(prefix + '/npsi',     npsi)
    ctod.add_field(prefix + '/nside',    nside)
    # Polarization angles
    if instrument in ('LFT'):
        polang = [0 for it in det_labels]
        """
        polang = [(0 + np.pi/4 * int(it.split('_')[3][0] == 'U') + np.pi/2 *
                   int(it.split('_')[5][0] == 'B')) 
                  * (1-0*int(it.split('_')[3][1] == 'B'))  for it in det_labels]
        # change 0 to 2 to change angle
        """
    elif instrument in ('HFT'):
        polang = [0 for it in det_labels]
        #polang = [0 + np.pi/4 * int(it.split('_')[3][0] == 'U') + np.pi/2 *
        #          int(it.split('_')[5][0] == 'B') for it in det_labels]
    elif instrument in ('MFT'):
        polang = [0 for it in det_labels]
        """
        polang = [int(it.split('_')[3][:2])*np.pi/180 + np.pi/2 *
                  int(it.split('_')[5][0] == 'B')
                  * (1-0*int(it.split('_')[3][2] == 'B')) for it in det_labels]
    else:
        print("Error: wrong instrument")
    """
    #print(polang)
    ctod.add_field(prefix + '/polang',   polang)
    # Main Beam
    mbang  =  np.array([0.0] * ndets)
    ctod.add_field(prefix + '/mbang',    mbang)
    # Sampling Rate
    ctod.add_field(prefix + '/fsamp',    fsamp)

    # Closing the file
    ctod.finalize_file()
    ctod.outFile.close()
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




