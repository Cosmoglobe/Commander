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
# to see progress bar
from tqdm import tqdm
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


def get_data(nside, k, dfile):

    print(f"Test statement {k}")
    with h5py.File(dfile, 'r') as readin_file:
        tod_cmb  = np.array(readin_file.get("tod_cmb"))
        tod_dip  = np.array(readin_file.get("tod_dip")) 
        tod_fg   = np.array(readin_file.get("tod_fg")) 
        tod_wn   = np.array(readin_file.get("tod_wn")) 
        tod_cadd = tod_cmb + tod_dip + tod_fg + tod_wn
        # TODO: Add calculation of pixels <= pointings
        pointings = np.array(readin_file.get("pointings"))
        theta     = pointings[:,:,0]
        phi       = pointings[:,:,1]
        pixels    = hp.ang2pix(nside, theta, phi)


    # Keeping track of current index
    #return {k: tod_cadd}
    return (tod_cadd, pixels)




def main():
    """
    Main method of the script 
    """
    nprocs = 64 #joblib.cpu_count() - 28 #16
    #lbdata_dir = pathlib.Path("/Users/maksymb/Desktop/litebird_db/data") 
    lbdata_dir = pathlib.Path(
            "/mn/stornext/d22/cmbco/litebird/e2e_ns512/sim0000/detectors_LFT_L1-060_T+B/tods"
            )
    #output_dir = pathlib.Path(__file__).parent.joinpath("test_litebird_sims")
    output_dir = pathlib.Path(
            "/mn/stornext/d5/data/maksymb/litebird"
            )
    if not pathlib.Path.is_dir(output_dir):
        pathlib.Path.mkdir(output_dir)

    # Getting file names inside specified directory and removing the path component
    lbdata_files = sorted(lbdata_dir.rglob('*_rank????.hdf5')) 
    #lbdata_files = [data_file.name for data_file in lbdata_files]

    # Getting all the data sizes so later will split everything into chunks
    with h5py.File(lbdata_files[0], 'r') as readin_file:
        tod_cmb  = np.array(readin_file.get("tod_cmb"))

    todlen = len(tod_cmb[0])

    nfiles = len(lbdata_files)

    # These parameters are specified by the user (!)
    # The size of one scan in a new file (to be created/output) 
    scan_size = 2**16
    # The amount of scans to have in a given (to be created/output) file
    scan_num  = 20 
    # The amount of time (in sec) the given scan will have 
    # (this value is unused and it just for us to see)
    scan_time = scan_size / 19.0

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
    print(f"The number ratio of new scan to old scan: {scan_num // N_scans_from_file:.2f}")
    print(f"Number of files to open (per CPU process): {nfiles_to_open:.2f}")
    print("#------------------")


    nfiles_to_open = int(nfiles_to_open)
    # Splitting it the workload (number of files) into equal batches
    batches = [ int(nfiles_to_open * nprocs)
            for i in range( nfiles // int(nfiles_to_open * nprocs) ) ]
    # Appending the remnant batch to the list of batches
    batches.append(nfiles - np.sum(batches))
    # 
    file_ranges = [0]
    for i in range(len(batches)):
        if i == len(batches)-1:
            file_ranges.append(nfiles)
        else:
            file_ranges.append((i+1)*batches[i])

    print(f"Batches are: {batches}")
    print(f"File ranges to process: {file_ranges}")
    print("#------------------")
    exit()

    #exit()
    # Adding unequal part
    #for i in range( nfiles % (nfiles_to_open * nprocs) ):
    #    workloads[i] += 1
    #print(workloads)
    #print(len(workloads))
    #print(np.sum(workloads))

    #print(88768128 / 192)
    #print(f"{lbdata_files}")

    # TODO: do the loop in batches (just split the total data files into
    # batches nfiles_tot/4) and do the processing in batches
    #remnants = np.zeros((48, something here))
    nside = 512
    #print(len(lbdata_files[:batches[0]]))
    for i in range(1, 2):#len(file_ranges)):
    #for i in tqdm(range(1, 2)):#len(batches)):
        print(f"Working with i = {i}:  {file_ranges[i-1]} -- {file_ranges[i]}")

        #with joblib.parallel_backend(backend="multiprocessing", n_jobs=nprocs):
        # This will pass e.g. 192 files (1 batch) if scan_size = 2**16,
        # scan_num = 20, nprocs = 64. In such a way each core will get 3 files
        # to work with. The resulting array will be of size (192, 48, 462334) 
        superTOD = joblib.Parallel(n_jobs=nprocs,
                verbose=2)(joblib.delayed(get_data)
                (nside, k, dfile) 
                for k, dfile in enumerate(lbdata_files[file_ranges[i-1]:file_ranges[i]])) 
                #(nside, k, dfile) for k, dfile in enumerate(lbdata_files))

        superTOD = list(map(list, zip(*superTOD)))
        superPix = superTOD[1]
        superTOD = superTOD[0]
        print(np.concatenate(superTOD, axis=1).shape)
        superTOD = np.concatenate(superTOD, axis=1)
        superPix = np.concatenate(superPix, axis=1)
        # TODO: Split these into chunks of equal length using the end number of scans
        print(f"TOD[0]:\n{superTOD[0]}")
        print(superTOD.shape)
        print(f"Pix[0]:\n{superPix[0]}")
        print(f"Pix[4]:\n{superPix[4]}")
        print(superPix.shape)

        # TODO: 
        # Save the superTOD chunk into a file 

        #make_od()


    exit()

    #toddir = "/Users/maksymb/Desktop/litebird_db/maksymb_lb/commander3/todscripts/litebird/test_litebird_sims"
    #todfile = "LB_L1-060_000001.h5"
    #outtod = h5py.File(f"{toddir}/{todfile}")
    nside = 512
    outmap = np.zeros((48, 12 * nside ** 2))
    nhits = np.zeros((48, 12 * nside ** 2))

    #scans = [scan for scan in list(outtod) if scan != 'common']

    #print(list(outtod[scans[0]]))

    #det_names = [det for det in list(outtod[scans[0]]) if det != 'common'] 

    #print(scans)
    #for scan in scans:
    #    print(scan)
    for i in range(48):
        #currscan_pixs = outtod[scan][detector]['pix']
        #print(currscan_pixs)
        #currscan_tod  = outtod[scan][detector]['tod']
        #print(currscan_tod)
        #nhits[i, currscan_pixs] += 1
        #outmap[i, currscan_pixs] += currscan_tod
        nhits[i, superPix[i]] += 1
        outmap[i, superPix[i]] += superTOD[i]
    print("Test message")
    outmap[np.where(nhits != 0)] /= nhits[np.where(nhits != 0)]
    outmap[np.where(nhits == 0)] = hp.UNSEEN

    for i in range(48):
        print(f"Working with sim_{i:03}.png")
        hp.mollview(outmap[i])
        plt.savefig(
                #f"/Users/maksymb/Desktop/litebird_db/maksymb_lb/commander3/todscripts/litebird/test_litebird_sims/sim_{i:03}.png", dpi=800
                f"/mn/stornext/u3/maksymb/commander/maksymb_lb/commander3/todscripts/litebird/test_litebird_sims/sim_{i:03}.png", dpi=800
                )
                #'/home/eirik/temp/sim_{i:03}.png', dpi=800)
        plt.clf()
    exit()




    #print(superTOD[3])
    #print(superTOD[5])
    #(joblib.delayed(concat_data)(superTOD[:][k]) for k in range(48))
    #print(np.shape(superTOD))
    #print(np.transpose(superTOD, axes=(48,len(superTOD[0][0])*10)).shape)
    # Convert list of dictionaries into a dictionary 
    # (in this way we ensure that each data array has its own index)
    #superTOD = {k: v for d in superTOD for k, v in d.items()}
    # Populating the new array with correct indeces (in order)
    #superTOD = [v for k, v in superTOD.items()]
    #for i in range(len(superTOD)):
    #    superTOD2 = np.concatenate()
    #superTOD = np.concatenate( ) for i in range(len(superTOD))
    #print(superTOD.shape)



    #print(f"Number of files {len(lbdata_files)}")
    ##tod_tot = np.array(h5py.File(lbdata_dir/lbdata_files[0],"r").get)
    for i, dfile in enumerate(lbdata_files):
        with h5py.File(lbdata_dir / dfile, 'r') as readin_file:
            tod_cmb  = np.array(readin_file.get("tod_cmb"))
            tod_dip  = np.array(readin_file.get("tod_dip")) 
            tod_fg   = np.array(readin_file.get("tod_fg")) 
            tod_wn   = np.array(readin_file.get("tod_wn")) 
            tod_cadd = tod_cmb + tod_dip + tod_fg + tod_wn

            pointings = np.array(readin_file.get("pointings"))
            theta     = pointings[:,:,0]
            phi       = pointings[:,:,1]
            pixels    = hp.ang2pix(nside, theta, phi)
            if i == 0:
                tod_tot = tod_cadd 
                pix_tot = pixels
            else: 
                tod_tot = np.concatenate((tod_tot, tod_cadd),axis=1)
                pix_tot = np.concatenate((pix_tot, pixels),axis=1)
            #print(tod_tot.shape)
            #print(tod_cadd.shape)
        #print(dfile)
    print("Test")
    print(tod_tot[0])
    print(tod_tot.shape)
    print(pix_tot[0])
    print(pix_tot[4])
    print(pix_tot.shape)
    #print(tod_tot[5])
    #arrs = np.array_split(tod_tot, 2**15)
    #for i in range(3):
    #    print(arrs.shape)


    exit()



    version = np.string_('0.0.1')
    # To access the database etc.
    freqs = ['L1-060']

    # 1296 files (rank: 0 - 1295), 1296/16 = 81
    # Test to see how to split into chunks
    #testdata1 = [[1,2,3], [4,5,6]]
    #testdata2 = [[7,8,9], [10,11,12]]
    #points1  = [["meat1", "meat2", "meat3"], ["juice1", "juice2", "juice3"]]
    # end result: det1: [1,2,3,7,8,9], det2: [4,5,6,10,11,12]

    # Initialising commander tod object
    manager = mp.Manager()
    dicts = {freqs[0]:manager.dict()}#, 44:manager.dict(), 70:manager.dict()}
    ctod = comm_tod.commander_tod(output_dir, 'LB', version, dicts=dicts, overwrite=True)

    # TODO: extend this to the entire LFT L1-060 files 
    # and add parallel implementation
    make_od(lbdata_dir, lbdata_files, ctod, version, freqs)
    # Making filelists.txt after working with all the data was finalized
    ctod.make_filelists()




def make_od(lbdata_dir, lbdata_files, ctod, version, freqs):#, 
        #level3_ces_nums, ctod, version, freqs, dicts, k):
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
    k=0
    readin_file = h5py.File(lbdata_dir / lbdata_files[k], 'r')
    print(f"Working with file: {lbdata_files[k]}")
    # Things to include per detector
    psi     = np.array(readin_file.get("psi"))
    tod_cmb = np.array(readin_file.get("tod_cmb"))
    tod_dip = np.array(readin_file.get("tod_dip")) 
    tod_fg  = np.array(readin_file.get("tod_fg")) 
    tod_wn  = np.array(readin_file.get("tod_wn")) 
    # Coadding tods from a single file
    tod_tot = tod_cmb + tod_dip + tod_fg + tod_wn
    # Pointings <= should be compressed
    nside = 512
    pointings = np.array(readin_file.get("pointings"))
    theta     = pointings[:,:,0]
    phi       = pointings[:,:,1]
    pixels    = hp.ang2pix(nside, theta, phi)
    print(pixels.shape)
    #---------------------------------------------
    # Writing data to a file
    #---------------------------------------------
    # Huffmann compression
    huffman = ['huffman', {'dictNum':1}]
    # Creating new file
    od   = 1
    scan = 1
    imo_db_interface = lbs.Imo()
    imo_version = 'v1.3'
    instrument = 'LFT'
    imo_db_datapath = f"/releases/{imo_version}/satellite"

    # Getting Data From LB Database
    scan_params = imo_db_interface.query(
            f"{imo_db_datapath}/scanning_parameters"
            )
    metadata = scan_params.metadata
    # {'spin_sun_angle_deg': 45.0, 'precession_period_min': 192.348,
    # 'spin_rate_rpm': 0.05, 'mission_duration_year': 3.0,
    # 'observation_duty_cycle': 0.85, 'cosmic_ray_loss': 0.95, 'margin': 0.95,
    # 'detector_yield': 0.8}
    #print("Test statement")
    #print(metadata["spin_sun_angle_deg"])
    #print(metadata)
    channel_info = imo_db_interface.query(
            f"{imo_db_datapath}/{instrument}/{freqs[0]}/channel_info"
            )
    metadata = channel_info.metadata
    #print(metadata)
    # Number of Detectos for a given Channel
    ndets      = metadata["number_of_detectors"]
    # Detector Lables for a given Instrument 
    det_labels = metadata["detector_names"]
    # Sampling rate in Hz
    fsamp      = metadata["sampling_rate_hz"]
    # Knee frequency in MHz 
    fknee      = metadata["fknee_mhz"]
    # Alpha 
    alpha      = metadata["alpha"]
    # FWHM in arcmin 
    fwhm       = metadata["fwhm_arcmin"]
    #print(metadata["detector_names"])

    #detector_info = imo_db_interface.query(
    #        f"{imo_db_datapath}/{instrument}/{freqs[0]}/{det_labels[0]}/detector_info"
    #        )
    #metadata = detector_info.metadata
    #print(metadata)

    #method_list = [func for func in dir(imo_db_interface) if
    #        callable(getattr(imo_db_interface, func))]
    # ['__class__', '__delattr__', '__dir__', '__eq__', '__format__', '__ge__',
    # '__getattribute__', '__gt__', '__hash__', '__init__',
    # '__init_subclass__', '__le__', '__lt__', '__ne__', '__new__',
    # '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__',
    # '__str__', '__subclasshook__', 'get_list_of_data_files',
    # 'get_queried_data_files', 'get_queried_entities',
    # 'get_queried_quantities', 'query', 'query_data_file', 'query_entity',
    # 'query_quantity']
    #print(method_list)

    imo = LitebirdImo(imo_db_interface, imo_version, instrument) 
    det_labels = imo.get_channel_dets(freqs[0])
    #print(det_labels)
    # Pythonic way to get the methods defined within a given class
    #method_list = [func for func in dir(imo) if callable(getattr(imo, func))]
    # ['__class__', '__delattr__', '__dir__', '__eq__', '__format__', '__ge__',
    # '__getattribute__', '__gt__', '__hash__', '__init__',
    # '__init_subclass__', '__le__', '__lt__', '__ne__', '__new__',
    # '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__',
    # '__str__', '__subclasshook__', 'get_channel_dets', 'get_channel_names',
    # 'get_detector_bandwidth', 'get_detector_frequency', 'get_detector_fwhm',
    # 'get_detector_property']
    #print(method_list)
    det_bandwidth = imo.get_detector_bandwidth(freqs[0])
    # 14.0
    #print(det_bandwidth)
    det_names  = imo.get_channel_names()
    # outputs:
    # ['L1-040', 'L2-050', 'L1-060', 'L3-068', 'L2-068', 'L4-078', 'L1-078',
    # 'L3-089', 'L2-089', 'L4-100', 'L3-119', 'L4-140']
    #print(det_names)
    det_fwhm = imo.get_detector_fwhm(freqs[0])
    # 51.1
    #print(det_fwhm)

    fsamps = [imo.get_detector_property(freqs[0], det, 'sampling_rate_hz') for det
            in det_labels]
    #gain   = [imo.get_detector_property(freqs[0], det, 'gain') for det in
    #        det_labels]
 
    
    #print(det_labels)


    ctod.init_file(freqs[0], od, mode='w')

    for det_idx, det_label in enumerate(det_labels):
        #prefix = f'{d}'.zfill(6) + '/' + label 
        prefix = f"{scan}".zfill(6) + "/" + det_label
        #---------------------------------------------
        # Adding fields
        #---------------------------------------------
        #ctod.add_field(f'{prefix}/theta', chunk[field_idx][i, :, 0])
        #ctod.add_field(f'{prefix}/phi', chunk[field_idx][i, :, 1])
        # TODO: To test mapmaking remove huffman compression
        ctod.add_field(prefix + "/pix", pixels[det_idx], huffman)
        ctod.add_field(prefix + "/tod", tod_tot[det_idx])
        ctod.add_field(prefix + "/psi", psi[det_idx], huffman)
        # Getting Scalars, namely:
        # gain, sigma0, fknee, alpha
        detector_info = imo_db_interface.query(
                f"{imo_db_datapath}/{instrument}/{freqs[0]}/{det_label}/detector_info"
                )
        metadata = detector_info.metadata
        fknee      = metadata["fknee_mhz"]
        # Alpha 
        alpha      = metadata["alpha"]
        # Making up some values
        gain   = 1.0 
        sigma0 = 1.0
        scalars = np.array([gain, sigma0, fknee, alpha])
        ctod.add_field(prefix + '/scalars', scalars)
    #---------------------------------------------
    # Things common for each ces scan
    prefix = 'common'
    # Detector Labels
    ctod.add_field(prefix + '/det',      np.string_(det_labels))
    #ctod.add_field(prefix + '/datatype', np.string_(datatype))
    #ctod.add_field(prefix + '/npsi',     npsi)
    ctod.add_field(prefix + '/nside',    nside)
    # Polarization angles
    if instrument in ('LFT', 'HFT'):
        polang = [0 + 45 * int(it.split('_')[3][0] == 'U') + 90 *
                int(it.split('_')[3][1] == 'B') for it in det_labels] 
    else:
        polang = [int(it.split('_')[3][:1]) + 90 *
                int(it.split('_')[3][2] == 'B') for it in det_labels]
    print(polang)
    #polang = np.array([0, 0, 0, 0])
    ctod.add_field(prefix + '/polang',   polang)
    # Main Beam
    mbang  = np.array([0.0] * ndets)
    ctod.add_field(prefix + '/mbang',    mbang)
    # Sampling Rate
    ctod.add_field(prefix + '/fsamp',    fsamp)
    #---------------------------------------------
    # TODO: finish by adding the Huffmann compression
    #---------------------------------------------
    print(f"Running finalize_chunk on file: {lbdata_files[k]}")
    ctod.finalize_chunk(f'{od}'.zfill(6))
    print("finalize_chunk has finished")
    ctod.finalize_file()






if __name__ == '__main__':
    start_time = time.time()
    print("Script has started!")
    main()
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Script run time: {total_time:.2f} s")




