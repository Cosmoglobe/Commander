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


def main():
    """
    Main method of the script 
    """
    lbdata_dir = pathlib.Path("/Users/maksymb/Desktop/litebird_db/data") 
    output_dir = pathlib.Path(__file__).parent.joinpath("test_litebird_sims")
    if not pathlib.Path.is_dir(output_dir):
        pathlib.Path.mkdir(output_dir)

    # Getting file names inside specified directory and removing the path component
    lbdata_files = sorted(lbdata_dir.rglob('*.hdf5')) 
    lbdata_files = [data_file.name for data_file in lbdata_files]

    print(f"{lbdata_files}")


    print(2**22)
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




