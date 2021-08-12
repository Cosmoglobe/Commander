#---------------------------------------------
import sys
import os
import numpy as np
import healpy as hp        
import h5py
import math
import time
import argparse
import re
from tqdm import tqdm
import multiprocessing as mp
from joblib import cpu_count, Parallel, delayed, parallel_backend
from pathlib import Path
# Getting full path to Mathew's library as an object
commander_tools_path = Path(__file__).absolute().parents[2].joinpath('python','commander_tools').resolve()
# Appending the path to `PYTHONPATH`, so no need to modify it externally (in your `.bashrc` etc.)
sys.path.append(str(commander_tools_path))
# Importing necessary modules from Mathew's library 
from tod_tools import commander_tod as comm_tod
from tod_tools import huffman
#---------------------------------------------

def main():
    """
    Main method of the script
    """
    #nprocs = mp.cpu_count()
    #nprocs = joblib.cpu_count()
    nprocs = cpu_count()
    level3_dir = Path('/mn/stornext/d16/cmbco/bp/maksym/quiet/data/Q/ces/patch_gc') 
    output_dir = Path('/mn/stornext/d16/cmbco/bp/maksym/quiet/data/Q/ces/patch_gc/output')
    if not Path.is_dir(output_dir):
        Path.mkdir(output_dir)
    version = np.string_('0.0.2')
    freqs = ['Q']
    #---------------------------------------------
    # Retrieving data
    #---------------------------------------------
    # Getting file names inside specified directory and removing the path component
    level3_data_files = sorted(level3_dir.rglob('*.hdf')) 
    level3_data_files = [data_file.name for data_file in level3_data_files]
    # Retrieving CES values from the file names 
    compiled_pattern = re.compile('[\d]')
    level3_ces_nums = [int("".join(compiled_pattern.findall(data_file))) 
            for data_file in level3_data_files] 
    #---------------------------------------------
    with parallel_backend(backend="multiprocessing", n_jobs=nprocs):
        manager = mp.Manager()
        # Initialising tod object
        dicts = {freqs[0]:manager.dict()}#, 44:manager.dict(), 70:manager.dict()}
        ctod = comm_tod.commander_tod(output_dir, 'QUIET', version, dicts=dicts, overwrite=True)
        #
        x = Parallel(verbose=2)(delayed(make_od)
                (level3_dir, level3_data_files, level3_ces_nums, ctod, 
                    version, freqs, dicts, k) for k in range(len(level3_ces_nums)))
        #make_od(level3_dir, output_dir, version, freqs, dicts)
        # making filelist
        ctod.make_filelists()


def make_od(level3_dir, level3_data_files, 
        level3_ces_nums, ctod, version, freqs, dicts, k):
    """
    Method to process one file/CES
    """
    # Working with all the files for a given patch
    #for i in tqdm(range(3)):#len(level3_ces_nums)):
    # Retrieving data from old Level3 files 
    readin_file = h5py.File(level3_dir / level3_data_files[k], 'r')
    print(f"Working with file: {level3_data_files[k]}")
    # Things to include per detector
    alpha        = np.array(readin_file.get('alpha'))
    fknee        = np.array(readin_file.get('fknee'))
    gain         = np.array(readin_file.get('gain'))
    sigma0       = np.array(readin_file.get('sigma0'))
    tods         = np.array(readin_file.get('tod'))
    tp           = np.array(readin_file.get('tp'))
    # Things to include into common group 
    coord_sys    = np.array(readin_file.get('coord_sys'))
    nside        = np.array(readin_file.get('nside'))
    samprate     = np.array(readin_file.get('samprate'))
    scanfreq     = np.array(readin_file.get('scanfreq'))
    time_vals    = np.array(readin_file.get('time'))
    time_gain    = np.array(readin_file.get('time_gain'))
    # Retrieving pointings which will be compressed
    pointing     = np.array(readin_file.get('point'))
    point_objrel = np.array(readin_file.get('point_objrel'))
    # Converting pointings to pixels
    phi          = pointing[:,:,0]
    theta        = pointing[:,:,1]
    psi          = pointing[:,:,2]
    pixels       = hp.ang2pix(nside, theta, phi)
    #---------------------------------------------
    # Writing data to a file
    #---------------------------------------------
    # Huffmann compression
    huffman = ['huffman', {'dictNum':1}]
    # Digitization values for \psi 
    npsi = 4096
    psiBins = np.linspace(0, 2*np.pi, npsi)
    datatype = 'QUIET'
    det_list = []
    # Creating new file
    ces = level3_ces_nums[k]
    ctod.init_file(freqs[0], ces, mode='w')
    #---------------------------------------------
    # Looping through 19 amplifiers (1 ampl has 4 diodes)
    # and adding new fields to a file
    i = 0
    #for det in tqdm(range(0, 19, 1)):
    for det in range(0, 19, 1):
        label  = str(det+1).zfill(2) #+ f'{diode_labels[diode]}'
        prefix = f'{ces}'.zfill(6) + '/' + label 
        # Digitizing \psi
        if(len(psi[det]) > 0):
            psiArray = np.where(psi[det] < 0, 2*np.pi + psi[det], psi[det])
            psiArray = np.where(psi[det] >= 2*np.pi, psi[det] - 2*np.pi, psi[det])
            psiIndices = np.digitize(psiArray, psiBins)
        # This field is concatination of 4 tods (i.e. 4 diodes to form amplifier)
        diodes_tods    = []
        diodes_scalars = []
        diodes_names   = []
        for diode in range(0, 4, 1):
            diodes_tods.append(tods[i])
            diodes_scalars.append(np.array([gain[i][0], sigma0[i], fknee[i], alpha[i]]))
            diodes_names.append(f'ref{diode}')
            i = i + 1
        diodes_tods    = np.array(diodes_tods)
        diodes_scalars = np.array(diodes_scalars)
        diodes_flags   = np.zeros_like(diodes_tods)
        #---------------------------------------------
        # Adding fields
        ctod.add_field(prefix + '/psi',    psiIndices,  huffman)#, [psiDigitize, huffman])
        ctod.add_field(prefix + '/pix',    pixels[det], huffman)
        # TODs should be 4xN_tod per amplifier
        ctod.add_field(prefix + '/diodes',   diodes_tods)
        # Scalars should be 4x4 per amplifier
        ctod.add_field(prefix + '/dscalars', diodes_scalars)
        # Flags for accepted/rejected vals in Commander3
        ctod.add_matrix(prefix + '/dflag',   diodes_flags, diodes_names, huffman)
        #
        det_list.append(label)
    #---------------------------------------------
    # Things common for each ces scan
    prefix = 'common'
    ctod.add_field(prefix + '/det',      np.string_(det_list))
    ctod.add_field(prefix + '/datatype', np.string_(datatype))
    ctod.add_field(prefix + '/npsi',     npsi)
    ctod.add_field(prefix + '/nside',    nside)
    polang = np.array([0, 0, 0, 0])
    nbang  = np.array([0, 0, 0, 0])
    ctod.add_field(prefix + '/polang',   polang)
    ctod.add_field(prefix + '/nbang',    nbang)
    ctod.add_field(prefix + '/fsamp',    samprate)
    #---------------------------------------------
    '''
    Thing commong for each detector

    [x] vsun <= velocity of the Earth wrt Sun in Galactic Coordinates, x,y,z, put np.array([0,0,0]) for now
    [x] time <= take np.array([time[0], 0, 0]) 
    [x] huffsym <= for flags, pix and the psi generated huffman symbols (for decoding)
    [x] hufftree <= for flags, pix and psi generated huffman tree (encoded values)
    '''
    prefix = f'{ces}'.zfill(6) + '/common'
    vsun = np.array([0, 0, 0])
    ctod.add_field(prefix + '/vsun', vsun)
    ctod.add_field(prefix + '/time', np.array([time_vals[0], 0, 0]))
    #---------------------------------------------
    print("Running finalize_chunk on file: {level3_data_files[k]}")
    ctod.finalize_chunk(f'{ces}'.zfill(6))
    print("finalize_chunk has finished")
    ctod.finalize_file()



if __name__ == '__main__':
    start_time = time.time()
    print("Script has started!")
    main()
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Script run time: {total_time:.2f} s")

