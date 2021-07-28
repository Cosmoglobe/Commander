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
from pathlib import Path
# Getting full path to Mathew's library as an object
commander_tools_path = Path(__file__).absolute().parents[2].joinpath('python','commander_tools').resolve()
# Appending the path to `PYTHONPATH`, so no need to modify it externally (in your `.bashrc` etc.)
sys.path.append(str(commander_tools_path))
# Importing necessary modules from Mathew's library 
from tod_tools import commander_tod as comm_tod
from tod_tools import huffman
#---------------------------------------------

"""
The file structure suitable for Commander3 (after deciphering `devel` branch)
file/
    scan/
        det/
                [x] scalars = np.array([gain[0], sigma0, fknees, alpha])
                [x] pix = hp.ang2pix() <= the compressed pixels
                [ ] psi <= huffman compressed array of values (psiDigitize, look down) 
                    comm_tod.add_field(prefix + '/psi', psiArray, compArray)
                    compArray = [lfi.psiDigitize, lfi.huffman]
                [x] flag <= compressed accepted/rejected data point, all zeros for len(tod)
                [x] tod <= just tods
        common/
                [x] vsun <= velocity of the Earth wrt Sun in Galactic Coordinates, x,y,z, put np.array([0,0,0]) for now
                [x] time <= take np.array([time[0], 0, 0]) 
                [x] huffsym <= for flags, pix and the psi generated huffman symbols (for decoding)
                [x] hufftree <= for flags, pix and psi generated huffman tree (encoded values)
    common/
        [x] det      <= list of detectors, e.g. '27M, 27S, 28M, 28S'
        [x] datatype <= write 'QUIET'
        [x] nside    <= single number
        [x] npsi     <= the number of bins you descretize the polarization angle to, e.g. 2048 or 4096 
                (the psi angle should be descretized somehow, lfi.psiDigitize variable is a clue):
                Inside file the lfi.py: 
                https://github.com/Cosmoglobe/Commander/blob/devel/commander3/python/commander_tools/tod_tools/lfi.py
                npsi = 4096
                psiDigitize = ['digitize', {'min': 0, 'max': 2*np.pi, 'nbins': npsi}]
        [x] fsamp <= sampling frequency, one value <= I have `samprate` instead, assume it is the same thing
        [x] polang <= planck only thing? np.array([0, 0, 0, 0])
        [x] nbang <= planck only thing? np.array([0, 0, 0, 0])
        
"""

'''
The `commander_tod(self, outPath, name, version=None, dicts=None, overwrite=False)` needs to be instantiated.
The `init_file(self, freq, od, mode='r')` is what you need to create a file. Instead of ODs we will use CES
The `add_field(self, fieldName, data, compression=None)` is what you need to create a new field inside the file.
'''
def main():
    print("Script has started")
    parser = argparse.ArgumentParser(description='Adjust Level3 QUIET data for Commander3')
    # 
    make_od()

def make_od():
    # 
    level3_dir = Path('/mn/stornext/d16/cmbco/bp/maksym/quiet/data/Q/ces/patch_gc') 
    output_dir = Path('/mn/stornext/d16/cmbco/bp/maksym/quiet/data/Q/ces/patch_gc/output')
    if not Path.is_dir(output_dir):
        Path.mkdir(output_dir)
    version = np.string_('0.0.1')
    # Getting file names inside specified directory and removing the path component
    level3_data_files = sorted(level3_dir.rglob('*.hdf')) 
    level3_data_files = [data_file.name for data_file in level3_data_files]
    #print(level3_data_files[0])
    # Retrieving CES values from the file names 
    compiled_pattern = re.compile('[\d]')
    level3_ces_nums = [int("".join(compiled_pattern.findall(data_file))) for data_file in level3_data_files] 
    #print(level3_ces[0])

    # TODO: 
    # [x] Retrieve polarization angle from `pointings` 
    # [x] Store these (together with `pixels_idx`?) in Huffmann compression form. 
    #     Look into lfy.py and/or others, because it is not just string, but a list 
    #     containing a dictionary
    #     Note: Needed to change a couple of lines inside commander_tod.py to make it work
    # [x] Change the `commander_tod.py` so it will now write not only numeric 
    #     frequencies in the name of the file
    # [ ] Create only fields as described above (and/or below) ignoring others for now
    # [ ] Scale this to the entire `patch_gc`, i.e. include looping in parallel
    # [ ] Scale this even further to include all patches for the same frequency 
    #     (just dump all CES into one dir? Each CES seems to be unique anyway)
    # [ ] Include W band as well
    # [ ] Include `filelist`s production
    # Retrieving data from old Level3 files 
    readin_file = h5py.File(level3_dir / level3_data_files[0], 'r')
    # Things to include per detector
    alpha        = np.array(readin_file.get('alpha'))
    fknee        = np.array(readin_file.get('fknee'))
    gain         = np.array(readin_file.get('gain'))
    sigma0       = np.array(readin_file.get('sigma0'))
    tods         = np.array(readin_file.get('tod'))
    tp           = np.array(readin_file.get('tp'))
    # Retrieving pointings which will be compressed
    pointing     = np.array(readin_file.get('point'))
    point_objrel = np.array(readin_file.get('point_objrel'))
    # Things to include into common group 
    apex_time    = np.array(readin_file.get('apex/time'))
    apex_value   = np.array(readin_file.get('apex/value'))
    bias_time    = np.array(readin_file.get('bias/time'))
    bias_value   = np.array(readin_file.get('bias/value'))
    coord_sys    = np.array(readin_file.get('coord_sys'))
    corr         = np.array(readin_file.get('corr'))
    corr_freqs   = np.array(readin_file.get('corr_freqs'))
    cryo_time    = np.array(readin_file.get('cryo/time'))
    cryo_value   = np.array(readin_file.get('cryo/value'))
    decimation   = np.array(readin_file.get('decimation'))
    diode_stats  = np.transpose(np.array(readin_file.get('diode_stats')))
    encl_time    = np.array(readin_file.get('encl/time'))
    encl_value   = np.array(readin_file.get('encl/value'))
    filter_par   = np.transpose(np.array(readin_file.get('filter_par')))
    nside        = np.array(readin_file.get('nside'))
    orig_point   = np.array(readin_file.get('orig_point'))
    peri_time    = np.array(readin_file.get('peri/time'))
    peri_value   = np.array(readin_file.get('peri/value'))
    pixels_idx   = np.array(readin_file.get('pixels'))
    samprate     = np.array(readin_file.get('samprate'))
    scanfreq     = np.array(readin_file.get('scanfreq'))
    stats        = np.array(readin_file.get('stats'))
    time1        = np.array(readin_file.get('time'))
    time_gain    = np.array(readin_file.get('time_gain'))
    # 
    phi          = pointing[:,:,0]
    theta        = pointing[:,:,1]
    psi          = pointing[:,:,2]
    pixels       = hp.ang2pix(nside, theta, phi)
    # 
    npsi = 4096
    psiDigitize = ['digitize', {'min': 0, 'max': 2*np.pi, 'nbins': npsi}]
    datatype = 'QUIET'
    # Writing data into new file(s)
    # Initialising tod object
    ctod = comm_tod.commander_tod(output_dir, 'QUIET', version, dicts=None, overwrite=False)
    # Huffmann compression
    huffman = ['huffman', {'dictNum':1}]
    # Creating new file
    #for ces in level3_ces:
    #    comm_tod.init_file(1, ces, mode='w')
    ctod.init_file('Q', level3_ces_nums[0], mode='w')
    # TODO: figure out whether this is the correct interpretation/nomenclature
    diode_labels = ['Qp', 'Up', 'Um', 'Qm']
    #detector_labels = ''
    det_list = []
    # Adding new fields to a file
    for det in tqdm(range(0, 19, 1)):
        #prefix = f'{level3_ces[0]}'.zfill(6) + '/' + str(det+1).zfill(2)
        #ctod.add_field(prefix + '/psi', psi[det])
        #ctod.add_field(prefix + '/pix', pixels[det])#, huffman)
        # In level3 files the `pixels_idx` were called `pixels`
        #ctod.add_field(prefix + '/pixels_idx', pixels_idx[det])
        for diode in range(0, 4, 1):
            label  = str(det+1).zfill(2) + f'{diode_labels[diode]}'
            prefix = f'{level3_ces_nums[0]}'.zfill(6) + '/' + label 
            '''
            [x] scalars = np.array([gain[0], sigma0, fknees, alpha])
            [x] pix = hp.ang2pix() <= the compressed pixels
            [ ] psi <= huffman compressed array of values (psiDigitize). relevant pieces of code: 
                    -- From `.../commander3/python/commander_tools/tod_tools/lfi.py`:
                    psiDigitize = ['digitize', {'min':0, 'max':2*np.pi,'nbins':npsi}]
                    -- From `.../commander3/todscripts/lfi/lfitohdf5.py`:
                    psiArray = fileName[str(horn) + hornType + '/PSI'][pid_start:pid_end] + r.angle_ref(fileName[str(horn) + hornType + '/THETA'][pid_start:pid_end], fileName[str(horn) + hornType + '/PHI'][pid_start:pid_end]) + math.radians(rimo[1].data.field('psi_pol')[rimo_i])
                    if(len(psiArray) > 0):
                        psiArray = np.where(psiArray < 0, 2*np.pi + psiArray, psiArray)
                        psiArray = np.where(psiArray >= 2*np.pi, psiArray - 2*np.pi, psiArray)
                        compArray = None
                        if not args.no_compress:
                            compArray = [lfi.psiDigitize, lfi.huffman]
                        comm_tod.add_field(prefix + '/psi', psiArray, compArray)
            [x] flag <= compressed accepted/rejected data point, all zeros for len(tod)
            [x] tod <= just tods
            '''
            scalars = np.array([gain[det+diode][0], sigma0[det+diode], fknee[det+diode], alpha[det+diode]])
            ctod.add_field(prefix + '/scalars', scalars)
            ctod.add_field(prefix + '/pix',     pixels[det], huffman)
            #ctod.add_field(prefix + '/psi',     psi[det])
            flag = np.zeros_like(tods[det+diode])
            ctod.add_field(prefix + '/flag',    flag, huffman)
            ctod.add_field(prefix + '/tod',     tods[det+diode])
            #
            '''
            ctod.add_field(prefix + '/psi',         psi[det])
            ctod.add_field(prefix + '/pixels_idx',  pixels_idx[det])
            ctod.add_field(prefix + '/alpha',       alpha[det+diode])
            ctod.add_field(prefix + '/fknee',       fknee[det+diode])
            ctod.add_field(prefix + '/gain',        gain[det+diode])
            ctod.add_field(prefix + '/sigma0',      sigma0[det+diode])
            ctod.add_field(prefix + '/tp',          tp[det+diode])
            ctod.add_field(prefix + '/diode_stats', diode_stats[det+diode])
            ctod.add_field(prefix + '/filter_par',  filter_par[det+diode])
            '''
            det_list.append(label)
            #if (det == 18 and diode == 3):
            #    detector_labels += label 
            #else:
            #    detector_labels += label + ", "

    # Things common for each scan
    # There are `common` and `scan/common`
    # For common:
    prefix = 'common'
    #print(detector_labels)
    #print(det_list)
    #print(np.string_(det_list))
    # Converting string to the byte string?
    ctod.add_field(prefix + '/det',      np.string_(det_list))
    ctod.add_field(prefix + '/datatype', np.string_(datatype))
    ctod.add_field(prefix + '/npsi',     npsi)
    ctod.add_field(prefix + '/nside',    nside)
    polang = np.array([0, 0, 0, 0])
    nbang  = np.array([0, 0, 0, 0])
    ctod.add_field(prefix + '/polang',   polang)
    ctod.add_field(prefix + '/nbang',    nbang)
    ctod.add_field(prefix + '/fsamp',    samprate)
    # For `scan/common`
    '''
    [x] vsun <= velocity of the Earth wrt Sun in Galactic Coordinates, x,y,z, put np.array([0,0,0]) for now
    [x] time <= take np.array([time[0], 0, 0]) 
    [x] huffsym <= for flags, pix and the psi generated huffman symbols (for decoding)
    [x] hufftree <= for flags, pix and psi generated huffman tree (encoded values)
    '''
    prefix = f'{level3_ces_nums[0]}'.zfill(6) + '/common'
    vsun = np.array([0, 0, 0])
    ctod.add_field(prefix + '/vsun', vsun)
    ctod.add_field(prefix + '/time', np.array([time1[0], 0, 0]))
    '''
    ctod.add_field(prefix + '/apex/time',  apex_time)
    ctod.add_field(prefix + '/apex/value', apex_value)
    ctod.add_field(prefix + '/bias/time',  bias_time)
    ctod.add_field(prefix + '/bias/value', bias_value)
    ctod.add_field(prefix + '/coord_sys',  coord_sys)
    ctod.add_field(prefix + '/corr',       corr)
    ctod.add_field(prefix + '/corr_freqs', corr_freqs)
    ctod.add_field(prefix + '/cryo/time',  cryo_time)
    ctod.add_field(prefix + '/cryo/value', cryo_value)
    ctod.add_field(prefix + '/decimation', decimation)
    ctod.add_field(prefix + '/encl/time',  encl_time)
    ctod.add_field(prefix + '/encl/value', encl_value)
    ctod.add_field(prefix + '/orig_point', orig_point)
    ctod.add_field(prefix + '/peri/time',  peri_time)
    ctod.add_field(prefix + '/peri/value', peri_value)
    ctod.add_field(prefix + '/samprate',   samprate)
    ctod.add_field(prefix + '/scanfreq',   scanfreq)
    ctod.add_field(prefix + '/stats',      stats)
    ctod.add_field(prefix + '/time_gain',  time_gain)
    '''
    # [[None, {'Dict': 1}] ...]

    #writeout_file = h5py.File(output_dir.joinpath("QUIET_001_000001.h5")  , 'r')
    #print(np.array(writeout_file.get('common/time')))
    print("Running finalize_chunk")
    ctod.finalize_chunk(f'{level3_ces_nums[0]}'.zfill(6))
    print("finalize_chunk has finished")
    ctod.finalize_file()

if __name__ == '__main__':
    start_time = time.time()
    main()
    end_time = time.time()
    run_time = end_time - start_time
    print(f"Script run time: {run_time}s")



