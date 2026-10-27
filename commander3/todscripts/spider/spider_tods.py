#================================================================================
#
# Copyright (C) 2024 Institute of Theoretical Astrophysics, University of Oslo.
#
# This file is part of Commander3.
#
# Commander3 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Commander3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Commander3. If not, see <https://www.gnu.org/licenses/>.
#
#================================================================================

from commander_tools.tod_tools.spider import spider
from commander_tools.tod_tools import commander_tod as tod
import argparse
import multiprocessing as mp
import os
import numpy as np
import math
from astropy.io import fits
import healpy as hp
import sys
import random
import h5py
import scipy.signal as sig
import scipy.stats as stats

import spider_analysis as sa

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--out-dir', type=str, action='store', default=os.getcwd(), help='path to output data structure you want to generate')

    parser.add_argument('--num-procs', type=int, action='store', default=1, help='number of processes to use')

    parser.add_argument('--freqs', type=int, nargs='+', default=spider.freqs, help='which spider frequencies to operate on')

    parser.add_argument('--chunks', type=int, nargs=2, default=[255, 2088], help='the data chunks to operate on')

    parser.add_argument('--default_scalars', type=float, nargs=3, default=[0.001, 0.15, -2], help="default scalars for Commander 'oof' model: sigma0, fknee, alpha. Gain is read in by --gain-file")

    parser.add_argument('--gain-file', type=str, default='/mn/stornext/u3/sureng/spider_tools/analysis/config/bolotables/cal_fwhm_29.txt', help='SPIDER-style bolotable with a cal column (ADU / uK_CMB per detector)') 

    parser.add_argument('--no-compress', action='store_true', default=False, help='Produce uncompressed data output')

    parser.add_argument('--restart', action='store_true', default=False, help="restart from a previous run that didn't finish")

    parser.add_argument('--produce-filelist', action='store_true', default=False, help='force the production of a filelist even if only some files are present')

    parser.add_argument('--downsample', action='store', type=int, default=5, help='downsampling factor to be used in the TODs')

    args = parser.parse_args()

    args.version = 5

    U = sa.Unifile(default_latest=True, flight=1)

    args.fsamp = U.freq

    random.seed()

    os.environ['OMP_NUM_THREADS'] = '1'

    pool = mp.Pool(processes=args.num_procs)
    manager = mp.Manager()

    U = sa.Unifile(default_latest=True, flight=1)    

    # Grab every chunk from the start (inclusive) to the end (exclusive)
    chunks = range(args.chunks[0], args.chunks[1])
    dicts = {90:manager.dict(), 150:manager.dict()} # ?

    # Grab all the good SPIDER detectors in lists per frequency 
    args.dets = {90:list(U.good_channels(['X2', 'X4', 'X6'])), 150:list(U.good_channels(['X1', 'X3', 'X5']))}

    # Generate the default gains per detector from the SPIDER bolotable
    # We save this in args because that is already being broadcasted to all processes
    args.gains = read_bolotable_gains(args.gain_file)

    # Check for missing detectors. Very pythonic lol
    missing = [det for freq in args.freqs for det in args.dets[freq] if det not in args.gains]

    # If any are missing, report it to the user. This shouldn't happen in practice
    if missing:
        raise ValueError(str(len(missing)) + ' detectors have no gain in ' \
                         + args.gain_file)

    args.comp_array = ['huffman']

    #write detlist
    for freq in args.freqs:
        f = open(os.path.join(args.out_dir, spider.detlist_name(freq)), 'w')
        for det in args.dets[freq]:
            f.write(det + '\n')
        f.close()

    comm_tod = tod.commander_tod(args.out_dir, 'spider', args.version, dicts, not args.restart)

    x = [[pool.apply_async(make_chunk, args=[comm_tod, freq, chunk, args]) for freq in args.freqs] for chunk in chunks]


    for res1 in np.array(x):
        for res in res1:
            res.get()

    pool.close()
    pool.join()

    # Partial runs shouldn't automatically make a filelist
    if((args.chunks[0] == 255 and args.chunks[1] == 2088) or args.produce_filelist):
        comm_tod.make_filelists()

def block_reshape(x, ds, nsamps):
    # Drop any remainder samples and reshape by adding a downsample dimension
    x = np.asarray(x)
    return x[:nsamps*ds].reshape((nsamps, ds) + x.shape[1:])

def read_bolotable_gains(filename):
    with open(filename) as f:
        header = f.readline()

    # grab the columns index from the file
    columns = [c.strip() for c in header.split(':', 1)[1].split(',')]
    cal_col = columns.index('cal')

    # Iterate over the elements and save them in a gains dict, with the commander gain convention and units
    gains = {}
    for det, cal in np.genfromtxt(filename, dtype=str, usecols=(0, cal_col), comments='#', ndmin=2):
        cal = float(cal)
        if np.isfinite(cal) and cal != 0:
            # convert the gains from uK / ADU to -ADU / K
            # This makes gains a positive quantity, to match the sign flipped TOD
            gains[det] = -1e6 / cal
            

    return gains

def make_chunk(comm_tod, freq, chunk, args):

    print(freq, chunk)

    comm_tod.init_file(freq, chunk, mode='w')

    if(args.restart and comm_tod.exists):
        comm_tod.finalize_file()
        print('Skipping existing file ' + comm_tod.outName)
        return

    U = sa.Unifile(default_latest=True, flight=1)
    Um = sa.UnifileMap(default_latest=True, flight=1) # ? Probably just Um can do all of this
    chunkdict = U.event_partition('time_chunk_10min')[chunk]

    prefix = 'common'
    
    fsamp =  args.fsamp / args.downsample
    comm_tod.add_field(prefix + '/fsamp', fsamp)

    comm_tod.add_field(prefix + '/nside', spider.nside)

    polangs = []
    mbangs = []

    for det in args.dets[freq]:
        polangs.append(0)
        mbangs.append(0)

    comm_tod.add_field(prefix + '/det', spider.detlist_name(freq))
    comm_tod.add_field(prefix + '/mbang', mbangs) #? what is this
    comm_tod.add_field(prefix + '/polang', polangs)

    # now do the parts in the per chunk subdirectory
    # even though each file only has one chunk for spider
    prefix = str(chunk).zfill(6) + '/common'

    comm_tod.add_field(prefix + '/satpos', [0,0,0])
    #harald had numbers here for vsun, I'm not sure where they would get used 
    #or where he got them from
    comm_tod.add_field(prefix + '/vsun', [0,0,0])

    time = U.get_time(**chunkdict)
    comm_tod.add_field(prefix + '/time', time[0])

    #ntod 
    nsamps = int(U.get_sample_count(**chunkdict) / args.downsample)
    comm_tod.add_field(prefix + '/ntod', nsamps)   

    if args.no_compress:
        compArray = []
        psiArray = []
    else:
        compArray = [spider.huffman]
        psiArray = [spider.psiDigitize, spider.huffman]


    for det in args.dets[freq]:

        prefix = str(chunk).zfill(6) + '/' + det

        # Read in the detector gain to add to the scalar list
        gain = args.gains[det]
        scalars = [gain, *args.default_scalars]
        
        # Read in detector pointing for this chunk
        try:
            lon, lat, psi= Um.get_det_coords(det, coord='G', **chunkdict)
        except ValueError as e:
            print('No det pointing + psi in: ', freq, chunk, det)
            # add fake, completely flagged data to fill the chunk
            fake_data = np.ones(nsamps)
            comm_tod.add_field(prefix + '/pix', fake_data, compArray)
            comm_tod.add_field(prefix + '/psi', fake_data, compArray)
            comm_tod.add_field(prefix + '/tod', fake_data)
            comm_tod.add_field(prefix + '/flag', fake_data, compArray)
            comm_tod.add_field(prefix + '/scalars', scalars)

            continue


        # Need to restimate effective beam ? FIX
        ########################################################################
        # psi: can't naively average. Need to handle blocks straddling ±180°
        # First, make each angle psi a point on the unit circle: cos(psi) + i sin(psi)
        psi_ucirc = np.exp(1j*np.deg2rad(psi))

        # Block reshape the psi's to get ready to average
        psi_ucirc_block = block_reshape(psi_ucirc, args.downsample, nsamps)

        # Average down, convert back to degrees. np.angle runs arctan2(Im z, Re z)
        psi_mean = np.rad2deg(np.angle(psi_ucirc_block.mean(axis=1)))

        # SPIDER TOD seems to use a different convention for PSI than Commander. Possibly IAU instead of COSMO
        psi_dec = -1* psi_mean
        psi_rad = np.deg2rad(psi_dec)

        # rewrap from 0 to 2 pi
        psi_rad[psi_rad<=0] += 2*np.pi 

        # There is a chance that psi is zero, which when wrapped will be 2*np.pi, which may be problematic if digitized
        psi_rad = np.clip(psi_rad, 0, 2*np.pi*(1 - 1e-12))
        comm_tod.add_field(prefix + '/psi', psi_rad, psiArray)

        ########################################################################
        # pix: average the pointing unit vectors. First generate them
        vec = hp.ang2vec(lon, lat, lonlat=True)

        # Block reshape and average
        vec_mean = block_reshape(vec, args.downsample, nsamps).mean(axis=1)

        # convert to pix using healpy and the default SPIDER nside of 512
        pix = hp.vec2pix(spider.nside, vec_mean[:, 0], vec_mean[:, 1], vec_mean[:, 2])

        comm_tod.add_field(prefix + '/pix', pix, compArray)

        ########################################################################
        # TOD: First load the TOD and subtract stepstitch
        tod = U.getdata(det, product='dcclean08', **chunkdict)
        stepstitch = U.getdata(det, product='stepstitch07', **chunkdict)
        tod = tod - stepstitch

        # Block reshape and average
        tod = block_reshape(tod, args.downsample, nsamps).mean(axis=1)

        # Mean subtract to remove chunk-to-chunk offsets. Add a negative sign because
        # SPIDER detector data is inversely related to the sky
        tod = -1.*(tod - np.mean(tod))
        comm_tod.add_field(prefix + '/tod', tod)

        ########################################################################
        # flag: Same flags and masks as SPIDER's own mapmaking (get_tot_flag with
        # dcclean08 + stepstitch07). These already contain the stepstitch07 
        # flags, so that product's flag is not read separately.
        flag_comb = U.getdata(det, product='flag_comb08', **chunkdict)
        flag_extra = U.getdata(det, product='flag_extra02', **chunkdict)        

        # Generate binary flags 

        # Use all 0-28 bits but the 23rd bit, which is xjump_reset. 
        # This bit marks where a correction was made earlier in the analysis
        flag_comb_binary = (np.asarray(flag_comb).astype(np.uint64) & np.uint64(0x1F7FFFFF)) != 0
        
        # Use all 0-12 bits 
        flag_extra_binary = (np.asarray(flag_extra).astype(np.uint64) & np.uint64(0x1FFF)) != 0
        flag_tot = (flag_comb_binary | flag_extra_binary).astype(np.uint8)

        # Downsample by reshaping and taking the max per group of args.downsample
        flag_downsampled = block_reshape(flag_tot, args.downsample, nsamps).max(axis=1)
        comm_tod.add_field(prefix + '/flag', flag_downsampled, compArray)

        ########################################################################
        # scalars: Add the default scalar for both frequencies, for all chunks 
        # and detectors

        comm_tod.add_field(prefix + '/scalars', scalars)

        U.raw_close()

    comm_tod.finalize_chunk(chunk)
    comm_tod.finalize_file()

if __name__ == '__main__':
    main()
