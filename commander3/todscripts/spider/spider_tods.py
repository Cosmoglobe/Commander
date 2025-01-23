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

    parser.add_argument('--chunks', type=int, nargs=2, default=[255, 2084], help='the data chunks to operate on')

    parser.add_argument('--no-compress', action='store_true', default=False, help='Produce uncompressed data output')

    parser.add_argument('--restart', action='store_true', default=False, help="restart from a previous run that didn't finish")

    parser.add_argument('--produce-filelist', action='store_true', default=False, help='force the production of a filelist even if only some files are present')

    parser.add_argument('--downsample', action='store', type=int, default=5, help='downsampling factor to be used in the TODs')

    args = parser.parse_args()

    args.version = 4

    U = sa.Unifile(default_latest=True, flight=1)

    args.fsamp = U.freq

    random.seed()

    os.environ['OMP_NUM_THREADS'] = '1'

    pool = mp.Pool(processes=args.num_procs)
    manager = mp.Manager()

    U = sa.Unifile(default_latest=True, flight=1)    

    chunks = range(args.chunks[0], args.chunks[1])
    dicts = {90:manager.dict(), 150:manager.dict()}

    args.dets = {90:list(U.good_channels(['X2', 'X4', 'X6'])), 150:list(U.good_channels(['X1', 'X3', 'X5']))}

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

    if((args.chunks[0] == 255 and args.chunks[1] == 2084) or args.produce_filelist):
        comm_tod.make_filelists()

def make_chunk(comm_tod, freq, chunk, args):

    print(freq, chunk)

    comm_tod.init_file(freq, chunk, mode='w')

    if(args.restart and comm_tod.exists):
        comm_tod.finalize_file()
        print('Skipping existing file ' + comm_tod.outName)
        return

    U = sa.Unifile(default_latest=True, flight=1)
    Um = sa.UnifileMap(default_latest=True, flight=1)
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
    comm_tod.add_field(prefix + '/mbang', mbangs)
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
        
        # pix and psi
        try:
            lon, lat, psi= Um.get_det_coords(det, coord='G', **chunkdict)
        except ValueError as e:
            print('Value error in ', freq, chunk, det)
            # add fake, completely flagged data to fill the chunk
            fake_data = np.ones(nsamps)
            comm_tod.add_field(prefix + '/pix', fake_data, compArray)
            comm_tod.add_field(prefix + '/psi', fake_data, compArray)
            comm_tod.add_field(prefix + '/tod', fake_data)
            comm_tod.add_field(prefix + '/flag', fake_data, compArray)
            comm_tod.add_field(prefix + '/scalars', [0,0,0,0])

            continue


        #psi_dec = sig.decimate(psi, args.downsample)
        psi_dec = psi[::args.downsample]

        psi_rad = np.deg2rad(psi_dec)
        psi_rad[psi_rad<=0] += 2*np.pi 

        comm_tod.add_field(prefix + '/psi', psi_rad, psiArray)

        #lat_down = sig.decimate(lat_gal, args.downsample)
        #lon_down = stats.circmean(np.split(lon_gal, args.downsample), high=360.0, axis=0)
        lat_down = lat[::args.downsample]
        lon_down = lon[::args.downsample]

        pix = hp.ang2pix(spider.nside, lon_down, lat_down, lonlat=True)
        comm_tod.add_field(prefix + '/pix', pix, compArray)


        #TOD
        #tod = sig.decimate(U.getdata(det, product='dcclean08', **chunkdict), args.downsample)
        tod = U.getdata(det, product='dcclean08', **chunkdict)
        stepstitch = U.getdata(det, product='stepstitch07', **chunkdict)

        tod = tod - stepstitch
        tod = tod[::args.downsample]
        tod = tod - np.mean(tod)

        comm_tod.add_field(prefix + '/tod', tod)

        # flag
        flag = U.getdata(det, product='flag_comb08', **chunkdict)
        flag_extra = U.getdata(det, product='flag_extra02', **chunkdict)        
        flag_stepstitch = U.getdata(det + '_flag', product='stepstitch07', **chunkdict)
        
        flag_tot = np.where(np.bitwise_and(np.uint64(flag + 2**32*flag_extra) ,35180601409535) == 0, 0, 1)
               
 
        #downsample flag by taking bitwise or of each n entries
        #nsplit = np.split(flag_tot, args.downsample)
        #flag_downsampled = np.bitwise_or.reduce(nsplit, 0)
        flag_downsampled = flag_tot[::args.downsample]
        comm_tod.add_field(prefix + '/flag', flag_downsampled, compArray)

        # scalars
        scalars = [1, 1, 0.1, -2] # gain, sigma0, fknee, alpha
        comm_tod.add_field(prefix + '/scalars', scalars)

        U.raw_close()

    comm_tod.finalize_chunk(chunk)
    comm_tod.finalize_file()

if __name__ == '__main__':
    main()
