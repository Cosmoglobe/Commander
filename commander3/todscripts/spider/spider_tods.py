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
import scipy.signal.decimate

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

    os.environ['OMP_NUM_THREADS'] = 1

    pool = mp.Pool(processes=args.num_procs)
    manager = mp.Manager()

    U = sa.Unifile(default_latest=True, flight=1)    

    chunks = range(args.chunks[0], args.chunks[1])
    dicts = {90:manager.dict(), 150:manager.dict()}

    args.dets = {90:U.good_channels(['X2', 'X4', 'X6']), 150:U.good_channels(['X1', 'X3', 'X5'])}

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

    if((args.chunk[0] == 255 and args.chunk[1] == 2084) or args.produce_filelist):
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
    chunkdict = U.event_partition()[chunk]

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

    comm_tod.add_field(prefix + '/satpos', [0,0])
    #harald had numbers here for vsun, I'm not sure where they would get used 
    #or where he got them from
    comm_tod.add_field(prefix + '/vsun', [0,0,0])

    time = U.get_time(star=chunkdict['start'], end=chunkdict['end'], index_mode=chunkdict['index_mode'])
    comm_tod.add_field(prefix + '/time', time[0])

    #ntod 
    nsamps = U.get_sample_count(**chunkdict) / args.downsample
    comm_tod.add_field(prefix + '/ntod', nsamps)   

    # get telescope boreshight quaternion
    bore_quats, pflag = Um.get_bore_quat('point07', return_flag=True, hack_bore_sign=True, **chunkdict)

    offsets = Um.get_offset_quats(args.dets[freq]) 

    for det, offset in zip(args.dets[freq], offsets):
        prefix = str(chunk).zfill(6) + '/' + det
        #TOD
        tod = decimate(U.getdata(det, product='dcclean08', **chunkdict), params.downsample)

        # pix and psi
        pix, psi = Um.bore2pix(offset, None, bore_quats, nside=512)

        # flag

        # scalars
        scalars = [1, 1, 0.1, -2] # gain, sigma0, fknee, alpha
        comm_tod.add_field(prefix + '/scalars', scalars)



    comm_tod.finalize_chunk(chunk, loadBalance=1)
    comm_tod.finalize_file()



if __name__ == '__main__':
    main()
