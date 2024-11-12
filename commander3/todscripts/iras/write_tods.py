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

import argparse
import multiprocessing as mp
import os
import numpy as np
import math
from astropy.io import fits
from astropy import units as u
import healpy as hp
import sys
import random
import h5py
import scipy.signal as sig
import scipy.stats as stats
from astropy.time import Time
from pathlib import Path

from cosmoglobe.tod_tools import TODLoader

from astropy.coordinates import SkyCoord, FK4

from astropy.coordinates import (
    HeliocentricMeanEcliptic,
    get_body,
    solar_system_ephemeris,
)


NSIDE = 2048
t0 = Time('1981-01-01', scale='utc')
IRAS_DATA_PATH = Path("/mn/stornext/d5/data/duncanwa/IRAS/sopobs_data")

FSAMP = {'12':16, '25':16, '60':8, '100':4}

BANDPASS = {12: 25e12*u.Hz, 25:12e12*u.Hz, 60:5e12*u.Hz, 100:3e12*u.Hz}

gc = SkyCoord(l=0*u.degree, b=0*u.degree, frame='galactic')
Omegas = {1: 14.5,
        2: 12.7,
        3: 13.0,
        4: 11.53,
        5: 12.0,
        6: 12.4,
        7:12.6,
        8:7.2,
        9:6.7,
        10:6.6,
        11:2.8,
        12:4.3,
        13:6.6,
        14:6.1,
        15:6.2,
        16:3.5,
        18:3.6,
        19:2.8,
        21:2.8,
        22:3.1,
        23:2.9,
        24:3.0,
        25:3.2,
        26:1.2,
        27:2.0,
        28:3.1,
        29:2.5,
        30:2.8,
        55:7.1,
        56:14.0,
        57:13.2,
        58:11.2,
        59:11.7,
        60:13.3,
        61:13.5,
        62:10.6,
        31:2.1,
        32:6.4,
        33:5.9,
        34:6.5,
        35:6.3,
        37:6.6,
        38:3.9,
        39:1.4,
        40:3.1,
        41:3.1,
        42:3.4,
        43:3.2,
        44:3.2,
        45:3.2,
        46:2.4,
        47:0.77,
        48:3.1,
        49:2.9,
        50:3.0,
        51:2.7,
        52:2.5,
        53:2.8,
        54:2.0}

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--out-dir', type=str, action='store', default=os.getcwd(), help='path to output data structure you want to generate')

    parser.add_argument('--num-procs', type=int, action='store', default=1, help='number of processes to use')

    parser.add_argument('--freqs', type=int, nargs='+',
            default=np.array([12,25,60,100]), help='which IRAS bands to operate on')

    parser.add_argument('--chunks', type=int, nargs=2, default=[1, 5787], help='the data chunks to operate on')

    parser.add_argument('--no-compress', action='store_true', default=False, help='Produce uncompressed data output')

    parser.add_argument('--restart', action='store_true', default=False, help="restart from a previous run that didn't finish")

    parser.add_argument('--produce-filelist', action='store_true', default=False, help='force the production of a filelist even if only some files are present')

    args = parser.parse_args()

    args.version = 0


    os.environ['OMP_NUM_THREADS'] = '1'

    pool = mp.Pool(processes=args.num_procs)
    manager = mp.Manager()

    chunks = range(args.chunks[0]//10+1, args.chunks[1]//10+1)
    dicts = {12:manager.dict(), 25:manager.dict(), 60:manager.dict(), 100:manager.dict()}

    args.dets = {
        12:np.concatenate((np.arange(23, 31), np.arange(48, 54))),
        25: np.concatenate((np.array([16,18,19,21, 22]), np.arange(40,46))),
        60: np.concatenate((np.array([8,9,10,13,14,15]), np.array([32,33,34,35,37]))),
        100: np.concatenate((np.arange(1, 8), np.arange(56, 62)))
        }

    args.comp_array = ['huffman']

    #write detlist
    for freq in args.freqs:
        f = open(os.path.join(args.out_dir, f'detlist_{str(freq)}.txt'), 'w')
        for det in args.dets[freq]:
            f.write(f'{det:02}\n')
        f.close()

    comm_tod = TODLoader(args.out_dir, 'iras', args.version, dicts, not args.restart)


    if args.num_procs == 1:
        for chunk in chunks:
            for freq in args.freqs:
                make_chunk(comm_tod, freq, chunk, args)
    else:

        x = [[pool.apply_async(make_chunk, args=[comm_tod, freq, chunk, args]) for freq in args.freqs] for chunk in chunks]


        for res1 in np.array(x):
            for res in res1:
                res.get()

    pool.close()
    pool.join()

    if((args.chunks[0] == 255 and args.chunks[1] == 2084) or args.produce_filelist):
        comm_tod.make_filelists()

def make_chunk(comm_tod, freq, chunk, args):

    #if chunk == 273:
    #    return

    #print(freq, chunk)

    comm_tod.init_file(freq, chunk, mode='w')

    if(args.restart and comm_tod.exists):
        comm_tod.finalize_file()
        print('Skipping existing file ' + comm_tod.outName)
        return


    prefix = 'common'
    
    fsamp = FSAMP[str(freq)]
    comm_tod.add_field(prefix + '/fsamp', fsamp)

    comm_tod.add_field(prefix + '/nside', NSIDE)

    polangs = []
    mbangs = []

    for det in args.dets[freq]:
        polangs.append(0)
        mbangs.append(0)



    comm_tod.add_field(prefix + '/det', 'detlist_' + str(freq) + '.txt')
    comm_tod.add_field(prefix + '/mbang', mbangs)
    comm_tod.add_field(prefix + '/polang', polangs)




    compArray = [['huffman', {'dictNum':1}]]
    psiDigitize = ['digitize', {'min':0, 'max':2*np.pi,'nbins':8, 'offset':1}] 
    psiArray = [psiDigitize, ['huffman', {'dictNum':1}]]


    for scan in range((chunk-1)*10 + 1, chunk*10 + 1):

        if scan == 273:
            continue


        found = False

        for det in args.dets[freq]:

            prefix = f'{scan:06}/{det:02}'
            f = f'{IRAS_DATA_PATH}/det_{det:02}/sopobs_{scan:04}.npy'
            if os.path.exists(f):
                t, lon, lat, tod = np.load(f)
                found = True
            else:
                continue


            sorted_inds = np.argsort(t)
            t = t[sorted_inds]
            lon = lon[sorted_inds]
            lat = lat[sorted_inds]

            tod = tod/(Omegas[det]*1e-7)
            tod *= u.W/u.m**2/u.sr
            tod /= BANDPASS[freq]
            tod = tod.to('MJy/sr').value

            

            #print(chunk, det, freq, t.shape)


            # flag
            flag_tot = np.zeros(len(t))
            flag_tot[~np.isfinite(lon)] += 2**0
            flag_tot[~np.isfinite(tod)] += 2**1
             
            tod[flag_tot != 0] = 0
            lon[flag_tot != 0] = 0
            lat[flag_tot != 0] = 0
            comm_tod.add_field(prefix + '/tod', tod)

            # pix and psi
            psi = np.zeros_like(tod)

            comm_tod.add_field(prefix + '/psi', psi, psiArray)


            good_data = flag_tot == 0

            #lon[good_data], lat[good_data] = r(lon[good_data], lat[good_data], lonlat=True)
            sc = SkyCoord(ra=lon[good_data], dec=lat[good_data], unit='deg',
                    equinox='B1950.0', obstime='J1983.5', frame=FK4)
            coords = sc.transform_to(gc)
            lon[good_data] = coords.l.value
            lat[good_data] = coords.b.value
            pix = hp.ang2pix(NSIDE, lon, lat, lonlat=True) 
            comm_tod.add_field(prefix + '/pix', pix, compArray)
 
            comm_tod.add_field(prefix + '/flag', flag_tot, compArray)

            # scalars
            scalars = [1, 1, 0.1, -2] # gain, sigma0, fknee, alpha
            comm_tod.add_field(prefix + '/scalars', scalars)


        if found:
            prefix = str(scan).zfill(6) + '/common'
            time = t[0]*u.s + t0
            time_end = t[-1]*u.s + t0

            earth_pos = get_body("earth", Time(time.mjd, format="mjd")).transform_to(
        HeliocentricMeanEcliptic
    )
            earth_pos = earth_pos.cartesian.xyz.to(u.AU).transpose()

            earth_pos_end = get_body("earth", Time(time_end.mjd, format="mjd")).transform_to(
        HeliocentricMeanEcliptic
    )
            earth_pos_end = earth_pos_end.cartesian.xyz.to(u.AU).transpose()

            comm_tod.add_field(f'{scan:06}/common/satpos', earth_pos)
            comm_tod.add_field(f'{scan:06}/common/earthpos', earth_pos)
            comm_tod.add_field(f'{scan:06}/common/satpos_end', earth_pos_end)
            comm_tod.add_field(f'{scan:06}/common/earthpos_end', earth_pos_end)
            comm_tod.add_field(f'{scan:06}/common/vsun', [0,0,0])
            comm_tod.add_field(prefix + '/time', [time.mjd,0,0])
            comm_tod.add_field(prefix + '/time_end', [time_end.mjd,0,0])

            nsamps = len(t)
            comm_tod.add_field(prefix + '/ntod', nsamps)   

            comm_tod.finalize_chunk(scan)
    comm_tod.finalize_file()

if __name__ == '__main__':
    main()
