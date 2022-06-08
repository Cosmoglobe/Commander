#================================================================================
#
# Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
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

import sys
sys.path.append('/mn/stornext/u3/duncanwa/Commander/commander3/python/')
from commander_tools.tod_tools import commander_tod
import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
import h5py
import healpy as hp

from glob import glob
import multiprocessing as mp
from multiprocessing import Pool
from joblib import Parallel, delayed

import huffman


from scipy.interpolate import interp1d
from joblib import Parallel, delayed
import os, sys

from tqdm import tqdm

from astroquery.jplhorizons import Horizons
from astropy.time import Time
from datetime import datetime
from get_gain_model import get_gain

all_band_labels = np.array([
  'K113', 'K114', 'K123', 'K124',
  'Ka113', 'Ka114', 'Ka123', 'Ka124',
  'Q113', 'Q114', 'Q123', 'Q124',
  'Q213', 'Q214', 'Q223', 'Q224',
  'V113', 'V114', 'V123', 'V124',
  'V213', 'V214', 'V223', 'V224',
  'W113', 'W114', 'W123', 'W124',
  'W213', 'W214', 'W223', 'W224',
  'W313', 'W314', 'W323', 'W324',
  'W413', 'W414', 'W423', 'W424'])


Nobs_array = np.array([12, 12, 15, 15, 20, 20, 30, 30, 30, 30])
fknees = np.array([0.7,  0.5, 1.2, 0.6,
                   1.1, 0.6, 3.0, 4.0,
                   1.3, 2.2, 1.5, 4.0,
                  16.17,15.05, 9.02, 7.47,
                   1.84, 2.39, 46.5, 26.0]) # mHz
fknees *= 1e-3

alphas = np.array([-0.8, -0.7, -1.0, -0.9,
                   -1.0, -1.0, -1.1, -1.1,
                   -1.0, -0.9, -1.1, -1.0,
                   -1.0, -1.0, -1.0, -1.0,
                   -1.0, -1.0, -1.0, -1.0])


# All of the events listed in the explanatory supplement. Each is listed with a
# "beginning" and an ending".
data = np.loadtxt('events2.txt', dtype=str, usecols=(0,1))
t = Time(data)
t.format = 'mjd'
events = t.value
events_flags = np.loadtxt('events2.txt', usecols=(2,3,4,5,6,7,8,9,10,11))
print(events.shape)
print(events_flags[0])
print(events_flags[0].shape)
print(events_flags.shape)


# From Table 2 of Jarosik et al. 2003, "On-orbit radiometer
# characterization", using the higher of the in-flight versus GSFC measurements.


# version 15 uses the pre-calibrated data, gain = 1
# version 16 is an attempt to fix the quaternions
# version 17 is using center = False
# version 18 computes the planet exclusion flags according to Bennett et al.
# version 19 uses the new planet exclusion flags, also uses center = True
# version 20 splits up the data into 25 chunks
# version 21 splits up the data into 24 chunks
# version 22 splits up the data into 24 chunks and uses precalibrated data
# version 24 uses a more accurate gain model
# version 25 is same as version 24, but uses precalibrated data.
# version 26 is same as version 24, but does not recompute planet flags
# version 27 is same as version 26, center = False
# version 28 converts the velocity vector to Galactic coordinates
# version 29 converts the velocity to meters, corrects format of mbang, polang.
# version 30 adds the genflags to the daflags
# version 31 has center=True
# version 'cal' uses the WMAP precalibrated data
# version 34 has fixed a bug in the flag timing
# version 35 makes the time in MJD.
# version 36 reverts to the original planet flagging, since events like solar
# storms are also included in there (see Table 8 of the ExSupp)
# version 37 uses the planet flag based on radii
# version 38 compresses the TODs, adjusts the TODs and gains so that they are all positive
# version 39 uses weeklong scans
# version 40 computes the polarization angles in a different way
# version 41 fixes the unit vectors in the polarization angle calculation
# version 42 uses the radius-based planet flag, changes zipped TOD to ztod
# version 43 uses default flags, changes zipped TOD to ztod
# version 44 changes todtree, todsymb to hufftree2, huffsymb2
# version 45 uses commmander_tools package, and fixes sign flip of timestream to be constant over tod
# version 46 uses JPL Horizons ephemerides directly
# version 47 uses same as above, just uses the precalibrated data
# version 48 reverts to the center=False parameter
# version 49 back to center False

from time import sleep
from time import time as timer

# from https://towardsdatascience.com/how-to-profile-your-code-in-python-e70c834fad89
import cProfile
import pstats
from functools import wraps


def profile(output_file=None, sort_by='cumulative', lines_to_print=None, strip_dirs=False):
    """A time profiler decorator.
    Inspired by and modified the profile decorator of Giampaolo Rodola:
    http://code.activestate.com/recipes/577817-profile-decorator/
    Args:
        output_file: str or None. Default is None
            Path of the output file. If only name of the file is given, it's
            saved in the current directory.
            If it's None, the name of the decorated function is used.
        sort_by: str or SortKey enum or tuple/list of str/SortKey enum
            Sorting criteria for the Stats object.
            For a list of valid string and SortKey refer to:
            https://docs.python.org/3/library/profile.html#pstats.Stats.sort_stats
        lines_to_print: int or None
            Number of lines to print. Default (None) is for all the lines.
            This is useful in reducing the size of the printout, especially
            that sorting by 'cumulative', the time consuming operations
            are printed toward the top of the file.
        strip_dirs: bool
            Whether to remove the leading path info from file names.
            This is also useful in reducing the size of the printout
    Returns:
        Profile of the decorated function
    """

    def inner(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            _output_file = output_file or func.__name__ + '.prof'
            pr = cProfile.Profile()
            pr.enable()
            retval = func(*args, **kwargs)
            pr.disable()
            pr.dump_stats(_output_file)

            with open(_output_file, 'w') as f:
                ps = pstats.Stats(pr, stream=f)
                if strip_dirs:
                    ps.strip_dirs()
                if isinstance(sort_by, (tuple, list)):
                    ps.sort_stats(*sort_by)
                else:
                    ps.sort_stats(sort_by)
                ps.print_stats(lines_to_print)
            return retval

        return wrapper

    return inner



def get_all_ephems(files):
    planets = ['mars', 'jupiter', 'saturn', 'neptune', 'uranus']

    first_file = files[0]
    last_file = files[-1]

    data_first = fits.open(first_file)
    data_last  = fits.open(last_file)
    t2jd = data_first[1].header['time2jd']

    t0 = data_first[2].data['TIME'][0] + t2jd
    tl = data_last[2].data['TIME'][-1] + t2jd
    #tl = 2454928.5     # April 7, 2009, when the JPL ephemerides run out for WMAP
    times = np.array([t0, tl])

    for p in planets:
        t_jds, ras, decs = get_planet_radec(times, p)

        np.savetxt(f'{p}_ephem.txt', np.array([t_jds, ras, decs]).T)

    return

TARGET_ALIASES = {"mars": 499,  # (4, 499)
        "jupiter": 599, # (5, 599)
        "saturn": 699, # (6, 699)
        "uranus": 799, # (7, 799)
        "neptune": 899} # (8, 899)


def get_planet_radec(times, planet):
    '''
    Uses the JPL Horizons API, takes input planet as string, and time as JD
    500@-165 is the label for the WMAP object.
    It is also not entirely clear if we should use the barycenter (n) or the
    "regular" object, (n99).
    '''

    RAs = []
    Decs = []

    mydates_jd = Time(times, format='jd')

    print(mydates_jd[0].iso, mydates_jd[-1].iso)


    query = Horizons(
        id=TARGET_ALIASES[planet],
        id_type="majorbody",
        #location="500@-165",
        location="500@32", # at L2
        #location="500", # geocentric
        epochs={'start': mydates_jd[0].iso,
                'stop': mydates_jd[-1].iso,
                'step': '1h'}
    )
    ephemerides = query.ephemerides(quantities='1,2')
    RAs = ephemerides['RA'].value
    Decs = ephemerides['DEC'].value
    t_jd = ephemerides['datetime_jd'].value

    t_jd = np.array(t_jd)
    RAs = np.array(RAs)
    Decs = np.array(Decs)

    return t_jd, RAs, Decs



def get_ephem(time, planet):
    # Obtained these from https://ssd.jpl.nasa.gov/horizons.cgi
    t_p, ra_p, dec_p = np.loadtxt(f'{planet}_ephem.txt').T

    f_ra = interp1d(t_p, ra_p, fill_value='extrapolate')
    f_dec = interp1d(t_p, dec_p, fill_value='extrapolate')

    return f_ra(time), f_dec(time)


def get_flags(data, test=False, center=True, bands=np.arange(10)):

    Nobs_array = np.array([12, 12, 15, 15, 20, 20, 30, 30, 30, 30])

    t2jd = data[1].header['TIME2JD']

    quat = data[1].data['QUATERN']
    Nobs_arr = Nobs_array[bands]

    ll_A, ll_B, p_A, p_B, t_list = quat_to_sky_coords(quat, lonlat=True,
        center=center, ret_times=True,
            coord_out='C', Nobs_array = Nobs_arr)

    time_majorframe = data[2].data['TIME'] + t2jd

    daflags = data[2].data['daflags']

    planets = ['mars', 'jupiter', 'saturn', 'uranus', 'neptune']
    radii   = np.array([
              [2.0,     3.0,      2.0,      2.0,       2.0], #K (yr!=2)
              [1.5,     2.5,      1.5,      1.5,       1.5], #Ka
              [1.5,     2.5,      1.5,      1.5,       1.5], #Q1
              [1.5,     2.5,      1.5,      1.5,       1.5], #Q2
              [1.5,     2.2,      1.5,      1.5,       1.5], #V1
              [1.5,     2.2,      1.5,      1.5,       1.5], #V2
              [1.5,     2.0,      1.5,      1.5,       1.5], #W1
              [1.5,     2.0,      1.5,      1.5,       1.5], #W2
              [1.5,     2.0,      1.5,      1.5,       1.5], #W3
              [1.5,     2.0,      1.5,      1.5,       1.5], #W4
              ])
    #radii   = np.array([
    #          [7,   7,      7,      7,       7], #K (yr!=2)
    #          [7,   7,      7,      7,       7], #Ka
    #          [7,   7,      7,      7,       7], #Q1
    #          [7,   7,      7,      7,       7], #Q2
    #          [7,   7,      7,      7,       7], #V1
    #          [7,   7,      7,      7,       7], #V2
    #          [7,   7,      7,      7,       7], #W1
    #          [7,   7,      7,      7,       7], #W2
    #          [7,   7,      7,      7,       7], #W3
    #          [7,   7,      7,      7,       7], #W4
    #          ])


    radii = radii*(np.pi/180)
    # I should be calculating the distances using each individual observation,
    # not the daflags shape, since that's going to be the same size as the major
    # frames, not the individual ones. Additionally, I am using the time for
    # each major frame, not the actual time for each individual observation.
    myflags = []
    daflags_copy = []
    for band in bands:
        myflags.append([])
        daflags_copy.append([])

    for b,band in enumerate(bands):
        myflags[b] = np.zeros(len(t_list[b]))
        daflags_copy[b] = np.zeros(len(t_list[b]))
        for i,p in enumerate(planets):
            t = t_list[b] + time_majorframe[0]
            ra_p, dec_p = get_ephem(t, p)
            ll_p = np.array([ra_p, dec_p])

            d_A = hp.rotator.angdist(ll_A[b].T, ll_p, lonlat=True)
            inds = (d_A <= radii[band][i])
            myflags[b][inds] += 2**(2*i+1)

            d_B = hp.rotator.angdist(ll_B[b].T, ll_p, lonlat=True)
            inds = (d_B <= radii[band][i])
            myflags[b][inds] += 2**(2*i+1+1)


    for b, band in enumerate(bands):
        Nobs = Nobs_array[band]
        for i in range(Nobs):
            daflags_copy[b][i::Nobs] = daflags[:,band]
        ind1 = (daflags_copy[b] % 2 == 1)
        myflags[b] = np.where(ind1, daflags_copy[b], myflags[b])

    if test:
        return daflags_copy, myflags
    else:
        return myflags


def test_flags(band=2):
    prefix = '/mn/stornext/d16/cmbco/ola/wmap/tods/'
    outdir = '/mn/stornext/d16/cmbco/bp/wmap/data/'
    files = glob(prefix + 'uncalibrated/*.fits')
    files.sort()
    #get_all_ephems(files)
    '''
    Given a series of pointings in RA and Dec, first, for a single horn,
    determine the distance from each planet, and if it is within a certain
    radius, you need to bit-encode it.

    You need to take into account that there are horns A and B, and you need the
    distance from the planet for both of the horns, so you need 2 n_planets x
    n_tod x n_bands to account for the distances.

    The final product will be n_tod x n_bands, and each element will be a flag
    2**(i+1) where i is the planet index from 0 to nplanets - 1.
    '''
    Nobs_arr = np.array([12, 12, 15, 15, 20, 20, 30, 30, 30, 30])
    Nobs = Nobs_arr[band]
    prefix = '/mn/stornext/d16/cmbco/ola/wmap/tods/'
    files = glob(prefix + 'uncalibrated/*.fits')
    files.sort()
    file_input = files[87]

    file_input = np.random.choice(files[-365:])
    #file_input = files[-765]
    print(file_input)

    #file_input = files[3197]

    data = fits.open(file_input)
    daflags, myflags = get_flags(data, test=True)



    plt.figure()

    tod_tile = data[2].data['Q113']
    tod = np.zeros(tod_tile.size)
    #daflag_all = np.zeros(tod_tile.size)
    myflag_all = np.zeros(tod_tile.size)
    for i in range(Nobs):
        tod[i::Nobs] = tod_tile[:,i]
        #daflag_all[i::Nobs] = daflags[:,band]
    daflag_all = daflags[band]
    myflag_all[:] = myflags[band]
    ind1 = (daflag_all % 2 == 1)
    myflag_all = np.where(ind1, daflag_all, myflag_all)

    inds = (daflag_all != 0) & ((daflag_all % 2) != 1)
    t = np.arange(len(inds))*1.536/Nobs
    plt.plot(t[tod > 1e3], tod[tod > 1e3], 'k.', alpha=0.5, ms=1)
    plt.plot(t[inds], tod[inds], 'g.', ms=10)
    inds = (myflag_all != 0) & ((myflag_all % 2) != 1)
    plt.plot(t[inds], tod[inds], 'r.', ms=7)
    #plt.xlim([7230, 7250])
    #plt.xlim([27420, 27440])
    #plt.xlim([0,10000])


    plt.figure()
    plt.plot(t, daflag_all, '.', label='DAflags', ms=3)
    plt.yscale('log')

    plt.semilogy(t, myflag_all, '.', label='My flags', ms=2)
    #plt.xlim([7230, 7250])
    #plt.xlim([27420, 27440])
    #plt.xlim([0,10000])
    plt.legend(loc='best')

    ax = plt.gca()

    ax.set_yticks([2**i for i in range(11)])
    ax.set_yticks([], minor=True)

    ax.yaxis.set_major_formatter(lambda x, pos: f"{int(x):011b}")


    plt.show()


    return

def write_file_parallel(comm_tod, i, obsid, obs_ind, daflags, TODs, gain_guesses,
        baseline_guesses,
        band_labels, band, psi_A, psi_B, pix_A, pix_B, alpha, n_per_day,
         npsi, psiBins, nside, fsamp, pos, vel, time, version,
        compress=False, precal=False):

    for ni in range(len(events)):
        if any(time > events[ni][0]) or (any(time < events[ni][1])):
            print(events[ni])

    for t in TODs_all:
      print(len(t), np.log2(len(t)), len(t) - 2**(int(np.log2(len(t)))))

    dt0 = np.diff(time).mean()
    det_list = []
    for j in range(len(band_labels)):
        label = band_labels[j]
        if label[:-2] == band.upper():
            todi = TODs[j]
            gain = gain_guesses[j][i]
            if gain < 0:
              todi = -todi
              gain = -gain


            sigma_0 = np.diff(todi).std()/2**0.5 # Using Eqn 18 of BP06
            scalars = np.array([gain, sigma_0, fknees[j//2], alphas[j//2]])
            if precal:
                baseline = 0
            else:
                baseline = np.median(todi)

            pixA = np.array_split(pix_A[j//4], n_per_day)[i]

            pixB = np.array_split(pix_B[j//4], n_per_day)[i]


            psiA = np.array_split(psi_A[j//4], n_per_day)[i]
            psiA = np.where(psiA < 0,           2*np.pi+psiA,   psiA)
            psiA = np.where(psiA >= 2*np.pi,    psiA - 2*np.pi, psiA)

            psiB = np.array_split(psi_B[j//4], n_per_day)[i]
            psiB = np.where(psiB < 0,           2*np.pi+psiB,   psiB)
            psiB = np.where(psiB >= 2*np.pi,    psiB - 2*np.pi, psiB)

            N = Nobs_array[j//4]
            flags_new = daflags[j//4]
            flags = np.array_split(flags_new, n_per_day)[i]

            if label[-2:] == '13':
                comm_tod.init_file(label.replace('KA', 'Ka')[:-2], obsid, mode='w')
                if compress:
                    huffman = ['huffman', {'dictNum':1}]
                    compArr = [huffman]
                    comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2] + '/flag', 
                             flags, compArr)
                    
                    comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/pixA',
                            pixA, compArr)
                    comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/pixB',
                            pixB, compArr)

                    psiDigitize = ['digitize', {'min':0, 'max':2*np.pi,'nbins':npsi}]
                    compArray = [psiDigitize, huffman]
                    comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/psiA',
                            psiA, compArray)
                    comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/psiB',
                            psiB, compArray)
                else:
                    comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2] + '/flag', 
                             flags)
                    
                    comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/pixA',
                            pixA)
                    comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/pixB',
                            pixB)

                    comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/psiA',
                            psiA)
                    comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/psiB',
                            psiB)
            if precal:
                comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')+ '/tod',
                        data=todi)
            else:
              if compress:
                huffTod = ['huffman', {'dictNum':2}]
                compArr = [huffTod]
                comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')+ '/ztod',
                        todi, compArr)
              else:
                comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')+ '/tod',
                        data=todi)
            # Link to the pointing and flag information
            comm_tod.add_softlink(obsid + '/' + label.replace('KA','Ka') + '/flag',
                        '/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/flag')
            comm_tod.add_softlink(obsid + '/' + label.replace('KA','Ka') + '/pixA',
                        '/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/pixA')
            comm_tod.add_softlink(obsid + '/' + label.replace('KA','Ka') + '/pixB',
                        '/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/pixB')
            comm_tod.add_softlink(obsid + '/' + label.replace('KA','Ka') + '/psiA',
                        '/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/psiA')
            comm_tod.add_softlink(obsid + '/' + label.replace('KA','Ka') + '/psiB',
                        '/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/psiB')



            det_list.append(label.replace('KA','Ka'))

            comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')+'/scalars', 
                data=scalars)
            comm_tod.add_attribute(obsid + '/' + label.replace('KA', 'Ka') + '/scalars',
                'index','gain, sigma0, fknee, alpha')

            comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')+ '/baseline',
                    data=baseline)
            comm_tod.add_attribute(obsid + '/' + label.replace('KA','Ka') + '/baseline', 
                'index', 'baseline')
            # filler 
            comm_tod.add_field(obsid +'/'+label.replace('KA', 'Ka') + '/outP',
                data=np.array([0,0]))




    #satelite position
    comm_tod.add_field(obsid +  '/common/satpos',  np.array_split(pos,n_per_day)[i][0])
    comm_tod.add_attribute(obsid + '/common/satpos','index','X, Y, Z')
    comm_tod.add_attribute(obsid + '/common/satpos','coords','heliocentric')


    comm_tod.add_field(obsid + '/common/vsun', np.array_split(vel,n_per_day)[i][0])
    comm_tod.add_attribute(obsid + '/common/vsun','index', '[x, y, z]')
    comm_tod.add_attribute(obsid + '/common/vsun','coords','galactic')





    dt = dt0/len(TODs[0])
    time_band = np.arange(time.min(), time.min() + dt*len(TOD), dt)

    #time field
    comm_tod.add_field(obsid + '/common/time',[np.array_split(time_band, n_per_day)[i][0],0,0])
    comm_tod.add_attribute(obsid + '/common/time','index','MJD, null, null')

    comm_tod.add_field(obsid + '/common/ntod',
            data=[len(np.array_split(TOD,n_per_day)[i])])

    comm_tod.add_field('/common/det', data=np.string_(', '.join(det_list)))
    comm_tod.add_field('/common/fsamp', fsamp)
    comm_tod.add_field('/common/nside', [nside])

    # fillers

    comm_tod.add_field('/common/polang', np.array([0,0,0,0]))
    comm_tod.add_attribute('/common/polang', 'index', ', '.join(det_list))

    comm_tod.add_field('/common/mbang', np.array([0,0,0,0]))
    comm_tod.add_attribute('/common/mbang', 'index', ', '.join(det_list))

    comm_tod.finalize_chunk(int(obsid), loadBalance=np.array([0,0]))
    comm_tod.finalize_file()

    return

def write_file_serial(comm_tod, i, obsid, obs_ind, daflags, TODs, gain_guesses,
        labels, band, psiA, psiB, pixA, pixB, alpha, n_per_day,
         npsi, psiBins, nside, fsamp, pos, vel, time, version, Nobs,
        compress=False, precal=False):

    dt0 = np.diff(time).mean()
    det_list = []
    for j in range(len(labels)):
        label = labels[j]
        todi = TODs[j]
        gain = gain_guesses[label == all_band_labels][0]
        if gain < 0:
          todi = -todi
          gain = -gain


        sigma_0 = np.diff(todi).std()/2**0.5 # Using Eqn 18 of BP06
        scalars = np.array([gain, sigma_0, fknees[j//2], alphas[j//2]])
        if precal:
            baseline = 0
        else:
            baseline = np.median(todi)

        psiA = np.where(psiA < 0,           2*np.pi+psiA,   psiA)
        psiA = np.where(psiA >= 2*np.pi,    psiA - 2*np.pi, psiA)

        psiB = np.where(psiB < 0,           2*np.pi+psiB,   psiB)
        psiB = np.where(psiB >= 2*np.pi,    psiB - 2*np.pi, psiB)

        N = Nobs
        flags = daflags

        if j == 0:
            comm_tod.init_file(label.replace('KA', 'Ka')[:-2], obsid, mode='w')
            if compress:
                huffman = ['huffman', {'dictNum':1}]
                compArr = [huffman]
                comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2] + '/flag', 
                         flags, compArr)
                
                comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/pixA',
                        pixA, compArr)
                comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/pixB',
                        pixB, compArr)

                psiDigitize = ['digitize', {'min':0, 'max':2*np.pi,'nbins':npsi}]
                compArray = [psiDigitize, huffman]
                comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/psiA',
                        psiA, compArray)
                comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/psiB',
                        psiB, compArray)
            else:
                comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2] + '/flag', 
                         flags)
                
                comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/pixA',
                        pixA)
                comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/pixB',
                        pixB)

                comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/psiA',
                        psiA)
                comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/psiB',
                        psiB)
        if precal:
            comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')+ '/tod',
                    data=todi)
        else:
          if compress:
            huffTod = ['huffman', {'dictNum':2}]
            compArr = [huffTod]
            comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')+ '/ztod',
                    todi, compArr)
          else:
            comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')+ '/tod',
                    data=todi)
        # Link to the pointing and flag information
        comm_tod.add_softlink(obsid + '/' + label.replace('KA','Ka') + '/flag',
                    '/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/flag')
        comm_tod.add_softlink(obsid + '/' + label.replace('KA','Ka') + '/pixA',
                    '/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/pixA')
        comm_tod.add_softlink(obsid + '/' + label.replace('KA','Ka') + '/pixB',
                    '/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/pixB')
        comm_tod.add_softlink(obsid + '/' + label.replace('KA','Ka') + '/psiA',
                    '/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/psiA')
        comm_tod.add_softlink(obsid + '/' + label.replace('KA','Ka') + '/psiB',
                    '/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/psiB')



        det_list.append(label.replace('KA','Ka'))

        comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')+'/scalars', 
            data=scalars)
        comm_tod.add_attribute(obsid + '/' + label.replace('KA', 'Ka') + '/scalars',
            'index','gain, sigma0, fknee, alpha')

        comm_tod.add_field(obsid + '/' + label.replace('KA','Ka')+ '/baseline',
                data=baseline)
        comm_tod.add_attribute(obsid + '/' + label.replace('KA','Ka') + '/baseline', 
            'index', 'baseline')
        # filler 
        comm_tod.add_field(obsid +'/'+label.replace('KA', 'Ka') + '/outP',
            data=np.array([0,0]))




    #satelite position
    comm_tod.add_field(obsid +  '/common/satpos',  pos)
    comm_tod.add_attribute(obsid + '/common/satpos','index','X, Y, Z')
    comm_tod.add_attribute(obsid + '/common/satpos','coords','heliocentric')


    comm_tod.add_field(obsid + '/common/vsun', vel)
    comm_tod.add_attribute(obsid + '/common/vsun','index', '[x, y, z]')
    comm_tod.add_attribute(obsid + '/common/vsun','coords','galactic')





    #time field
    comm_tod.add_field(obsid + '/common/time',[time[0],0,0])
    comm_tod.add_attribute(obsid + '/common/time','index','MJD, null, null')

    comm_tod.add_field(obsid + '/common/ntod',
            data=[len(TODs[0])])

    comm_tod.add_field('/common/det', data=np.string_(', '.join(det_list)))
    comm_tod.add_field('/common/fsamp', fsamp)
    comm_tod.add_field('/common/nside', [nside])

    # fillers

    comm_tod.add_field('/common/polang', np.array([0,0,0,0]))
    comm_tod.add_attribute('/common/polang', 'index', ', '.join(det_list))

    comm_tod.add_field('/common/mbang', np.array([0,0,0,0]))
    comm_tod.add_attribute('/common/mbang', 'index', ', '.join(det_list))

    comm_tod.finalize_chunk(int(obsid), loadBalance=np.array([0,0]))
    comm_tod.finalize_file()

    return

def coord_trans(pos_in, coord_in, coord_out, lonlat=False):
    if coord_in == coord_out:
        return pos_in
    r = hp.rotator.Rotator(coord=[coord_in, coord_out])
    pos_out = r(pos_in.T).T

    if lonlat:
        if pos_out.shape[1] == 2:
            return pos_out
        elif pos_out.shape[1] == 3:
            return hp.vec2dir(pos_out.T, lonlat=lonlat).T
    else:
        return pos_out



def Q2M(Q):
    '''
    PURPOSE:
        Converts quaternions to rotation matrices.
        
    CALLING SEQUENCE:
        M =Q2M(Q)
    
    INPUTS:
        Q - Quaternions.  May be a 2-D array dimensioned 4xN or
            simply a vector dimensioned 4.
    
    OUTPUTS:  
        M - Cube of attitude rotation matrices, 3x3xN (or 3x3
            if only one input quaternion).
    '''
    q1=-Q[0,:]
    q2=-Q[1,:]
    q3=-Q[2,:]
    q4= Q[3,:]
    
    q11=q1**2
    q22=q2**2
    q33=q3**2
    q44=q4**2
    s=q11+q22+q33+q44
    w = (abs(s-1.0) > 1e-5)
    if sum(w) > 0:
        s=np.sqrt(s)
        q1=q1/s
        q2=q2/s
        q3=q3/s
        q4=q4/s
    
    q12=q1*q2
    q13=q1*q3
    q14=q1*q4
    q23=q2*q3
    q24=q2*q4
    q34=q3*q4

    
    M = np.zeros((len(q1), 3,3))
    
    M[:,0,0] = q11 - q22 - q33 + q44
    M[:,0,1] = 2. * ( q12 + q34 )
    M[:,0,2] = 2. * ( q13 - q24 )
    M[:,1,0] = 2. * ( q12 - q34 )
    M[:,1,1] = -q11 + q22 - q33 + q44
    M[:,1,2] = 2. * ( q23 + q14 )
    M[:,2,0] = 2. * ( q13 + q24 )
    M[:,2,1] = 2. * ( q23 - q14 )
    M[:,2,2] = -q11 - q22 + q33 + q44

    M = np.transpose(M, [1,2,0])



    return M

def gamma_from_pol(gal, pol):
    # gal and pol are galactic lonlat vectors
    dir_gal = hp.ang2vec(gal[:,0]%np.pi,gal[:,1]%(2*np.pi), lonlat=False)
    dir_pol = hp.ang2vec(pol[:,0]%np.pi,pol[:,1]%(2*np.pi), lonlat=False)

    dir_Z = np.array([0,0,1])

    sin_theta = np.sqrt(dir_gal[:,0]**2 + dir_gal[:,1]**2)

    dir_west_x = dir_gal[:,1]/sin_theta
    dir_west_y = -dir_gal[:,0]/sin_theta
    dir_west_z = dir_gal[:,1]*0
    dir_west = np.array([dir_west_x, dir_west_y, dir_west_z]).T
    dir_north = (dir_Z - dir_gal[:,2][:,np.newaxis]*dir_gal)/sin_theta[:,np.newaxis]

    # If north Galactic pole is observed
    ind = (sin_theta == 0)
    dir_west[ind] = np.array([1,0,0])
    dir_north[ind] = np.array([0,1,0])


    sin_gamma = dir_pol[:,0]*dir_west[:,0] + dir_pol[:,1]*dir_west[:,1] + dir_pol[:,2]*dir_west[:,2]
    cos_gamma = dir_pol[:,0]*dir_north[:,0] + dir_pol[:,1]*dir_north[:,1] + dir_pol[:,2]*dir_north[:,2]

    cos_2_gamma = 2*cos_gamma**2 - 1
    sin_2_gamma = 2*sin_gamma*cos_gamma

    return sin_gamma, cos_gamma



def quat_to_sky_coords(quat, center=True, lonlat=False, nointerp=False,
    ret_times=False,
        coord_out='G',Nobs_array = np.array([12, 12, 15, 15, 20, 20, 30, 30, 30,
          30])):
    '''
    Quaternion is of form (N_frames, 30, 4), with one redundant frame at the
    beginning and two redundant ones at the end, that match the adjacent frames.
    '''
    nt = len(quat)
    Q = np.zeros( (4, 33, nt))
    q0 = quat[:,0::4]
    q1 = quat[:,1::4]
    q2 = quat[:,2::4]
    q3 = quat[:,3::4]
    q0 = np.array([q0[0,0]] + q0[:,1:-2].flatten().tolist() +
            q0[-1,-2:].tolist())
    q1 = np.array([q1[0,0]] + q1[:,1:-2].flatten().tolist() +
            q1[-1,-2:].tolist())
    q2 = np.array([q2[0,0]] + q2[:,1:-2].flatten().tolist() +
            q2[-1,-2:].tolist())
    q3 = np.array([q3[0,0]] + q3[:,1:-2].flatten().tolist() +
            q3[-1,-2:].tolist())
    Q = np.zeros((4, 30*nt + 3))
    Q[0] = q0
    Q[1] = q1
    Q[2] = q2
    Q[3] = q3
    t0 = np.arange(0, 30*nt + 3)

    dir_A_los = np.array([
                [  0.03993743194318,  0.92448267167832, -0.37912635267982],
                [ -0.03836350153280,  0.92543717887494, -0.37695393578810],
                [ -0.03157188095163,  0.95219265474988, -0.30386241059657],
                [  0.03193385161530,  0.95220162163922, -0.30379647935526],
                [ -0.03317333754910,  0.94156429439011, -0.33519577742792],
                [  0.03337676771235,  0.94149468374332, -0.33537106592570],
                [ -0.00918939185649,  0.93943847522010, -0.34259437583453],
                [ -0.00950701394255,  0.94586439605663, -0.32442281201900],
                [  0.00980040822398,  0.94576779947882, -0.32469558276581],
                [  0.00980808738477,  0.93934799994236, -0.34282522723123]])
    dir_B_los = np.array([
                [  0.03794083653062, -0.92391755783762, -0.38070571212253],
                [ -0.04002167684949, -0.92463440201100, -0.37874726137612],
                [ -0.03340297596219, -0.95176877819247, -0.30499251475222],
                [  0.03014337784306, -0.95192770480751, -0.30483605690947],
                [ -0.03503633693827, -0.94094544143324, -0.33674045100040],
                [  0.03144454385558, -0.94113854675448, -0.33655530968115],
                [ -0.01147317267740, -0.93883247845653, -0.34418300902847],
                [ -0.01159000320270, -0.94535005109668, -0.32585112047876],
                [  0.00768184749607, -0.94540702221088, -0.32580139897397],
                [  0.00751408106677, -0.93889226303920, -0.34412912836731  ]])

    dir_A_pol = np.array([  
                [  0.69487757242271, -0.29835139515692, -0.65431766318192, ],
                [ -0.69545992357813, -0.29560553030986, -0.65494493291187, ],
                [  0.71383872060219, -0.19131247543171, -0.67367189173456, ],
                [ -0.71390969181845, -0.19099503229669, -0.67368675923286, ],
                [ -0.69832280289930, -0.26176968417604, -0.66619959126169, ],
                [  0.69826122350352, -0.26204606404493, -0.66615548040223, ],
                [  0.70944248806767, -0.23532277684296, -0.66431509603747, ],
                [ -0.70476543555624, -0.23649685267332, -0.66886091193973, ],
                [  0.70468980214241, -0.23690904054153, -0.66879472879665, ],
                [ -0.70959923775957, -0.23501806310177, -0.66425554705017]])
    dir_B_pol = np.array([  
                [  0.69546590081501,  0.29798590641998, -0.65385899120425,],
                [ -0.69486414021667,  0.29814186328140, -0.65442742607568, ],
                [  0.71423586688235,  0.19072845484161, -0.67341650037147, ],
                [ -0.71357469183546,  0.19306390125546, -0.67345192048426, ],
                [ -0.69775710213559,  0.26425762446771, -0.66580998365151, ],
                [  0.69876566230957,  0.26145991550208, -0.66585678772745, ],
                [  0.71002796142313,  0.23471528678222, -0.66390438178103, ],
                [ -0.70422900931886,  0.23906270891214, -0.66851366750529, ],
                [  0.70521159225086,  0.23611413753036, -0.66852578425466, ],
                [ -0.70903152581832,  0.23766935833457, -0.66391834701609]])

    M = Q2M(Q)
    M = np.transpose(M, [2,0,1])

    gal_A = []
    pol_A = []
    gal_B = []
    pol_B = []
    t_list = []
    for n, Nobs in enumerate(Nobs_array):
        # for each group from 0--4, the interpolation is valid between 1.5--2.5,
        # which is equivalent to cutting out the first 1.5 time units from the
        # beginning of the total array and the final set of quaternions does not
        # need the last half of the time interval.
        # for i in range(nt):
        #     for j in range(30):
        #          for k in range(Nobs):
        #               offset = 1 + (k+0.5)/Nobs
        #               or
        #               offset = k/Nobs + 1 + 0.5/Nobs
        #               interp(qt, offset, qout)
        if nointerp:
            M2 = M[1:-2]
            Npts = 30*nt
        else:
            Npts = 30*nt*Nobs
            k = np.arange(Npts)
            if center:
                t = 1+(k+0.5)/Nobs
                #t = np.arange(0, 30*nt, 1/Nobs) + 0.5
            else:
                t = 1+k/Nobs
                #t = np.arange(0, 30*nt, 1/Nobs) 

            M2 = np.zeros((len(t), 3, 3))
            for i in range(3):
                for j in range(3):
                    inds = np.isfinite(M[:,i,j])
                    f = interp1d(t0[inds], M[:,i,j][inds], kind='cubic',
                            fill_value='extrapolate')
                    M2[:,i,j] = f(t)
            # T is the sample number here, but I want the actual time elapsed
            t *= 1.536/Nobs
            # This is in seconds, need it in days
            t /= 3600*24
            t_list.append(t)


        dir_A_los_cel = []
        dir_B_los_cel = []
        dir_A_los_cel = np.sum(M2*np.tile(dir_A_los[n, np.newaxis, np.newaxis,:], (Npts,3,1)),axis=2)
        dir_B_los_cel = np.sum(M2*np.tile(dir_B_los[n, np.newaxis, np.newaxis,:], (Npts,3,1)),axis=2)

        dir_A_los_gal = coord_trans(dir_A_los_cel, 'C', coord_out)
        Pll_A = np.array(hp.vec2ang(dir_A_los_gal, lonlat=lonlat))

        dir_B_los_gal = coord_trans(dir_B_los_cel, 'C', coord_out)
        Pll_B = np.array(hp.vec2ang(dir_B_los_gal, lonlat=lonlat))
        gal_A.append(Pll_A.T)
        gal_B.append(Pll_B.T)

        dir_A_pol_cel = np.sum(M2*np.tile(dir_A_pol[n, np.newaxis, np.newaxis,:], (Npts,3,1)),axis=2)
        dir_B_pol_cel = np.sum(M2*np.tile(dir_B_pol[n, np.newaxis, np.newaxis,:], (Npts,3,1)),axis=2)

        dir_A_pol_gal = coord_trans(dir_A_pol_cel, 'C', coord_out)
        Pll_A = np.array(hp.vec2ang(dir_A_pol_gal, lonlat=lonlat))

        dir_B_pol_gal = coord_trans(dir_B_pol_cel, 'C', coord_out)
        Pll_B = np.array(hp.vec2ang(dir_B_pol_gal, lonlat=lonlat))
        pol_A.append(Pll_A.T)
        pol_B.append(Pll_B.T)




    if ret_times:
        return gal_A, gal_B, pol_A, pol_B, t_list
    else:
        return gal_A, gal_B, pol_A, pol_B


def get_psi_band(gal, pol):
    psi = []
    sing, cosg = gamma_from_pol(gal, pol)
    return np.arctan2(sing, cosg)
def get_psi(gal, pol, band_labels):
    psi = []
    for band in range(len(band_labels)):
        sing, cosg = gamma_from_pol(gal[band], pol[band])
        psi.append(np.arctan2(sing, cosg))
    return psi

def ang2pix_multiprocessing(nside, theta, phi):
    return hp.ang2pix(nside, theta, phi)

#@profile(sort_by='cumulative', strip_dirs=True)
def fits_to_h5(comm_tod, file_input, file_ind, compress, plot, version, center,
    precal, simulate):
    t0 = timer()

    # from programs.pars
    #       Gains and baselines given signs and values from the first
    #       attempt at cal_map for Pass 1.  These are median values
    #       from the hourly calibration files.
    #
    gain_guesses0=np.array([ -0.9700,  0.9938,  1.1745, -1.1200, 
                              0.8668, -0.8753, -1.0914,  1.0033, 
                              1.0530, -0.9834,  0.4914, -0.5365, 
                             -0.9882,  1.0173, -0.8135,  0.7896, 
                              0.4896, -0.5380, -0.5840,  0.5840, 
                             -0.4948,  0.4872,  0.4096, -0.3802, 
                              0.3888, -0.4139,  0.3290, -0.3003, 
                             -0.3587,  0.3701,  0.3655, -0.3666, 
                             -0.3255,  0.3517, -0.3291,  0.3225, 
                              0.2841, -0.2918,  0.3796, -0.3591 ])
    baseline = [ 32136.98,  31764.96,  31718.19,  32239.29, 
                 31489.19,  32356.00,  32168.49,  31634.28, 
                 25621.62,  25502.45,  25500.11,  25667.74, 
                 26636.99,  24355.67,  26751.75,  24240.62, 
                 19050.87,  19129.09,  19380.14,  19081.39, 
                 19291.37,  18966.04,  18730.91,  19505.05, 
                 12428.31,  13567.44,  13049.69,  12930.21, 
                 13516.42,  12477.95,  12229.22,  13363.53, 
                 12678.23,  12934.65,  12730.91,  12692.85, 
                 11759.38,  13704.71,  11537.42,  13956.94  ]

    Nobs_array = np.array([12, 12, 15, 15, 20, 20, 30, 30, 30, 30])


    # From Jarosik et al. 2003, Figure 3.
    # Get the alphas and fknees from the mean of each of the Commander runs
    alpha = -1

    nside = 512
    npsi = 2048
    psiBins = np.linspace(0, 2*np.pi, npsi)
    fsamp = 1/1.536 # A single TOD record contains 30 1.536 second major science frames
    n_per_day = 1


    bands = ['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']

    t2jd = 2.45e6

    TODs_all = [
        [],[],[],[],[],[],[],[],[],[],
        [],[],[],[],[],[],[],[],[],[],
        [],[],[],[],[],[],[],[],[],[],
        [],[],[],[],[],[],[],[],[],[]]
    flags_all = [[],[],[],[],[],[],[],[],[],[]]
    psi_A_all = [[],[],[],[],[],[],[],[],[],[]]
    psi_B_all = [[],[],[],[],[],[],[],[],[],[]]
    pix_A_all = [[],[],[],[],[],[],[],[],[],[]]
    pix_B_all = [[],[],[],[],[],[],[],[],[],[]]
    pos_all = []
    vel_all = []
    time_all = []
    for fi in file_input:
        data = fits.open(fi, memmap=False)

        band_labels = data[2].columns.names[1:-6]

        # Returns the gain model estimate at the start of each frame.
        gain_guesses = []
        for i in range(len(band_labels)):
            t_jd, gain = get_gain(data, band_labels[i])
            gain_split = np.array_split(gain, n_per_day)
            gain_guesses.append([gain_guesses0[i] for g in gain_split])
            #gain_guesses.append([g.mean() for g in gain_split])
        gain_guesses = np.array(gain_guesses)
        if precal or simulate:
            gain_guesses = gain_guesses*0 + 1



        TODs = []
        for index, key in enumerate(band_labels):
            tod = data[2].data[key].flatten()
            #tod = np.zeros(TOD.size)
            #for n in range(len(TOD[0])):
            #    tod[n::len(TOD[0])] = TOD[:,n]
            TODs.append(tod)
            if len(TODs_all[index]) == 0:
              TODs_all[index] = tod
            else:
              TODs_all[index] = np.concatenate((TODs_all[index], tod))

        
        
   
        # position (and velocity) in km(/s) in Sun-centered coordinates
        pos = data[1].data['POSITION']*1e3
        vel = data[1].data['VELOCITY']*1e3


        vel = coord_trans(vel, 'C', 'G', lonlat=False)

        # The raw time that is recorded by WMAP is modified reduced julian day.
        # "To convert a modified reduced Julian day to a Julian day, add
        # 2,450,000 to its value."
        # I subtracted 2.4e6 + 0.5 to convert to modified Julian date, the
        # preferred LFI format.
        time = data[2].data['TIME'] + t2jd - 2_400_000.5
        
        if len(pos_all) == 0:
          pos_all = pos
          vel_all = vel
          time_all = time
        else:
          pos_all = np.concatenate((pos_all, pos))
          vel_all = np.concatenate((vel_all, pos))
          time_all = np.concatenate((time_all, time))

        dt0 = np.median(np.diff(time))

        quat = data[1].data['QUATERN']
        gal_A, gal_B, pol_A, pol_B = quat_to_sky_coords(quat, center=center)

        genflags = data[2].data['genflags']*2**11
        daflags = data[2].data['daflags']
        daflags = get_flags(data, center=center)
        for i in range(10):
            Nobs = Nobs_array[i]
            for j in range(Nobs):
                daflags[i][j::Nobs] += genflags
            if len(flags_all[i]) == 0:
                flags_all[i] = daflags[i]
            else:
                flags_all[i] = np.concatenate((flags_all[i], daflags[i]))
        data.close()


        psi_A = get_psi(gal_A, pol_A, band_labels[::4])
        psi_B = get_psi(gal_B, pol_B, band_labels[1::4])


        args_A = [(nside, gal_A[i][:,0], gal_A[i][:,1]) for i in range(len(gal_A))]
        args_B = [(nside, gal_B[i][:,0], gal_B[i][:,1]) for i in range(len(gal_B))]
        pix_A = []
        pix_B = []
        for i in range(len(args_A)):
            pix_A.append(ang2pix_multiprocessing(*args_A[i]))
            pix_B.append(ang2pix_multiprocessing(*args_B[i]))
        for b in range(10):
          if len(psi_A_all[b]) == 0:
            psi_A_all[b] = psi_A[b]
            psi_B_all[b] = psi_B[b]
            pix_A_all[b] = pix_A[b]
            pix_B_all[b] = pix_B[b]
          else:
            psi_A_all[b] = np.concatenate((psi_A_all[b], psi_A[b]))
            psi_B_all[b] = np.concatenate((psi_B_all[b], psi_B[b]))
            pix_A_all[b] = np.concatenate((pix_A_all[b], pix_A[b]))
            pix_B_all[b] = np.concatenate((pix_B_all[b], pix_B[b]))


    if simulate:
        #TODs_all[index] = np.concatenate((TODs_all[index], tod))
        # Overwrite all bands with the simple dipole + CMB fluctuations
        T, Q, U = hp.read_map('cmb_plus_dip.fits', field=(0,1,2))
        for b in range(40):
            pixA = pix_A_all[b//4]
            pixB = pix_B_all[b//4]
            psiA = psi_A_all[b//4]
            psiB = psi_B_all[b//4]
            if (b % 4) < 2:
                TODs_all[b] = np.rint((T[pixA] + Q[pixA]*np.cos(2*psiA) + U[pixA]*np.sin(2*psiA)) \
                                    - (T[pixB] + Q[pixB]*np.cos(2*psiB) + U[pixB]*np.sin(2*psiB)))
            else:
                TODs_all[b] = np.rint((T[pixA] - Q[pixA]*np.cos(2*psiA) - U[pixA]*np.sin(2*psiA)) \
                                    - (T[pixB] - Q[pixB]*np.cos(2*psiB) - U[pixB]*np.sin(2*psiB)))


    for e in events:
      if any(np.searchsorted(e, time) == 1):
        print('Found jump at ', e, file_input)

    #asdfsad
    # Do we just have two obs_inds per file?
    #obs_inds = np.arange(n_per_day) + n_per_day*file_ind + 1
    #obsids = [str(obs_ind).zfill(6) for obs_ind in obs_inds]
    #pos_all = np.array(pos_all)
    #vel_all = np.array(vel_all)
    #time_all = np.array(time_all)
    #for ind, band in enumerate(bands):
    #    args = [(comm_tod, i, obsids[i], obs_inds[i], flags_all, TODs_all, gain_guesses, baseline,
    #            band_labels, band, psi_A_all, psi_B_all, pix_A_all, pix_B_all, 
    #            alpha, n_per_day, npsi, psiBins, nside,
    #            fsamp*Nobs_array[ind], pos_all, vel_all, time_all, version,
    #            compress, precal) for i in range(len(obs_inds))]
    #    for i in range(n_per_day):
    #        write_file_parallel(*args[i])



    return

def split_pow2(comm_tod, band='K1', band_ind=0,
        par=True, plot=False, compress=True, nfiles=sys.maxsize, version=18,
        center=True, precal=False, simulate=False):
    prefix = '/mn/stornext/d16/cmbco/ola/wmap/tods/'
    if (precal):
        files = glob(prefix + 'calibrated/*.fits')
        outdir = '/mn/stornext/d16/cmbco/bp/wmap/data_precal/'
    else:
        files = glob(prefix + 'uncalibrated/*.fits')
        outdir = '/mn/stornext/d16/cmbco/bp/wmap/data_2n_temp/'
    files.sort()

    if (simulate):
        outdir = '/mn/stornext/d16/cmbco/bp/wmap/data_sim/'
        
    gain_guesses0=np.array([ -0.9700,  0.9938,  1.1745, -1.1200, 
                              0.8668, -0.8753, -1.0914,  1.0033, 
                              1.0530, -0.9834,  0.4914, -0.5365, 
                             -0.9882,  1.0173, -0.8135,  0.7896, 
                              0.4896, -0.5380, -0.5840,  0.5840, 
                             -0.4948,  0.4872,  0.4096, -0.3802, 
                              0.3888, -0.4139,  0.3290, -0.3003, 
                             -0.3587,  0.3701,  0.3655, -0.3666, 
                             -0.3255,  0.3517, -0.3291,  0.3225, 
                              0.2841, -0.2918,  0.3796, -0.3591 ])

    print(band)
    i = band_ind
    labels = ['13', '14', '23', '24']
    Nobs_array = np.array([12, 12, 15, 15, 20, 20, 30, 30, 30, 30])

    npsi = 2048
    psiBins = np.linspace(0, 2*np.pi, npsi)
    fsamp = 1/1.536 
    alpha = -1
    obs_ind = 0
    n_per_day = 1
    band_labels = [f'{band}{labels[i]}' for i in range(4)]

    #files = files[:10]
    inds = np.arange(len(files))

    if ('V' in band) or ('W' in band):
      N = 2**22
    else:
      N = 2**21
    nside = 512

    TOD_all = [[],[],[],[]]
    time_all = []
    pos_all = []
    vel_all = []
    TODs = [[],[],[],[]]
    times = []
    flags_all = []
    psi_A_all = []
    psi_B_all = []
    pix_A_all = []
    pix_B_all = []
    plot = False
    for f_ind, f in enumerate(files):
        data = fits.open(f, memmap=False)
        t2jd = data[1].header['time2jd']


        pos = data[1].data['POSITION']*1e3
        vel = data[1].data['VELOCITY']*1e3
        vel = coord_trans(vel, 'C', 'G', lonlat=False)





        genflags = data[2].data['genflags']*2**11
        daflags = data[2].data['daflags']
        daflags = get_flags(data, center=center, bands=np.array([i]))
        daflags = np.array(daflags).flatten()
        Nobs = Nobs_array[i]
        for j in range(Nobs):
            daflags[j::Nobs] += genflags




        quat = data[1].data['QUATERN']
        gal_A, gal_B, pol_A, pol_B = quat_to_sky_coords(quat, center=center,
            Nobs_array=[Nobs])
        psi_A = get_psi_band(gal_A[0], pol_A[0])
        psi_B = get_psi_band(gal_B[0], pol_B[0])

        args_A = (nside, gal_A[0][:,0], gal_A[0][:,1])
        args_B = (nside, gal_B[0][:,0], gal_B[0][:,1])
        pix_A = ang2pix_multiprocessing(*args_A)
        pix_B = ang2pix_multiprocessing(*args_B)




        if len(TODs[0]) == 0:
            for n, lab in enumerate(labels):
               TODs[n] = data[2].data[f'{band}{lab}'].flatten()

            times = data[2].data['TIME'] + t2jd - 2_400_000.5
            dt = np.diff(times)[0]*len(times)/len(TODs[0])
            times = np.arange(times[0], times[0] + len(TODs[0])*dt, dt)

            t_lores = data[1].data['time']
            t_hires = times - t2jd + 2_400_000.5

            pos_arr = np.array([interp1d(t_lores, pos[:,i], fill_value='extrapolate')(t_hires) 
              for i in range(3)])
            vel_arr = np.array([interp1d(t_lores, vel[:,i], fill_value='extrapolate')(t_hires) 
              for i in range(3)])


            pos_arr   = pos_arr.T
            vel_arr   = vel_arr.T
            pix_A_arr = pix_A
            pix_B_arr = pix_B
            psi_A_arr = psi_A
            psi_B_arr = psi_B
            flags_arr = daflags
        else:
            for n, lab in enumerate(labels):
               TODs_ = data[2].data[f'{band}{lab}'].flatten()
               TODs[n] = np.concatenate((TODs[n], TODs_))
            times_ = data[2].data['TIME'] + t2jd - 2_400_000.5

            dt = np.diff(times_)[0]*len(times_)/len(TODs_)
            times_ = np.arange(times_[0], times_[0] + len(TODs_)*dt, dt)
            times = np.concatenate((times, times_))

            t_lores = data[1].data['time']
            t_hires = times_ - t2jd + 2_400_000.5

            pos_ = np.array([interp1d(t_lores, pos[:,i], fill_value='extrapolate')(t_hires) 
              for i in range(3)]).T
            vel_ = np.array([interp1d(t_lores, vel[:,i], fill_value='extrapolate')(t_hires) 
              for i in range(3)]).T


            pos_arr = np.concatenate((pos_arr, pos_))
            vel_arr = np.concatenate((vel_arr, vel_))
            pix_A_arr = np.concatenate((pix_A_arr, pix_A))
            pix_B_arr = np.concatenate((pix_B_arr, pix_B))
            psi_A_arr = np.concatenate((psi_A_arr, psi_A))
            psi_B_arr = np.concatenate((psi_B_arr, psi_B))
            flags_arr = np.concatenate((flags_arr, daflags))


        for e, ef in zip(events, events_flags):
          if any(np.searchsorted(e, times) == 1) and (ef[i] == 1):
            print(band, e, ef)
            i0, i1 = np.searchsorted(times, e)
            TODs_old = [[],[],[],[]]
            TODs_ev  = [[],[],[],[]]
            if (i1 < len(times)) & (i0 > 0):
              times_old = times[:i0]
              for n in range(4):
                  TODs_old[n] = TODs[n][:i0]
              pos_old =     pos_arr[:i0]
              vel_old =     vel_arr[:i0]
              pix_A_old = pix_A_arr[:i0]
              pix_B_old = pix_B_arr[:i0]
              psi_A_old = psi_A_arr[:i0]
              psi_B_old = psi_B_arr[:i0]
              flags_old = flags_arr[:i0]
              if len(times_old) >= N:
                time_all.append(times_old[:N])
                for n in range(4):
                    TOD_all[n].append(TODs_old[n][:N])
                pos_all.append(    pos_old[0])
                vel_all.append(    vel_old[0])
                pix_A_all.append(pix_A_old[:N])
                pix_B_all.append(pix_B_old[:N])
                psi_A_all.append(psi_A_old[:N])
                psi_B_all.append(psi_B_old[:N])
                flags_all.append(flags_old[:N])

                times_old = times_old[N:]
                for n in range(4):
                    TODs_old[n] = TODs_old[n][N:]
                pos_old =     pos_old[N:]
                vel_old =     vel_old[N:]
                pix_A_old = pix_A_old[N:]
                pix_B_old = pix_B_old[N:]
                psi_A_old = psi_A_old[N:]
                psi_B_old = psi_B_old[N:]
                flags_old = flags_old[N:]

              #N0 = min(20, int(np.log2(len(TODs_old[-1]))))
              #while ((len(TODs_old[-1]) > 2**19)
              #    and (not np.log2(len(TODs_old[-1])).is_integer())):
              #  time_all.append(times_old[:2**N0])
              #  for n in range(4):
              #      TOD_all[n].append(TODs_old[n][:2**N0])
              #  pos_all.append(    pos_old[0])
              #  vel_all.append(    vel_old[0])
              #  pix_A_all.append(pix_A_old[:2**N0])
              #  pix_B_all.append(pix_B_old[:2**N0])
              #  psi_A_all.append(psi_A_old[:2**N0])
              #  psi_B_all.append(psi_B_old[:2**N0])
              #  flags_all.append(flags_old[:2**N0])

              #  times_old = times_old[2**N0:]
              #  pos_old =     pos_old[2**N0:]
              #  vel_old =     vel_old[2**N0:]
              #  pix_A_old = pix_A_old[2**N0:]
              #  pix_B_old = pix_B_old[2**N0:]
              #  psi_A_old = psi_A_old[2**N0:]
              #  psi_B_old = psi_B_old[2**N0:]
              #  flags_old = flags_old[2**N0:]
              #  for n in range(4):
              #      TODs_old[n] = TODs_old[n][2**N0:]

              #  N0 -= 1


              times_ev = times[i0:i1]
              for n in range(4):
                  TODs_ev[n] = TODs[n][i0:i1]
              pos_ev =     pos_arr[i0:i1]
              vel_ev =     vel_arr[i0:i1]
              pix_A_ev = pix_A_arr[i0:i1]
              pix_B_ev = pix_B_arr[i0:i1]
              psi_A_ev = psi_A_arr[i0:i1]
              psi_B_ev = psi_B_arr[i0:i1]
              flags_ev = flags_arr[i0:i1]
              # Explicitly flags all data in the event
              flags_ev[(flags_ev % 2 == 0)] += 1

              times = times[i1:]
              for n in range(4):
                  TODs[n] = TODs[n][i1:]
              pos_arr =     pos_arr[i1:]
              vel_arr =     vel_arr[i1:]
              pix_A_arr = pix_A_arr[i1:]
              pix_B_arr = pix_B_arr[i1:]
              psi_A_arr = psi_A_arr[i1:]
              psi_B_arr = psi_B_arr[i1:]
              flags_arr = flags_arr[i1:]


              time_all.append(times_old)
              for n in range(4):
                  TOD_all[n].append(TODs_old[n])
              pos_all.append(pos_old[0])
              vel_all.append(vel_old[0])
              pix_A_all.append(pix_A_old)
              pix_B_all.append(pix_B_old)
              psi_A_all.append(psi_A_old)
              psi_B_all.append(psi_B_old)
              flags_all.append(flags_old)

              time_all.append(times_ev)
              for n in range(4):
                  TOD_all[n].append(TODs_ev[n])
              pos_all.append(pos_ev[0])
              vel_all.append(vel_ev[0])
              pix_A_all.append(pix_A_ev)
              pix_B_all.append(pix_B_ev)
              psi_A_all.append(psi_A_ev)
              psi_B_all.append(psi_B_ev)
              flags_all.append(flags_ev)
            elif (i0 > 0):
              times_old = times[:i0]
              for n in range(4):
                  TODs_old[n] = TODs[n][:i0]
              pos_old =     pos_arr[:i0]
              vel_old =     vel_arr[:i0]
              pix_A_old = pix_A_arr[:i0]
              pix_B_old = pix_B_arr[:i0]
              psi_A_old = psi_A_arr[:i0]
              psi_B_old = psi_B_arr[:i0]
              flags_old = flags_arr[:i0]

              time_all.append(times_old)
              for n in range(4):
                  TOD_all[n].append(TODs_old[n])
              pos_all.append(pos_old[0])
              vel_all.append(vel_old[0])
              pix_A_all.append(pix_A_old)
              pix_B_all.append(pix_B_old)
              psi_A_all.append(psi_A_old)
              psi_B_all.append(psi_B_old)
              flags_all.append(flags_old)

              times = times[i1:]
              for n in range(4):
                  TODs[n] = TODs[n][i1:]
              pos_arr =     pos_arr[i1:]
              vel_arr =     vel_arr[i1:]
              pix_A_arr = pix_A_arr[i1:]
              pix_B_arr = pix_B_arr[i1:]
              psi_A_arr = psi_A_arr[i1:]
              psi_B_arr = psi_B_arr[i1:]
              flags_arr = flags_arr[i1:]


        TODs_old = [[],[],[],[]]
        if len(TODs[0]) >= N:
          times_old = times[:N]
          for n in range(4):
              TODs_old[n] = TODs[n][:N]
          pos_old =     pos_arr[:N]
          vel_old =     vel_arr[:N]
          pix_A_old = pix_A_arr[:N]
          pix_B_old = pix_B_arr[:N]
          psi_A_old = psi_A_arr[:N]
          psi_B_old = psi_B_arr[:N]
          flags_old = flags_arr[:N]

          times = times[N:]
          for n in range(4):
              TODs[n] = TODs[n][N:]
          pos_arr   =   pos_arr[N:]
          vel_arr   =   vel_arr[N:]
          pix_A_arr = pix_A_arr[N:]
          pix_B_arr = pix_B_arr[N:]
          psi_A_arr = psi_A_arr[N:]
          psi_B_arr = psi_B_arr[N:]
          flags_arr = flags_arr[N:]

          time_all.append(times_old)
          for n in range(4):
              TOD_all[n].append(TODs_old[n])
          pos_all.append(pos_old[0])
          vel_all.append(vel_old[0])
          pix_A_all.append(pix_A_old)
          pix_B_all.append(pix_B_old)
          psi_A_all.append(psi_A_old)
          psi_B_all.append(psi_B_old)
          flags_all.append(flags_old)
        #if (times[-1] > 52140) and plot:
        #  fig, axes = plt.subplots(nrows=4, sharex=True)
        #  for n in range(4):
        #      for k in range(len(TOD_all[0])):
        #        axes[n].plot(time_all[k], TOD_all[n][k])
        #  plt.savefig(f'{band}.png')
        #  #plt.show()
        #  for k in range(len(TOD_all[0])):
        #    plt.plot(time_all[k], psi_A_all[k])
        #  #plt.show()
        #  plot = False
        #  plt.close('all')
        #  break


        if (len(time_all) > 0) & ((len(time_all) % 5) == 0):
            print(f_ind, band)
            for ii in range(len(time_all)):
              TOD_in = np.array([TOD_all[n][ii] for n in range(4)])
              obs_ind += 1
              obsid = str(obs_ind).zfill(6)
              write_file_serial(comm_tod, ii, obsid, obs_ind, flags_all[ii],
                  TOD_in, gain_guesses0, band_labels, band, psi_A_all[ii], psi_B_all[ii],
                  pix_A_all[ii], pix_B_all[ii], alpha, n_per_day, npsi, psiBins, nside,
                  fsamp*Nobs, pos_all[ii], vel_all[ii], time_all[ii], version, Nobs, compress=compress,
                  precal=precal)
              #print(f'   {ii} of {band}')

            TOD_all = [[],[],[],[]]
            time_all = []
            pos_all = []
            vel_all = []
            flags_all = []
            psi_A_all = []
            psi_B_all = []
            pix_A_all = []
            pix_B_all = []

    time_all.append(times)
    for n in range(4):
        TOD_all[n].append(TODs[n])
    pos_all.append(pos_arr[0])
    vel_all.append(vel_arr[0])
    pix_A_all.append(pix_A_arr)
    pix_B_all.append(pix_B_arr)
    psi_A_all.append(psi_A_arr)
    psi_B_all.append(psi_B_arr)
    flags_all.append(flags_arr)

    print(f"On {band}'s final iteration")

    for ii in range(len(time_all)):
      TOD_in = np.array([TOD_all[n][ii] for n in range(4)])
      obs_ind += 1
      obsid = str(obs_ind).zfill(6)
      write_file_serial(comm_tod, ii, obsid, obs_ind, flags_all[ii],
          TOD_in, gain_guesses0, band_labels, band, psi_A_all[ii], psi_B_all[ii],
          pix_A_all[ii], pix_B_all[ii], alpha, n_per_day, npsi, psiBins, nside,
          fsamp*Nobs, pos_all[ii], vel_all[ii], time_all[ii], version, Nobs, compress=compress,
          precal=precal)
      #print(f'   {ii} of {band}')
    print(f'finished {band}')



def main(par=True, plot=False, compress=True, nfiles=sys.maxsize, version=18,
        center=True, precal=False, simulate=False):


    prefix = '/mn/stornext/d16/cmbco/ola/wmap/tods/'
    if (precal):
        files = glob(prefix + 'calibrated/*.fits')
        outdir = '/mn/stornext/d16/cmbco/bp/wmap/data_precal/'
    else:
        files = glob(prefix + 'uncalibrated/*.fits')
        outdir = '/mn/stornext/d16/cmbco/bp/wmap/data_2n_temp/'
    files.sort()

    if (simulate):
        outdir = '/mn/stornext/d16/cmbco/bp/wmap/data_sim/'
        

    get_all_ephems(files)



    inds = np.arange(len(files))
    #inds = inds[:8]
    #files = np.array(files)[:8]
    #inds = inds[:len(files)//4]
    #files = np.array(files)[:len(files)//4]
    #inds = inds[len(files)//4:2*len(files)//4]
    #files = np.array(files)[len(files)//4:2*len(files)//4]
    #inds = inds[2*len(files)//4:3*len(files)//4]
    #files = np.array(files)[2*len(files)//4:3*len(files)//4]
    #inds = inds[3*len(files)//4:]
    #files = np.array(files)[3*len(files)//4:]

    #files = np.array_split(np.array(files), 3280//7)
    files = np.array_split(np.array(files), 3280//3)
    #files = np.array_split(np.array(files), 3280)
    inds = np.arange(len(files))



    manager = mp.Manager()
    dicts = {'K1':manager.dict(), 'Ka1':manager.dict(), 'Q1':manager.dict(),
             'Q2':manager.dict(), 'V1':manager.dict(), 'V2':manager.dict(),
             'W1':manager.dict(), 'W2':manager.dict(), 'W3':manager.dict(),
             'W4':manager.dict(),}
    comm_tod = commander_tod.commander_tod(outdir, 'wmap', version, dicts,
        overwrite=True)

    if par:
        nprocs = 128
        #nprocs = 72
        nprocs = 64
        nprocs = 47
        #nprocs = 32
        #nprocs = 2
        os.environ['OMP_NUM_THREADS'] = '1'


        pool = Pool(processes=nprocs)
        print('pool set up')
        x = [pool.apply_async(fits_to_h5, args=[comm_tod, f, i, compress, plot,
          version, center, precal, simulate]) for i, f in zip(inds, files)]
        for i in tqdm(range(len(x)), smoothing=0):
            x[i].get()
            #res.wait()
        pool.close()
        pool.join()

        comm_tod.make_filelists()
    else:
        for i, f in zip(inds, files):
            fits_to_h5(comm_tod, f,i,compress, plot, version, center, precal,
                simulate)

def main2():
    manager = mp.Manager()
    dicts = {'K1':manager.dict(), 'Ka1':manager.dict(), 'Q1':manager.dict(),
             'Q2':manager.dict(), 'V1':manager.dict(), 'V2':manager.dict(),
             'W1':manager.dict(), 'W2':manager.dict(), 'W3':manager.dict(),
             'W4':manager.dict(),}
    outdir = '/mn/stornext/d16/cmbco/bp/wmap/data_2n_temp/'
    version = 50
    comm_tod = commander_tod.commander_tod(outdir, 'wmap', version, dicts,
        overwrite=True)
    bands = ['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']
    pool = Pool(processes=10)
    x = [pool.apply_async(split_pow2, args=[comm_tod, band, ind])
           for ind, band in enumerate(bands)]
    for i in tqdm(range(len(x)), smoothing=0):
        x[i].get()
        #res.wait()
    pool.close()
    pool.join()
    comm_tod.make_filelists()

if __name__ == '__main__':
    #main(version=49, precal=False, compress=True, center=True, par=False)
    #main(version=49, precal=True, compress=True, center=True, simulate=False)
    main2()
