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

import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
import h5py
import healpy as hp

from pytest import approx
from glob import glob
import multiprocessing as mp
from multiprocessing import Pool
from joblib import Parallel, delayed

import huffman


from scipy.interpolate import interp1d
from joblib import Parallel, delayed
import os, sys

from tqdm import tqdm

Nobs_array = np.array([12, 12, 15, 15, 20, 20, 30, 30, 30, 30])
fknees = np.array([6.13, 5.37, 1.66, 1.29,
                   3.21, 3.13, 5.76, 8.62,
                   2.56, 4.49, 2.43, 8.35,
                  16.17,15.05, 9.02, 7.47,
                   1.84, 2.39, 46.5, 26.0]) # mHz
fknees *= 1e-3
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
# version 42 uses the radius-based planet flag, changes tipped TOD to ztod

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

from get_gain_model import get_gain


def get_ephem(time, planet):
    # Obtained these from https://ssd.jpl.nasa.gov/horizons.cgi
    t_p, ra_p, dec_p = np.loadtxt(f'{planet.lower()}_ephem.txt').T

    f_ra = interp1d(t_p, ra_p, fill_value='extrapolate')
    f_dec = interp1d(t_p, dec_p, fill_value='extrapolate')

    return f_ra(time), f_dec(time)


def get_flags(data, test=False):

    t2jd = 2.45e6

    quat = data[1].data['QUATERN']

    ll_A, ll_B, p_A, p_B = quat_to_sky_coords(quat, lonlat=True, nointerp=True,
            coord_out='C')

    time = data[2].data['TIME'] + t2jd

    daflags = data[2].data['daflags']

    bands = np.arange(10)
    planets = ['Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
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

    dists = np.zeros((2*len(planets), daflags.shape[0], daflags.shape[1]))
    for i,p in enumerate(planets):
        ra_p, dec_p = get_ephem(time, p)
        ll_p = np.array([ra_p, dec_p])
        for band in bands:
            d_A = hp.rotator.angdist(ll_A[band].T, ll_p, lonlat=True)
            d_B = hp.rotator.angdist(ll_B[band].T, ll_p, lonlat=True)
            dists[2*i  ,:,band] = d_A*180/np.pi
            dists[2*i+1,:,band] = d_B*180/np.pi

    myflags = np.zeros(daflags.shape, dtype=int)
    for i in range(2*len(planets)):
        for band in bands:
            inds = (dists[i,:,band] < radii[band][i//2]*2)
            #inds = (dists[i,:,band] < radii[band][i//2])
            #inds = (dists[i,:,band] < 7)
            myflags[inds,band]= 2**(i+1)


    ind1 = (daflags % 2 == 1)
    myflags = np.where(ind1, daflags, myflags)

    if test:
        return daflags, myflags
    else:
        return myflags


def test_flags(band=0):
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
    prefix = '/mn/stornext/d16/cmbco/ola/wmap/tods/'
    files = glob(prefix + 'uncalibrated/*.fits')
    #file_input = np.random.choice(files)
    files.sort()
    file_input = files[87]
    #file_input = files[3197]
    print(file_input)

    data = fits.open(file_input)
    daflags, myflags = get_flags(data, test=True)

    band = 0
    plt.figure()
    plt.plot(daflags[:,band], '.', label='DAflags', ms=10)
    plt.yscale('log')

    plt.semilogy(myflags[:,band], '.', label='My flags', ms=7)
    plt.legend(loc='best')

    plt.show()

    return

def write_file_parallel(file_ind, i, obsid, obs_ind, daflags, TODs, gain_guesses,
        baseline_guesses,
        band_labels, band, psi_A, psi_B, pix_A, pix_B, alpha, n_per_day,
        ntodsigma, npsi, psiBins, nside, fsamp, pos, vel, time, version, compress=False):
    prefix = '/mn/stornext/d16/cmbco/bp/wmap/'
    file_out =  prefix + f'data/wmap_{band}_{str(file_ind+1).zfill(6)}_v{version}.h5'
    dt0 = np.diff(time).mean()
    det_list = []
    # make huffman code tables
    # Pixel, Psi, Flag
    pixArray = [[], [], []]
    todArray = []
    for j in range(len(band_labels)):
        label = band_labels[j]
        if label[:-2] == band.upper():
            TOD = TODs[j]
            gain = gain_guesses[j][i]
            if gain < 0:
              TOD = -TOD
              gain = -gain
                    
            todi = np.array_split(TOD, n_per_day)[i].astype('int')
            sigma_0 = np.diff(todi).std()/2**0.5 # Using Eqn 18 of BP06
            scalars = np.array([gain, sigma_0, fknees[j//2], alpha])

            delta = np.diff(todi)
            delta = np.insert(delta, 0, todi[0])
            todArray.append(delta)


            pix = np.array_split(pix_A[j//4], n_per_day)[i]
            delta = np.diff(pix)
            delta = np.insert(delta, 0, pix[0])
            pixArray[0].append(delta)

            pix = np.array_split(pix_B[j//4], n_per_day)[i]
            delta = np.diff(pix)
            delta = np.insert(delta, 0, pix[0])
            pixArray[0].append(delta)


            psi = np.array_split(psi_A[j//4], n_per_day)[i]
            psi = np.where(psi < 0,         2*np.pi + psi,  psi)
            psi = np.where(psi >= 2*np.pi,  psi - 2*np.pi,  psi)
            psiIndexes = np.digitize(psi, psiBins)
            delta = np.diff(psiIndexes)
            delta = np.insert(delta, 0, psiIndexes[0])
            pixArray[1].append(delta)

            psi = np.array_split(psi_B[j//4], n_per_day)[i]
            psi = np.where(psi < 0,         2*np.pi + psi,  psi)
            psi = np.where(psi >= 2*np.pi,  psi - 2*np.pi,  psi)
            psiIndexes = np.digitize(psi, psiBins)
            delta = np.diff(psiIndexes)
            delta = np.insert(delta, 0, psiIndexes[0])
            pixArray[1].append(delta)


            N = Nobs_array[j//4]
            flags = daflags[:,j//4]
            flags_new = np.zeros(len(flags)*N)
            for n in range(N):
              if len(flags_new[n::N]) == len(flags):
                flags_new[n::N] = flags
              else:
                flags_new[n::N] = flags[:len(flags_new[n::N])-len(flags)]
            flags = np.array_split(flags_new, n_per_day)[i]
            delta = np.diff(flags)
            delta = np.insert(delta, 0, flags[0])
            # just necessary to make the array have the correct shape. Redundant
            # info.
            pixArray[2].append(delta)
            pixArray[2].append(delta)


    if compress:
      h = huffman.Huffman("", nside)
      h.GenerateCode(pixArray)
      h_Tod = huffman.Huffman("", nside)
      h_Tod.GenerateCode(todArray)

      huffarray = np.append(np.append(np.array(h.node_max), h.left_nodes), h.right_nodes)
      huffarray_Tod = np.append(np.append(np.array(h_Tod.node_max), h_Tod.left_nodes), h_Tod.right_nodes)


    f = h5py.File(file_out, 'a')
    for j in range(len(band_labels)):
        label = band_labels[j]
        if label[:-2] == band.upper():
            TOD = TODs[j]
            gain = gain_guesses[j][i]
            if gain < 0:
              TOD = -TOD
              gain = -gain


            todi = np.array_split(TOD, n_per_day)[i].astype('int')
            deltatod = np.diff(todi)
            deltatod = np.insert(deltatod, 0, todi[0])

            sigma_0 = np.diff(todi).std()/2**0.5 # Using Eqn 18 of BP06
            scalars = np.array([gain, sigma_0, fknees[j//2], alpha])
            if (version == 'cal'):
                baseline = 0
            else:
                baseline = np.median(todi)

            pixA = np.array_split(pix_A[j//4], n_per_day)[i]
            deltapixA = np.diff(pixA)
            deltapixA = np.insert(deltapixA, 0, pixA[0])


            pixB = np.array_split(pix_B[j//4], n_per_day)[i]
            deltapixB = np.diff(pixB)
            deltapixB = np.insert(deltapixB, 0, pixB[0])


            psiA = np.array_split(psi_A[j//4], n_per_day)[i]
            psiA = np.where(psiA < 0,           2*np.pi+psiA,   psiA)
            psiA = np.where(psiA >= 2*np.pi,    psiA - 2*np.pi, psiA)
            psiIndexesA = np.digitize(psiA, psiBins)
            deltapsiA = np.diff(psiIndexesA)
            deltapsiA = np.insert(deltapsiA, 0, psiIndexesA[0])

            psiB = np.array_split(psi_B[j//4], n_per_day)[i]
            psiB = np.where(psiB < 0,           2*np.pi+psiB,   psiB)
            psiB = np.where(psiB >= 2*np.pi,    psiB - 2*np.pi, psiB)
            psiIndexesB = np.digitize(psiB, psiBins)
            deltapsiB = np.diff(psiIndexesB)
            deltapsiB = np.insert(deltapsiB, 0, psiIndexesB[0])

            N = Nobs_array[j//4]
            flags = daflags[:,j//4]
            flags_new = np.zeros(len(flags)*N)
            for n in range(N):
              if len(flags_new[n::N]) == len(flags):
                flags_new[n::N] = flags
              else:
                flags_new[n::N] = flags[:len(flags_new[n::N])-len(flags)]
            flags = np.array_split(flags_new, n_per_day)[i]


            deltaflag = np.diff(flags)
            deltaflag = np.insert(deltaflag, 0, flags[0])


            if version != 'cal':
              if compress:
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/ztod',
                        data=np.void(bytes(h_Tod.byteCode(deltatod))))
              else:
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/tod',
                        data=todi)
            else:
              if compress:
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/ztod',
                        data=np.void(bytes(h_Tod.byteCode(deltatod))))
              else:
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/tod',
                        data=todi)
            if label[-2:] == '13':
                if compress:
                    f.create_dataset(obsid + '/' + label.replace('KA','Ka')[:-2] + '/flag',
                            data=np.void(bytes(h.byteCode(deltaflag))))
                    
                    f.create_dataset(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/pixA',
                            data=np.void(bytes(h.byteCode(deltapixA))))
                    f.create_dataset(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/pixB',
                            data=np.void(bytes(h.byteCode(deltapixB))))

                    f.create_dataset(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/psiA',
                            data=np.void(bytes(h.byteCode(deltapsiA))))
                    f.create_dataset(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/psiB',
                            data=np.void(bytes(h.byteCode(deltapsiB))))
                else:
                    f.create_dataset(obsid + '/' + label.replace('KA','Ka')[:-2] + '/flag',
                            data=flags)
                    
                    f.create_dataset(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/pixA',
                            data=pixA)
                    f.create_dataset(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/pixB',
                            data=pixB)

                    f.create_dataset(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/psiA',
                            data=psiA)
                    f.create_dataset(obsid + '/' + label.replace('KA','Ka')[:-2]+ '/psiB',
                            data=psiB)
            # Link to the pointing and flag information
            f[obsid + '/' + label.replace('KA','Ka') + '/flag'] = h5py.SoftLink('/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/flag')
            
            f[obsid + '/' + label.replace('KA','Ka') + '/pixA'] = h5py.SoftLink('/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/pixA')
            f[obsid + '/' + label.replace('KA','Ka') + '/pixB'] = h5py.SoftLink('/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/pixB')
            f[obsid + '/' + label.replace('KA','Ka') + '/psiA'] = h5py.SoftLink('/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/psiA')
            f[obsid + '/' + label.replace('KA','Ka') + '/psiB'] = h5py.SoftLink('/' + obsid + '/' + label.replace('KA','Ka')[:-2] + '/psiB')



            det_list.append(label.replace('KA','Ka'))

            f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/scalars',
                    data=scalars)
            f[obsid + '/' + label.replace('KA','Ka') + '/scalars'].attrs['legend'] = 'gain, sigma0, fknee, alpha'
            f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/baseline',
                    data=np.array([baseline]))
            f[obsid + '/' + label.replace('KA','Ka') + '/baseline'].attrs['legend'] = 'baseline'
            # filler 
            f.create_dataset(obsid + '/' + label.replace('KA','Ka') + '/outP',
                    data=np.array([0,0]))



    if compress:
        f.create_dataset(obsid + '/common/todtree', data=huffarray_Tod)
        f.create_dataset(obsid + '/common/todsymb', data=h_Tod.symbols)

        f.create_dataset(obsid + '/common/hufftree', data=huffarray)
        f.create_dataset(obsid + '/common/huffsymb', data=h.symbols)
    
    f.create_dataset(obsid + '/common/satpos',
            data=np.array_split(pos,n_per_day)[i][0])
    f[obsid + '/common/satpos'].attrs['info'] = '[x, y, z]'
    f[obsid + '/common/satpos'].attrs['coords'] = 'galactic'

    f.create_dataset(obsid + '/common/vsun',
            data=np.array_split(vel,n_per_day)[i][0])
    f[obsid + '/common/vsun'].attrs['info'] = '[x, y, z]'
    f[obsid + '/common/vsun'].attrs['coords'] = 'galactic'

    dt = dt0/len(TODs[0])
    time_band = np.arange(time.min(), time.min() + dt*len(TOD), dt)
    f.create_dataset(obsid + '/common/time',
            data=[np.array_split(time_band, n_per_day)[i][0],0,0])
    f[obsid + '/common/time'].attrs['type'] = 'MJD, null, null'

    f.create_dataset(obsid + '/common/ntod',
            data=len(np.array_split(TOD,n_per_day)[i]))

    if "/common/fsamp" not in f:
        f.create_dataset('/common/fsamp', data=fsamp)
        f.create_dataset('/common/nside', data=nside)
        f.create_dataset('/common/npsi', data=npsi)
        f.create_dataset('/common/det', data=np.string_(', '.join(det_list)))
        f.create_dataset('/common/datatype', data='WMAP')
        # fillers
        f.create_dataset('/common/mbang', data=np.array([0,0,0,0]))
        f.create_dataset('/common/ntodsigma', data=100)
        f.create_dataset('/common/polang', data=np.array([0,0,0,0]))
    with open(prefix + f'data/filelist_{band}_v{version}.txt', 'a') as file_list: 
        file_list.write(f'{str(obs_ind).zfill(6)}\t"{file_out}"\t1\t0\t0\n')
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

def q_interp(q_arr, t):
    '''
    Copied from interpolate_quaternions.pro

    This is an implementation of Lagrange polynomials, equation 3.2.1 of
    numerical recipes 3rd edition.


    ;   input_q  - Set of 4 evenly-spaced quaternions (in a 4x4 array).
    ;          See the COMMENTS section for how this array should
    ;          be arranged.
    ;   offset   - Dimensionless time offset relative to the first quaternion.
    ;   This routine expects a unifomly sampled set of quaternions Q1,Q2,Q3,Q4.
    ;   It interpolate a quaternion for any time between Q1 and Q4, inclusive.
    ;   The output is calculated at a time T_Out, expressed in terms of the
    ;   sampling of the input quaternions:
    ;
    ;                   T_Out - T(Q1)
    ;       Offset = -----------------
    ;                   T(Q2) - T(Q1)
    ;
    ;   where T(Q1) is the time at quaternion Q1, and so forth.  That is,
    ;   the time for the output quaternion (variable OFFSET) should be
    ;   a number in the range -1.000 to 4.000 inclusive.  Input values outside
    ;   that range result in an error.  Input values outside 0.0 to 3.0 result
    ;   in extrapolation instead of interpolation.
    ;
    ;       In other words, Offset is essentially a floating point subscript,
    ;       similar to the those used by the IDL intrinsic routine INTERPOLATE.
    ;
    ;   For optimal results, OFFSET should be in the range [1.0, 2.0] -- that
    ;   is, the input quaternions Q1...Q4 should be arranged such that 2 come
    ;   before the desired output and 2 come after.

    '''
    xp0 = t-1
    xn0 = -xp0
    xp1 = xp0 + 1
    xn1 = xp0 - 1
    xn2 = xp0 - 2
    w = np.array([xn0*xn1*xn2/6, xp1*xn1*xn2/2, xp1*xn0*xn2/2, xp1*xp0*xn1/6])
    Qi = q_arr.dot(w)
    Qi = Qi/np.sum(Qi**2, axis=0)**0.5
    return Qi


def quat_to_sky_coords(quat, center=True, lonlat=False, nointerp=False,
        coord_out='G'):
    Nobs_array = np.array([12, 12, 15, 15, 20, 20, 30, 30, 30, 30])
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
                [ 0.69487757242271, -0.29835139515692, -0.65431766318192, ],
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
                [ 0.69546590081501,  0.29798590641998, -0.65385899120425,],
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




    return gal_A, gal_B, pol_A, pol_B


def get_psi(gal, pol, band_labels):
    psi = []
    for band in range(len(band_labels)):
        sing, cosg = gamma_from_pol(gal[band], pol[band])
        psi.append(np.arctan2(sing, cosg))
    return psi

def ang2pix_multiprocessing(nside, theta, phi):
    return hp.ang2pix(nside, theta, phi)

#@profile(sort_by='cumulative', strip_dirs=True)
def fits_to_h5(file_input, file_ind, compress, plot, version, center):
    prefix = '/mn/stornext/d16/cmbco/bp/wmap/'
    file_out = prefix + f'data/wmap_K1_{str(file_ind+1).zfill(6)}_v{version}.h5'
    if (os.path.exists(file_out) and file_ind != 1):
        return
    t0 = timer()

    # from programs.pars
    #       Gains and baselines given signs and values from the first
    #       attempt at cal_map for Pass 1.  These are median values
    #       from the hourly calibration files.
    #
    gain_guesses =np.array([ -0.9700,  0.9938,  1.1745, -1.1200, 
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


    alpha = -1

    nside = 512
    ntodsigma = 100
    npsi = 2048
    psiBins = np.linspace(0, 2*np.pi, npsi)
    fsamp = 1/1.536 # A single TOD record contains 30 1.536 second major science frames
    chunk_size = 1875
    nsamp = chunk_size*fsamp
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
            gain_guesses.append([g.mean() for g in gain_split])
        gain_guesses = np.array(gain_guesses)
        if (version == 'cal'):
            gain_guesses = gain_guesses*0 + 1



        TODs = []
        for index, key in enumerate(band_labels):
            TOD = data[2].data[key]
            tod = np.zeros(TOD.size)
            for n in range(len(TOD[0])):
                tod[n::len(TOD[0])] = TOD[:,n]
            TODs.append(tod)
            if len(TODs_all[index]) == 0:
              TODs_all[index] = tod
            else:
              TODs_all[index] = np.concatenate((TODs_all[index], tod))

        
        
   
        # position (and velocity) in km(/s) in Sun-centered coordinates
        pos = data[1].data['POSITION']*1e3
        vel = data[1].data['VELOCITY']*1e3


        pos = coord_trans(pos, 'C', 'G', lonlat=False)
        vel = coord_trans(vel, 'C', 'G', lonlat=False)

        time_aihk = data[1].data['TIME'] + t2jd
        time = data[2].data['TIME'] + t2jd - 24000000.5
        
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

        # If genflags == 1, there is an issue with the spacecraft attitude. This
        # appears to be the quaternion problem.
        # v14: add genflags somehow.
        genflags = data[2].data['genflags']*2**11
        daflags = data[2].data['daflags']
        daflags = get_flags(data)
        for i in range(10):
            daflags[:,i] += genflags
            if len(flags_all[i]) == 0:
              flags_all[i] = daflags[:,i]
            else:
              flags_all[i] = np.concatenate((flags_all[i], daflags[:,i]))
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


    obs_inds = np.arange(n_per_day) + n_per_day*file_ind + 1
    obsids = [str(obs_ind).zfill(6) for obs_ind in obs_inds]
    pos_all = np.array(pos_all)
    vel_all = np.array(vel_all)
    time_all = np.array(time_all)
    flags_all = np.array(flags_all).T
    #print('TODs[1].shape')
    #print(TODs[1].shape)
    #print('TODs_all[1].shape')
    #print(TODs_all[1].shape)
    #print('psi_A[0].shape')
    #print(psi_A[0].shape)
    #print('len(psi_A[0])')
    #print(len(psi_A[0]))
    #print('len(psi_A_all[0])')
    #print(len(psi_A_all[0]))
    #print('pos.shape')
    #print(pos.shape)
    #print('pos_all.shape')
    #print(pos_all.shape)
    #print('time_all.shape')
    #print(time_all.shape)
    #print('time.shape')
    #print(time.shape)
    #print('daflags.shape')
    #print(daflags.shape)
    #print('flags_all.shape')
    #print(flags_all.shape)
    for ind, band in enumerate(bands):
        #args = [(file_ind, i, obsids[i], obs_inds[i], daflags, TODs, gain_guesses, baseline,
        #            band_labels, band, psi_A, psi_B, pix_A, pix_B, 
        #            alpha, n_per_day, ntodsigma, npsi, psiBins, nside,
        #            fsamp*Nobs_array[ind], pos, vel, time, version, compress) for i in range(len(obs_inds))]
        args = [(file_ind, i, obsids[i], obs_inds[i], flags_all, TODs_all, gain_guesses, baseline,
                band_labels, band, psi_A_all, psi_B_all, pix_A_all, pix_B_all, 
                alpha, n_per_day, ntodsigma, npsi, psiBins, nside,
                fsamp*Nobs_array[ind], pos_all, vel_all, time_all, version, compress) for i in range(len(obs_inds))]
        for i in range(n_per_day):
            write_file_parallel(*args[i])



    return

def main(par=True, plot=False, compress=False, nfiles=sys.maxsize, version=18,
        center=False):
    '''
    Make 1 hdf5 file for every 10 fits files
    # Actually, just doing 1 hdf5 file for every fits file. Too much clashing is
    # happening right now.
    '''

    prefix = '/mn/stornext/d16/cmbco/ola/wmap/tods/'
    if (version == 'cal'):
        files = glob(prefix + 'calibrated/*.fits')
    else:
        files = glob(prefix + 'uncalibrated/*.fits')
    files.sort()
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

    files = np.array_split(np.array(files), 3280//7)
    inds = np.arange(len(files))



    if par:
        nprocs = 64
        nprocs = 128
        os.environ['OMP_NUM_THREADS'] = '1'

        pool = Pool(processes=nprocs)
        print('pool set up')
        x = [pool.apply_async(fits_to_h5, args=[f, i, compress, plot, version, center]) for i, f in zip(inds, files)]
        for i in tqdm(range(len(x)), smoothing=0):
            x[i].get()
            #res.wait()
        pool.close()
        pool.join()
    else:
        for i, f in zip(inds, files):
            fits_to_h5(f,i,compress, plot, version, center)


if __name__ == '__main__':
    #main(par=True, plot=False, compress=True, version=31, center=True)
    #main(par=True, plot=False, compress=True, version='cal', center=True)
    #main(par=True, plot=False, compress=True, version=34, center=True)
    #main(par=True, plot=False, compress=True, version=35, center=True)
    #main(par=True, plot=False, compress=True, version=36, center=True)
    #main(par=True, plot=False, compress=True, version=37, center=True)
    #main(par=True, plot=False, compress=True, version=38, center=True)
    #main(par=True, plot=False, compress=True, version=39, center=True)
    #main(par=True, plot=False, compress=True, version=40, center=True)
    #main(par=True, plot=False, compress=True, version=41, center=True)
    main(par=True, plot=False, compress=True, version=42, center=True)
    #test_flags()
