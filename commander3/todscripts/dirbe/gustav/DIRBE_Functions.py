import healpy as hp
import numpy as np
#import matplotlib.pyplot as plt
#import astropy.time as at
#import h5py
#import cosmoglobe
import pixfunc as pf
from astropy.io import fits
#from astropy import units as u
from tqdm import tqdm
#from numba import jit
#from timeit import default_timer as timer
#from datetime import datetime

# Function converting band number string into wavelength string
def wavelength(band):
    bands = np.char.mod('%02d', np.arange(1,11))
    waves = np.array(('1.25', '2.2', '3.5', '4.9', '12', '25', '60', \
            '100', '140', '240'))
    wave = waves[np.where(bands==band)[0]][0]
    return wave

# Function converting an array of QuadCube pixels into HEALpix
# pixels for two choices of nside resolution
def csc2healpix(pixn, subpos, sbsbpos, nside):

    if nside == 128:
        coords = pf.pix2coord(pixn, res=9, coord='E')
        return hp.ang2pix(128, coords.lon.deg, coords.lat.deg, lonlat=True)

    # Currently not functioning as intended
    elif nside == 256:
        sub_pix = (4096*pixn + 16*subpos + sbsbpos).astype('int')
        coords = pf.pix2coord(sub_pix, res=15, coord='E')
        return hp.ang2pix(256, coords.lon.deg, coords.lat.deg, lonlat=True)

# Function calculating the polarisation angle Theta
def polarisation_theta(fits_data):
    skymap_info = fits.open('/mn/stornext/d16/cmbco/ola/dirbe/auxdata/\
    DIRBE_SKYMAP_INFO.FITS')

    # Reading data from fits_data and ordering by time
    time_order = np.argsort(fits_data['Time'])
    pix_num = fits_data['Pixel_No'][time_order]
    attackv = fits_data['AttackV'][time_order]

    # Reading longitude and latitude from skymap_info
    lamda = skymap_info[1].data['ELON-CSC'][pix_num]*np.pi/180
    beta = skymap_info[1].data['ELAT-CSC'][pix_num]*np.pi/180

    # Calculating Theta following recipe by original DIRBE paper
    att_x = attackv[:,0]*np.cos(lamda) + attackv[:,1]*np.sin(lamda)
    att_y = attackv[:,1]*np.cos(lamda) - attackv[:,0]*np.sin(lamda)

    colatitude = np.pi/2 - beta

    x_prime = att_x*np.cos(colatitude) - attackv[:,2]*np.sin(colatitude)

    theta = np.arctan(att_y/x_prime)

    return theta

# Function making a full-mission binned map for a single wavelength
def full_data_binmap(file, band, detector, save=True):
    chunks = np.char.mod('%06d', np.arange(1,286))
    tot_int = np.zeros(hp.nside2npix(nside))
    hits = np.zeros(hp.nside2npix(nside))

    # The flagging bits deemed relevant to the current analysis
    rel_flags = 1 + 2 + 2**2 + 2**6 + 2**7 + 2**10 + 2**11 + 2**12

    # Counters
    unused = 0
    used   = 0

    for chunk in tqdm(chunks):
        # Loading datasets
        tod   = file['/' + chunk + '/' + detector + '/tod'][()]
        pix   = file['/' + chunk + '/' + detector + '/pix'][()]
        flags = file['/' + chunk + '/' + detector + '/flag'][()]

        # Removing flagged data
        tod[np.where(np.bitwise_and(flags, rel_flags) != 0)[0]] = 0

        # Removing sentinel values
        tod[tod <= -16375] = 0

        unused += len(tod[tod==0])
        used   += len(tod[tod!=0])

        if np.any(tod > 0):
            pix = pix[tod != 0]
            tod = tod[tod != 0]

            for i, p in enumerate(pix):
                tot_int[p] += tod[i]
                hits[p] += 1

    #print('%g percent of the data is discarded'\
    # % (unused/(unused + used)*100))

    if (save==True):
        np.save('npy_files/binmap' + band + detector + '.npy', tot_int)
        np.save('npy_files/hitsmap' + band + detector + '.npy', hits)

# Function creating binned maps for all wavelengths using the full_
# data_binmap function
def binmap_all_bands():
    pol_bands = np.array(('01', '02', '03'))
    bands = np.char.mod('%02d', np.arange(1,11))
    detectors = np.array(('A', 'B', 'C'))
    for band in tqdm(bands):
        file = h5py.File('dirbe_hdf5_files/Phot' + band + '_' + \
               str(nside) + '_.hdf5', 'r')
        if (band in pol_bands):
            for detector in detectors:
                print('\n\n', band, detector)
                full_data_binmap(file, band, detector)
        else:
            print('\n\n', band, 'A')
            full_data_binmap(file, band, 'A')

# Function writing all binned maps to FITS file format
def write_maps():
    detectors = np.array(('01A', '01B', '01C', '02A', '02B', '02C', \
    '03A', '03B', '03C', '04A', '05A', '06A', '07A', '08A', \
    '09A', '10A'))
    for detector in detectors:
        binmap = np.load('npy_files/binmap' + detector + '.npy')
        hitsmap = np.load('npy_files/hitsmap' + detector + '.npy')
        hitsmap[hitsmap == 0] = 1
        map = binmap/hitsmap
        hp.fitsfunc.write_map('fits_maps/binmap' + detector + \
        '_gen1.fits', map)
