import os
import sys
import time

import astropy.units as u
import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from funcs import dust
from globals_mapmaker import IFG_SIZE, SPEC_SIZE

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
import globals as g
import my_utils as mu


def sim_dust():

    dust_map_downgraded_mjy = fits.open("test_output/dust_map_downgraded.fits")
    # get map data from the fits file
    dust_map_downgraded_mjy = dust_map_downgraded_mjy[0].data

    nu0_dust = 545 * u.GHz # Planck 2015
    A_d = 163 * u.uK
    T_d = 21 * u.K
    beta_d = 1.53

    frequencies = mu.generate_frequencies("ll", "ss", 257) * u.GHz

    signal = dust(frequencies, A_d, nu0_dust, beta_d, T_d).value
    # check for invalid value encountered in divide
    signal = np.nan_to_num(signal)
    
    plt.plot(frequencies, signal)
    # plt.show()
    plt.savefig("test_output/signal.png")
    plt.close()

    return dust_map_downgraded_mjy, frequencies, signal

def white_noise(ntod, sigma_min=0.001, sigma_max=0.1):
    sigmarand = np.random.uniform(sigma_min, sigma_max, ntod)
    noise = np.random.normal(0, sigmarand[:, np.newaxis], (ntod, IFG_SIZE))

    # save noise in a npz file
    np.savez("test_output/white_noise.npz", noise=sigmarand)

    return noise

def scanning_strategy():
    """
    Get the pixels for each of the data points in the original sky data in order to generate a realistic scanning strategy.
    """
    sky_data = h5py.File(
        g.PREPROCESSED_DATA_PATH_SKY,
        "r",
    )["df_data"]
    
    return np.array(sky_data["pix_gal"][:], dtype=int)

if __name__ == "__main__":
    dust_map_downgraded_mjy, frequencies, sed = sim_dust()
    sed = np.nan_to_num(sed)

    spec = dust_map_downgraded_mjy[:, np.newaxis] * sed[np.newaxis, :]

    print("Calculating and plotting IFGs")

    # time ifg making
    time_start = time.time()

    ifg = np.fft.irfft(spec, axis=1)

    fig, ax = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    # add phase to ifg
    ifg = np.roll(ifg, 360, axis=1)
    # turn ifg into real signal
    ifg = ifg.real

    # introduce scanning strategy
    pix_gal = scanning_strategy()
    ifg_scanning = np.zeros((len(pix_gal), IFG_SIZE))
    for i, pix in enumerate(pix_gal):
        ifg_scanning[i] = ifg[pix]

    # add noise to ifg
    ifg_scanning = ifg_scanning + white_noise(ifg_scanning.shape[0])

    # save ifg products in a npz file
    np.savez("test_output/ifgs.npz", ifg=ifg_scanning, pix=pix_gal)

    time_end = time.time()
    print(f"Time elapsed for IFGs: {(time_end - time_start)/60} minutes")
