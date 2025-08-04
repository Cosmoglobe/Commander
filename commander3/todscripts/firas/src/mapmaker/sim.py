import os
import sys
import time

import astropy.units as u
import h5py
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from funcs import dust
from globals_mapmaker import IFG_SIZE, SPEC_SIZE

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
import globals as g


def sim_dust():

    dust_map_downgraded_mjy = fits.open("test_output/dust_map_downgraded.fits")
    # get map data from the fits file
    dust_map_downgraded_mjy = dust_map_downgraded_mjy[0].data

    nu0_dust = 545 * u.GHz # Planck 2015
    A_d = 163 * u.uK
    T_d = 21 * u.K
    beta_d = 1.53

    dnu = 13.604162
    frequencies = np.linspace(1e-5, dnu * SPEC_SIZE, SPEC_SIZE) * u.GHz
    # print(f"frequencies: {frequencies}")

    signal = dust(frequencies, A_d, nu0_dust, beta_d, T_d).value
    # check for invalid value encountered in divide
    signal = np.nan_to_num(signal)
    
    plt.plot(frequencies, signal)
    # plt.show()
    plt.savefig("test_output/signal.png")
    plt.clf()

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
    dust_map_downgraded_mjy, frequencies, signal = sim_dust()
    # check signal for nans
    print(f"Number of nans in signal: {np.isnan(signal).sum()}")
    print(f"shape of signal: {signal.shape}")
    # print(f"signal: {signal}")
    signal = np.nan_to_num(signal)
    print(f"sizes: {dust_map_downgraded_mjy.shape} and {signal.shape}")
    spec = dust_map_downgraded_mjy[:, np.newaxis] * signal[np.newaxis, :]
    spec512 = np.zeros((spec.shape[0], IFG_SIZE))
    print(f"spec shape: {spec.shape}")
    spec512[:, :SPEC_SIZE] = spec
    for i in range(SPEC_SIZE, IFG_SIZE - 1, 1):
        spec512[:, i] = spec[:, SPEC_SIZE - i]

    # visualise spec_complex
    plt.imshow(np.abs(spec512), aspect="auto")
    plt.colorbar()
    # plt.show()
    plt.savefig("test_output/spec512.png")
    plt.clf()

    # plot some spec_complex
    for i in range(0, len(spec), 100):  
        plt.plot(np.abs(spec512[i]), color="black", alpha=0.5)
        plt.plot(spec512[i], color="red", alpha=0.5)
    # plt.show()
    plt.savefig("test_output/spec512_plot.png")
    plt.clf()

    # plot real and imag parts of spec_complex
    fig, ax = plt.subplots(2, 1)
    for i in range(0, len(spec), 100):
        ax[0].plot(np.real(spec512[i]), color="black", alpha=0.5)
        ax[1].plot(np.imag(spec512[i]), color="red", alpha=0.5)
    # plt.show()
    plt.savefig("test_output/spec512_real_imag_plot.png")
    plt.clf()

    # check ifg for nans
    x_cm = np.linspace(
        0, 1.76, IFG_SIZE
    )  # cm - does this make sense? sky frequencies are cropped but how does that relate to space?

    print("Calculating and plotting IFGs")

    # time ifg making
    time_start = time.time()

    frequencies_icm = (frequencies).to(1 / u.cm, equivalencies=u.spectral()).value

    # dft matrix
    # IW = np.zeros((IFG_SIZE, IFG_SIZE), dtype=complex)
    # IW[0, :] = 1
    # IW[:, 0] = 1
    # omega = np.exp(2j * np.pi / IFG_SIZE)
    # for xi in range(1, IFG_SIZE):
    #     for nui in range(1, IFG_SIZE):
    #         IW[xi, nui] = omega ** ((xi * nui) % IFG_SIZE) # the mod operator just avoids calculating high exponents
    # IW = IW / IFG_SIZE

    # ifg = np.dot(IW, spec512.T).T
    ifg = np.fft.irfft(spec, axis=1)
    print(f"ifg shape: {ifg.shape}")

    # add phase to ifg
    # ifg = ifg * np.exp(1j * np.pi * (x_cm - 1.22))
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

    # check for nans
    print(f"Number of nans in IFGs: {np.isnan(ifg_scanning).sum()}")

    # save ifg products in a npz file
    np.savez("test_output/ifgs.npz", ifg=ifg_scanning, pix=pix_gal)

    time_end = time.time()
    print(f"Time elapsed for IFGs: {(time_end - time_start)/60} minutes")
