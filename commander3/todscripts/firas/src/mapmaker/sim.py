import time

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from funcs import dust
from globals import IFG_SIZE, SPEC_SIZE


def sim_dust():

    dust_map_downgraded_mjy = fits.open("tests/dust_map_downgraded.fits")
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

    return dust_map_downgraded_mjy, frequencies, signal

def white_noise(ntod):
    noise = np.random.normal(0, 0.01, (ntod, IFG_SIZE))

    # save noise in a npz file
    np.savez("tests/white_noise.npz", noise=noise)
    print(f"Noise: {noise}")

    return noise

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
    plt.show()

    # plot some spec_complex
    for i in range(0, len(spec), 100):  
        plt.plot(np.abs(spec512[i]), color="black", alpha=0.5)
        plt.plot(spec512[i], color="red", alpha=0.5)
    plt.show()

    # plot real and imag parts of spec_complex
    fig, ax = plt.subplots(2, 1)
    for i in range(0, len(spec), 100):
        ax[0].plot(np.real(spec512[i]), color="black", alpha=0.5)
        ax[1].plot(np.imag(spec512[i]), color="red", alpha=0.5)
    plt.show()

    # check ifg for nans
    x_cm = np.linspace(
        0, 1.76, IFG_SIZE
    )  # cm - does this make sense? sky frequencies are cropped but how does that relate to space?

    print("Calculating and plotting IFGs")

    # time ifg making
    time_start = time.time()

    frequencies_icm = (frequencies).to(1 / u.cm, equivalencies=u.spectral()).value

    # dft matrix
    IW = np.zeros((IFG_SIZE, IFG_SIZE), dtype=complex)
    IW[0, :] = 1
    IW[:, 0] = 1
    omega = np.exp(2j * np.pi / IFG_SIZE)
    for xi in range(1, IFG_SIZE):
        for nui in range(1, IFG_SIZE):
            IW[xi, nui] = omega ** ((xi * nui) % IFG_SIZE) # the mod operator just avoids calculating high exponents
    IW = IW / IFG_SIZE

    ifg = np.dot(IW, spec512.T).T
    print(f"ifg shape: {ifg.shape}")

    # add phase to ifg
    # ifg = ifg * np.exp(1j * np.pi * (x_cm - 1.22))
    ifg = np.roll(ifg, 360, axis=1)

    ifgnp = np.fft.ifft(spec, n=IFG_SIZE, axis = 1)
    ifgnp = np.roll(ifgnp, 360, axis=1)

    # plot real and imaginary parts of ifg
    fig, ax = plt.subplots(2, 1)
    for i in range(0, len(ifg), 100):
        ax[0].plot(np.real(ifg[i]), color="black", alpha=0.5)
        # ax[1].plot(np.imag(ifg[i]), color="red", alpha=0.5)
        ax[1].plot(np.real(ifgnp[i]), color="blue", alpha=0.5)
    plt.show()

    # compare ifg from IW and ifft for a spec with 257 elements
    # IW = np.zeros((IFG_SIZE, SPEC_SIZE), dtype=complex)
    # IW[0, :] = 1
    # IW[:, 0] = 1
    # omega = np.exp(2j * np.pi / IFG_SIZE)
    # for xi in range(1, IFG_SIZE):
    #     for nui in range(1, SPEC_SIZE):
    #         IW[xi, nui] = omega ** ((xi * nui) % IFG_SIZE) # the mod operator just avoids calculating high exponents
    # IW = IW / IFG_SIZE

    # ifg2 = np.dot(IW, spec.T).T

    # npifg2 = np.fft.ifft(spec, n=IFG_SIZE, axis = 1)

    # fig, ax = plt.subplots(2, 1)
    # for i in range(0, len(ifg), 100):
    #     ax[0].plot(np.real(ifg2[i]), color="black", alpha=0.5)
    #     ax[1].plot(np.real(npifg2[i]), color="blue", alpha=0.5)
    # plt.show()

    # turn ifg into real signal
    ifg = ifg.real

    # add noise to ifg
    # ifg = ifg + white_noise(ifg.shape[0])

    # check for nans
    print(f"Number of nans in IFGs: {np.isnan(ifg).sum()}")

    # save ifg products in a npz file
    np.savez("tests/ifgs.npz", ifg=ifg)

    time_end = time.time()
    print(f"Time elapsed for IFGs: {(time_end - time_start)/60} minutes")
