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

    # frequency mapping
    # nu0 = {"ss": 68.020812, "lf": 23.807283}
    # dnu = {"ss": 13.604162, "lf": 3.4010405}
    # nf = {"lh_ss": 210, "ll_lf": 182, "ll_ss": 43, "rh_ss": 210, "rl_lf": 182, "rl_ss": 43}

    # f_ghz = {}
    # # let's say we use rhss for now
    # for channel in ["rh"]:
    #     for mode in ["ss"]:
    #         if not (mode == "lf" and (channel == "lh" or channel == "rh")):
    #             f_ghz[f"{channel}_{mode}"] = np.linspace(
    #                 nu0[mode],
    #                 nu0[mode] + dnu[mode] * (nf[f"{channel}_{mode}"] - 1),
    #                 nf[f"{channel}_{mode}"],
    #             )

    nu0_dust = 545 * u.GHz # Planck 2015
    # beta_dust = 1.51
    # t_dust = 19.6
    # tau_dust = 9.6e-7
    A_d = 163 * u.uK
    T_d = 21 * u.K
    beta_d = 1.53

    # calculate dust at each frequency
    # for freq in f_ghz["rh_ss"]:
    #     dust_map = dust_map_downgraded_mjy * dust(freq, nu0_dust, beta_dust, t_dust)
    #     hp.mollview(dust_map, title=f"{int(freq):04d} GHz", unit="MJy/sr", min=0, max=200)
    #     plt.savefig(f"tests/dust_maps/{int(freq):04d}.png")
    #     plt.clf()
    #     plt.close()

    dnu = 13.604162
    frequencies = np.linspace(1e-5, dnu * SPEC_SIZE, SPEC_SIZE) * u.GHz
    # print(f"frequencies: {frequencies}")

    # signal = dust(frequencies, tau_dust, nu0_dust, beta_dust, t_dust)
    signal = dust(frequencies, A_d, nu0_dust, beta_d, T_d)
    # check for invalid value encountered in divide
    signal = np.nan_to_num(signal)
    
    plt.plot(frequencies, signal)
    # plt.show()

    return dust_map_downgraded_mjy, frequencies, signal

# pad f_ghz to 512
# frequencies = np.zeros(257)
# frequencies[0:5] = 0
# frequencies[5 : len(f_ghz["rh_ss"]) + 5] = f_ghz["rh_ss"]
# frequencies[len(f_ghz["rh_ss"]) + 5 :] = 0

if __name__ == "__main__":

    dust_map_downgraded_mjy, frequencies, signal = sim_dust()
    # check signal for nans
    print(f"Number of nans in signal: {np.isnan(signal).sum()}")
    print(f"shape of signal: {signal.shape}")
    # print(f"signal: {signal}")
    signal = np.nan_to_num(signal)
    print(f"sizes: {dust_map_downgraded_mjy.shape} and {signal.shape}")
    spec = dust_map_downgraded_mjy[:, np.newaxis] * signal[np.newaxis, :]

    # visualise spec_complex
    plt.imshow(np.abs(spec.value), aspect="auto")
    plt.colorbar()
    plt.show()

    # plot some spec_complex
    for i in range(0, len(spec), 100):  
        plt.plot(np.abs(spec[i]), color="black", alpha=0.5)
        plt.plot(spec[i], color="red", alpha=0.5)
    plt.show()

    # plot real and imag parts of spec_complex
    fig, ax = plt.subplots(2, 1)
    for i in range(0, len(spec), 100):
        ax[0].plot(np.real(spec[i]), color="black", alpha=0.5)
        ax[1].plot(np.imag(spec[i]), color="red", alpha=0.5)
    plt.show()

    # ifg = np.zeros((len(dust_map_downgraded_mjy), IFG_SIZE))
    # check ifg for nans
    x_cm = np.linspace(
        0, 1.76, IFG_SIZE
    )  # cm - does this make sense? sky frequencies are cropped but how does that relate to space?
    # nu_icm = (f_ghz["rh_ss"] * u.GHz).to(1 / u.cm, equivalencies=u.spectral()).value  # cm-1
    # get the IFGs from the spectra

    print("Calculating and plotting IFGs")

    # time ifg making
    time_start = time.time()

    frequencies_icm = (frequencies).to(1 / u.cm, equivalencies=u.spectral()).value
    # x = np.linspace(0, 1.76, 512)
    # for i in range(len(dust_map_downgraded_mjy)):
    #     for j in range(IFG_SIZE):
    #         ifg[i, j] = np.sum(
    #             # dust_map_downgraded_mjy[i] *
    #             signal
    #             * np.cos(
    #                 2
    #                 * np.pi
    #                 * frequencies_icm
    #                 * (x[j] - 1.22)
    #                 # * x[j]
    #             )
    #         )

    #     # check ifg[i] for nans
    #     # print(f"Number of nans in IFGs 2: {np.isnan(ifg[i]).sum()}")
    #     # ifg[i] = np.fft.irfft(dust_map_downgraded_mjy[i] * signal)
    #     ifg[i] = dust_map_downgraded_mjy[i] * ifg[i]

    # dft matrix
    IW = np.zeros((IFG_SIZE, SPEC_SIZE), dtype=complex)
    IW[0, :] = 1
    IW[:, 0] = 1
    omega = np.exp(2j * np.pi / IFG_SIZE)
    for xi in range(1, IFG_SIZE):
        for nui in range(1, SPEC_SIZE):
            IW[xi, nui] = omega ** ((xi * nui) % IFG_SIZE) # the mod operator just avoids calculating high exponents
    IW = IW / np.sqrt(IFG_SIZE)

    # IF = np.zeros((IFG_SIZE, SPEC_SIZE), dtype=complex)
    
    # for i in range(SPEC_SIZE):
    #     x = np.zeros(SPEC_SIZE)
    #     x[i] = 1
    #     y = np.fft.irfft(x, n = IFG_SIZE)
    #     print(f"y: {y.shape}")
    #     IF[:, i] = y

    # ifg = np.dot(F, spec_complex.T).T
    ifg = np.dot(IW, spec.T).T
    print(f"ifg shape: {ifg.shape}")

    # add phase to ifg
    # ifg = ifg * np.exp(1j * np.pi * (x_cm - 1.22))
    # ifg = np.roll(ifg, 360, axis=1)

    # check for nans
    print(f"Number of nans in IFGs: {np.isnan(ifg).sum()}")

    # save ifg products in a npz file
    np.savez("tests/ifgs.npz", ifg=ifg.value)

    time_end = time.time()
    print(f"Time elapsed for IFGs: {(time_end - time_start)/60} minutes")
