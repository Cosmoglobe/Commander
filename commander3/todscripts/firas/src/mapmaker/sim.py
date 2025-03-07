import time

import astropy.units as u
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

    nu0_dust = 545  # Planck
    beta_dust = 1.55
    t_dust = 23

    # calculate dust at each frequency
    # for freq in f_ghz["rh_ss"]:
    #     dust_map = dust_map_downgraded_mjy * dust(freq, nu0_dust, beta_dust, t_dust)
    #     hp.mollview(dust_map, title=f"{int(freq):04d} GHz", unit="MJy/sr", min=0, max=200)
    #     plt.savefig(f"tests/dust_maps/{int(freq):04d}.png")
    #     plt.clf()
    #     plt.close()

    frequencies = np.linspace(0, 13.604162/2 * SPEC_SIZE, SPEC_SIZE)

    signal = dust(frequencies, nu0_dust, beta_dust, t_dust)

    return dust_map_downgraded_mjy, frequencies, signal

# pad f_ghz to 512
# frequencies = np.zeros(257)
# frequencies[0:5] = 0
# frequencies[5 : len(f_ghz["rh_ss"]) + 5] = f_ghz["rh_ss"]
# frequencies[len(f_ghz["rh_ss"]) + 5 :] = 0

if __name__ == "__main__":

    dust_map_downgraded_mjy, frequencies, signal = sim_dust()

    ifg = np.zeros((len(dust_map_downgraded_mjy), IFG_SIZE))
    x = np.linspace(
        0, 1.76, IFG_SIZE
    )  # cm - does this make sense? sky frequencies are cropped but how does that relate to space?
    # nu_icm = (f_ghz["rh_ss"] * u.GHz).to(1 / u.cm, equivalencies=u.spectral()).value  # cm-1
    # get the IFGs from the spectra

    print("Calculating and plotting IFGs")

    # time ifg making
    time_start = time.time()

    # x = np.linspace(0, 1.76, 512)
    for i in range(len(dust_map_downgraded_mjy)):
        for j, freq in enumerate(frequencies):
            ifg[i, j] = np.sum(
                # dust_map_downgraded_mjy[i] *
                signal[j]
                * np.cos(
                    2
                    * np.pi
                    * (frequencies * u.GHz).to(1 / u.cm, equivalencies=u.spectral()).value
                    * (x[j] - 1.22)
                    # * x[j]
                )
            )

        # ifg[i] = np.fft.irfft(dust_map_downgraded_mjy[i] * signal)
        ifg[i] = dust_map_downgraded_mjy[i] * ifg[i]

    # save ifg products in a npz file
    np.savez("tests/ifgs.npz", ifg=ifg)

    time_end = time.time()
    print(f"Time elapsed for IFGs: {(time_end - time_start)/60} minutes")
