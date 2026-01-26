"""
This script is meant to plot a curve of the intensity of the galaxy and outside the galaxy over the frequencies.
"""

import astropy.units as u
import cosmoglobe as cg
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import globals as g
import utils.my_utils as utils

# data = np.load(g.PROCESSED_DATA_PATH)

modes = {"ss": 0, "lf": 3}
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}

path = "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/fits_files/"
firas_path = "maps/frequency_maps/"

spectral_lines = [
    115,
    230,
    345,
    425,
    460,
    492,
    575,
    690,
    805,
    810,
    920,
    1035,
    1114,
    1150,
    1462,
    1897,
    2053,
    2306,
    2457,
    2584,
]

high_lat_mask = hp.read_map("/mn/stornext/d16/cmbco/ola/masks/HI_mask_4e20_n1024.fits")
mask_alm = hp.sphtfunc.map2alm(high_lat_mask, pol=False)
high_lat_mask = hp.alm2map(mask_alm, g.NSIDE, pol=False)
high_lat_mask = np.where(high_lat_mask < 0.5, 1, 0)
low_lat_mask = hp.read_map(
    "/mn/stornext/d16/cmbco/ola/masks/mask_common_dx12_firas_n0016.fits"
)
mask_alm = hp.sphtfunc.map2alm(low_lat_mask, pol=False)
low_lat_mask = hp.alm2map(mask_alm, g.NSIDE, pol=False)
low_lat_mask = np.where(low_lat_mask < 0.5, 0, 1)

for channel in g.CHANNELS_PLOT:
    for mode in g.MODES_PLOT:
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            print(f"Processing {channel} {mode}")
            f_ghz = utils.generate_frequencies(channel, mode)

            monopole = utils.planck(f_ghz, np.array(g.T_CMB))
            high_lat_amp = np.zeros(len(f_ghz))
            low_lat_amp = np.zeros(len(f_ghz))

            m = np.zeros((g.NPIX, len(f_ghz)))
            for fi, freq in enumerate(f_ghz):
                m[:, fi] = fits.open(
                    f"{path}{firas_path}{channel}_{mode}/galactic/{int(freq):04d}_nside32.fits"
                )[0].data

            # high_lat_map = np.where(
            #     high_lat_mask[:, np.newaxis] == 1, np.nan, m
            # )  # + monopole
            high_lat_map = m
            # print(high_lat_map)
            # check for nans
            if np.isnan(high_lat_map).any():
                print("NaNs found in high latitude map")
            low_lat_map = np.where(low_lat_mask[:, np.newaxis] == 1, np.nan, m)

            # get rid of 5sigma outliers
            # std_high = np.nanstd(high_lat_map, axis=0)
            # median_high = np.nanmedian(high_lat_map, axis=0)
            # dist = np.abs(high_lat_map - median_high) / std_high

            # high_lat_map = (
            #     np.where(dist > 5, np.nan, high_lat_map) + monopole
            # )  # adding monopole back in since the input maps are monopole subtracted
            # std_high = np.nanstd(high_lat_map, axis=0)

            # std_low = np.nanstd(low_lat_map, axis=0)
            # median_low = np.nanmedian(low_lat_map, axis=0)
            # dist = np.abs(low_lat_map - median_low) / std_low

            # low_lat_map = np.where(dist > 5, np.nan, low_lat_map)  # - monopole
            # std_low = np.nanstd(low_lat_map, axis=0)

            high_lat_amp = np.nanmedian(high_lat_map, axis=0)
            print(high_lat_amp)
            low_lat_amp = np.nanmedian(low_lat_map, axis=0)

            for freqi, freq in enumerate(f_ghz):
                fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(12, 12))
                plt.subplots_adjust(wspace=0.3)
                plt.axes(ax1)

                # ax1.fill_between(
                #     f_ghz,
                #     high_lat_amp - std_high,
                #     high_lat_amp + std_high,
                #     alpha=0.2,
                #     color="blue",
                #     label="1sigma",
                # )
                ax1.plot(
                    f_ghz,
                    high_lat_amp,
                    markevery=[freqi],
                    marker="o",
                    markersize=5,
                )
                ax1.set_title(f"Frequency Curve of {channel} {mode} Data")
                ax1.set_xlabel("Frequency (GHz)")
                ax1.set_ylabel("Intensity (MJy/sr)")
                ax1.legend(["High Latitude"])
                # ax1.set_ylim(0, 600)

                if channel == "ll":
                    max = 10
                    min = -10
                else:
                    max = 400
                    min = 0

                plt.axes(ax2)
                hp.mollview(
                    high_lat_map[:, freqi],
                    title=f"High Latitude {channel} {mode} Data @ {int(f_ghz[freqi]):04d} GHz",
                    unit="Intensity (MJy/sr)",
                    hold=True,
                    min=min,
                    max=max,
                    cmap="RdBu_r",
                )
                plt.savefig(
                    f"{g.SAVE_PATH}/plots/frequency_curve/{channel}_{mode}/high_lat/{int(f_ghz[freqi]):04d}.png",
                    bbox_inches="tight",
                )
                plt.close(fig)

                fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(12, 12))
                plt.subplots_adjust(wspace=0.3)
                plt.axes(ax1)
                # ax1.fill_between(
                #     f_ghz,
                #     low_lat_amp - std_low,
                #     low_lat_amp + std_low,
                #     alpha=0.2,
                #     color="blue",
                #     label="1sigma",
                # )
                ax1.plot(
                    f_ghz,
                    low_lat_amp,
                    markevery=[freqi],
                    marker="o",
                    markersize=5,
                )
                if channel[1] == "l":
                    max_freq = 7
                else:
                    max_freq = -1
                ax1.vlines(
                    spectral_lines[:max_freq],
                    np.nanmin(low_lat_amp),
                    np.nanmax(low_lat_amp),
                    linestyles="dashed",
                    colors="red",
                )
                ax1.set_title(
                    f"Frequency Curve of {channel.upper()}{mode.upper()} Data (monopole subtracted)"
                )
                ax1.set_xlabel("Frequency (GHz)")
                ax1.set_ylabel("Intensity (MJy/sr)")
                ax1.legend(["Low Latitude"])
                # ax1.set_ylim(-10, 75)
                plt.axes(ax2)
                if channel == "ll":
                    max = 20
                elif channel == "rh":
                    max = 200
                hp.mollview(
                    low_lat_map[:, freqi],
                    title=f"Low Latitude {channel.upper()}{mode.upper()} Data @ {int(f_ghz[freqi]):04d} GHz",
                    unit="Intensity (MJy/sr)",
                    min=0,
                    max=max,
                    # cmap="jet",
                    hold=True,
                    # norm="symlog2",
                    # cbar=False,
                )
                # cg.plot(
                #     low_lat_map[:, freqi],
                #     title=f"{int(freq):04d} GHz as seen by {channel.upper()}{mode.upper()}",
                #     unit="MJy/sr",
                #     min=0,
                #     max=200,
                #     comp="freqmap",
                #     freq=int(freq) * u.GHz,
                #     cmap="planck",
                #     hold=True,
                # )
                # plt.show()
                plt.savefig(
                    f"{g.SAVE_PATH}/plots/frequency_curve/{channel}_{mode}/low_lat/{int(f_ghz[freqi]):04d}.png",
                    bbox_inches="tight",
                )
                plt.close(fig)
