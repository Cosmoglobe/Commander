"""
This script is meant to plot a curve of the intensity of the galaxy and outside the galaxy over the frequencies.
"""

import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import globals as g
import my_utils as mu

data = np.load(g.PROCESSED_DATA_PATH)

modes = {"ss": 0, "lf": 3}
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}

spectral_lines = [115, 230, 345, 425, 460, 492, 575, 690, 805, 810, 920, 1035, 1114, 1150, 1462, 1897, 2053, 2306, 2457, 2584]

high_lat_mask = hp.read_map("/mn/stornext/d16/cmbco/ola/masks/HI_mask_4e20_n1024.fits")
mask_alm = hp.sphtfunc.map2alm(high_lat_mask, pol=False)
high_lat_mask = hp.alm2map(mask_alm, g.NSIDE, pol=False)
high_lat_mask = np.where(high_lat_mask < 0.5, 1, 0)
low_lat_mask = hp.read_map("/mn/stornext/d16/cmbco/ola/masks/mask_common_dx12_firas_n0016.fits")
mask_alm = hp.sphtfunc.map2alm(low_lat_mask, pol=False)
low_lat_mask = hp.alm2map(mask_alm, g.NSIDE, pol=False)
low_lat_mask = np.where(low_lat_mask < 0.5, 0, 1)

for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            f_ghz = mu.generate_frequencies(channel, mode)

            pix_gal = data[f"pix_gal_{mode}"]
            sky = np.abs(data[f"{channel}_{mode}"])
            monopole = mu.planck(f_ghz, np.array(g.T_CMB))
            high_lat_amp = np.zeros(len(f_ghz))
            low_lat_amp = np.zeros(len(f_ghz))

            hpxmap = np.zeros((g.NPIX, len(f_ghz)))
            data_density = np.zeros(g.NPIX)

            for todi in range(len(pix_gal)):
                hpxmap[pix_gal[todi]] += sky[todi]
                data_density[pix_gal[todi]] += 1

            mask = data_density == 0
            high_lat_mask = high_lat_mask | mask
            low_lat_mask = low_lat_mask | mask  

            m = hpxmap / data_density[:, np.newaxis]  # divide by data density

            high_lat_map = np.where(high_lat_mask[:, np.newaxis] == 1, np.nan, m)
            low_lat_map = np.where(low_lat_mask[:, np.newaxis] == 1, np.nan, m)

            # get rid of 5sigma outliers
            std_high = np.nanstd(high_lat_map, axis=0)
            median_high = np.nanmedian(high_lat_map, axis=0)
            dist = np.abs(high_lat_map - median_high)/std_high
            
            high_lat_map = np.where(dist > 5, np.nan, high_lat_map)
            std_high = np.nanstd(high_lat_map, axis=0)

            std_low = np.nanstd(low_lat_map, axis=0)
            median_low = np.nanmedian(low_lat_map, axis=0)
            dist = np.abs(low_lat_map - median_low)/std_low

            low_lat_map = np.where(dist > 5, np.nan, low_lat_map) - monopole
            std_low = np.nanstd(low_lat_map, axis=0)

            high_lat_amp = np.nanmean(high_lat_map, axis=0)
            low_lat_amp = np.nanmean(low_lat_map, axis=0)


            for freqi in range(len(f_ghz)):
                fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(12, 12))
                plt.subplots_adjust(wspace=0.3)
                plt.axes(ax1)

                ax1.fill_between(f_ghz, high_lat_amp - std_high, high_lat_amp + std_high, alpha=0.2, color="blue", label="1sigma")
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

                plt.axes(ax2)
                hp.mollview(
                    high_lat_map[:, freqi],
                    title=f"High Latitude {channel} {mode} Data",
                    unit="Intensity (MJy/sr)",
                    cmap="jet",
                    hold=True,
                    min=np.nanmin(high_lat_amp),
                    max=np.nanmax(high_lat_amp),
                )
                # plt.show()
                plt.savefig(f"{g.SAVE_PATH}/plots/frequency_curve/{channel}_{mode}/high_lat/{int(f_ghz[freqi]):04d}.png")
                plt.clf()

                fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(12, 12))
                plt.subplots_adjust(wspace=0.3)
                plt.axes(ax1)
                ax1.fill_between(f_ghz, low_lat_amp - std_low, low_lat_amp + std_low, alpha=0.2, color="blue", label="1sigma")
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
                ax1.vlines(spectral_lines[:max_freq], np.nanmin(low_lat_amp), np.nanmax(low_lat_amp), linestyles="dashed", colors="red")
                ax1.set_title(f"Frequency Curve of {channel.upper()}{mode.upper()} Data (monopole subtracted)")
                ax1.set_xlabel("Frequency (GHz)")
                ax1.set_ylabel("Intensity (MJy/sr)")
                ax1.legend(["Low Latitude"])
                # ax1.set_ylim(-10, 75)
                plt.axes(ax2)
                hp.mollview(
                    low_lat_map[:, freqi],
                    title=f"Low Latitude {channel.upper()}{mode.upper()} Data",
                    unit="Intensity (MJy/sr)",
                    min=np.nanmin(low_lat_amp),
                    max=np.nanmax(low_lat_amp),
                    cmap="jet",
                    hold=True
                )
                # plt.show()
                plt.savefig(f"{g.SAVE_PATH}/plots/frequency_curve/{channel}_{mode}/low_lat/{int(f_ghz[freqi]):04d}.png")
                plt.clf()