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

            high_lat_amp = np.nanmean(high_lat_map, axis=0)
            low_lat_amp = np.nanmean(low_lat_map, axis=0)
            plt.plot(
                f_ghz,
                high_lat_amp,
                label=f"{channel}_{mode} high lat",
                marker="o",
                markersize=2,
            )
            plt.show()
            plt.plot(
                f_ghz,
                low_lat_amp,
                label=f"{channel}_{mode} low lat",
                marker="x",
                markersize=2,
            )
            plt.xlabel("Frequency (GHz)")
            plt.ylabel("Intensity (MJy/sr)")
            plt.title("Frequency Curve of Firas Data")
            plt.legend()
            plt.grid()
            plt.show()