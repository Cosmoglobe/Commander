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
high_lat_mask = np.where(high_lat_mask < 0.5, 0, 1)
low_lat_mask = hp.read_map("/mn/stornext/d16/cmbco/ola/masks/mask_common_dx12_firas_n0016.fits")
mask_alm = hp.sphtfunc.map2alm(low_lat_mask, pol=False)
low_lat_mask = hp.alm2map(mask_alm, g.NSIDE, pol=False)
low_lat_mask = np.logical_not(low_lat_mask)
low_lat_mask = np.where(low_lat_mask < 0.5, 0, 1)

for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            f_ghz = mu.generate_frequencies(channel, mode)

            pix_gal = data[f"pix_gal_{mode}"]
            sky = np.abs(data[f"{channel}_{mode}"])

            hpxmap = np.zeros((g.NPIX, len(sky[0])))
            data_density = np.zeros(g.NPIX)

            for i in range(len(pix_gal)):
                hpxmap[pix_gal[i]] += sky[i]
                data_density[pix_gal[i]] += 1
            
            monopole = mu.planck(f_ghz, np.array(g.T_CMB))
            mask = data_density == 0
            print(f"Shape of hpxmap: {hpxmap.shape}, data_density: {data_density.shape}, monopole: {monopole.shape}")
            hpxmap = hpxmap / data_density[:, np.newaxis] - monopole
            hpxmap[mask] = np.nan

            high_lat = np.where((high_lat_mask == 0)[:, np.newaxis], np.nan, hpxmap)
            low_lat = np.where((low_lat_mask == 0)[:, np.newaxis], np.nan, hpxmap)
            for i in range(len(f_ghz)):
                hp.mollview(
                    high_lat[:, i],
                    title=f"{channel} {mode} high lat",
                    unit="MJy/sr",
                    min=0,
                    max=200,
                )
                hp.graticule()
                plt.show()