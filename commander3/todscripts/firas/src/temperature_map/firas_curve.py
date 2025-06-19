import os
import sys

import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
import globals as g
import my_utils as mu
from utils.config import gen_nyquistl

T_CMB = 2.72548  # Fixsen 2009
modes = {"ss": 0, "lf": 3}
channels = {"rl": 1, "ll": 3}

mask_gal = hp.read_map("/mn/stornext/d16/cmbco/ola/masks/HI_mask_4e20_n1024.fits")
mask_alm = hp.sphtfunc.map2alm(mask_gal, pol=False)
mask_gal = hp.alm2map(mask_alm, g.NSIDE, pol=False)
mask_gal = np.where(mask_gal < 0.5, 1, 0)

data = np.load(g.PROCESSED_DATA_PATH)

# subtract dust emission
nu0_dust = 353  # using Planck values
beta_dust = 1.55
t_dust = np.array(20.8)  # u.K
optical_depth_nu0 = 9.6 * 10 ** (-7)

for mode in modes.keys():
    pix_gal = data[f"pix_gal_{mode}"]
    for channel in channels.keys():
        sky = np.abs(data[f"{channel}_{mode}"])
        f_ghz = mu.generate_frequencies(channel, mode)

        print("Calculating BB curve")
        
        dust = (
            optical_depth_nu0
            * mu.planck(f_ghz, t_dust)
            * (f_ghz / nu0_dust) ** beta_dust
        )

        hpxmap = np.zeros((g.NPIX, len(f_ghz)))
        data_density = np.zeros(g.NPIX)

        for todi in range(len(sky)):
            hpxmap[pix_gal[todi]] += sky[todi]
            data_density[pix_gal[todi]] += 1

        mask_nodata = data_density == 0
        
        m = (hpxmap - dust) / data_density[:, np.newaxis]
        
        mask_high = np.where(m > 1000, 1, 0)
        mask_low = np.where(m < 0, 1, 0)

        mask_total = mask_nodata[:, np.newaxis] | mask_gal[:, np.newaxis] | mask_high | mask_low

        m = np.where(mask_total == 1, np.nan, m)

        # average all points at this frequency to get each frequency point in the bb curve
        bb_curve = np.nanmean(m, axis=0)

        t0 = np.array(g.T_CMB)
        fit = minimize(mu.residuals, t0, args=(f_ghz, bb_curve))

        print(f"Plotting BB curve: {fit.x[0]}")

        plt.plot(f_ghz, bb_curve, label="Data")
        plt.plot(
            f_ghz,
            mu.planck(f_ghz, fit.x[0]),
            label="Fit",
        )
        plt.plot(
            f_ghz,
            mu.planck(f_ghz, t0),
            label="Original",
        )
        plt.plot(f_ghz, dust, label="Dust")
        plt.xlabel("Frequency [GHz]")
        plt.ylabel("Brightness [MJy/sr]")
        plt.title("BB curve")
        plt.legend()
        plt.savefig(f"{g.SAVE_PATH}/plots/bb_curve/{f"{channel}_{mode}"}.png")
        plt.clf()
