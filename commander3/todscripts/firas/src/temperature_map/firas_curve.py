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

        print(f"Calculating BB curve for {channel.upper()}{mode.upper()}...")
        
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
        mask_total = mask_nodata[:, np.newaxis] | mask_gal[:, np.newaxis]
        m = np.where(mask_total == 1, np.nan, m)

        # bb_curve = np.nanmean(m, axis=0)

        t0 = np.array(g.T_CMB)
        temps = np.zeros(g.NPIX)
        for pixi in range(g.NPIX):
            if not np.isnan(m[pixi, :]).all():
                temps[pixi] = minimize(mu.residuals, t0, args=(f_ghz, m[pixi, :])).x[0]
            else:
                temps[pixi] = np.nan

        hp.mollview(
            (temps - g.T_CMB)*1e3,
            title=f"{channel.upper()}{mode.upper()} Deviation from Fixsen 2009 CMB temperature",
            unit="uK",
            min=-100,
            max=100,
            cmap="seismic",
            norm="linear",
            fig=1,
        )
        plt.savefig(f"{g.SAVE_PATH}/maps/bb_temp/{channel}_{mode}.png")
        plt.close()

        # get rid of 5sigma outliers
        std = np.nanstd(temps)
        median = np.nanmedian(temps)
        print(f"Median temperature: {median:.5f} K, Standard deviation: {std:.5f} K")
        dist = np.abs(temps - median) / std
        temps = np.where(dist > 5, np.nan, temps)
        bb_temp = np.nanmean(temps)

        m_excl_outliers = np.where(dist[:, np.newaxis] > 5, np.nan, m)
        bb_curve_data = np.nanmean(m_excl_outliers, axis=0)

        # plt.plot(f_ghz, bb_curve, label="Data")
        tf = minimize(mu.residuals, t0, args=(f_ghz, bb_curve_data)).x[0]
        print(f"BB temperature according to {channel.upper()}{mode.upper()} average temperature after fitting per pixel and fitting after averaging over sky: {bb_temp:.5f} K and {tf:.5f} K")
        plt.plot(
            f_ghz,
            mu.planck(f_ghz, bb_temp),
            label=f"Averaged after fitting: {bb_temp:.5f} K",
        )
        plt.plot(
            f_ghz,
            mu.planck(f_ghz, tf),
            label=f"Fitting after averaging: {tf:.5f} K",
        )
        plt.plot(f_ghz, bb_curve_data, label="Data")
        plt.plot(
            f_ghz,
            mu.planck(f_ghz, t0),
            label="Original",
        )
        plt.plot(f_ghz, dust, label="Dust")
        plt.xlabel("Frequency [GHz]")
        plt.ylabel("Brightness [MJy/sr]")
        plt.title(f"{channel.upper()}{mode.upper()}")
        plt.grid()
        plt.legend()
        plt.savefig(f"{g.SAVE_PATH}/plots/bb_curve/{f"{channel}_{mode}"}.png")
        plt.clf()
