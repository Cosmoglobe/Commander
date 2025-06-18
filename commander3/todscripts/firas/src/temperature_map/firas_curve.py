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
from my_utils import planck, residuals
from utils.config import gen_nyquistl

T_CMB = 2.72548  # Fixsen 2009
modes = {"ss": 0, "lf": 3}
channels = {"rl": 1, "ll": 3}

mask = fits.open("BP_CMB_I_analysis_mask_n1024_v2.fits")
mask = mask[1].data.astype(int)

data = np.load(g.PROCESSED_DATA_PATH)

sky = {}
pix_gal = {}
for channel in channels.keys():
    for mode in modes.keys():
        sky[f"{channel}_{mode}"] = np.abs(data[f"{channel}_{mode}"])
        pix_gal[mode] = data[f"pix_gal_{mode}"]

# frequency mapping
nu0 = {"ss": 68.020812, "lf": 23.807283}
dnu = {"ss": 13.604162, "lf": 3.4010405}
nf = {"lh_ss": 210, "ll_lf": 182, "ll_ss": 43, "rh_ss": 210, "rl_lf": 182, "rl_ss": 43}

f_ghz = {}
for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            f_ghz[f"{channel}_{mode}"] = np.linspace(
                nu0[mode],
                nu0[mode] + dnu[mode] * (nf[f"{channel}_{mode}"] - 1),
                nf[f"{channel}_{mode}"],
            )

# print(f"sky shape: {sky.shape} and pix_gal shape: {pix_gal.shape}")

print("Calculating BB curve")

for channel in channels.keys():
    for mode in modes.keys():
        bb_curve = np.zeros(len(f_ghz[f"{channel}_{mode}"]))

        # subtract dust emission
        nu0_dust = 353  # using Planck values
        beta_dust = 1.55
        t_dust = np.array(20.8)  # u.K
        optical_depth_nu0 = 9.6 * 10 ** (-7)
        dust = (
            optical_depth_nu0
            * planck(f_ghz[f"{channel}_{mode}"], t_dust)
            * (f_ghz[f"{channel}_{mode}"] / nu0_dust) ** beta_dust
        )

        for freq in range(len(f_ghz[f"{channel}_{mode}"])):
            hpxmap = np.zeros(g.NPIX)
            data_density = np.zeros(g.NPIX)

            for i in range(len(pix_gal)):
                hpxmap[pix_gal[mode][i]] += sky[f"{channel}_{mode}"][i][freq]
                data_density[pix_gal[mode][i]] += 1

            m = np.zeros(g.NPIX)
            mask2 = data_density == 0
            m[~mask2] = (hpxmap[~mask2] - dust[freq]) / data_density[~mask2]
            # m[~mask2] = hpxmap[~mask2] / data_density[~mask2]

            # mask the map with the points of no data
            m[mask2] = hp.UNSEEN
            # mask the galaxy
            m[mask] = hp.UNSEEN

            # average all points at this frequency to get each frequency point in the bb curve
            bb_curve[freq] = np.mean(m[m != hp.UNSEEN])

        # print(f"Fitting BB curve: {bb_curve}")

        t0 = np.array(2.726)
        fit = minimize(residuals, t0, args=(f_ghz[f"{channel}_{mode}"], bb_curve))

        print(f"Plotting BB curve: {fit.x[0]}")

        plt.plot(f_ghz[f"{channel}_{mode}"], bb_curve, label="Data")
        plt.plot(
            f_ghz[f"{channel}_{mode}"],
            planck(f_ghz[f"{channel}_{mode}"], fit.x[0]),
            label="Fit",
        )
        plt.plot(
            f_ghz[f"{channel}_{mode}"],
            planck(f_ghz[f"{channel}_{mode}"], t0),
            label="Original",
        )
        plt.plot(f_ghz[f"{channel}_{mode}"], dust, label="Dust")
        plt.xlabel("Frequency [GHz]")
        plt.ylabel("Brightness [MJy/sr]")
        plt.title("BB curve")
        plt.legend()
        plt.savefig(f"{g.SAVE_PATH}/plots/bb_curve/{f"{channel}_{mode}"}.png")
        plt.clf()
