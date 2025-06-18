import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
import globals as g
from my_utils import planck

channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}
modes = {"ss": 0, "lf": 3}

T_CMB = 2.72548  # Fixsen 2009

data = np.load("../../output/data/lambda_sky.npz", allow_pickle=True)
# print(data.files)
sky = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            sky[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]
            # print(f"{channel}_{mode}", sky[f"{channel}_{mode}"])

            plt.imshow(
                np.abs(sky[f"{channel}_{mode}"]).T,
                aspect="auto",
                extent=[
                    0,
                    len(sky[f"{channel}_{mode}"]),
                    0,
                    len(sky[f"{channel}_{mode}"][0]),
                ],
                vmin=0,
                vmax=500,
            )
            plt.savefig(
                f"{g.SAVE_PATH}plots/sky_over_time_lambda_{f"{channel}_{mode}"}.png"
            )
            plt.clf()


# frequency mapping
nu0 = {"ss": 68.020812, "lf": 23.807283}
dnu = {"ss": 13.604162, "lf": 3.4010405}
nf = {"lh_ss": 210, "ll_lf": 182, "ll_ss": 43, "rh_ss": 210, "rl_lf": 182, "rl_ss": 43}

f_ghz = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            f_ghz[f"{channel}_{mode}"] = np.linspace(
                nu0[mode],
                nu0[mode] + dnu[mode] * (nf[f"{channel}_{mode}"] - 1),
                nf[f"{channel}_{mode}"],
            )

# get the position on the sky
data_path = "/mn/stornext/d16/cmbco/ola/firas/coadded_interferograms/"
pix_gal = {}
for channel in channels.keys():
    for mode in modes.keys():
        data = fits.open(
            data_path
            + f"FIRAS_COADDED_SKY_INTERFEROGRAMS_{channel.upper()}{mode.upper()}.FITS"
        )[1].data

        gal_lon = data["GAL_LON"]
        gal_lat = data["GAL_LAT"]
        pix_gal[f"{channel}_{mode}"] = hp.ang2pix(
            g.NSIDE, gal_lon, gal_lat, lonlat=True
        ).astype(int)

hpxmap = {}
data_density = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            hpxmap[f"{channel}_{mode}"] = np.zeros(
                (g.NPIX, len(f_ghz[f"{channel}_{mode}"]))
            )
            data_density[f"{channel}_{mode}"] = np.zeros(g.NPIX)
            for i in range(len(pix_gal[f"{channel}_{mode}"])):
                hpxmap[f"{channel}_{mode}"][pix_gal[f"{channel}_{mode}"][i]] += np.abs(
                    sky[f"{channel}_{mode}"][i]
                )
                data_density[f"{channel}_{mode}"][pix_gal[f"{channel}_{mode}"][i]] += 1

m = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            m[f"{channel}_{mode}"] = np.zeros((g.NPIX, len(sky[f"{channel}_{mode}"][0])))
            mask = data_density[f"{channel}_{mode}"] == 0
            monopole = planck(f_ghz[f"{channel}_{mode}"], np.array(T_CMB))
            m[f"{channel}_{mode}"][~mask] = (
                hpxmap[f"{channel}_{mode}"][~mask]
                / data_density[f"{channel}_{mode}"][~mask, np.newaxis]
            ) - monopole
            m[f"{channel}_{mode}"][mask] = hp.UNSEEN

for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            for freq in range(len(f_ghz[f"{channel}_{mode}"])):
                hp.mollview(
                    m[f"{channel}_{mode}"][:, freq],
                    title=f"{f_ghz[f"{channel}_{mode}"][freq]:.2f} GHz",
                    unit="MJy/sr",
                    # norm="hist",
                    min=0,
                    max=200,
                )
                plt.savefig(
                    f"{g.SAVE_PATH}maps/lambda_ifg/{f"{channel}_{mode}"}/{int(f_ghz[f"{channel}_{mode}"][freq]):04d}.png"
                )
                plt.close()
