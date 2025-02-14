"""
Script to take the previously generated cal spectra (cal.npy) and plot them
"""

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from my_utils import planck
import os

T_CMB = 2.72548  # Fixsen 2009
modes = {"ss": 0, "lf": 3}
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}

data = np.load("../../data/processed_cal.npz")
user = os.environ['USER']
# print(data.files)

cal = {}
pix_gal = {}
scan = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            cal[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]
            cal[f"xcal_{mode}"] = data[f"xcal_{mode}"]
            cal[f"ical_{mode}"] = data[f"ical_{mode}"]


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
            # print(f"{channel}_{mode}: {f_ghz[f"{channel}_{mode}"]}")

print("plotting cal")


from pathlib import Path
Path(f"/mn/stornext/d16/www_cmb/{user}/firas/plots/cal_over_diff").mkdir(parents=True, exist_ok=True)
for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            # plot the cal
            print(f"plotting {channel}_{mode}")
            inds = np.argsort(cal[f"xcal_{mode}"] - cal[f"ical_{mode}"])
            plt.imshow(
                cal[f"{channel}_{mode}"][inds].T.real,
                aspect="auto",
                extent=[
                    0,
                    len(cal[f"{channel}_{mode}"]),
                    0,
                    len(cal[f"{channel}_{mode}"][0]),
                ],
                vmax=200,
                vmin=-200,
                cmap='RdBu_r',
            )
            plt.title(f"{channel}_{mode}")
            plt.colorbar(label=r'$\mathrm{MJy/sr}$')
            plt.savefig(
                f"/mn/stornext/d16/www_cmb/{user}/firas/plots/cal_over_diff/{f"{channel}_{mode}"}_real.png"
            )
            plt.clf()

for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            # plot the cal
            print(f"plotting {channel}_{mode}")
            inds = np.argsort(cal[f"xcal_{mode}"] - cal[f"ical_{mode}"])
            plt.imshow(
                cal[f"{channel}_{mode}"][inds].T.imag,
                aspect="auto",
                extent=[
                    0,
                    len(cal[f"{channel}_{mode}"]),
                    0,
                    len(cal[f"{channel}_{mode}"][0]),
                ],
                vmax=100,
                vmin=-100,
                cmap='RdBu_r',
            )
            plt.title(f"{channel}_{mode}")
            plt.colorbar(label=r'$\mathrm{MJy/sr}$')
            plt.savefig(
                f"/mn/stornext/d16/www_cmb/{user}/firas/plots/cal_over_diff/{f"{channel}_{mode}"}_imag.png"
            )
            plt.clf()

for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            errors = (abs(cal[f'{channel}_{mode}'])).mean(axis=1)
            plt.scatter(cal[f"xcal_{mode}"],cal[f"ical_{mode}"], 
                    s=1, c=errors, vmin=0, vmax=200)
            plt.xlabel(r'$T_\mathrm{xcal}$')
            plt.ylabel(r'$T_\mathrm{ical}$')
            plt.colorbar(label='Mean of residuals')
            plt.savefig(f"/mn/stornext/d16/www_cmb/{user}/firas/plots/cal_scatter/{f"{channel}_{mode}"}.png",
                    bbox_inches='tight')
            plt.close()


for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            # plot the cal
            print(f"plotting {channel}_{mode}")
            plt.imshow(
                cal[f"{channel}_{mode}"].T.real,
                aspect="auto",
                extent=[
                    0,
                    len(cal[f"{channel}_{mode}"]),
                    0,
                    len(cal[f"{channel}_{mode}"][0]),
                ],
                vmax=200,
                vmin=-200,
                cmap='RdBu_r',
            )
            plt.title(f"{channel}_{mode}")
            plt.colorbar(label=r'$\mathrm{MJy/sr}$')
            plt.savefig(
                f"/mn/stornext/d16/www_cmb/{user}/firas/plots/cal_over_time/{f"{channel}_{mode}"}_real.png"
            )
            plt.clf()

for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            # plot the cal
            print(f"plotting {channel}_{mode}")
            plt.imshow(
                cal[f"{channel}_{mode}"].T.imag,
                aspect="auto",
                extent=[
                    0,
                    len(cal[f"{channel}_{mode}"]),
                    0,
                    len(cal[f"{channel}_{mode}"][0]),
                ],
                vmax=100,
                vmin=-100,
                cmap='RdBu_r',
            )
            plt.title(f"{channel}_{mode}")
            plt.colorbar(label=r'$\mathrm{MJy/sr}$')
            plt.savefig(
                f"/mn/stornext/d16/www_cmb/{user}/firas/plots/cal_over_time/{f"{channel}_{mode}"}_imag.png"
            )
            plt.clf()

for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            # plot the cal
            print(f"plotting {channel}_{mode}")
            inds = np.argsort(cal[f"ical_{mode}"])
            plt.imshow(
                cal[f"{channel}_{mode}"][inds].T.real,
                aspect="auto",
                extent=[
                    0,
                    len(cal[f"{channel}_{mode}"]),
                    0,
                    len(cal[f"{channel}_{mode}"][0]),
                ],
                vmax=200,
                vmin=-200,
                cmap='RdBu_r',
            )
            plt.title(f"{channel}_{mode}")
            plt.colorbar(label=r'$\mathrm{MJy/sr}$')
            plt.savefig(
                f"/mn/stornext/d16/www_cmb/{user}/firas/plots/cal_over_ical/{f"{channel}_{mode}"}_real.png"
            )
            plt.clf()

for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            # plot the cal
            print(f"plotting {channel}_{mode}")
            inds = np.argsort(cal[f"ical_{mode}"])
            plt.imshow(
                cal[f"{channel}_{mode}"][inds].T.imag,
                aspect="auto",
                extent=[
                    0,
                    len(cal[f"{channel}_{mode}"]),
                    0,
                    len(cal[f"{channel}_{mode}"][0]),
                ],
                vmax=100,
                vmin=-100,
                cmap='RdBu_r',
            )
            plt.title(f"{channel}_{mode}")
            plt.colorbar(label=r'$\mathrm{MJy/sr}$')
            plt.savefig(
                f"/mn/stornext/d16/www_cmb/{user}/firas/plots/cal_over_ical/{f"{channel}_{mode}"}_imag.png"
            )
            plt.clf()

for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            # plot the cal
            print(f"plotting {channel}_{mode}")
            inds = np.argsort(cal[f"xcal_{mode}"])
            plt.imshow(
                cal[f"{channel}_{mode}"][inds].T.real,
                aspect="auto",
                extent=[
                    0,
                    len(cal[f"{channel}_{mode}"]),
                    0,
                    len(cal[f"{channel}_{mode}"][0]),
                ],
                vmax=200,
                vmin=-200,
                cmap='RdBu_r',
            )
            plt.title(f"{channel}_{mode}")
            plt.colorbar(label=r'$\mathrm{MJy/sr}$')
            plt.savefig(
                f"/mn/stornext/d16/www_cmb/{user}/firas/plots/cal_over_xcal/{f"{channel}_{mode}"}_real.png"
            )
            plt.clf()

for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            # plot the cal
            print(f"plotting {channel}_{mode}")
            inds = np.argsort(cal[f"xcal_{mode}"])
            plt.imshow(
                cal[f"{channel}_{mode}"][inds].T.imag,
                aspect="auto",
                extent=[
                    0,
                    len(cal[f"{channel}_{mode}"]),
                    0,
                    len(cal[f"{channel}_{mode}"][0]),
                ],
                vmax=100,
                vmin=-100,
                cmap='RdBu_r',
            )
            plt.title(f"{channel}_{mode}")
            plt.colorbar(label=r'$\mathrm{MJy/sr}$')
            plt.savefig(
                f"/mn/stornext/d16/www_cmb/{user}/firas/plots/cal_over_xcal/{f"{channel}_{mode}"}_imag.png"
            )
            plt.clf()

