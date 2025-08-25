"""
This script is meant to plot the residuals of the calibration data.
"""

import os
import sys

import matplotlib.gridspec as grd
import matplotlib.pyplot as plt
import numpy as np

import globals as g
import utils.my_utils as utils

data = np.load(g.PROCESSED_DATA_PATH_CAL)

modes = {"ss": 0, "lf": 3}
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}


for mode in modes.keys():
    xcal = data[f"xcal_{mode}"]
    sorted_indices = np.argsort(xcal)
    xcal = xcal[sorted_indices]

    ical = data[f"ical_{mode}"]
    ical = ical[sorted_indices]

    for channel in channels.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            cal = data[f"{channel}_{mode}"]

            frequency = utils.generate_frequencies(channel, mode)

            # sort cal data by increasing xcal temp
            cal = cal[sorted_indices]

            ifig = plt.figure(figsize=(8, 8))
            # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 16))
            gs = grd.GridSpec(
                3, 2, height_ratios=[1, 1, 1], width_ratios=[16, 1], wspace=0.1
            )

            ax1 = plt.subplot(gs[0])
            im = ax1.imshow(
                cal.T.real,
                aspect="auto",
                extent=[
                    0,
                    len(cal),
                    frequency[0],
                    frequency[-1],
                ],
                vmax=200,
                vmin=-200,
                cmap="RdBu_r",
                origin="lower",
            )

            ax1.set_title(f"{channel.upper()}{mode.upper()} Spectra Residuals")
            ax1.set_ylabel("Frequency (GHz)")

            colorAx = plt.subplot(gs[1])
            plt.colorbar(im, cax=colorAx, label="Residuals (MJy/sr)")

            ax2 = plt.subplot(gs[2])
            ax2.plot(xcal)
            ax2.set_xlabel("Index")
            ax2.set_ylabel("XCAL Temperature (K)")
            ax2.set_xlim(0, len(cal))

            ax3 = plt.subplot(gs[4])
            ax3.plot(ical)
            ax3.set_xlabel("Index")
            ax3.set_ylabel("ICAL Temperature (K)")
            ax3.set_xlim(0, len(cal))

            plt.tight_layout()
            # plt.show()
            plt.savefig(
                f"/mn/stornext/d16/www_cmb/{os.environ['USER']}/firas/plots/residuals_spectra/{channel}{mode}.png",
                bbox_inches="tight",
            )
            plt.clf()
