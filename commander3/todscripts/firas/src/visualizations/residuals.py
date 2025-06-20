"""
This script is meant to plot the residuals of the calibration data.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import globals as g
import my_utils as mu

data = np.load(g.PROCESSED_DATA_PATH_CAL)

modes = {"ss": 0, "lf": 3}
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}

for channel in channels.keys():
    for mode in modes.keys():
        if not(mode == "lf" and (channel == "lh" or channel == "rh")):
            cal = data[f"{channel}_{mode}"]
            frequency = mu.generate_frequencies(channel, mode)

            plt.imshow(
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
                cmap='RdBu_r',
                origin='lower',
            )
            plt.colorbar(label="Residuals (MJy/sr)")
            plt.title(f"{channel.upper()}{mode.upper()} Spectra Residuals")
            plt.ylabel("Frequency (GHz)")
            # plt.show()
            plt.savefig(
                f"/mn/stornext/d16/www_cmb/{os.environ['USER']}/firas/plots/residuals_spectra/{channel}{mode}.png"
            )
            plt.clf()

