"""
Estimates the PSD by subtracting the consecutive IFGs to get rid of the signal.
"""

import os
import sys

import numpy as np

import utils.config as config
from noise import plots

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
grandparent = os.path.dirname(parent)
sys.path.append(grandparent)

import globals as g
from utils.config import gen_nyquistl

# sky_data = h5py.File(
#     g.PREPROCESSED_DATA_PATH_SKY,
#     "r",
# )

save_path = "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/src/noise/output/sub/"

for channel in g.CHANNELS.keys():
    cal_data = np.load(
        f"{g.PREPROCESSED_DATA_PATH}/cal_{channel}.npz",
    )

    mtm_length = cal_data["mtm_length"]
    mtm_speed = cal_data["mtm_speed"]

    for mode in g.MODES.keys():
        if mode[0] == "s":
            length_filter = mtm_length == 0
        else:
            length_filter = mtm_length == 1
        if mode[1] == "s":
            speed_filter = mtm_speed == 0
        else:
            speed_filter = mtm_speed == 1

        mode_filter = length_filter & speed_filter

        ifgs = cal_data["ifg"][mode_filter]

        ifgs = ifgs - np.median(ifgs, axis=1, keepdims=True)
        ifgs_sub = ifgs[:-1] - ifgs[1:]

        psd = np.abs(np.fft.rfft(ifgs_sub, axis=1)) ** 2

        fnyq_hz = config.gen_nyquistl(
            "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
        )["hz"][4 * (g.CHANNELS[channel] % 2) + g.MODES[mode]]
        print(f"fnyq_hz for channel {channel}, mode {mode}: {fnyq_hz}")
        fnyq_icm = config.gen_nyquistl(
            "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
        )["icm"][4 * (g.CHANNELS[channel] % 2) + g.MODES[mode]]
        freqs_hz = np.fft.rfftfreq(ifgs_sub.shape[1], d=1 / (2 * fnyq_hz))
        freqs_icm = np.fft.rfftfreq(ifgs_sub.shape[1], d=1 / (2 * fnyq_icm))

        plots.plot_psd(
            freqs_hz,
            freqs_icm,
            psd,
            channel,
            mode=mode,
            path=save_path,
        )
        plots.plot_mean_psd(
            psd,
            freqs_hz,
            freqs_icm,
            channel,
            mode=mode,
            path=save_path,
        )
