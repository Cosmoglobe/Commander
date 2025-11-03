"""
Estimates the PSD by subtracting the consecutive IFGs to get rid of the signal.
"""

import os
import sys

import h5py
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
cal_data = h5py.File(
    g.PREPROCESSED_DATA_PATH,
    "r",
)

print(cal_data["df_data/ifg_ll"].shape)

ifgs_ll = cal_data["df_data/ifg_ll"]
ifgs_rl = cal_data["df_data/ifg_rl"]

ifgs_ll = ifgs_ll - np.median(ifgs_ll, axis=1, keepdims=True)
ifgs_rl = ifgs_rl - np.median(ifgs_rl, axis=1, keepdims=True)

ifgs_sub = ifgs_ll[:-1] - ifgs_ll[1:]

psd = np.abs(np.fft.rfft(ifgs_sub, axis=1)) ** 2

fnyq_ll_hz = config.gen_nyquistl(
    "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
)["hz"][4 * (g.CHANNELS["ll"] % 2) + g.MODES["ss"]]
fnyq_ll_icm = config.gen_nyquistl(
    "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
)["icm"][4 * (g.CHANNELS["ll"] % 2) + g.MODES["ss"]]
freqs_ll_hz = np.fft.rfftfreq(ifgs_sub.shape[1], d=1 / (2 * fnyq_ll_hz))
freqs_ll_icm = np.fft.rfftfreq(ifgs_sub.shape[1], d=1 / (2 * fnyq_ll_icm))

plots.plot_psd(freqs_ll_hz, psd, "ll")
