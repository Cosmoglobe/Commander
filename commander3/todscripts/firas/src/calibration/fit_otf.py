"""
Script to fit an optical transfer function to the calibration data. For now, it follows the model:
S_sky/XCAL = 1/OTF * 1/ETF * 1/Bol * Y - R,
where Y is the Fourier transformed interferogram and R is defined by
R = 1/OTF * sum over all emitters i (excluding XCAL) of E_i * P(T_i).
"""

import os
import sys

import numpy as np

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import globals as g
import my_utils as mu

data = np.load(g.PREPROCESSED_DATA_PATH_CAL)

channels = {"ll": 3}
modes = {"ss": 0}

for mode in modes:
    xcal = data[f"xcal_{mode}"][:]
    ical = data[f"ical_{mode}"][:]
    dihedral = data[f"dihedral_{mode}"][:]
    refhorn = data[f"refhorn_{mode}"][:]
    skyhorn = data[f"skyhorn_{mode}"][:]
    collimator = data[f"collimator_{mode}"][:]
    bolometer_ll = data[f"bolometer_ll_{mode}"][:]
    bolometer_lh = data[f"bolometer_lh_{mode}"][:]
    bolometer_rl = data[f"bolometer_rl_{mode}"][:]
    bolometer_rh = data[f"bolometer_rh_{mode}"][:]

    temps = {
        "xcal": xcal,
        "ical": ical,
        "dihedral": dihedral,
        "refhorn": refhorn,
        "skyhorn": skyhorn,
        "collimator": collimator,
        "bolometer_ll": bolometer_ll,
        "bolometer_lh": bolometer_lh,
        "bolometer_rl": bolometer_rl,
        "bolometer_rh": bolometer_rh,
    }

    adds_per_group = data[f"adds_per_group_{mode}"][:]
    sweeps = data[f"sweeps_{mode}"][:]

    for channel in channels:
        bol_cmd_bias = data[f"bol_cmd_bias_{channel}_{mode}"][:]
        bol_volt = data[f"bol_volt_{channel}_{mode}"][:]
        gain = data[f"gain_{channel}_{mode}"][:]

        # Calculate the OTF for the channel
        frequencies = mu.generate_frequencies(channel, mode, 257)
        otf = np.zeros(257, dtype=np.complex128)