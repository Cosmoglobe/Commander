"""
This script is meant to plot a curve of the intensity of the galaxy and outside the galaxy over the frequencies.
"""

import os
import sys

import numpy as np

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import globals as g
from my_utils import generate_frequencies

data = np.load(g.PROCESSED_DATA_PATH)

modes = {"ss": 0, "lf": 3}
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}

f_ghz = {}
pix_gal = {}
sky = {}
hpxmap = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            f_ghz[f"{channel}_{mode}"] = generate_frequencies(channel, mode)

            pix_gal[mode] = data[f"pix_gal_{mode}"]
            sky[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]

            hpxmap[f"{channel}_{mode}"] = np.zeros(g.NPIX)

            for i in range(len(pix_gal[mode])):
                hpxmap[f"{channel}_{mode}"] = 