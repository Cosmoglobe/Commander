
import os
import sys

import healpy as hp
import numpy as np

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
import globals as g

# modes = {"ss": 0, "lf": 3}
# modes = {"ss": 0}
mode = "ss"

data = np.load(
    g.PREPROCESSED_DATA_PATH_CAL,
    allow_pickle=True,
)

# print(data.files)
variables_names = [
    "stat_word_1",
    "stat_word_12",
    "stat_word_13",
    "stat_word_16",
    "stat_word_5",
    "stat_word_9",
    "lvdt_stat_a",
    "lvdt_stat_b",
]

mask = hp.read_map("/mn/stornext/d16/cmbco/ola/masks/HI_mask_4e20_n1024.fits")
mask = hp.sphtfunc.map2alm(mask, pol=False)
mask = hp.alm2map(mask, g.NSIDE, pol=False)
mask = np.where(mask < 0.5, 0, 1)

variables = {}
print(f"mode: {mode}")
pix_gal = data[f"pix_gal_{mode}"]
for name in variables_names:
    variables[name] = data[f"{name}_{mode}"][mask[pix_gal] == 1]

# print shapes
for name in variables_names:
    print(f"{name}: {variables[name].shape}")

ind = 10

print(f"stat_word_1: {variables["stat_word_1"][ind]}")
print(f"stat_word_12: {variables["stat_word_12"][ind]}")
print(f"stat_word_13: {variables["stat_word_13"][ind]}")
print(f"stat_word_16: {variables["stat_word_16"][ind]}")
print(f"stat_word_5: {variables["stat_word_5"][ind]}")
print(f"stat_word_9: {variables["stat_word_9"][ind]}")
print(f"lvdt_stat_a: {variables["lvdt_stat_a"][ind]}")
print(f"lvdt_stat_b: {variables["lvdt_stat_b"][ind]}")
