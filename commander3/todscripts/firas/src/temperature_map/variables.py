import globals as g
import numpy as np
from astropy.io import fits

# modes = {"ss": 0, "lf": 3}
# modes = {"ss": 0}
mode = "ss"

data = np.load(
    g.PROCESSED_DATA_PATH,
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

mask = fits.open("BP_CMB_I_analysis_mask_n1024_v2.fits")
mask = mask[1].data.astype(int)

variables = {}
print(f"mode: {mode}")
pix_gal = data[f"pix_gal_{mode}"]
for name in variables_names:
    variables[name] = data[f"{name}_{mode}"][mask[pix_gal] == 1]

# print shapes
for name in variables_names:
    print(f"{name}: {variables[name].shape}")

ind = 1223

print(f"stat_word_1: {variables["stat_word_1"][ind]}")
print(f"stat_word_12: {variables["stat_word_12"][ind]}")
print(f"stat_word_13: {variables["stat_word_13"][ind]}")
print(f"stat_word_16: {variables["stat_word_16"][ind]}")
print(f"stat_word_5: {variables["stat_word_5"][ind]}")
print(f"stat_word_9: {variables["stat_word_9"][ind]}")
print(f"lvdt_stat_a: {variables["lvdt_stat_a"][ind]}")
print(f"lvdt_stat_b: {variables["lvdt_stat_b"][ind]}")
