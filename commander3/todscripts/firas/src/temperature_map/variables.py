import numpy as np
from astropy.io import fits

# modes = {"ss": 0, "lf": 3}
# modes = {"ss": 0}
mode = "lf"

data = np.load(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/data/sky.npz",
    allow_pickle=True,
)

# print(data.files)
variables_names = [
    "stat_word_16",
    "stat_word_13",
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

print(f"stat_word_16: {variables["stat_word_16"][76700]}")
print(f"stat_word_13: {variables["stat_word_13"][83700]}")
print(f"stat_word_5: {variables["stat_word_5"][83600]}")
print(f"stat_word_9: {variables["stat_word_9"][1067]}")
print(f"lvdt_stat_a: {variables["lvdt_stat_a"][51491]}")
print(f"lvdt_stat_b: {variables["lvdt_stat_b"][51491]}")
