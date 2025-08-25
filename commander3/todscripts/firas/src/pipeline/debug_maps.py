import os
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
import globals as g
from utils.config import gen_nyquistl

modes = {"ss": 0, "lf": 3}
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}

reference_maps = "/mn/stornext/d16/cmbco/ola/firas/healpix_maps/"
# data = np.load(
#     # g.PROCESSED_DATA_PATH,
#     # g.PROCESSED_DATA_PATH_CAL,
#     allow_pickle=True,
# )

cal_data = h5py.File(
    g.PREPROCESSED_DATA_PATH_CAL,
    "r",
)

print(cal_data["df_data"].keys())

variable_names = [
    # "gmt",
    "a_ical",
    # # "b_ical",
    "a_xcal",
    # # "b_xcal",
    # "a_dihedral",
    # # "b_dihedral",
    # "a_refhorn",
    # # "b_refhorn",
    # "a_skyhorn",
    # # "b_skyhorn",
    # "a_bol_assem_rh",
    "b_bol_assem_rh",
    "a_bol_assem_rl",
    "b_bol_assem_rl",
    "a_bol_assem_lh",
    "b_bol_assem_lh",
    "a_bol_assem_ll",
    "b_bol_assem_ll",
    # # "mtm_length",
    # # "mtm_speed",
    # # "pix_gal",
    # "adds_per_group",
    "bol_cmd_bias_ll",
    "bol_cmd_bias_lh",
    "bol_cmd_bias_rl",
    # "bol_cmd_bias_rh",
    # "dwell_stat_a",
    # # "dwell_stat_b",
    # "ext_cal_temp_a",
    # # "ext_cal_temp_b",
    # "stat_word_1",
    # "stat_word_12",
    # "stat_word_13",
    # "stat_word_16",
    # "stat_word_4",
    # "stat_word_5",
    # "stat_word_8",
    # "stat_word_9",
    # "lvdt_stat_a",
    # "lvdt_stat_b",
    # "power_a_status_a",
    # # "power_a_status_b",
    # # "power_b_status_a",
    # # "power_b_status_b",
    # "ref_hrn_temp_a",
    # # "ref_hrn_temp_b",
    # "micro_stat_bus_1",
    # "micro_stat_bus_2",
    # "micro_stat_bus_3",
    # "micro_stat_bus_4",
    # "grt_addr_a",
    # # "grt_addr_b",
    # "hot_spot_cmd_a",
    # # "hot_spot_cmd_b",
    # "int_ref_temp_a",
    # # "int_ref_temp_b",
    # "sky_hrn_temp_a",
    # # "sky_hrn_temp_b",
    # # "bol_volt_ll",
    # # "bol_volt_lh",
    # # "bol_volt_rl",
    # "bol_volt_rh",
]

fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

n_cols = 4
n_rows = len(variable_names) // n_cols + 1


variables = {}
masks = {}
mtm_speed = cal_data[f"df_data/mtm_speed"][()]
mtm_length = cal_data[f"df_data/mtm_length"][()]
masks["ss"] = (mtm_speed == 0) & (mtm_length == 0)

ifg_ll = cal_data[f"df_data/ifg_ll"][masks["ss"]]
ifg_ll = ifg_ll - np.median(ifg_ll, axis=1)[:, np.newaxis]

for variable_name in variable_names:
    variables[variable_name] = cal_data[f"df_data/{variable_name}"][masks["ss"]]

mask_grt = (
    (variables["a_bol_assem_rh"] > 0)
    & (variables["a_bol_assem_rl"] >0)
    & (variables["a_bol_assem_lh"] > 0)
    & (variables["b_bol_assem_lh"] > 0)
    & (variables["a_bol_assem_ll"] > 0)
    & (variables["b_bol_assem_ll"] > 0)
)

fig, ax = plt.subplots(n_rows, n_cols, figsize=(20, 20), sharex=True)

for i, name in enumerate(variable_names):
    ax[i // n_cols, i % n_cols].plot(variables[name])
    ax[i // n_cols, i % n_cols].set_title(name, fontsize=10)
# ax[1, 1].set_ylim(18536-1, 18536+1)

ax[-1, -1].imshow(
    ifg_ll.T,
    aspect="auto",
    extent=[
        0,
        len(ifg_ll),
        0,
        len(ifg_ll[0]),
    ],
    vmax=500,
    vmin=0,
)
ax[-1, -1].set_title("ifg", fontsize=10)
plt.show()
