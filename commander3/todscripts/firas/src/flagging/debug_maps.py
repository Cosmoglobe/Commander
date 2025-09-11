import os
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np

import globals as g
from flagging import filter
from utils.config import gen_nyquistl

modes = {"ss": 0, "lf": 3}
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}

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
    # "b_bol_assem_rh",
    # "a_bol_assem_rl",
    # "b_bol_assem_rl",
    # "a_bol_assem_lh",
    # "b_bol_assem_lh",
    # "a_bol_assem_ll",
    # "b_bol_assem_ll",
    # # "mtm_length",
    # # "mtm_speed",
    # # "pix_gal",
    # "adds_per_group",
    # "bol_cmd_bias_ll",
    # "bol_cmd_bias_lh",
    # "bol_cmd_bias_rl",
    # "bol_cmd_bias_rh",
    "dwell_stat_a",
    "dwell_stat_b",
    # "ext_cal_temp_a",
    # # "ext_cal_temp_b",
    "stat_word_1",
    "stat_word_12",
    "stat_word_13",
    "stat_word_16",
    "stat_word_4",
    "stat_word_5",
    "stat_word_8",
    "stat_word_9",
    "lvdt_stat_a",
    "lvdt_stat_b",
    "power_a_status_a",
    "power_a_status_b",
    "power_b_status_a",
    "power_b_status_b",
    # "ref_hrn_temp_a",
    # # "ref_hrn_temp_b",
    # "micro_stat_bus_1",
    # "micro_stat_bus_2",
    # "micro_stat_bus_3",
    # "micro_stat_bus_4",
    "grt_addr_a",
    "grt_addr_b",
    "hot_spot_cmd_a",
    "hot_spot_cmd_b",
    # "int_ref_temp_a",
    # # "int_ref_temp_b",
    # "sky_hrn_temp_a",
    # # "sky_hrn_temp_b",
    # # "bol_volt_ll",
    # # "bol_volt_lh",
    # # "bol_volt_rl",
    # "bol_volt_rh",
]


mtm_speed = cal_data[f"df_data/mtm_speed"][()]
mtm_length = cal_data[f"df_data/mtm_length"][()]

masks = {}
masks["ss"] = (mtm_speed == 0) & (mtm_length == 0)

ifg_ll = cal_data[f"df_data/ifg_ll"][masks["ss"]]
ifg_ll = ifg_ll - np.median(ifg_ll, axis=1)[:, np.newaxis]

stat_word_9 = cal_data[f"df_data/stat_word_9"][masks["ss"]]
lvdt_stat_b = cal_data[f"df_data/lvdt_stat_b"][masks["ss"]]

flags = filter.flag(stat_word_9, lvdt_stat_b)

a_bol_assem_rh = cal_data[f"df_data/a_bol_assem_rh"][masks["ss"]]
b_bol_assem_rh = cal_data[f"df_data/b_bol_assem_rh"][masks["ss"]]
a_bol_assem_rl = cal_data[f"df_data/a_bol_assem_rl"][masks["ss"]]
b_bol_assem_rl = cal_data[f"df_data/b_bol_assem_rl"][masks["ss"]]
a_bol_assem_lh = cal_data[f"df_data/a_bol_assem_lh"][masks["ss"]]
b_bol_assem_lh = cal_data[f"df_data/b_bol_assem_lh"][masks["ss"]]
a_bol_assem_ll = cal_data[f"df_data/a_bol_assem_ll"][masks["ss"]]
b_bol_assem_ll = cal_data[f"df_data/b_bol_assem_ll"][masks["ss"]]
bol_cmd_bias_rh = cal_data[f"df_data/bol_cmd_bias_rh"][masks["ss"]]
bol_cmd_bias_rl = cal_data[f"df_data/bol_cmd_bias_rl"][masks["ss"]]
bol_cmd_bias_lh = cal_data[f"df_data/bol_cmd_bias_lh"][masks["ss"]]
bol_cmd_bias_ll = cal_data[f"df_data/bol_cmd_bias_ll"][masks["ss"]]

filter_bol = filter.filter_bol(
    a_bol_assem_rh,
    a_bol_assem_rl,
    a_bol_assem_lh,
    a_bol_assem_ll,
    b_bol_assem_rh,
    b_bol_assem_rl,
    b_bol_assem_lh,
    b_bol_assem_ll,
    bol_cmd_bias_rh,
    bol_cmd_bias_rl,
    bol_cmd_bias_lh,
    bol_cmd_bias_ll,
)

index = cal_data["df_data/id"][masks["ss"]]

idx = filter.flag_bad_ifgs()
# mask the corresponding indices
flag_idx = ~np.isin(index, idx)

ifg_ll = ifg_ll[flags & filter_bol & flag_idx]

variables_plot = {}
for name in variable_names:
    variables_plot[name] = cal_data[f"df_data/{name}"][masks["ss"]]
    variables_plot[name] = variables_plot[name][flags & filter_bol & flag_idx]
    print(f"{name}: {np.min(variables_plot[name])} to {np.max(variables_plot[name])}")

n_cols = 4
n_rows = (len(variable_names)) // n_cols + 1
fig, ax = plt.subplots(n_rows, n_cols, figsize=(20, 20), sharex=True)

i = 0
for key in variables_plot.keys():
    ax[i // n_cols, i % n_cols].plot(variables_plot[key])
    ax[i // n_cols, i % n_cols].set_title(key, fontsize=10)
    i += 1

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
    interpolation="none",
)
ax[-1, -1].set_title("ifg", fontsize=10)
plt.show()
