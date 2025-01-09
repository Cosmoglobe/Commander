import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

# m = hp.fitsfunc.read_map(
#     "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/maps/fits/0179.fits",
# )

# hp.mollview(m, coord="G", title="179", unit="MJy/sr", min=0, max=500)
# plt.show()


sky = np.load("../../output/data/sky.npy")
data = h5py.File(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v4.1.h5",
    "r",
)

print(sky.shape)

spec = np.mean(np.abs(sky[:, 1:]), axis=1)

variable_names = [
    "a_dihedral",
    "a_ical",
    "b_dihedral",
    "b_ical",
    "ext_cal_temp_a",
    "ext_cal_temp_b",
    "hot_spot_cmd_a",
    "hot_spot_cmd_b",
    "int_ref_temp_a",
    "int_ref_temp_b",
    "lvdt_stat_a",
    "lvdt_stat_b",
    "mtm_length",
    "mtm_speed",
    "power_a_status_a",
    "power_a_status_b",
    "power_b_status_a",
    "power_b_status_b",
    "ref_hrn_temp_a",
    "ref_hrn_temp_b",
    "stat_word_1",
    "stat_word_12",
    "stat_word_13",
    "stat_word_16",
    "stat_word_4",
    "stat_word_5",
    "stat_word_8",
    "stat_word_9",
]
channel_dependent = ["adds_per_group", "bol_cmd_bias", "bol_volt", "gain", "sweeps"]

variables = {}

for name in variable_names:
    variables[name] = np.array(data["df_data/" + name][()])

for name in channel_dependent:
    variables[name] = np.array(data["df_data/" + name + "_ll"][()])

print(len(variables.keys()))

fig, ax = plt.subplots(6, 6, figsize=(20, 20), sharex=True)

for i, name in enumerate(variables.keys()):
    ax[i // 6, i % 6].plot(variables[name])
    ax[i // 6, i % 6].set_title(name)

ax[-1, -1].plot(spec)
ax[-1, -1].set_title("spec")

plt.tight_layout()
plt.show()
