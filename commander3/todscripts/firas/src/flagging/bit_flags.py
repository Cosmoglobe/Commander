"""
This script looks at the different existing flags to see if looking at them in binary points out any more obvious bit flags. The relevant flags should be (at least for calibration):
- stat_word_5
- stat_word_9
- stat_word_12
- stat_word_13
- stat_word_16
- lvdt_stat_a
- lvdt_stat_b
- hot_spot_cmd_a
- hot_spot_cmd_b
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np

import globals as g

cal_data = h5py.File(
    g.PREPROCESSED_DATA_PATH_CAL,
    "r",
)

stat_word_5 = cal_data["df_data/stat_word_5"][:].astype(int)
stat_word_9 = cal_data["df_data/stat_word_9"][:].astype(int)
stat_word_12 = cal_data["df_data/stat_word_12"][:].astype(int)
stat_word_13 = cal_data["df_data/stat_word_13"][:].astype(int)
stat_word_16 = cal_data["df_data/stat_word_16"][:].astype(int)

ifg_ll = cal_data["df_data/ifg_ll"][:]
ifg_lh = cal_data["df_data/ifg_lh"][:]
ifg_rl = cal_data["df_data/ifg_rl"][:]
ifg_rh = cal_data["df_data/ifg_rh"][:]

mtm_speed = cal_data["df_data/mtm_speed"][:]
mtm_length = cal_data["df_data/mtm_length"][:]

idx = cal_data["df_data/index"][:]
print(f"idx: {idx}")

sw5 = np.zeros((len(stat_word_5), 16), dtype=int)
sw9 = np.zeros((len(stat_word_9), 16), dtype=int)
sw12 = np.zeros((len(stat_word_12), 16), dtype=int)
sw13 = np.zeros((len(stat_word_13), 16), dtype=int)
sw16 = np.zeros((len(stat_word_16), 16), dtype=int)
for i in range(len(stat_word_12)):
    bitstr5 = np.binary_repr(stat_word_5[i], width=16)
    bitstr9 = np.binary_repr(stat_word_9[i], width=16)
    bitstr12 = np.binary_repr(stat_word_12[i], width=16)
    bitstr13 = np.binary_repr(stat_word_13[i], width=16)
    bitstr16 = np.binary_repr(stat_word_16[i], width=16)
    sw5[i, :] = np.array(list(bitstr5), dtype=int)
    sw9[i, :] = np.array(list(bitstr9), dtype=int)
    sw12[i, :] = np.array(list(bitstr12), dtype=int)
    sw13[i, :] = np.array(list(bitstr13), dtype=int)
    sw16[i, :] = np.array(list(bitstr16), dtype=int)

print(f"sw5: {sw5}")

plt.imshow(sw5, aspect="auto", cmap="gray_r", interpolation="none")
plt.colorbar(label="Bit value")
plt.xlabel("Bit number")
plt.ylabel("Interferogram index")
plt.title("Bit flags for stat_word_5")
plt.savefig("./flagging/output/bit_flags_stat_word_5.png")
plt.clf()

# percentage of times each bit is set
fig, ax = plt.subplots(figsize=(10, 6))
ax.bar(range(1, 17), np.sum(sw5, axis=0) / len(stat_word_5) * 100)
ax.set_xticks(range(1, 17))
ax.set_xlabel("Bit number")
ax.set_ylabel("Percentage (%)")
ax.set_title("Bit flag percentages for stat_word_5")
plt.savefig("./flagging/output/bit_flag_percentages_stat_word_5.png")

print("stat_word_5")
for i in range(16):
    print(f"Bit {i+1}: {(np.sum(sw5, axis=0)[i] / len(stat_word_5)) * 100} %")

# indices of the ifgs that have bit 13 set
sw5_bit13_idx = np.where(sw5[:, 12] == 1)[0]

for i in range(0, len(sw5_bit13_idx), 100):
    id = sw5_bit13_idx[i]
    fig, ax = plt.subplots(2, 2, figsize=(10, 6))
    ax[0, 0].plot(ifg_ll[id, :])
    ax[0, 0].set_title("ifg_ll")
    ax[0, 1].plot(ifg_lh[id, :])
    ax[0, 1].set_title("ifg_lh")
    ax[1, 0].plot(ifg_rl[id, :])
    ax[1, 0].set_title("ifg_rl")
    ax[1, 1].plot(ifg_rh[id, :])
    ax[1, 1].set_title("ifg_rh")
    fig.suptitle(
        f"Interferograms with stat_word_5 bit 13 set, index {idx[id]}, mtm_speed {mtm_speed[id]}, mtm_length {mtm_length[id]}"
    )
    plt.tight_layout()
    plt.savefig(f"./flagging/output/ifgs_stat_word_5_bit13_{id}.png")
    plt.clf()

# TODO: fraction analysis
# TODO: correlation between flags
# TODO: add all flags
# TODO: see if some flags are subset of others
