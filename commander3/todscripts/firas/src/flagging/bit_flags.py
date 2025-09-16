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


def get_bin_matrix(stat_word_number):
    status_word = cal_data[f"df_data/stat_word_{stat_word_number}"][:].astype(int)
    sw = np.zeros((len(status_word), 16), dtype=int)
    for i in range(len(status_word)):
        bitstr = np.binary_repr(status_word[i], width=16)
        sw[i, :] = np.array(list(bitstr), dtype=int)[::-1]
    return sw


def print_percentages(sw, num):
    print(f"status_word_{num}")
    for i in range(16):
        print(f"Bit {i+1}: {(np.sum(sw, axis=0)[i] / len(sw)) * 100} %")


def plot_matrix(sw, num):
    plt.imshow(sw, aspect="auto", cmap="gray_r", interpolation="none")
    plt.colorbar(label="Bit value")
    plt.xlabel("Bit number")
    plt.ylabel("Interferogram index")
    plt.title(f"Bit flags for stat_word_{num}")
    plt.savefig(f"./flagging/output/bit_flags_stat_word_{num}.png")
    plt.clf()


def plot_hist(sw, num):
    # percentage of times each bit is set
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(range(1, 17), np.sum(sw, axis=0) / len(sw) * 100)
    ax.set_xticks(range(1, 17))
    ax.set_xlabel("Bit number")
    ax.set_ylabel("Percentage (%)")
    ax.set_title(f"Bit flag percentages for stat_word_{num}")
    plt.savefig(f"./flagging/output/bit_flag_percentages_stat_word_{num}.png")
    plt.clf()


def plot_ifgs(sw, sw_num, bit_num, set_to=1):
    # indices of the ifgs that have the specified bit set
    idx_bits = np.where(sw[:, bit_num - 1] == set_to)[0]

    num_idx = len(idx_bits)

    for i in range(0, num_idx, num_idx // 10):
        id = idx_bits[i]

        mode = ""
        mode = "s" if mtm_length[id] == 0 else "l"
        mode = mode + ("s" if mtm_speed[id] == 0 else "f")

        highfreq_lim = np.max([np.abs(ifg_rh[id, :]), np.abs(ifg_lh[id, :])])
        lowfreq_lim = np.max([np.abs(ifg_ll[id, :]), np.abs(ifg_rl[id, :])])

        fig, ax = plt.subplots(3, 2, figsize=(10, 6))
        ax[0, 0].plot(ifg_rh[id, :])
        ax[0, 0].set_title("ifg_rh")
        ax[0, 0].axvline(x=g.PEAK_POSITIONS[f"rh_{mode}"], color="r", linestyle="--")
        ax[0, 0].set_ylim(-1.1 * highfreq_lim, 1.1 * highfreq_lim)
        ax[1, 0].plot(ifg_lh[id, :])
        ax[1, 0].set_title("ifg_lh")
        ax[1, 0].axvline(x=g.PEAK_POSITIONS[f"lh_{mode}"], color="r", linestyle="--")
        ax[1, 0].set_ylim(-1.1 * highfreq_lim, 1.1 * highfreq_lim)
        ax[2, 0].plot(
            ifg_lh[id, :]
            + ifg_rh[id, :]
            * np.max(np.abs(ifg_lh), axis=1)[:, np.newaxis]
            / np.max(np.abs(ifg_rh), axis=1)[:, np.newaxis]
        )
        ax[2, 0].plot(
            ifg_rh[id, :]
            + ifg_lh[id, :]
            * np.max(np.abs(ifg_rh), axis=1)[:, np.newaxis]
            / np.max(np.abs(ifg_lh), axis=1)[:, np.newaxis]
        )
        ax[2, 0].set_title("ifg_lh + ifg_rh")
        ax[2, 0].axvline(x=g.PEAK_POSITIONS[f"lh_{mode}"], color="r", linestyle="--")
        ax[2, 0].set_ylim(-1.1 * highfreq_lim, 1.1 * highfreq_lim)
        ax[0, 1].plot(ifg_rl[id, :])
        ax[0, 1].set_title("ifg_rl")
        ax[0, 1].axvline(x=g.PEAK_POSITIONS[f"rl_{mode}"], color="r", linestyle="--")
        ax[0, 1].set_ylim(-1.1 * lowfreq_lim, 1.1 * lowfreq_lim)
        ax[1, 1].plot(ifg_ll[id, :])
        ax[1, 1].set_title("ifg_ll")
        ax[1, 1].axvline(x=g.PEAK_POSITIONS[f"ll_{mode}"], color="r", linestyle="--")
        ax[1, 1].set_ylim(-1.1 * lowfreq_lim, 1.1 * lowfreq_lim)
        ax[2, 1].plot(
            ifg_ll[id, :]
            + ifg_rl[id, :]
            * np.max(np.abs(ifg_ll), axis=1)[:, np.newaxis]
            / np.max(np.abs(ifg_rl), axis=1)[:, np.newaxis]
        )
        ax[2, 1].plot(
            ifg_rl[id, :]
            + ifg_ll[id, :]
            * np.max(np.abs(ifg_rl), axis=1)[:, np.newaxis]
            / np.max(np.abs(ifg_ll), axis=1)[:, np.newaxis]
        )
        ax[2, 1].set_title("ifg_ll + ifg_rl")
        ax[2, 1].axvline(x=g.PEAK_POSITIONS[f"ll_{mode}"], color="r", linestyle="--")
        ax[2, 1].set_ylim(-1.1 * lowfreq_lim, 1.1 * lowfreq_lim)
        plt.subplots_adjust(hspace=0.5)
        fig.suptitle(
            f"Interferograms with stat_word_{sw_num} bit {bit_num} set to {set_to}\nIndex {idx[id]}\nMode {mode.upper()}\nICAL: {ical[id]:.2f}, XCAL: {xcal[id]:.2f}",
        )
        plt.tight_layout()
        plt.savefig(
            f"./flagging/output/ifgs_stat_word_{sw_num}/bit{bit_num}_{idx[id]}.png"
        )
        plt.close(fig)


cal_data = h5py.File(
    g.PREPROCESSED_DATA_PATH_CAL,
    "r",
)

ifg_ll = cal_data["df_data/ifg_ll"][:]
ifg_lh = cal_data["df_data/ifg_lh"][:]
ifg_rl = cal_data["df_data/ifg_rl"][:]
ifg_rh = cal_data["df_data/ifg_rh"][:]
ifg_ll = ifg_ll - np.median(ifg_ll, axis=1)[:, None]
ifg_lh = ifg_lh - np.median(ifg_lh, axis=1)[:, None]
ifg_rl = ifg_rl - np.median(ifg_rl, axis=1)[:, None]
ifg_rh = ifg_rh - np.median(ifg_rh, axis=1)[:, None]

mtm_speed = cal_data["df_data/mtm_speed"][:]
mtm_length = cal_data["df_data/mtm_length"][:]

a_xcal = cal_data["df_data/a_xcal"][:]
b_xcal = cal_data["df_data/b_xcal"][:]
xcal = (a_xcal + b_xcal) / 2
a_ical = cal_data["df_data/a_ical"][:]
b_ical = cal_data["df_data/b_ical"][:]
ical = (a_ical + b_ical) / 2

idx = cal_data["df_data/index"][:]

stat_word_nums = [5, 9, 12, 13, 16]
for num in stat_word_nums:
    sw = get_bin_matrix(num)
    print_percentages(sw, num)
    plot_matrix(sw, num)
    plot_hist(sw, num)

sw5 = get_bin_matrix(5)
plot_ifgs(sw5, 5, 4, set_to=1)
plot_ifgs(sw5, 5, 7, set_to=1)
plot_ifgs(sw5, 5, 8, set_to=0)
plot_ifgs(sw5, 5, 9, set_to=1)
plot_ifgs(sw5, 5, 10, set_to=1)
plot_ifgs(sw5, 5, 11, set_to=0)
plot_ifgs(sw5, 5, 12, set_to=0)

# sw9 = get_bin_matrix(9)
# plot_ifgs(sw9, 9, 16, set_to=0)
# plot_ifgs(sw9, 9, 13, set_to=0)
# plot_ifgs(sw9, 9, 14, set_to=0)

# sw12 = get_bin_matrix(12)
# plot_ifgs(sw12, 12, 1, set_to=1)
# plot_ifgs(sw12, 12, 3, set_to=1)
# plot_ifgs(sw12, 12, 4, set_to=0)
# plot_ifgs(sw12, 12, 6, set_to=0)
# plot_ifgs(sw12, 12, 8, set_to=1)
# plot_ifgs(sw12, 12, 9, set_to=1)
# plot_ifgs(sw12, 12, 14, set_to=1)

# TODO: fraction analysis
# TODO: correlation between flags
# TODO: add all flags
# TODO: see if some flags are subset of others
