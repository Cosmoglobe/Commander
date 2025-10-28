import h5py
import matplotlib.pyplot as plt
import numpy as np

import globals as g
from flagging import filter

# modes = {"ss": 0, "lf": 3}
# modes = {"ss": 0}
if __name__ == "__main__":
    mode = "ss"

    cal_data = h5py.File(g.PREPROCESSED_DATA_PATH_CAL, "r")

    mtm_speed = cal_data["df_data/mtm_speed"][:]
    mtm_length = cal_data["df_data/mtm_length"][:]

    ss_filter = (mtm_speed == 0) & (mtm_length == 0)

    stat_word_9 = cal_data["df_data/stat_word_9"][ss_filter]
    lvdt_stat_b = cal_data["df_data/lvdt_stat_b"][ss_filter]

    flags = filter.flag(stat_word_9, lvdt_stat_b)

    a_bol_assem_rh = cal_data[f"df_data/a_bol_assem_rh"][ss_filter]
    b_bol_assem_rh = cal_data[f"df_data/b_bol_assem_rh"][ss_filter]
    a_bol_assem_rl = cal_data[f"df_data/a_bol_assem_rl"][ss_filter]
    b_bol_assem_rl = cal_data[f"df_data/b_bol_assem_rl"][ss_filter]
    a_bol_assem_lh = cal_data[f"df_data/a_bol_assem_lh"][ss_filter]
    b_bol_assem_lh = cal_data[f"df_data/b_bol_assem_lh"][ss_filter]
    a_bol_assem_ll = cal_data[f"df_data/a_bol_assem_ll"][ss_filter]
    b_bol_assem_ll = cal_data[f"df_data/b_bol_assem_ll"][ss_filter]
    bol_cmd_bias_rh = cal_data[f"df_data/bol_cmd_bias_rh"][ss_filter]
    bol_cmd_bias_rl = cal_data[f"df_data/bol_cmd_bias_rl"][ss_filter]
    bol_cmd_bias_lh = cal_data[f"df_data/bol_cmd_bias_lh"][ss_filter]
    bol_cmd_bias_ll = cal_data[f"df_data/bol_cmd_bias_ll"][ss_filter]

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

    ind = 1063

    index = cal_data["df_data/id"][ss_filter]

    idx = filter.flag_bad_ifgs()
    flag_idx = ~np.isin(index, idx)

    index = index[flags & filter_bol & flag_idx]
    stat_word_9 = stat_word_9[flags & filter_bol & flag_idx]
    lvdt_stat_b = lvdt_stat_b[flags & filter_bol & flag_idx]

    ifg_rh = cal_data["df_data/ifg_rh"][ss_filter]
    ifg_rh = ifg_rh[flags & filter_bol & flag_idx]
    ifg_rl = cal_data["df_data/ifg_rl"][ss_filter]
    ifg_rl = ifg_rl[flags & filter_bol & flag_idx]
    ifg_lh = cal_data["df_data/ifg_lh"][ss_filter]
    ifg_lh = ifg_lh[flags & filter_bol & flag_idx]
    ifg_ll = cal_data["df_data/ifg_ll"][ss_filter]
    ifg_ll = ifg_ll[flags & filter_bol & flag_idx]

    print(f"index: {index[ind]}")
    print(f"stat_word_9: {stat_word_9[ind]}")
    print(f"lvdt_stat_b: {lvdt_stat_b[ind]}")

    fig, ax = plt.subplots(2, 2, figsize=(10, 6))
    ax[0, 0].plot(ifg_rh[ind])
    ax[0, 0].set_title("ifg_rh")
    ax[0, 0].axvline(357, color="r", linestyle="--")
    ax[0, 1].plot(ifg_rl[ind])
    ax[0, 1].set_title("ifg_rl")
    ax[0, 1].axvline(360, color="r", linestyle="--")
    ax[1, 0].plot(ifg_lh[ind])
    ax[1, 0].set_title("ifg_lh")
    ax[1, 0].axvline(357, color="r", linestyle="--")
    ax[1, 1].plot(ifg_ll[ind])
    ax[1, 1].set_title("ifg_ll")
    ax[1, 1].axvline(360, color="r", linestyle="--")

    fig.suptitle(f"ID: {index[ind]}")

    plt.tight_layout()
    plt.show()
