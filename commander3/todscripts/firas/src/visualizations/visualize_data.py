import matplotlib.pyplot as plt
import h5py
import numpy as np

data = h5py.File("../../data/df_v11.h5", "r")

peak_positions = {
    "lh ss": 357,
    "rh ss": 357,
    "lh sf": 359,
    "rh sf": 359,
    "lh lf": 355,
    "rh lf": 355,
    # ----------------
    "ll ss": 360,
    "rl ss": 360,
    "ll fs": 90,
    "rl fs": 90,
    "ll fl": 90,
    "rl fl": 90,
    "ll lf": 90,
    "rl lf": 90,
}

# for i in range(0, len(data["id"]), 100):
for i in range(0, len(data["df_data/gmt"]), 1):
    # if str(data["df_data/gmt"][i])[7:9] == "99":
    speed = data["df_data/mtm_speed"][i]
    speed = "s" if speed == 0 else "f"
    length = data["df_data/mtm_length"][i]
    length = "s" if length == 0 else "l"

    # xcal position
    if data["df_data/xcal_pos"][i] == 1:
        xcal_pos = "in"
    elif data["df_data/xcal_pos"][i] == 2:
        xcal_pos = "out"

    if (
        xcal_pos == "in"
        and speed == "s"
        and length == "s"
        and data["df_data/ical"][i] > 10
    ):  # conditions for plotting
        fig, ax = plt.subplots(sharex=True, nrows=4)
        ax[0].plot(data["df_data/ifg_lh"][i])
        try:
            ax[0].axvline(peak_positions[f"lh {length}{speed}"], color="r", ls="--")
        except KeyError:
            pass
        ax[0].set_ylabel("LH")

        ax[1].plot(data["df_data/ifg_ll"][i])
        try:
            ax[1].axvline(peak_positions[f"ll {speed}{length}"], color="r", ls="--")
        except KeyError:
            pass
        ax[1].set_ylabel("LL")

        ax[2].plot(data["df_data/ifg_rh"][i])
        try:
            ax[2].axvline(peak_positions[f"rh {length}{speed}"], color="r", ls="--")
        except KeyError:
            pass
        ax[2].set_ylabel("RH")

        ax[3].plot(data["df_data/ifg_rl"][i])
        try:
            ax[3].axvline(peak_positions[f"rl {speed}{length}"], color="r", ls="--")
        except KeyError:
            pass
        ax[3].set_ylabel("RL")

        fig.suptitle(
            f"{length.upper()}{speed.upper()}, XCAL {xcal_pos}, ICAL: {data['df_data/ical'][i]:.2f} K, XCAL: {data['df_data/xcal'][i]:.2f} K"
        )
        plt.tight_layout()
        plt.savefig(f"../../plots/{data['df_data/gmt'][i]}.png")
        plt.clf()
        plt.close()
