import numpy as np
import matplotlib.pyplot as plt

data = np.load("./data/data_ifgs_mtm.npz")

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

for i in range(0, len(data["id"]), 100):
    speed = data["mtm_speed"][i]
    speed = "s" if speed == 0 else "f"
    length = data["mtm_length"][i]
    length = "s" if length == 0 else "l"

    # for i in range(10):
    fig, ax = plt.subplots(sharex=True, nrows=4)
    ax[0].plot(data["ifg_lh"][i])
    try:
        ax[0].axvline(peak_positions[f"lh {length}{speed}"], color="r", ls="--")
    except KeyError:
        pass
    ax[0].set_ylabel("LH")

    ax[1].plot(data["ifg_ll"][i])
    try:
        ax[1].axvline(peak_positions[f"ll {speed}{length}"], color="r", ls="--")
    except KeyError:
        pass
    ax[1].set_ylabel("LL")

    ax[2].plot(data["ifg_rh"][i])
    try:
        ax[2].axvline(peak_positions[f"rh {length}{speed}"], color="r", ls="--")
    except KeyError:
        pass
    ax[2].set_ylabel("RH")

    ax[3].plot(data["ifg_rl"][i])
    try:
        ax[3].axvline(peak_positions[f"rl {speed}{length}"], color="r", ls="--")
    except KeyError:
        pass
    ax[3].set_ylabel("RL")

    fig.suptitle(
        f"ID: {data['id'][i]}, time: {data['time'][i]}, length: {data['mtm_length'][i]}, speed: {data['mtm_speed'][i]}"
    )
    plt.tight_layout()
    plt.savefig(f"./plots/{data['id'][i]}.png")
    plt.clf()
    plt.close()
