"""
This script estimates the noise PSD and covariance matrix by adding the normalized left and right IFGs.
The way I'm doing it now doesn't estimate the white noise level properly because I have to normalize them in order to add them.
"""

import random

import h5py
import matplotlib.pyplot as plt
import numpy as np

import globals as g

cal_data = h5py.File(
    g.PREPROCESSED_DATA_PATH_CAL,
    "r",
)

# TODO: only looking at ss for now, generalize later

mtm_speed = cal_data["df_data/mtm_speed"][:]
mtm_length = cal_data["df_data/mtm_length"][:]

ss_filter = (mtm_speed == 0) & (mtm_length == 0)

ifg_rh = cal_data["df_data/ifg_rh"][ss_filter]
ifg_rl = cal_data["df_data/ifg_rl"][ss_filter]
ifg_lh = cal_data["df_data/ifg_lh"][ss_filter]
ifg_ll = cal_data["df_data/ifg_ll"][ss_filter]
ifg_rh = ifg_rh - np.median(ifg_rh, axis=1, keepdims=True)
ifg_rl = ifg_rl - np.median(ifg_rl, axis=1, keepdims=True)
ifg_lh = ifg_lh - np.median(ifg_lh, axis=1, keepdims=True)
ifg_ll = ifg_ll - np.median(ifg_ll, axis=1, keepdims=True)

gain_rh = cal_data["df_data/gain_rh"][ss_filter]
gain_rl = cal_data["df_data/gain_rl"][ss_filter]
gain_lh = cal_data["df_data/gain_lh"][ss_filter]
gain_ll = cal_data["df_data/gain_ll"][ss_filter]
sweeps = cal_data["df_data/sweeps"][ss_filter]

n = random.randint(0, ifg_rh.shape[0])

print(f"Using IFG number {n} for noise estimation")
print(
    f"Gain RH: {gain_rh[n]}, Gain RL: {gain_rl[n]}, Gain LH: {gain_lh[n]}, Gain LL: {gain_ll[n]}"
)
print(f"Sweeps: {sweeps[n]}")

# i thought this would normalize them "naturally" but it seems both the sides on high and both the sides on low each have the same gain value, so i will need to normalize them separately
ifg_rh = ifg_rh[n] / gain_rh[n] / sweeps[n]
ifg_rl = ifg_rl[n] / gain_rl[n] / sweeps[n]
ifg_lh = ifg_lh[n] / gain_lh[n] / sweeps[n]
ifg_ll = ifg_ll[n] / gain_ll[n] / sweeps[n]

# ifg_rh = np.roll(ifg_rh, -g.PEAK_POSITIONS["rh_ss"], axis=1)
# ifg_rl = np.roll(ifg_rl, -g.PEAK_POSITIONS["rl_ss"], axis=1)
# ifg_lh = np.roll(ifg_lh, -g.PEAK_POSITIONS["lh_ss"], axis=1)
# ifg_ll = np.roll(ifg_ll, -g.PEAK_POSITIONS["ll_ss"], axis=1)

# normalize by amplitude at PEAK_POSITIONS
# ifg_rh_norm = ifg_rh / ifg_rh[g.PEAK_POSITIONS["rh_ss"]]
# ifg_rl_norm = ifg_rl / ifg_rl[g.PEAK_POSITIONS["rl_ss"]]
# ifg_lh_norm = ifg_lh / ifg_lh[g.PEAK_POSITIONS["lh_ss"]]
# ifg_ll_norm = ifg_ll / ifg_ll[g.PEAK_POSITIONS["ll_ss"]]

ifg_low_normr = (
    ifg_rl
    - ifg_ll * ifg_rl[g.PEAK_POSITIONS["rl_ss"]] / ifg_ll[g.PEAK_POSITIONS["ll_ss"]]
)
ifg_low_norml = (
    ifg_ll
    - ifg_rl * ifg_ll[g.PEAK_POSITIONS["ll_ss"]] / ifg_rl[g.PEAK_POSITIONS["rl_ss"]]
)
ifg_high_normr = (
    ifg_rh
    - ifg_lh * ifg_rh[g.PEAK_POSITIONS["rh_ss"]] / ifg_lh[g.PEAK_POSITIONS["lh_ss"]]
)
ifg_high_norml = (
    ifg_lh
    - ifg_rh * ifg_lh[g.PEAK_POSITIONS["lh_ss"]] / ifg_rh[g.PEAK_POSITIONS["rh_ss"]]
)


max_amp_high = np.max([np.abs(ifg_lh).max(), np.abs(ifg_rh).max()])
max_amp_low = np.max([np.abs(ifg_ll).max(), np.abs(ifg_rl).max()])

fig, ax = plt.subplots(4, 2, figsize=(10, 6))

ax[0, 0].plot(ifg_rh)
ax[0, 0].axvline(g.PEAK_POSITIONS["rh_ss"], color="r", linestyle="--")
ax[0, 0].set_title("IFG RH")
ax[0, 0].set_ylim([-1.1 * max_amp_high, 1.1 * max_amp_high])

ax[1, 0].plot(ifg_lh)
ax[1, 0].axvline(g.PEAK_POSITIONS["lh_ss"], color="r", linestyle="--")
ax[1, 0].set_title("IFG LH")
ax[1, 0].set_ylim([-1.1 * max_amp_high, 1.1 * max_amp_high])

ax[2, 0].plot(ifg_high_norml)
ax[2, 0].axvline(g.PEAK_POSITIONS["lh_ss"], color="r", linestyle="--")
ax[2, 0].set_title("IFG High (LH + RH) Normalized to Left")
ax[2, 0].set_ylim([-1.1 * max_amp_high, 1.1 * max_amp_high])

ax[3, 0].plot(ifg_high_normr)
ax[3, 0].axvline(g.PEAK_POSITIONS["rh_ss"], color="r", linestyle="--")
ax[3, 0].set_title("IFG High (LH + RH) Normalized to Right")
ax[3, 0].set_ylim([-1.1 * max_amp_high, 1.1 * max_amp_high])

ax[0, 1].plot(ifg_rl)
ax[0, 1].axvline(g.PEAK_POSITIONS["rl_ss"], color="r", linestyle="--")
ax[0, 1].set_title("IFG RL")
ax[0, 1].set_ylim([-1.1 * max_amp_low, 1.1 * max_amp_low])

ax[1, 1].plot(ifg_ll)
ax[1, 1].axvline(g.PEAK_POSITIONS["ll_ss"], color="r", linestyle="--")
ax[1, 1].set_title("IFG LL")
ax[1, 1].set_ylim([-1.1 * max_amp_low, 1.1 * max_amp_low])

ax[2, 1].plot(ifg_low_norml)
ax[2, 1].axvline(g.PEAK_POSITIONS["ll_ss"], color="r", linestyle="--")
ax[2, 1].set_title("IFG Low (LL + RL) Normalized to Left")
ax[2, 1].set_ylim([-1.1 * max_amp_low, 1.1 * max_amp_low])

ax[3, 1].plot(ifg_low_normr)
ax[3, 1].axvline(g.PEAK_POSITIONS["rl_ss"], color="r", linestyle="--")
ax[3, 1].set_title("IFG Low (LL + RL) Normalized to Right")
ax[3, 1].set_ylim([-1.1 * max_amp_low, 1.1 * max_amp_low])

plt.tight_layout()
plt.show()
