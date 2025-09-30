"""
This script estimates the noise PSD and covariance matrix by adding the normalized left and right IFGs.
The way I'm doing it now doesn't estimate the white noise level properly because I have to normalize them in order to add them.
"""

import random

import h5py
import matplotlib.pyplot as plt
import numpy as np

import globals as g
import utils.config as config

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

# i thought this would normalize them "naturally" but it seems both the sides on high and both the sides on low each have the same gain value, so i will need to normalize them separately
ifg_rh = ifg_rh / gain_rh[:, np.newaxis] / sweeps[:, np.newaxis]
ifg_rl = ifg_rl / gain_rl[:, np.newaxis] / sweeps[:, np.newaxis]
ifg_lh = ifg_lh / gain_lh[:, np.newaxis] / sweeps[:, np.newaxis]
ifg_ll = ifg_ll / gain_ll[:, np.newaxis] / sweeps[:, np.newaxis]

print(f"Shape of ifg_rh: {ifg_rh.shape}")

# ifg_rh = np.roll(ifg_rh, -g.PEAK_POSITIONS["rh_ss"], axis=1)
# ifg_rl = np.roll(ifg_rl, -g.PEAK_POSITIONS["rl_ss"], axis=1)
# ifg_lh = np.roll(ifg_lh, -g.PEAK_POSITIONS["lh_ss"], axis=1)
# ifg_ll = np.roll(ifg_ll, -g.PEAK_POSITIONS["ll_ss"], axis=1)

# normalize by amplitude at PEAK_POSITIONS
# ifg_rh_norm = ifg_rh / ifg_rh[g.PEAK_POSITIONS["rh_ss"]]
# ifg_rl_norm = ifg_rl / ifg_rl[g.PEAK_POSITIONS["rl_ss"]]
# ifg_lh_norm = ifg_lh / ifg_lh[g.PEAK_POSITIONS["lh_ss"]]
# ifg_ll_norm = ifg_ll / ifg_ll[g.PEAK_POSITIONS["ll_ss"]]

ifg_low_normr = ifg_rl - (
    ifg_ll
    * np.divide(
        ifg_rl[:, g.PEAK_POSITIONS["rl_ss"]],
        ifg_ll[:, g.PEAK_POSITIONS["ll_ss"]],
        out=np.zeros_like(ifg_rl[:, g.PEAK_POSITIONS["rl_ss"]]),
        where=ifg_ll[:, g.PEAK_POSITIONS["ll_ss"]] != 0,
    )[:, np.newaxis]
)
ifg_low_norml = ifg_ll - (
    ifg_rl
    * np.divide(
        ifg_ll[:, g.PEAK_POSITIONS["ll_ss"]],
        ifg_rl[:, g.PEAK_POSITIONS["rl_ss"]],
        out=np.zeros_like(ifg_ll[:, g.PEAK_POSITIONS["ll_ss"]]),
        where=ifg_rl[:, g.PEAK_POSITIONS["rl_ss"]] != 0,
    )[:, np.newaxis]
)
ifg_high_normr = ifg_rh - (
    ifg_lh
    * np.divide(
        ifg_rh[:, g.PEAK_POSITIONS["lh_ss"]],
        ifg_lh[:, g.PEAK_POSITIONS["rh_ss"]],
        out=np.zeros_like(ifg_lh[:, g.PEAK_POSITIONS["lh_ss"]]),
        where=ifg_lh[:, g.PEAK_POSITIONS["rh_ss"]] != 0,
    )[:, np.newaxis]
)
ifg_high_norml = ifg_lh - (
    ifg_rh
    * np.divide(
        ifg_lh[:, g.PEAK_POSITIONS["lh_ss"]],
        ifg_rh[:, g.PEAK_POSITIONS["rh_ss"]],
        out=np.zeros_like(ifg_lh[:, g.PEAK_POSITIONS["lh_ss"]]),
        where=ifg_rh[:, g.PEAK_POSITIONS["rh_ss"]] != 0,
    )[:, np.newaxis]
)

n = random.randint(0, ifg_rh.shape[0])

max_amp_high = np.max([np.abs(ifg_lh[n]), np.abs(ifg_rh[n])])
max_amp_low = np.max([np.abs(ifg_ll[n]), np.abs(ifg_rl[n])])
print(f"Shape of max_amp_high: {max_amp_high.shape}")

fig, ax = plt.subplots(4, 2, figsize=(10, 6))

ax[0, 0].plot(ifg_rh[n])
ax[0, 0].axvline(g.PEAK_POSITIONS["rh_ss"], color="r", linestyle="--")
ax[0, 0].set_title("IFG RH")
ax[0, 0].set_ylim([-1.1 * max_amp_high, 1.1 * max_amp_high])

ax[1, 0].plot(ifg_lh[n])
ax[1, 0].axvline(g.PEAK_POSITIONS["lh_ss"], color="r", linestyle="--")
ax[1, 0].set_title("IFG LH")
ax[1, 0].set_ylim([-1.1 * max_amp_high, 1.1 * max_amp_high])

ax[2, 0].plot(ifg_high_norml[n])
ax[2, 0].axvline(g.PEAK_POSITIONS["lh_ss"], color="r", linestyle="--")
ax[2, 0].set_title("IFG High (LH + RH) Normalized to Left")
ax[2, 0].set_ylim([-1.1 * max_amp_high, 1.1 * max_amp_high])

ax[3, 0].plot(ifg_high_normr[n])
ax[3, 0].axvline(g.PEAK_POSITIONS["rh_ss"], color="r", linestyle="--")
ax[3, 0].set_title("IFG High (LH + RH) Normalized to Right")
ax[3, 0].set_ylim([-1.1 * max_amp_high, 1.1 * max_amp_high])

ax[0, 1].plot(ifg_rl[n])
ax[0, 1].axvline(g.PEAK_POSITIONS["rl_ss"], color="r", linestyle="--")
ax[0, 1].set_title("IFG RL")
ax[0, 1].set_ylim([-1.1 * max_amp_low, 1.1 * max_amp_low])

ax[1, 1].plot(ifg_ll[n])
ax[1, 1].axvline(g.PEAK_POSITIONS["ll_ss"], color="r", linestyle="--")
ax[1, 1].set_title("IFG LL")
ax[1, 1].set_ylim([-1.1 * max_amp_low, 1.1 * max_amp_low])

ax[2, 1].plot(ifg_low_norml[n])
ax[2, 1].axvline(g.PEAK_POSITIONS["ll_ss"], color="r", linestyle="--")
ax[2, 1].set_title("IFG Low (LL + RL) Normalized to Left")
ax[2, 1].set_ylim([-1.1 * max_amp_low, 1.1 * max_amp_low])

ax[3, 1].plot(ifg_low_normr[n])
ax[3, 1].axvline(g.PEAK_POSITIONS["rl_ss"], color="r", linestyle="--")
ax[3, 1].set_title("IFG Low (LL + RL) Normalized to Right")
ax[3, 1].set_ylim([-1.1 * max_amp_low, 1.1 * max_amp_low])

plt.tight_layout()
plt.savefig(f"./noise/output/ifg_left_right_normalized_{n}.png")
plt.close()

cov_rh = np.corrcoef(ifg_high_normr, rowvar=False)
print(cov_rh.shape)
cov_lh = np.corrcoef(ifg_high_norml, rowvar=False)
cov_rl = np.corrcoef(ifg_low_normr, rowvar=False)
cov_ll = np.corrcoef(ifg_low_norml, rowvar=False)

plt.imshow(cov_rh, cmap="RdBu_r", vmax=1, vmin=-1)
plt.colorbar()
plt.title("Correlation Coefficient Matrix of IFG RH")
plt.xlabel("Sample Index")
plt.ylabel("Sample Index")
plt.savefig("./noise/output/cov_rh.png")
plt.close()

plt.imshow(cov_lh, cmap="RdBu_r", vmax=1, vmin=-1)
plt.colorbar()
plt.title("Correlation Coefficient Matrix of IFG LH")
plt.xlabel("Sample Index")
plt.ylabel("Sample Index")
plt.savefig("./noise/output/cov_lh.png")
plt.close()

plt.imshow(cov_rl, cmap="RdBu_r", vmax=1, vmin=-1)
plt.colorbar()
plt.title("Correlation Coefficient Matrix of IFG RL")
plt.xlabel("Sample Index")
plt.ylabel("Sample Index")
plt.savefig("./noise/output/cov_rl.png")
plt.close()

plt.imshow(cov_ll, cmap="RdBu_r", vmax=1, vmin=-1)
plt.colorbar()
plt.title("Correlation Coefficient Matrix of IFG LL")
plt.xlabel("Sample Index")
plt.ylabel("Sample Index")
plt.savefig("./noise/output/cov_ll.png")
plt.close()

# take the psd
psd_rh = np.abs(np.fft.rfft(ifg_high_normr, axis=0)) ** 2
psd_lh = np.abs(np.fft.rfft(ifg_high_norml, axis=0)) ** 2
psd_rl = np.abs(np.fft.rfft(ifg_low_normr, axis=0)) ** 2
psd_ll = np.abs(np.fft.rfft(ifg_low_norml, axis=0)) ** 2

# get nyquist frequency
fnyq_rh = config.gen_nyquistl(
    "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
)["hz"][4 * (g.CHANNELS["rh"] % 2) + g.MODES["ss"]]
fnyq_rl = config.gen_nyquistl(
    "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
)["hz"][4 * (g.CHANNELS["rl"] % 2) + g.MODES["ss"]]
fnyq_lh = config.gen_nyquistl(
    "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
)["hz"][4 * (g.CHANNELS["lh"] % 2) + g.MODES["ss"]]
fnyq_ll = config.gen_nyquistl(
    "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
)["hz"][4 * (g.CHANNELS["ll"] % 2) + g.MODES["ss"]]

# get the frequencies
freqs_rh = np.fft.rfftfreq(ifg_rh.shape[0], d=2 * fnyq_rh)
freqs_rl = np.fft.rfftfreq(ifg_rl.shape[0], d=2 * fnyq_rl)
freqs_lh = np.fft.rfftfreq(ifg_lh.shape[0], d=2 * fnyq_lh)
freqs_ll = np.fft.rfftfreq(ifg_ll.shape[0], d=2 * fnyq_ll)

plt.plot(freqs_rh, psd_rh)
plt.xscale("log")
plt.yscale("log")
plt.title("PSD of IFG RH")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power Spectral Density")
plt.show()
