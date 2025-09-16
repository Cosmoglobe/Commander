import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

import globals as g

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

# ifg_low = ifg_ll + ifg_rl
# ifg_high = ifg_lh + ifg_rh
ifg_low_right = (
    ifg_rl
    + ifg_ll
    * np.max(np.abs(ifg_rl), axis=1)[:, np.newaxis]
    / np.max(np.abs(ifg_ll), axis=1)[:, np.newaxis]
)
ifg_low_left = (
    ifg_ll
    + ifg_rl
    * np.max(np.abs(ifg_ll), axis=1)[:, np.newaxis]
    / np.max(np.abs(ifg_rl), axis=1)[:, np.newaxis]
)
ifg_high_right = (
    ifg_rh
    + ifg_lh
    * np.max(np.abs(ifg_rh), axis=1)[:, np.newaxis]
    / np.max(np.abs(ifg_lh), axis=1)[:, np.newaxis]
)
ifg_high_left = (
    ifg_lh
    + ifg_rh
    * np.max(np.abs(ifg_lh), axis=1)[:, np.newaxis]
    / np.max(np.abs(ifg_rh), axis=1)[:, np.newaxis]
)

std_low_right = np.std(ifg_low_right, axis=1)
std_low_left = np.std(ifg_low_left, axis=1)
std_high_right = np.std(ifg_high_right, axis=1)
std_high_left = np.std(ifg_high_left, axis=1)

bins = np.linspace(0, 1000, 100)
# bins[-1] = np.inf

# let's do the std of the stds and cut at 5 sigma
mean_std_low_right = np.mean(std_low_right)
mean_std_low_left = np.mean(std_low_left)
mean_std_high_right = np.mean(std_high_right)
mean_std_high_left = np.mean(std_high_left)

(mu, sigma) = norm.fit(std_low_right)

plt.hist(
    std_low_right,
    bins=bins,
    density=True,
)
y = norm.pdf(bins, mu, sigma)
plt.plot(bins, y, "r--")
plt.vlines(mean_std_low_right, 0, 0.03, colors="g", label="Mean of stds")
plt.vlines(sigma, 0, 0.03, colors="y", label="Sigma fit")
plt.xlabel("Standard Deviation of Interferogram")
plt.ylabel("Number of Interferograms")
plt.title("Distribution of Interferogram Standard Deviations - Low Frequency Right")
plt.legend()
plt.savefig("./flagging/output/std_low_right_hist.png")
plt.close()

(mu, sigma) = norm.fit(std_low_left)

plt.hist(std_low_left, bins=bins, density=True)
y = norm.pdf(bins, mu, sigma)
plt.vlines(mean_std_low_left, 0, 0.0175, colors="g", label="Mean of stds")
plt.vlines(sigma, 0, 0.0175, colors="y", label="Sigma fit")
plt.plot(bins, y, "r--")
plt.xlabel("Standard Deviation of Interferogram")
plt.ylabel("Number of Interferograms")
plt.title("Distribution of Interferogram Standard Deviations - Low Frequency Left")
plt.legend()
plt.savefig("./flagging/output/std_low_left_hist.png")
plt.close()

(mu, sigma) = norm.fit(std_high_right)
plt.hist(std_high_right, bins=bins, density=True)
y = norm.pdf(bins, mu, sigma)
plt.vlines(mean_std_high_right, 0, 0.03, colors="g", label="Mean of stds")
plt.vlines(sigma, 0, 0.03, colors="y", label="Sigma fit")
plt.plot(bins, y, "r--")
plt.xlabel("Standard Deviation of Interferogram")
plt.ylabel("Number of Interferograms")
plt.title("Distribution of Interferogram Standard Deviations - High Frequency Right")
plt.legend()
plt.savefig("./flagging/output/std_high_right_hist.png")
plt.close()

(mu, sigma) = norm.fit(std_high_left)
plt.hist(std_high_left, bins=bins, density=True)
y = norm.pdf(bins, mu, sigma)
plt.vlines(mean_std_high_left, 0, 0.0175, colors="g", label="Mean of stds")
plt.vlines(sigma, 0, 0.0175, colors="y", label="Sigma fit")
plt.plot(bins, y, "r--")
plt.xlabel("Standard Deviation of Interferogram")
plt.ylabel("Number of Interferograms")
plt.title("Distribution of Interferogram Standard Deviations - High Frequency Left")
plt.legend()
plt.savefig("./flagging/output/std_high_left_hist.png")
plt.close()
