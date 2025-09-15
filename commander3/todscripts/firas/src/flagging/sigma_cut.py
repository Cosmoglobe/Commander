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

ifg_low = ifg_ll + ifg_rl
ifg_high = ifg_lh + ifg_rh

std_low = np.std(ifg_low, axis=1)
std_high = np.std(ifg_high, axis=1)

bins = np.linspace(0, 1000, 100)
# bins[-1] = np.inf

# let's do the std of the stds and cut at 5 sigma
mean_std_low = np.mean(std_low)
mean_std_high = np.mean(std_high)

(mu, sigma) = norm.fit(std_low)

plt.hist(
    std_low,
    bins=bins,
    density=True,
)
y = norm.pdf(bins, mu, sigma)
plt.plot(bins, y, "r--")
plt.vlines(mean_std_low, 0, 0.03, colors="g", label="Mean of stds")
plt.vlines(sigma, 0, 0.03, colors="y", label="Sigma fit")
plt.xlabel("Standard Deviation of Interferogram")
plt.ylabel("Number of Interferograms")
plt.title("Distribution of Interferogram Standard Deviations - Low Frequency")
plt.legend()
plt.savefig("./flagging/output/std_low_hist.png")
plt.close()

(mu, sigma) = norm.fit(std_low)

plt.hist(std_high, bins=bins, density=True)
y = norm.pdf(bins, mu, sigma)
plt.vlines(mean_std_high, 0, 0.0175, colors="g", label="Mean of stds")
plt.vlines(sigma, 0, 0.0175, colors="y", label="Sigma fit")
plt.plot(bins, y, "r--")
plt.xlabel("Standard Deviation of Interferogram")
plt.ylabel("Number of Interferograms")
plt.title("Distribution of Interferogram Standard Deviations - High Frequency")
plt.legend()
plt.savefig("./flagging/output/std_high_hist.png")
plt.close()
