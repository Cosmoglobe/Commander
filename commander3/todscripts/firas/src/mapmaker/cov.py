import os
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import globals as g
from utils.config import gen_nyquistl

sky_data = h5py.File(
    g.PREPROCESSED_DATA_PATH_SKY,
    "r",
)

print(sky_data["df_data/ifg_ll"].shape)

ifgs_ll = sky_data["df_data/ifg_ll"]
ifgs_rl = sky_data["df_data/ifg_rl"]

ifgs_ll = ifgs_ll - np.median(ifgs_ll, axis=1, keepdims=True)
ifgs_rl = ifgs_rl - np.median(ifgs_rl, axis=1, keepdims=True)

# plt.plot(ifgs_ll[0])
# plt.plot(ifgs_rl[0])
# plt.title("IFGs LL and RL")
# plt.xlabel("Sample Index")
# plt.ylabel("Amplitude")
# plt.legend(["LL", "RL"])
# plt.show()

ifgs_sub = ifgs_ll[:-1] - ifgs_ll[1:]
# ifgs_sub = ifgs_ll + ifgs_rl
# TODO: subtract the whole model

cov = np.corrcoef(ifgs_sub, rowvar=False)
print(cov.shape)

print(cov)

plt.imshow(cov, cmap="RdBu_r", vmax=1, vmin=-1)
plt.colorbar()
plt.title("Correlation Coefficient Matrix of IFGs")
plt.xlabel("IFG Index")
plt.ylabel("IFG Index")
plt.show()

channel_value = 3
mode_value = 0

frec = 4 * (channel_value % 2) + mode_value

fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)["hz"][frec]

psd = np.abs(np.fft.rfft(ifgs_sub, axis=0))**2
freqs = np.fft.rfftfreq(ifgs_sub.shape[0], d=2*fnyq)

print(psd.shape)

plt.plot(freqs, psd)
plt.xscale("log")
plt.yscale("log")
plt.title("Power Spectral Density of IFGs")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power Spectral Density")
plt.savefig("./tests/psd.png")