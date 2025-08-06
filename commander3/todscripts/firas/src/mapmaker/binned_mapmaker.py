import os
import sys

import globals_mapmaker as g
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
import my_utils as mu

d = np.load("test_output/ifgs.npz")["ifg"]
pix = np.load("test_output/ifgs.npz")["pix"]

npix = hp.nside2npix(g.NSIDE)
frequencies = mu.generate_frequencies("ll", "ss", 257)

data_density = np.zeros(npix, dtype=int)
m = np.zeros((npix, 257), dtype=complex)
for i in range(d.shape[0]):
    m[pix[i]] += np.abs(np.fft.rfft(d[i]))
    data_density[pix[i]] += 1

mask = data_density == 0
m[~mask] = m[~mask] / data_density[~mask][:, np.newaxis]
m[mask] = np.nan

print("Finished generating map cube, saving to disk...")

# save m as maps
for nui in range(len(frequencies)):
    if g.FITS:
        hp.write_map(f"./test_output/binned_mapmaker/{int(frequencies[nui]):04d}.fits", np.abs(m[:, nui]), overwrite=True)
    if g.PNG:
        hp.mollview(
            m[:, nui],
            title=f"{int(frequencies[nui]):04d} GHz",
            unit="MJy/sr",
            min=0,
            max=200,
            xsize=2000,
        )
        plt.savefig(f"./test_output/binned_mapmaker/{int(frequencies[nui]):04d}.png")
        plt.close()
        plt.clf()