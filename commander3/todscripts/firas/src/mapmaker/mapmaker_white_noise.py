"""
Conjugate gradient mapmaker that solves the equation
    A x = b
or more explicitely
    (P^T M^T N^{-1} M P) m = P^T M^T N^{-1} d.
"""

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

if __name__ == "__main__":
    d = np.load("test_output/ifgs.npz")["ifg"]
    d = np.roll(d, -360, axis=1)

    ntod = d.shape[0]
    sigma = np.load("test_output/white_noise.npz")["noise"]

    npix = hp.nside2npix(g.NSIDE)
    # P = np.linspace(0, npix, npix, endpoint=False, dtype=int)
    P = np.load("test_output/ifgs.npz")["pix"]

    numerator = np.zeros((npix, d.shape[1]), dtype=float)
    for xi in range(d.shape[1]):
        numerator[:, xi] = np.bincount(P, weights=d[:, xi]/sigma**2, minlength=npix)
    denominator = np.bincount(P, weights=1/sigma**2, minlength=npix)
    m = numerator / denominator[:, np.newaxis]
    m = np.fft.rfft(m, axis=1)

    frequencies = mu.generate_frequencies("ll", "ss", 257)
    # save m as maps
    for nui in range(len(frequencies)):
        if g.FITS:
            hp.write_map(f"./test_output/mapmaker_white_noise/{int(frequencies[nui]):04d}.fits", np.abs(m[:, nui]), overwrite=True)
        if g.PNG:
            hp.mollview(
                np.abs(m[:, nui]),
                title=f"{int(frequencies[nui]):04d} GHz",
                unit="MJy/sr",
                min=0,
                max=200,
                xsize=2000,
                # norm='hist',
            )
            plt.savefig(f"./test_output/mapmaker_white_noise/{int(frequencies[nui]):04d}.png")
            plt.close()
            plt.clf()