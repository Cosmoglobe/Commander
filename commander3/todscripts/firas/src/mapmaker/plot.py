import os

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from sim import sim_dust


def plot_dust_maps():
    # plot map for each frequency
    for i in range(len(frequencies)):
        print(f"Plotting dust map for frequency {i}")
        dust_map = dust_map_downgraded_mjy * signal[i]
        hp.mollview(dust_map, title=f"{int(frequencies[i]):04d} GHz", unit="MJy/sr", min=0, max=200)
        plt.savefig(f"tests/dust_maps/{int(frequencies[i]):04d}.png")
        plt.close()
        plt.clf()

def plot_m_invert(frequencies):
    # clean previous maps
    for file in os.listdir("tests/m_invert"):
        os.remove(f"tests/m_invert/{file}")

    m = np.load("tests/m_invert.npz")['m']
    frequencies = frequencies[1:]
    for i in range(m.shape[1]):
        # print(f"Plotting m for frequency {i}")
        hp.mollview(m[:, i], title=f"{int(frequencies[i]):04d} GHz", min=0, max=200)
        plt.savefig(f"tests/m_invert/{int(frequencies[i]):04d}.png")
        plt.close()
        plt.clf()

if __name__ == "__main__":
    dust_map_downgraded_mjy, frequencies, signal = sim_dust()
    # plot_dust_maps()
    plot_m_invert(frequencies)