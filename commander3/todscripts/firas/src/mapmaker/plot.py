import os

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from globals_mapmaker import FITS, PNG
from sim import sim_dust


def plot_ifgs(ifg):
    # clean previous maps
    for file in os.listdir("./test_output/ifgs"):
        os.remove(f"./test_output/ifgs/{file}")

    for i in range(0, ifg.shape[0], 1000):
        # print(f"Plotting ifg {i}: {ifg[i]}")
        plt.plot(ifg[i])
        # plt.ylim(-300, 300)
        plt.savefig(f"./test_output/ifgs/{i}.png")
        plt.close()
        plt.clf()

def plot_dust_maps(dust_map_downgraded_mjy, frequencies, signal):
    # clean previous maps
    for file in os.listdir("./test_output/dust_maps"):
        os.remove(f"./test_output/dust_maps/{file}")
    # plot map for each frequency
    dust_map = dust_map_downgraded_mjy[:, np.newaxis] * signal[np.newaxis, :]
    for i in range(len(frequencies)):
        print(f"Plotting dust map for frequency {i}")
        # dust_map = dust_map_downgraded_mjy * signal[i]
        if PNG:
            hp.mollview(dust_map[:, i], title=f"{int(frequencies.value[i]):04d} GHz", unit="MJy/sr", min=0, max=200)
            plt.savefig(f"./test_output/dust_maps/{int(frequencies.value[i]):04d}.png")
            plt.close()
            plt.clf()
        if FITS:
            hp.write_map(f"./test_output/dust_maps/{int(frequencies.value[i]):04d}.fits", dust_map[:, i], overwrite=True)


def plot_m_invert(frequencies):
    # clean previous maps
    print("Cleaning previous maps")
    for file in os.listdir("./test_output/m_invert"):
        os.remove(f"./test_output/m_invert/{file}")

    m = np.load("./test_output/m_invert.npz")['m']
    # remove monopole

    if PNG:
        for i in range(m.shape[1]):
            # print(f"Plotting m for frequency {i}")
            hp.mollview(m[:, i].real, title=f"{int(frequencies.value[i]):04d} GHz", min=0, max=200, xsize=2000)
            plt.savefig(f"./test_output/m_invert/{int(frequencies.value[i]):04d}.png")
            plt.close()
            plt.clf()
    if FITS:
        for i in range(m.shape[1]):
            # print(f"Plotting m for frequency {i}")
            hp.write_map(f"./test_output/m_invert/{int(frequencies.value[i]):04d}.fits", m[:, i].real, overwrite=True)
            plt.close()
            plt.clf()

def plot_m_cg_per_tod(frequencies):
    # clean previous maps
    print("Cleaning previous maps")
    for file in os.listdir("./test_output/m_cg_per_tod"):
        os.remove(f"./test_output/m_cg_per_tod/{file}")

    m = np.load("./test_output/cg_per_tod.npz")['m']

    if PNG:
        for i in range(m.shape[1]):
            # print(f"Plotting m for frequency {i}")
            hp.mollview(m[:, i].real, title=f"{int(frequencies.value[i]):04d} GHz", min=0, max=200, xsize=2000)
            plt.savefig(f"./test_output/m_cg_per_tod/{int(frequencies.value[i]):04d}.png")
            plt.close()
            plt.clf()
    if FITS:
        for i in range(m.shape[1]):
            # print(f"Plotting m for frequency {i}")
            hp.write_map(f"./test_output/m_cg_per_tod/{int(frequencies.value[i]):04d}.fits", m[:, i].real, overwrite=True)
            plt.close()
            plt.clf()

if __name__ == "__main__":
    # open ifgs
    ifg = np.load("test_output/ifgs.npz")['ifg']
    plot_ifgs(ifg)

    dust_map_downgraded_mjy, frequencies, signal = sim_dust()

    # # # plot_dust_maps(dust_map_downgraded_mjy, frequencies, signal)
    # plot_m_invert(frequencies)
    # plot_m_cg_per_tod(frequencies)