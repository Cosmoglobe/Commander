"""
In this script, we want to solve for m per pixel directly by following the equation:
m = (P^T N^-1 P)^-1 P^T N^-1 d
"""

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from globals import IFG_SIZE, SPEC_SIZE


def solve_m(P, N, d):
    """
    Solve for m per pixel directly by following the equation:
    m = (P^T N^-1 P)^-1 P^T N^-1 d

    Parameters
    ----------
    P : np.ndarray
        The matrix P.
    N : np.ndarray
        The matrix N.
    d : np.ndarray
        The matrix d.

    Returns
    -------
    np.ndarray
        The matrix m.
    """
    N_inv = N # because for now we are assuming N = I
    P_T = np.transpose(P)
    P_T_N_inv = np.dot(P_T, N_inv)
    P_T_N_inv_P = np.dot(P_T_N_inv, P)
    P_T_N_inv_d = np.dot(P_T_N_inv, d)
    P_T_N_inv_P_inv = np.linalg.inv(P_T_N_inv_P)
    m = np.dot(P_T_N_inv_P_inv, P_T_N_inv_d)
    return m

if __name__ == "__main__":
    frequencies = np.linspace(0, 13.604162 * SPEC_SIZE, SPEC_SIZE) # GHz - for SS
    frequencies_icm = (frequencies * u.GHz).to(1 / u.cm, equivalencies=u.spectral()).value # cm-1
    x = np.linspace(
        0, 1.76, IFG_SIZE
    ) # cm

    d = np.load("tests/ifgs.npz")['ifg'][:, 1:]
    # find where nans are in d
    print(f"d: {d}")
    N = np.identity(d.shape[1]) # as we are doing this per pixel, this should be nfreq x nfreq
    P = np.zeros((d.shape[1], d.shape[1]))
    for i in range(d.shape[1]):
        for j in range(d.shape[1]):
            P[i, j] = np.abs(frequencies[i] - frequencies[j]) * np.sum(np.cos(2 * np.pi * frequencies_icm[i] * (x - 1.22))) # for SS - TODO: double check if this is correct

    # check for nans
    print(f"nans in P: {np.isnan(P).any()}")
    print(f"nans in N: {np.isnan(N).any()}")
    print(f"nans in d: {np.isnan(d).any()}")

    m = np.zeros((d.shape[0], d.shape[1]))
    for pixel in range(d.shape[0]):
        print(f"solving for pixel {pixel}")
        m[pixel] = solve_m(P, N, d[pixel])
        # print(f"m at pixel {pixel}: {m[pixel]}")

    # plot the result for a frequency
    hp.mollview(m[:, 63], title="m at frequency 857 GHz") # TODO: units?
    plt.savefig("tests/m.png")
        