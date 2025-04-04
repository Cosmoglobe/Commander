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
    # print(f"P: {P.shape}")
    # print(f"N: {N.shape}")
    N_inv = N # because for now we are assuming N = I
    # print(f"N_inv: {N_inv.shape}")
    P_T = np.transpose(P)
    # print(f"P_T: {P_T.shape}")
    P_T_N_inv = np.dot(P_T, N_inv)
    # print(f"P_T_N_inv: {P_T_N_inv.shape}")
    P_T_N_inv_P = np.dot(P_T_N_inv, P)
    # print(f"P_T_N_inv_P: {P_T_N_inv_P.shape}")
    
    assert np.allclose(P_T_N_inv_P, np.transpose(P_T_N_inv_P)), "P_T_N_inv_P is not symmetric"

    P_T_N_inv_d = np.dot(P_T_N_inv, d)
    P_T_N_inv_P_inv = np.linalg.inv(P_T_N_inv_P)
    m = np.dot(P_T_N_inv_P_inv, P_T_N_inv_d)
    return m

if __name__ == "__main__":
    dnu = 13.604162/2
    frequencies = np.linspace(0, dnu * SPEC_SIZE, SPEC_SIZE) # GHz - for SS
    frequencies_icm = (frequencies * u.GHz).to(1 / u.cm, equivalencies=u.spectral()).value # cm-1
    # print(f"frequencies: {frequencies}")
    # print(f"frequencies_icm: {frequencies_icm}")
    x = np.linspace(
        0, 1.76, IFG_SIZE
    ) # cm

    d = np.load("tests/ifgs.npz")['ifg']
    # find where nans are in d
    # print(f"d: {d}")
    N = np.identity(IFG_SIZE) # as we are doing this per pixel, this should be npoints x npoints
    # P = np.zeros((IFG_SIZE, SPEC_SIZE))
    # for xi in range(IFG_SIZE):
    #     for nui in range(SPEC_SIZE):
    #         P[xi, nui] = dnu * np.cos(2 * np.pi * frequencies_icm[nui] * (x[xi] - 1.22)) # for SS

    # alternatively for making the discrete fourier transform matrix apparently it should be done like this: https://en.wikipedia.org/wiki/DFT_matrix
    # W = np.zeros((IFG_SIZE, SPEC_SIZE), dtype=complex)
    # W[0, :] = 1
    # W[:, 0] = 1
    # omega = np.exp(-2j * np.pi / IFG_SIZE)
    # for xi in range(1, IFG_SIZE):
    #     for nui in range(1, SPEC_SIZE):
    #         W[xi, nui] = omega ** ((xi * nui) % IFG_SIZE) # the mod operator just avoids calculating high exponents
    # W = W / np.sqrt(IFG_SIZE)

    F = np.zeros((SPEC_SIZE, IFG_SIZE), dtype=complex)
    IF = np.zeros((IFG_SIZE, SPEC_SIZE), dtype=complex)

    # unit vector hammering method
    for i in range(IFG_SIZE):
        x = np.zeros(IFG_SIZE)
        x[i] = 1
        y = np.fft.rfft(x, n = IFG_SIZE)
        print(f"y: {y.shape}")
        F[:, i] = y
    
    for i in range(SPEC_SIZE):
        x = np.zeros(SPEC_SIZE)
        x[i] = 1
        y = np.fft.irfft(x, n = IFG_SIZE)
        print(f"y: {y.shape}")
        IF[:, i] = y

    # m = np.zeros((d.shape[0], SPEC_SIZE))
    # for pixel in range(d.shape[0]):
    #     print(f"solving for pixel {pixel}")
    #     # m[pixel] = solve_m(P, N, d[pixel])
    #     m[pixel] = solve_m(F, N, d[pixel])

    # the simplest way to get maps from d = Pm is using P-1d = m
    m = np.dot(F, d.T).T
    print(f"m: {m.shape}")
        
    # save map output
    np.savez("tests/m_invert.npz", m=m)