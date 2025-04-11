"""
In this script, we want to solve for m per pixel directly by following the equation:
m = (P^T N^-1 P)^-1 P^T N^-1 d
"""

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from globals import IFG_SIZE, SPEC_SIZE


def solve_m(P, P_T, N_inv, d):
    """
    Solve for m per pixel directly by following the equation:
    m = (P^T N^-1 P)^-1 P^T N^-1 d

    Parameters
    ----------
    P : np.ndarray
        The matrix P.
    P_T : np.ndarray
        The hermitian conjugate of the matrix P, because we are dealing with complex numbers. In the case of simply using a Fourier transform for the P operator, this is the inverse of the Fourier transform.
    N : np.ndarray
        The matrix N.
    d : np.ndarray
        The matrix d.

    Returns
    -------
    np.ndarray
        The matrix m.
    """
    # P_T = np.transpose(P)
    P_T_N_inv = np.dot(P_T, N_inv)
    P_T_N_inv_P = np.dot(P_T_N_inv, P)

    # print all sizes
    # print(f"P_T_N_inv_P: {P_T_N_inv_P.shape}")
    # print(f"P_T_N_inv: {P_T_N_inv.shape}")
    # print(f"N_inv: {N_inv.shape}")
    # print(f"P: {P.shape}")
    # print(f"d: {d.shape}")
    
    
    # print matrix to check symmetry
    assert np.allclose(P_T_N_inv_P, np.transpose(P_T_N_inv_P)), "P_T_N_inv_P is not symmetric"

    P_T_N_inv_d = np.dot(P_T_N_inv, d)
    P_T_N_inv_P_inv = np.linalg.inv(P_T_N_inv_P)
    m = np.dot(P_T_N_inv_P_inv, P_T_N_inv_d)

    # print the rest of the shapes
    # print(f"P_T_N_inv_d: {P_T_N_inv_d.shape}")
    # print(f"m: {m.shape}")
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
    ntod = d.shape[0]
    d = np.roll(d, -360, axis=1)

    # find where nans are in d
    # print(f"d: {d}")
    noise = np.load("tests/white_noise.npz")['noise']
    # N_inv = np.identity(IFG_SIZE)
    scale = 0.1
    N_inv = np.identity(IFG_SIZE) / scale **2
    print(f"N_inv: {N_inv.shape}")

    # alternatively for making the discrete fourier transform matrix apparently it should be done like this: https://en.wikipedia.org/wiki/DFT_matrix
    W = np.zeros((SPEC_SIZE, IFG_SIZE), dtype=complex)
    W[0, :] = 1
    W[:, 0] = 1
    omega = np.exp(-2j * np.pi / IFG_SIZE)
    for xi in range(1, IFG_SIZE):
        for nui in range(1, SPEC_SIZE):
            W[nui, xi] = omega ** ((xi * nui) % IFG_SIZE) # the mod operator just avoids calculating high exponents
    W = W / np.sqrt(IFG_SIZE)

    IW = np.zeros((IFG_SIZE, SPEC_SIZE), dtype=complex)
    IW[0, :] = 1
    IW[:, 0] = 1
    omega = np.exp(2j * np.pi / IFG_SIZE)
    for xi in range(1, IFG_SIZE):
        for nui in range(1, SPEC_SIZE):
            IW[xi, nui] = omega ** ((xi * nui) % IFG_SIZE) # the mod operator just avoids calculating high exponents
    IW = IW / np.sqrt(IFG_SIZE)

    # F = np.zeros((SPEC_SIZE, IFG_SIZE), dtype=complex)

    # # unit vector hammering method
    # for i in range(IFG_SIZE):
    #     x = np.zeros(IFG_SIZE)
    #     x[i] = 1
    #     y = np.fft.rfft(x, n = IFG_SIZE)
    #     print(f"y: {y.shape}")
    #     F[:, i] = y

    m = np.zeros((d.shape[0], SPEC_SIZE))
    for t in range(d.shape[0]):
        print(f"solving for ntod {t}")
        # m[pixel] = solve_m(P, N, d[pixel])

        m[t] = solve_m(IW, W, N_inv, d[t])

    # the simplest way to get maps from d = Pm is using P-1d = m
    # m = np.dot(W, d.T).T
        
    # save map output
    np.savez("tests/m_invert.npz", m=m)