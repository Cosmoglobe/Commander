"""
In this script, we want to solve for m per pixel directly by following the equation:
m = (P^T N^-1 P)^-1 P^T N^-1 d
"""

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from globals import IFG_SIZE, SPEC_SIZE


def solve_m(d):
    """
    Solve for m per pixel directly by following the equation:
    m = (P^T N^-1 P)^-1 P^T N^-1 d which simplifies to m = P^T d

    Parameters
    ----------
    N : np.ndarray
        The matrix N.
    d : np.ndarray
        The matrix d.

    Returns
    -------
    np.ndarray
        The matrix m.
    """
    m = np.fft.rfft(d, axis=1)
    
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
    # W = np.zeros((IFG_SIZE, IFG_SIZE), dtype=complex)
    # W[0, :] = 1
    # W[:, 0] = 1
    # omega = np.exp(-2j * np.pi / IFG_SIZE)
    # for xi in range(1, IFG_SIZE):
    #     for nui in range(1, IFG_SIZE):
    #         W[nui, xi] = omega ** ((xi * nui) % IFG_SIZE) # the mod operator just avoids calculating high exponents
    # W = W #/ np.sqrt(IFG_SIZE)

    # IW = np.zeros((IFG_SIZE, IFG_SIZE), dtype=complex)
    # IW[0, :] = 1
    # IW[:, 0] = 1
    # omega = np.exp(2j * np.pi / IFG_SIZE)
    # for xi in range(1, IFG_SIZE):
    #     for nui in range(1, IFG_SIZE):
    #         IW[xi, nui] = omega ** ((xi * nui) % IFG_SIZE) # the mod operator just avoids calculating high exponents
    # IW = IW / IFG_SIZE #IFG_SIZE#np.sqrt(IFG_SIZE)

    # F = np.zeros((SPEC_SIZE, IFG_SIZE), dtype=complex)
    # # unit vector hammering method
    # for i in range(IFG_SIZE):
    #     x = np.zeros(IFG_SIZE)
    #     x[i] = 1
    #     y = np.fft.rfft(x, n = IFG_SIZE)
    #     print(f"y: {y.shape}")
    #     F[:, i] = y

    # IF = np.zeros((IFG_SIZE, SPEC_SIZE), dtype=complex)
    # # unit vector hammering method
    # for i in range(SPEC_SIZE):
    #     x = np.zeros(SPEC_SIZE)
    #     x[i] = 1
    #     y = np.fft.irfft(x, n = IFG_SIZE)
    #     print(f"y: {y.shape}")
    #     IF[:, i] = y

    m = np.zeros((d.shape[0], SPEC_SIZE))
    m = solve_m(d)[:, :SPEC_SIZE]
    # m = np.dot(W, d.T).T[:, :SPEC_SIZE]
        
    # save map output
    np.savez("tests/m_invert.npz", m=m)