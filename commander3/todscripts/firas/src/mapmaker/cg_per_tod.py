"""
Script to run a conjugate gradient mapmaker that solves the equation
    A x = b
or more explicitely
    (P^t N^-1 P) S = P^T N^-1 d.
"""

import globals as g
import numpy as np

if __name__ == "__main__":
    # load data and define all matrices
    d = np.load("tests/ifgs.npz")['ifg']
    ntod = d.shape[0]
    noise = np.load("tests/white_noise.npz")['noise']

    W = np.zeros((g.SPEC_SIZE, g.IFG_SIZE), dtype=complex)
    W[0, :] = 1
    W[:, 0] = 1
    omega = np.exp(-2j * np.pi / g.IFG_SIZE)
    for xi in range(1, g.IFG_SIZE):
        for nui in range(1, g.SPEC_SIZE):
            W[nui, xi] = omega ** ((xi * nui) % g.IFG_SIZE) # the mod operator just avoids calculating high exponents
    W = W / np.sqrt(g.IFG_SIZE)

    IW = np.zeros((g.IFG_SIZE, g.SPEC_SIZE), dtype=complex)
    IW[0, :] = 1
    IW[:, 0] = 1
    omega = np.exp(2j * np.pi / g.IFG_SIZE)
    for xi in range(1, g.IFG_SIZE):
        for nui in range(1, g.SPEC_SIZE):
            IW[xi, nui] = omega ** ((xi * nui) % g.IFG_SIZE) # the mod operator just avoids calculating high exponents
    IW = IW / np.sqrt(g.IFG_SIZE)

    for t in range(ntod): # doing it per TOD
        N_inv = np.diag(1/noise[t]**2)