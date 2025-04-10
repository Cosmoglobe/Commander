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
    pixels = d.shape[0]

    N = np.ndarray((pixels, g.IFG_SIZE, g.IFG_SIZE))
    for pixel in pixels:
        N[pixel] = np.identity(g.IFG_SIZE) # for now we have no noise in the simulation
    # rearrange N to be in the shape (g.IFG_SIZE, G.IFG_SIZE, pixels)
    N = np.transpose(N, (1, 2, 0))

    