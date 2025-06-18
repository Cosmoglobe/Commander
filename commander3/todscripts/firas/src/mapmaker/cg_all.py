"""
Conjugate gradient mapmaker that solves the equation
    A x = b
or more explicitely
    (P^T M^T N^{-1} M P) m = P^T M^T N^{-1} d.
"""

import healpy as hp
import numpy as np
import scipy
import scipy.sparse

import globals as g


def calculate_b(P_T, N_inv, d):
    return d

if __name__ == "__main__":
    d = np.load("tests/ifgs.npz")["ifg"]
    d = np.roll(d, -360, axis=1)

    ntod = d.shape[0]

    d = np.reshape(d, (ntod * g.IFG_SIZE))
    print(f"d shape: {d.shape}")

    noise_scale = 0.1
    # N_inv = np.identity(g.IFG_SIZE * ntod) / noise_scale**2 # is this correct?
    N_inv = scipy.sparse.diags(
        np.ones(g.IFG_SIZE * ntod) / noise_scale**2, format="csr"
    )  # sparse matrix
    
    P = np.zeros((g.SPEC_SIZE, g.SPEC_SIZE*g.NPIX))

    m = np.zeros(g.SPEC_SIZE * ntod, dtype=complex)

    m = conjugate_gradient(N_inv, b)