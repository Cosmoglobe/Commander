"""
Script to run a conjugate gradient mapmaker that solves the equation
    A x = b
or more explicitely
    (P^t N^-1 P) S = P^T N^-1 d.
"""

import globals as g
import numpy as np


def calculate_A(P, P_T, N_inv):
    """
    Calculate the matrix A = P^T N^-1 P.

    Parameters
    ----------
    P : np.ndarray
        The matrix P.
    P_T : np.ndarray
        The hermitian conjugate of the matrix P, because we are dealing with complex numbers. In the case of simply using a Fourier transform for the P operator, this is the inverse of the Fourier transform.
    N : np.ndarray
        The matrix N.
    
    Returns
    -------
    np.ndarray
        The matrix A.
    """
    P_T_N_inv = np.dot(P_T, N_inv)
    A = np.dot(P_T_N_inv, P)

    return A

def calculate_b(P_T, N_inv, d):
    """
    Calculate the vector b = P^T N^-1 d.

    Parameters
    ----------
    P_T : np.ndarray
        The hermitian conjugate of the matrix P, because we are dealing with complex numbers. In the case of simply using a Fourier transform for the P operator, this is the inverse of the Fourier transform.
    N : np.ndarray
        The matrix N.
    d : np.ndarray
        The matrix d.

    Returns
    -------
    np.ndarray
        The vector b.
    """
    P_T_N_inv = np.dot(P_T, N_inv)
    b = np.dot(P_T_N_inv, d)

    return b

def conjugate_gradient(A, b, x0=None, tol=1e-10, maxiter=1000):
    """
    Solve the equation A x = b using the conjugate gradient method.

    Parameters
    ----------
    A : np.ndarray
        The matrix A.
    b : np.ndarray
        The vector b.
    x0 : np.ndarray, optional
        The initial guess for the solution. If None, a zero vector is used.
    tol : float, optional
        The tolerance for convergence. Default is 1e-10.
    maxiter : int, optional
        The maximum number of iterations. Default is 1000.

    Returns
    -------
    np.ndarray
        The solution vector x.
    """
    if x0 is None:
        x0 = np.zeros_like(b)

    x = x0

    r = b - np.dot(A, x0)
    delta = np.dot(r, r)
    delta0 = delta
    
    for i in range(maxiter):
        Ar = np.dot(A, r)
        alpha = delta / np.dot(r, Ar)

        x += alpha * r
        r -= alpha * Ar

        delta_new = np.dot(r, r)

        if delta_new < tol**2*delta0:
            break

        beta = delta_new / delta
        r = r + beta * r
        delta = delta_new

    return x

if __name__ == "__main__":
    # load data and define all matrices
    d = np.load("tests/ifgs.npz")['ifg']
    d = np.roll(d, -360, axis=1)
    ntod = d.shape[0]
    noise_scale = 0.1

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

    N_inv = np.identity(g.IFG_SIZE) / noise_scale ** 2

    m = np.zeros((ntod, g.SPEC_SIZE), dtype=complex)

    for t in range(ntod): # doing it per TOD
        A = calculate_A(IW, W, N_inv)
        b = calculate_b(W, N_inv, d[t])
        m[t] = conjugate_gradient(A, b)

    # save m in a npz file
    np.savez("tests/cg_per_tod.npz", m=m)