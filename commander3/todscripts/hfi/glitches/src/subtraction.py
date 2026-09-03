"""
In this script, we subtract all of the glitches at the same time using the equation
(T^T N^-1 T) x = T^T N^-1 d
and we solve for x, which is the amplitude of each glitch.
We then subtract the glitches from the data.
We can do this using a CG solver to solve the equation Ax = b
"""

import globals as g
import numpy as np
import templates
from scipy.sparse.linalg import cg


# let's build the T (template) matrix first. it should be a ntod x nglitch matrix
def build_template_matrix(glitch_idx, seconds, glitch_labels, glitch_amps):
    ntod = len(seconds)
    nglitch = len(glitch_idx)
    T = np.zeros((ntod, nglitch))

    oneminute = np.arange(0, g.NSECS, 1 / g.SAMPRATE) 

    print(f"Building template matrix with {nglitch} glitches and {ntod} samples. Shape of "
          f"glitch_idx: {glitch_idx.shape}, glitch_labels: {len(glitch_labels)}, glitch_amps: "
          f"{len(glitch_amps)}")
    for i, glitch_i in enumerate(glitch_idx):
        glitch_model = templates.glitch_model_func(oneminute, glitch_amps[i], band="143-2a",
                                                   glitch_type=glitch_labels[i])
        if glitch_i + len(glitch_model) < ntod:
            T[glitch_i:glitch_i + len(glitch_model), i] = glitch_model
        else:
            T[glitch_i:, i] = glitch_model[:ntod - glitch_i]
    return T

def calculate_b(res, T):
    N_invd = res / g.SIGMA**2
    b = T.T @ N_invd
    return b

def calculate_A(T):
    A = T.T @ T / g.SIGMA**2
    return A

def run_cg(glitch_idx, seconds, glitch_labels, glitch_amps, res):
    T = build_template_matrix(glitch_idx, seconds, glitch_labels, glitch_amps)
    b = calculate_b(res, T)
    A = calculate_A(T)
    x, info = cg(A, b)
    if info != 0:
        print(f"CG solver did not converge. Info: {info}")
    return x

def subtract_glitches_from_data(glitch_idx, seconds, glitch_labels, glitch_amps, res):
    # x = run_cg(glitch_idx, seconds, glitch_labels, glitch_amps, res)
    T = build_template_matrix(glitch_idx, seconds, glitch_labels, glitch_amps)
    b = calculate_b(res, T)
    A = calculate_A(T)

    # brute-force
    x = np.linalg.solve(A, b)

    result = res.copy()
    for (glitch_i, glitch_id), glitch_amp in zip(enumerate(glitch_idx), x):
        glitch_model = templates.glitch_model_func(np.arange(0, g.NSECS, 1 / g.SAMPRATE), glitch_amp,
                                                band="143-2a", glitch_type=glitch_labels[glitch_i])
        if glitch_id + len(glitch_model) < len(res):
            result[glitch_id:glitch_id + len(glitch_model)] -= glitch_model
        else:
            result[glitch_id:] -= glitch_model[:len(res) - glitch_id]
    return result