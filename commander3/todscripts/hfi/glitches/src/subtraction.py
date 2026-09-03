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
    for i, glitch_i in enumerate(glitch_idx):
        glitch_model = templates.glitch_model_func(oneminute, glitch_amps[i], band="143-2a",
                                                   glitch_type=glitch_labels[i])
        if glitch_i + len(glitch_model) < ntod:
            T[glitch_i:glitch_i + len(glitch_model), i] = glitch_model
        else:
            T[glitch_i:, i] = glitch_model[:ntod - glitch_i]
    return T

