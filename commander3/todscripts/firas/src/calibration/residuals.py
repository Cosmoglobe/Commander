import os
import sys

import h5py

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import globals as g
from simulations.main import generate_ifg

# get temperatures
data = h5py.File(g.PREPROCESSED_DATA_PATH_CAL, "r")
