import os
import sys

import numpy as np

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import globals as g

sky_data = np.load(g.PROCESSED_DATA_PATH_CAL)
print(sky_data.files)

ifgs_ll = sky_data["ifg_ll_ss"]

temp_xcal = sky_data["df_data"]