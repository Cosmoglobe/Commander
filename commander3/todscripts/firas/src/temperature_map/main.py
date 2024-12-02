"""
This script takes the interferograms into spectra, fits for the monopole temperature
and makes the temperature map.
"""

import h5py
import numpy as np

from my_utils import clean_ifg


data = h5py.File(
    "/mn/stornext/u3/duncanwa/Commander_debugging/commander3/todscripts/firas/data/df_v11.h5",
    "r",
)

# indices of the data where xcal is in
idx_xcal_in = np.where(data["df_data/xcal_pos"][()] == 1)

# getting the xcal data into their own variables
ical_xcal_in = data["df_data/ical"][idx_xcal_in]
xcal = data["df_data/xcal"][idx_xcal_in]
ifg_ll = data["df_data/ifg_ll"][idx_xcal_in]
ifg_lh = data["df_data/ifg_lh"][idx_xcal_in]
ifg_rl = data["df_data/ifg_rl"][idx_xcal_in]
ifg_rh = data["df_data/ifg_rh"][idx_xcal_in]
mtm_length = data["df_data/mtm_length"][idx_xcal_in]
mtm_speed = data["df_data/mtm_speed"][idx_xcal_in]

data.close()

# find the record with the smallest difference between the two temperatures
least_diff_idx = np.argmin(np.abs(ical_xcal_in - xcal))

clean_ifg(
    ifg_ll[least_diff_idx],
    mtm_length=mtm_length[least_diff_idx],
    mtm_speed=mtm_speed[least_diff_idx],
    channel=3,
    adds_per_group=1,
)
