"""
So for this version we will keep each channel separate and use the midpoint_time to interpolate to get the temperatures.
"""

import h5py
import numpy as np

import data.utils.my_utils as data_utils
import globals as g

# OPENING ORIGINAL DATA FILES
fdq_sdf = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5")
fdq_eng = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5")

# Longitudes and latitudes are stored in radians*1e4
fact = 180.0 / np.pi / 1e4

ct_head = fdq_eng["ct_head"]
eng_time = data_utils.binary_to_gmt(ct_head["time"][:])

en_stat = fdq_eng["en_stat"]
# side a
stat_word_1 = en_stat["stat_word_1"][:]
stat_word_4 = en_stat["stat_word_4"][:]
stat_word_5 = en_stat["stat_word_5"][:]
stat_word_8 = en_stat["stat_word_8"][:]
# side b
stat_word_9 = en_stat["stat_word_9"][:]
stat_word_12 = en_stat["stat_word_12"][:]
stat_word_13 = en_stat["stat_word_13"][:]
stat_word_16 = en_stat["stat_word_16"][:]

en_analog = fdq_eng["en_analog"]
grt = en_analog["grt"]
a_lo_xcal_tip = grt["a_lo_xcal_tip"][:]
a_lo_skyhorn = grt["a_lo_skyhorn"][:]
a_lo_refhorn = grt["a_lo_refhorn"][:]
a_lo_ical = grt["a_lo_ical"][:]
a_lo_dihedral = grt["a_lo_dihedral"][:]
a_lo_mirror = grt["a_lo_mirror"][:]
a_lo_xcal_cone = grt["a_lo_xcal_cone"][:]
a_lo_collimator = grt["a_lo_collimator"][:]
a_hi_xcal_tip = grt["a_hi_xcal_tip"][:]
a_hi_skyhorn = grt["a_hi_skyhorn"][:]
a_hi_refhorn = grt["a_hi_refhorn"][:]
a_hi_ical = grt["a_hi_ical"][:]
a_hi_dihedral = grt["a_hi_dihedral"][:]
a_hi_mirror = grt["a_hi_mirror"][:]
a_hi_xcal_cone = grt["a_hi_xcal_cone"][:]
a_hi_collimator = grt["a_hi_collimator"][:]
b_lo_xcal_tip = grt["b_lo_xcal_tip"][:]
b_lo_skyhorn = grt["b_lo_skyhorn"][:]
b_lo_refhorn = grt["b_lo_refhorn"][:]
b_lo_ical = grt["b_lo_ical"][:]
b_lo_dihedral = grt["b_lo_dihedral"][:]
b_lo_mirror = grt["b_lo_mirror"][:]
b_lo_xcal_cone = grt["b_lo_xcal_cone"][:]
b_lo_collimator = grt["b_lo_collimator"][:]
b_hi_xcal_tip = grt["b_hi_xcal_tip"][:]
b_hi_skyhorn = grt["b_hi_skyhorn"][:]
b_hi_refhorn = grt["b_hi_refhorn"][:]
b_hi_ical = grt["b_hi_ical"][:]
b_hi_dihedral = grt["b_hi_dihedral"][:]
b_hi_mirror = grt["b_hi_mirror"][:]
b_hi_xcal_cone = grt["b_hi_xcal_cone"][:]
b_hi_collimator = grt["b_hi_collimator"][:]

group1 = en_analog["group1"]

for channel, channel_i in g.CHANNELS.items():
    data = fdq_sdf[f"fdq_sdf_{channel}"]

    sci_head = data["sci_head"]
    gain = sci_head["gain"][:]
    mtm_speed = sci_head["mtm_speed"][:]
    mtm_length = sci_head["mtm_length"][:]
    upmode = sci_head["sc_head1a"][:]
    adds_per_group = sci_head["sc_head9"][:]
    sweeps = sci_head["sc_head11"][:]

    ifg_data = data["ifg_data"]
    ifg = ifg_data["ifg"][:]

    dq_data = data["dq_data"]
    fake = dq_data["fake"][:]
    xcal_pos = dq_data["xcal_pos"][:]

    collect_time = data["collect_time"]
    midpoint_time = data_utils.binary_to_gmt(collect_time["midpoint_time"][:])

    attitude = data["attitude"]
    cel_lon = attitude["ra"][:] * fact
    cel_lat = attitude["dec"][:] * fact
    terr_lat = attitude["terr_latitude"][:] * fact
    terr_lon = attitude["terr_longitude"][:] * fact
    earth_limb = attitude["earth_limb"][:] * fact
    sun_angle = attitude["sun_angle"][:] * fact
    moon_angle = attitude["moon_angle"][:] * fact
    gal_lon = attitude["galactic_longitude"][:] * fact
    gal_lat = attitude["galactic_latitude"][:] * fact
    ecl_lon = attitude["ecliptic_longitude"][:] * fact
    ecl_lat = attitude["ecliptic_latitude"][:] * fact
    solution = attitude["solution"][:]

    # engineering data based on channels
    bol_cmd_bias = en_stat["bol_cmd_bias"][:, channel_i]

    a_lo_bol_assem = grt["a_lo_bol_assem"][:, channel_i]
    a_hi_bol_assem = grt["a_hi_bol_assem"][:, channel_i]
    b_lo_bol_assem = grt["b_lo_bol_assem"][:, channel_i]
    b_hi_bol_assem = grt["b_hi_bol_assem"][:, channel_i]
    
    bol_volt = group1["bol_volt"][:, channel_i]