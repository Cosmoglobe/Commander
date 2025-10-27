"""
So for this version we will keep each channel separate and use the midpoint_time to interpolate to get the temperatures.
"""

import astropy.units as u
import h5py
import numpy as np

import data.utils.my_utils as data_utils
import globals as g
from data import stats
from engineering_timing.get_interpolated_times import get_data, get_interp

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
grts, grt_times = get_data()
# print(grt_times)
# quit()

group1 = en_analog["group1"]

chan = fdq_eng["chan"]

for channel, channel_i in g.CHANNELS.items():
    science_data = fdq_sdf[f"fdq_sdf_{channel}"]

    all_data = {}

    sci_head = science_data["sci_head"]
    all_data["gain"] = sci_head["gain"][:]
    all_data["mtm_speed"] = sci_head["mtm_speed"][:]
    all_data["mtm_length"] = sci_head["mtm_length"][:]
    all_data["upmode"] = sci_head["sc_head1a"][:]
    all_data["adds_per_group"] = sci_head["sc_head9"][:]
    all_data["sweeps"] = sci_head["sc_head11"][:]
    all_data["saturated"] = sci_head["sc_head20"][:]
    all_data["glitch"] = sci_head["sc_head21"][:]

    ifg_data = science_data["ifg_data"]
    all_data["ifg"] = ifg_data["ifg"][:]
    all_data["gltch"] = ifg_data["gltch"][:]

    dq_data = science_data["dq_data"]
    all_data["fake"] = dq_data["fake"][:]
    all_data["xcal_pos"] = dq_data["xcal_pos"][:]
    all_data["data_quality"] = dq_data["data_quality"][:]

    collect_time = science_data["collect_time"]
    # all_data["midpoint_time"] = data_utils.binary_to_gmt(
    #     collect_time["midpoint_time"][:]
    # )
    all_data["midpoint_time"] = collect_time["midpoint_time"][:]
    all_data["midpoint_time_s"] = (collect_time["midpoint_time"][:] * (100 * u.ns)).to(
        "s"
    )
    all_data["midpoint_time_gmt"] = data_utils.binary_to_gmt(
        collect_time["midpoint_time"][:]
    )

    attitude = science_data["attitude"]
    all_data["cel_lon"] = attitude["ra"][:] * fact
    all_data["cel_lat"] = attitude["dec"][:] * fact
    all_data["terr_lat"] = attitude["terr_latitude"][:] * fact
    all_data["terr_lon"] = attitude["terr_longitude"][:] * fact
    all_data["earth_limb"] = attitude["earth_limb"][:] * fact
    all_data["sun_angle"] = attitude["sun_angle"][:] * fact
    all_data["moon_angle"] = attitude["moon_angle"][:] * fact
    all_data["gal_lon"] = attitude["galactic_longitude"][:] * fact
    all_data["gal_lat"] = attitude["galactic_latitude"][:] * fact
    all_data["ecl_lon"] = attitude["ecliptic_longitude"][:] * fact
    all_data["ecl_lat"] = attitude["ecliptic_latitude"][:] * fact
    all_data["solution"] = attitude["solution"][:]

    print(f"\n\nChannel: {channel.upper()}")

    stats.table3_4(all_data["xcal_pos"], all_data["mtm_length"], all_data["mtm_speed"])

    fpp_fail, fakeit, xcal_transit = stats.table4_1(
        all_data["data_quality"], all_data["fake"], all_data["xcal_pos"]
    )

    for key in all_data:
        all_data[key] = all_data[key][fpp_fail == False]
        all_data[key] = all_data[key][fakeit == False]
        all_data[key] = all_data[key][xcal_transit == False]

    cal_data = {}
    sky_data = {}
    for key in all_data:
        cal_data[key] = all_data[key][all_data["xcal_pos"] == 1]
        sky_data[key] = all_data[key][all_data["xcal_pos"] == 2]

    (
        cal_saturated,
        cal_sw,
        cal_glitch_rate,
        sky_saturated,
        sky_sw,
        sky_glitch_rate,
        limb,
        no_solution,
    ) = stats.table4_2(
        cal_data["saturated"],
        cal_data["data_quality"],
        sky_data["saturated"],
        sky_data["data_quality"],
        sky_data["solution"],
    )

    for key in all_data:
        # TODO: unclear what they consider a glitch rate too high, check IFGs marked with that?
        cal_data[key] = cal_data[key][
            (cal_saturated == False) & (cal_sw == False) & (cal_glitch_rate == False)
        ]
        sky_data[key] = sky_data[key][
            (sky_saturated == False)
            & (sky_sw == False)
            & (sky_glitch_rate == False)
            & (limb == False)
            & (no_solution == False)
        ]

    # the next table to reproduce needs the ICAL temperatures so we need to match them now
    interpolators = get_interp(grts, grt_times)
    # Interpolate temperatures for sky data
    elements = ["ical", "dihedral"]
    sides = ["a", "b"]

    for element in elements:
        for side in sides:
            # Interpolate hi and lo temperatures
            sky_data[f"{side}_hi_{element}"] = interpolators[f"{side}_hi_{element}"](
                sky_data["midpoint_time_s"]
            )
            sky_data[f"{side}_lo_{element}"] = interpolators[f"{side}_lo_{element}"](
                sky_data["midpoint_time_s"]
            )
            # Vectorized temperature selection
            # TODO: we probably want to change this
            sky_data[f"{side}_{element}"] = data_utils.get_temperature_hl_vectorized(
                sky_data[f"{side}_lo_{element}"],
                sky_data[f"{side}_hi_{element}"],
                element,
                side,
            )

        # Average temperatures from both sides
        # TODO: we probably want to change this
        sky_data[element] = (sky_data[f"a_{element}"] + sky_data[f"b_{element}"]) / 2.0
        print(f"{element} temperatures matched for sky data.")

    earth_limb, wrong_ical_temp, sun_angle, wrong_sci_mode, dihedral_temp = (
        stats.table4_5(
            sky_data["earth_limb"],
            sky_data["midpoint_time_gmt"],
            sky_data["ical"],
            sky_data["sun_angle"],
            sky_data["upmode"],
            sky_data["dihedral"],
        )
    )

    for key in sky_data:
        sky_data[key] = sky_data[key][
            (earth_limb == False)
            & (wrong_ical_temp == False)
            & (sun_angle == False)
            & (wrong_sci_mode == False)
            & (dihedral_temp == False)
        ]

    # engineering data based on channels
    bol_cmd_bias = en_stat["bol_cmd_bias"][:, channel_i]

    a_lo_bol_assem = grt["a_lo_bol_assem"][:, channel_i]
    a_hi_bol_assem = grt["a_hi_bol_assem"][:, channel_i]
    b_lo_bol_assem = grt["b_lo_bol_assem"][:, channel_i]
    b_hi_bol_assem = grt["b_hi_bol_assem"][:, channel_i]

    bol_volt = group1["bol_volt"][:, channel_i]

    eng_upmode = chan["up_sci_mode"][:, channel_i]
    eng_fake = chan["fakeit"][:, channel_i]
    eng_adds_per_group = chan["up_adds_per_group"][:, channel_i]
    eng_sweeps = chan["up_swps_per_ifg"][:, channel_i]
    eng_mtm_speed = chan["xmit_mtm_speed"][:, channel_i]
    eng_mtm_length = chan["xmit_mtm_len"][:, channel_i]
    eng_gain = chan["sci_gain"][:, channel_i]
