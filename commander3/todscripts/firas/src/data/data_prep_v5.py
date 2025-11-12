"""
So for this version we will keep each channel separate and use the midpoint_time to interpolate to get the temperatures.
"""

import astropy.units as u
import h5py
import numpy as np
from scipy.interpolate import interp1d

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
eng_time = np.sort(ct_head["time"][:])
eng_time_gmt = data_utils.binary_to_gmt(eng_time)
eng_time_s = (eng_time * (100 * u.ns)).to("s")

en_stat = fdq_eng["en_stat"]
# side a
stat_word_1 = en_stat["stat_word_1"][:]
stat_word_4 = en_stat["stat_word_4"][:]
stat_word_5 = en_stat["stat_word_5"][:]
stat_word_8 = en_stat["stat_word_8"][:]
lvdt_stat_a = en_stat["lvdt_stat"][:, 0]
# side b
stat_word_9 = en_stat["stat_word_9"][:]
stat_word_12 = en_stat["stat_word_12"][:]
stat_word_13 = en_stat["stat_word_13"][:]
stat_word_16 = en_stat["stat_word_16"][:]
lvdt_stat_b = en_stat["lvdt_stat"][:, 1]

en_analog = fdq_eng["en_analog"]
grt = en_analog["grt"]
grts, grt_times = get_data()
# remove nans from grts and corresponding times
# valid_indices = ~np.isnan(grts).any(axis=1)
# grts = grts[valid_indices]
# grt_times = grt_times[valid_indices]

group1 = en_analog["group1"]

chan = fdq_eng["chan"]

for channel, channel_i in g.CHANNELS.items():
    science_data = fdq_sdf[f"fdq_sdf_{channel}"]

    all_data = {}

    sci_head = science_data["sci_head"]
    # all_data["gain"] = data_utils.convert_gain_array(sci_head["gain"][:])
    all_data["mtm_speed"] = sci_head["mtm_speed"][:]
    all_data["mtm_length"] = sci_head["mtm_length"][:]
    all_data["upmode"] = sci_head["sc_head1a"][:]
    all_data["adds_per_group"] = sci_head["sc_head9"][:]
    all_data["sweeps"] = sci_head["sc_head11"][:]
    all_data["saturated"] = sci_head["sc_head20"][:]

    ifg_data = science_data["ifg_data"]
    all_data["ifg"] = ifg_data["ifg"][:]

    dq_data = science_data["dq_data"]
    all_data["fake"] = dq_data["fake"][:]
    all_data["xcal_pos"] = dq_data["xcal_pos"][:]
    all_data["data_quality"] = dq_data["data_quality"][:]
    all_data["eng_time"] = np.sort(dq_data["eng_time"][:])

    collect_time = science_data["collect_time"]
    all_data["midpoint_time"] = np.sort(collect_time["midpoint_time"][:])
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

    # TODO: unclear what they consider a glitch rate too high, check IFGs marked with that?
    cal_cuts = cal_saturated | cal_sw | cal_glitch_rate
    sky_cuts = sky_saturated | sky_sw | sky_glitch_rate | limb | no_solution

    for key in all_data:
        cal_data[key] = cal_data[key][cal_cuts == False]
        sky_data[key] = sky_data[key][sky_cuts == False]

    # the next table to reproduce needs the ICAL temperatures so we need to match them now
    interpolators = get_interp(grts, grt_times)
    # Interpolate temperatures for sky data
    elements = ["ical", "dihedral", "refhorn", "skyhorn", "xcal_cone"]
    sides = ["a", "b"]

    print(f"Using reference file to decide between high and low currents")
    print(f"Taking the average of both sides")
    temps = {}
    for element in elements:
        for side in sides:
            # Interpolate hi and lo temperatures for cal data
            temps[f"{side}_hi_{element}"] = interpolators[f"{side}_hi_{element}"](
                cal_data["midpoint_time_s"]
            )
            temps[f"{side}_lo_{element}"] = interpolators[f"{side}_lo_{element}"](
                cal_data["midpoint_time_s"]
            )

            # Vectorized temperature selection
            temps[f"{side}_{element}"] = data_utils.get_temperature_hl_vectorized(
                temps[f"{side}_lo_{element}"],
                temps[f"{side}_hi_{element}"],
                element,
                side,
            )
        cal_data[element] = (temps[f"a_{element}"] + temps[f"b_{element}"]) / 2.0

        for side in sides:
            # Interpolate hi and lo temperatures
            temps[f"{side}_hi_{element}"] = interpolators[f"{side}_hi_{element}"](
                sky_data["midpoint_time_s"]
            )
            temps[f"{side}_lo_{element}"] = interpolators[f"{side}_lo_{element}"](
                sky_data["midpoint_time_s"]
            )

            # Vectorized temperature selection
            # TODO: we probably want to change this
            temps[f"{side}_{element}"] = data_utils.get_temperature_hl_vectorized(
                temps[f"{side}_lo_{element}"],
                temps[f"{side}_hi_{element}"],
                element,
                side,
            )
        # Average temperatures from both sides
        # TODO: we probably want to change this
        if "xcal" not in element:
            sky_data[element] = (temps[f"a_{element}"] + temps[f"b_{element}"]) / 2.0

    collimator_hi = interpolators["a_hi_collimator"](cal_data["midpoint_time_s"])
    collimator_lo = interpolators["a_lo_collimator"](cal_data["midpoint_time_s"])
    cal_data["collimator"] = (collimator_hi + collimator_lo) / 2.0

    collimator_hi = interpolators["a_hi_collimator"](sky_data["midpoint_time_s"])
    collimator_lo = interpolators["a_lo_collimator"](sky_data["midpoint_time_s"])
    sky_data["collimator"] = (collimator_hi + collimator_lo) / 2.0

    (
        earth_limb,
        wrong_ical_temp_cal,
        wrong_ical_temp_sky,
        sun_angle,
        wrong_sci_mode,
        dihedral_temp_cal,
        dihedral_temp_sky,
    ) = stats.table4_5(
        sky_data["earth_limb"],
        cal_data["midpoint_time_gmt"],
        sky_data["midpoint_time_gmt"],
        cal_data["ical"],
        sky_data["ical"],
        sky_data["sun_angle"],
        sky_data["upmode"],
        cal_data["dihedral"],
        sky_data["dihedral"],
    )

    cal_cuts = wrong_ical_temp_cal | dihedral_temp_cal
    for key in cal_data:
        cal_data[key] = cal_data[key][cal_cuts == False]
    sky_cuts = (
        earth_limb
        | wrong_ical_temp_sky
        | sun_angle
        | wrong_sci_mode
        | dihedral_temp_sky
    )
    for key in sky_data:
        sky_data[key] = sky_data[key][sky_cuts == False]

    # engineering data based on channels
    bol_cmd_bias = en_stat["bol_cmd_bias"][:, channel_i]
    bol_volt = group1["bol_volt"][:, channel_i]
    eng_upmode = chan["up_sci_mode"][:, channel_i]
    eng_fake = chan["fakeit"][:, channel_i]
    eng_adds_per_group = chan["up_adds_per_group"][:, channel_i]
    eng_sweeps = chan["up_swps_per_ifg"][:, channel_i]
    eng_mtm_speed = chan["xmit_mtm_speed"][:, channel_i]
    eng_mtm_length = chan["xmit_mtm_len"][:, channel_i]
    eng_gain = chan["sci_gain"][:, channel_i].astype(int)

    # Find two nearest engineering times for each science time using searchsorted
    # This assumes eng_time is sorted (which it should be)
    # Use searchsorted to find insertion indices

    # cal data
    indices = np.searchsorted(eng_time_gmt, cal_data["midpoint_time_gmt"])
    # Clip indices to valid range
    indices = np.clip(indices, 1, len(eng_time_gmt) - 1)

    # Get the two nearest neighbors (before and after)
    idx_before = indices - 1
    idx_after = indices

    # Calculate distances to both neighbors
    dist_before = np.abs(cal_data["midpoint_time_gmt"] - eng_time_gmt[idx_before])
    dist_after = np.abs(cal_data["midpoint_time_gmt"] - eng_time_gmt[idx_after])

    # Get indices of the two closest points
    idx0_cal = np.where(dist_before <= dist_after, idx_before, idx_after)
    idx1_cal = np.where(dist_before <= dist_after, idx_after, idx_before)

    # sky data
    indices = np.searchsorted(eng_time_gmt, sky_data["midpoint_time_gmt"])
    # Clip indices to valid range
    indices = np.clip(indices, 1, len(eng_time_gmt) - 1)

    # Get the two nearest neighbors (before and after)
    idx_before = indices - 1
    idx_after = indices

    # Calculate distances to both neighbors
    dist_before = np.abs(sky_data["midpoint_time_gmt"] - eng_time_gmt[idx_before])
    dist_after = np.abs(sky_data["midpoint_time_gmt"] - eng_time_gmt[idx_after])

    # Get indices of the two closest points
    idx0_sky = np.where(dist_before <= dist_after, idx_before, idx_after)
    idx1_sky = np.where(dist_before <= dist_after, idx_after, idx_before)

    print("\nOther data cuts")

    # cal data
    bol_cmd_bias_mismatch_cal = (bol_cmd_bias[idx0_cal] != bol_cmd_bias[idx1_cal]) | (
        bol_cmd_bias[idx0_cal] <= 0
    )
    cal_data["bol_cmd_bias"] = np.where(
        ~bol_cmd_bias_mismatch_cal,
        bol_cmd_bias[idx0_cal],
        0.0,  # or np.nan if you prefer
    )
    upmode_mismatch_cal = (eng_upmode[idx0_cal] != eng_upmode[idx1_cal]) | (
        eng_upmode[idx0_cal] != cal_data["upmode"]
    )
    fakeit_mismatch_cal = (eng_fake[idx0_cal] != eng_fake[idx1_cal]) | (
        eng_fake[idx0_cal] != cal_data["fake"]
    )
    adds_per_group_mismatch_cal = (
        eng_adds_per_group[idx0_cal] != eng_adds_per_group[idx1_cal]
    ) | (eng_adds_per_group[idx0_cal] != cal_data["adds_per_group"])
    sweeps_mismatch_cal = (
        (eng_sweeps[idx0_cal] != eng_sweeps[idx1_cal])
        | (eng_sweeps[idx0_cal] != cal_data["sweeps"])
        | (cal_data["sweeps"] == 1)
    )
    mtm_length_mismatch_cal = (eng_mtm_length[idx0_cal] != eng_mtm_length[idx1_cal]) | (
        eng_mtm_length[idx0_cal] != cal_data["mtm_length"]
    )
    mtm_speed_mismatch_cal = (eng_mtm_speed[idx0_cal] != eng_mtm_speed[idx1_cal]) | (
        eng_mtm_speed[idx0_cal] != cal_data["mtm_speed"]
    )
    gain_mismatch_cal = (eng_gain[idx0_cal] != eng_gain[idx1_cal]) | (
        eng_gain[idx0_cal] == 0
    )
    cal_data["gain"] = np.where(~gain_mismatch_cal, eng_gain[idx0_cal], np.nan)

    sw1_mismatch_cal = stat_word_1[idx0_cal] != stat_word_1[idx1_cal]
    sw4_mismatch_cal = stat_word_4[idx0_cal] != stat_word_4[idx1_cal]
    sw5_mismatch_cal = stat_word_5[idx0_cal] != stat_word_5[idx1_cal]
    sw8_mismatch_cal = stat_word_8[idx0_cal] != stat_word_8[idx1_cal]
    sw9_mismatch_cal = stat_word_9[idx0_cal] != stat_word_9[idx1_cal]
    sw12_mismatch_cal = stat_word_12[idx0_cal] != stat_word_12[idx1_cal]
    sw13_mismatch_cal = stat_word_13[idx0_cal] != stat_word_13[idx1_cal]
    sw16_mismatch_cal = stat_word_16[idx0_cal] != stat_word_16[idx1_cal]
    lvdt_stat_a_mismatch_cal = lvdt_stat_a[idx0_cal] != lvdt_stat_a[idx1_cal]
    lvdt_stat_b_mismatch_cal = lvdt_stat_b[idx0_cal] != lvdt_stat_b[idx1_cal]
    cal_data["stat_word_1"] = np.where(~sw1_mismatch_cal, stat_word_1[idx0_cal], 0)
    cal_data["stat_word_4"] = np.where(~sw4_mismatch_cal, stat_word_4[idx0_cal], 0)
    cal_data["stat_word_5"] = np.where(~sw5_mismatch_cal, stat_word_5[idx0_cal], 0)
    cal_data["stat_word_8"] = np.where(~sw8_mismatch_cal, stat_word_8[idx0_cal], 0)
    cal_data["stat_word_9"] = np.where(~sw9_mismatch_cal, stat_word_9[idx0_cal], 0)
    cal_data["stat_word_12"] = np.where(~sw12_mismatch_cal, stat_word_12[idx0_cal], 0)
    cal_data["stat_word_13"] = np.where(~sw13_mismatch_cal, stat_word_13[idx0_cal], 0)
    cal_data["stat_word_16"] = np.where(~sw16_mismatch_cal, stat_word_16[idx0_cal], 0)
    cal_data["lvdt_stat_a"] = np.where(
        ~lvdt_stat_a_mismatch_cal, lvdt_stat_a[idx0_cal], 0
    )
    cal_data["lvdt_stat_b"] = np.where(
        ~lvdt_stat_b_mismatch_cal, lvdt_stat_b[idx0_cal], 0
    )
    sw_mismatches_cal = (
        sw1_mismatch_cal
        | sw4_mismatch_cal
        | sw5_mismatch_cal
        | sw8_mismatch_cal
        | sw9_mismatch_cal
        | sw12_mismatch_cal
        | sw13_mismatch_cal
        | sw16_mismatch_cal
        | lvdt_stat_a_mismatch_cal
        | lvdt_stat_b_mismatch_cal
    )
    other_cuts_cal = (
        bol_cmd_bias_mismatch_cal
        | upmode_mismatch_cal
        | fakeit_mismatch_cal
        | adds_per_group_mismatch_cal
        | sweeps_mismatch_cal
        | mtm_length_mismatch_cal
        | mtm_speed_mismatch_cal
        | gain_mismatch_cal
        | sw_mismatches_cal
    )
    for key in cal_data:
        cal_data[key] = cal_data[key][other_cuts_cal == False]

    # sky data
    # Check if the two closest values match
    bol_cmd_bias_mismatch = (bol_cmd_bias[idx0_sky] != bol_cmd_bias[idx1_sky]) | (
        bol_cmd_bias[idx0_sky] <= 0
    )
    print(f"    bol_cmd_bias mismatch: {(bol_cmd_bias_mismatch).sum()}")
    # Assign values where they match
    sky_data["bol_cmd_bias"] = np.where(
        ~bol_cmd_bias_mismatch, bol_cmd_bias[idx0_sky], 0.0  # or np.nan if you prefer
    )

    # compare between id1 and idx2 but also with the one from the science data
    upmode_mismatch = (eng_upmode[idx0_sky] != eng_upmode[idx1_sky]) | (
        eng_upmode[idx0_sky] != sky_data["upmode"]
    )
    print(f"    upmode mismatch: {(upmode_mismatch).sum()}")

    fakeit_mismatch = (eng_fake[idx0_sky] != eng_fake[idx1_sky]) | (
        eng_fake[idx0_sky] != sky_data["fake"]
    )
    print(f"    fakeit mismatch: {(fakeit_mismatch).sum()}")

    adds_per_group_mismatch = (
        eng_adds_per_group[idx0_sky] != eng_adds_per_group[idx1_sky]
    ) | (eng_adds_per_group[idx0_sky] != sky_data["adds_per_group"])
    print(f"    adds_per_group mismatch: {(adds_per_group_mismatch).sum()}")

    sweeps_mismatch = (
        (eng_sweeps[idx0_sky] != eng_sweeps[idx1_sky])
        | (eng_sweeps[idx0_sky] != sky_data["sweeps"])
        | (sky_data["sweeps"] == 1)
    )
    print(f"    sweeps mismatch: {(sweeps_mismatch).sum()}")

    mtm_length_mismatch = (eng_mtm_length[idx0_sky] != eng_mtm_length[idx1_sky]) | (
        eng_mtm_length[idx0_sky] != sky_data["mtm_length"]
    )
    print(f"    mtm_length mismatch: {(mtm_length_mismatch).sum()}")

    mtm_speed_mismatch = (eng_mtm_speed[idx0_sky] != eng_mtm_speed[idx1_sky]) | (
        eng_mtm_speed[idx0_sky] != sky_data["mtm_speed"]
    )
    print(f"    mtm_speed mismatch: {(mtm_speed_mismatch).sum()}")

    # let's try to simply use the values from engineering for the gain
    # bad_gain = (sky_data["gain"] == np.nan) != (eng_gain[idx0] == eng_gain[idx1]) | (
    #     eng_gain[idx0] != sky_data["gain"]
    # )
    gain_mismatch = (eng_gain[idx0_sky] != eng_gain[idx1_sky]) | (
        eng_gain[idx0_sky] == 0
    )
    print(f"    Gain Mismatch: {gain_mismatch.sum()}")
    sky_data["gain"] = np.where(~gain_mismatch, eng_gain[idx0_sky], np.nan)

    moon_contamination = sky_data["moon_angle"] <= 22.0
    print(f"    Moon Angle <= 22.0: {moon_contamination.sum()}")

    # also the status words now
    sw1_mismatch = stat_word_1[idx0_sky] != stat_word_1[idx1_sky]
    sw4_mismatch = stat_word_4[idx0_sky] != stat_word_4[idx1_sky]
    sw5_mismatch = stat_word_5[idx0_sky] != stat_word_5[idx1_sky]
    sw8_mismatch = stat_word_8[idx0_sky] != stat_word_8[idx1_sky]
    sw9_mismatch = stat_word_9[idx0_sky] != stat_word_9[idx1_sky]
    sw12_mismatch = stat_word_12[idx0_sky] != stat_word_12[idx1_sky]
    sw13_mismatch = stat_word_13[idx0_sky] != stat_word_13[idx1_sky]
    sw16_mismatch = stat_word_16[idx0_sky] != stat_word_16[idx1_sky]
    lvdt_stat_a_mismatch = lvdt_stat_a[idx0_sky] != lvdt_stat_a[idx1_sky]
    lvdt_stat_b_mismatch = lvdt_stat_b[idx0_sky] != lvdt_stat_b[idx1_sky]
    sky_data["stat_word_1"] = np.where(~sw1_mismatch, stat_word_1[idx0_sky], 0)
    sky_data["stat_word_4"] = np.where(~sw4_mismatch, stat_word_4[idx0_sky], 0)
    sky_data["stat_word_5"] = np.where(~sw5_mismatch, stat_word_5[idx0_sky], 0)
    sky_data["stat_word_8"] = np.where(~sw8_mismatch, stat_word_8[idx0_sky], 0)
    sky_data["stat_word_9"] = np.where(~sw9_mismatch, stat_word_9[idx0_sky], 0)
    sky_data["stat_word_12"] = np.where(~sw12_mismatch, stat_word_12[idx0_sky], 0)
    sky_data["stat_word_13"] = np.where(~sw13_mismatch, stat_word_13[idx0_sky], 0)
    sky_data["stat_word_16"] = np.where(~sw16_mismatch, stat_word_16[idx0_sky], 0)
    sky_data["lvdt_stat_a"] = np.where(~lvdt_stat_a_mismatch, lvdt_stat_a[idx0_sky], 0)
    sky_data["lvdt_stat_b"] = np.where(~lvdt_stat_b_mismatch, lvdt_stat_b[idx0_sky], 0)
    sw_mismatches = (
        sw1_mismatch
        | sw4_mismatch
        | sw5_mismatch
        | sw8_mismatch
        | sw9_mismatch
        | sw12_mismatch
        | sw13_mismatch
        | sw16_mismatch
        | lvdt_stat_a_mismatch
        | lvdt_stat_b_mismatch
    )
    print(f"    Status Word Mismatches: {sw_mismatches.sum()}")

    other_cuts = (
        bol_cmd_bias_mismatch
        | upmode_mismatch
        | fakeit_mismatch
        | adds_per_group_mismatch
        | sweeps_mismatch
        | mtm_length_mismatch
        | mtm_speed_mismatch
        | gain_mismatch
        | moon_contamination
        | sw_mismatches
    )
    print(f"    Total Sky Records Failed Other Cuts: {other_cuts.sum()}")
    print(f"    Remaining Sky Records: {len(sky_data['ifg']) - other_cuts.sum()}")
    # only the bol_cmd_bias is actually cutting stuff, so i'm guessing these other cuts happen in the data quality flag but i will leave them here anyways
    for key in sky_data:
        sky_data[key] = sky_data[key][other_cuts == False]

    interp_func = interp1d(eng_time_s, bol_volt)
    cal_data["bol_volt"] = interp_func(cal_data["midpoint_time_s"])
    sky_data["bol_volt"] = interp_func(sky_data["midpoint_time_s"])

    # no high cuirrent readings for bolometers
    a_lo_bol_assem = grt["a_lo_bol_assem"][:, channel_i]
    b_lo_bol_assem = grt["b_lo_bol_assem"][:, channel_i]
    bol_assem = (a_lo_bol_assem + b_lo_bol_assem) / 2.0
    interp_func = interp1d(eng_time_s, bol_assem)
    cal_data["bolometer"] = interp_func(cal_data["midpoint_time_s"])
    bad_bolometer = cal_data["bolometer"] <= 0.0
    for key in cal_data:
        cal_data[key] = cal_data[key][bad_bolometer == False]
    sky_data["bolometer"] = interp_func(sky_data["midpoint_time_s"])
    bad_bolometer = sky_data["bolometer"] <= 0.0
    print(f"    Bad Bolometer Temperature Readings: {bad_bolometer.sum()}")
    for key in sky_data:
        sky_data[key] = sky_data[key][bad_bolometer == False]

    # save
    np.savez(
        f"{g.PREPROCESSED_DATA_PATH}/sky_{channel}.npz",
        **sky_data,
    )
    np.savez(
        f"{g.PREPROCESSED_DATA_PATH}/cal_{channel}.npz",
        **cal_data,
    )
