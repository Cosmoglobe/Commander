import matplotlib.pyplot as plt
import numpy as np

import data.utils.my_utils as data_utils


def table3_4(xcal_pos, mtm_length, mtm_speed):
    """Compute and print scan mode statistics as in table 3.4 of the FIRAS Explanatory Supplement."""

    science_data = xcal_pos == 2
    mtm_length = mtm_length[science_data]
    mtm_speed = mtm_speed[science_data]

    ss = (mtm_length == 0) & (mtm_speed == 0)
    sf = (mtm_length == 0) & (mtm_speed == 1)
    ls = (mtm_length == 1) & (mtm_speed == 0)
    lf = (mtm_length == 1) & (mtm_speed == 1)
    print("\nTable 3.4")
    print(
        f"    SS: {ss.sum()/len(mtm_length)*100:.1f}, SF: {sf.sum()/len(mtm_length)*100:.1f}, LS: {ls.sum()/len(mtm_length)*100:.1f}, LF: {lf.sum()/len(mtm_length)*100:.1f}"
    )


def table4_1(data_quality, fake, xcal_pos):
    print("\nTable 4.1")
    print(f"    Total Science Records: {len(data_quality)}")

    # so it seems that the FPP fails correspond to data quality flag 109 being 0 or flag 110 being 32 (always the same)
    fpp_fail = data_quality[:, -1] == 32
    print(f"    Records Failed by FPP: {fpp_fail.sum()}")

    fake = fake[fpp_fail == False]
    xcal_pos = xcal_pos[fpp_fail == False]
    fakeit = fake == 1
    print(f"    Fake-it Records: {fakeit.sum()}")

    xcal_pos = xcal_pos[fakeit == False]
    xcal_transit = xcal_pos == 3
    print(f"    XCAL Transit Records: {xcal_transit.sum()}")

    print(
        f"    Records Eliminated Before FDQ: {(fpp_fail.sum() + fakeit.sum() + xcal_transit.sum())}"
    )

    xcal_pos = xcal_pos[xcal_transit == False]
    cal = xcal_pos == 1
    sci = xcal_pos == 2
    print(f"    Cal Records Passed by FPP: {cal.sum()}")
    print(f"    Sky Records Passed by FPP: {sci.sum()}")

    return fpp_fail, fakeit, xcal_transit


def table4_2(
    cal_saturated, cal_data_quality, sky_saturated, sky_data_quality, solution
):
    # TODO: this is the most questionable data selection
    print("\nTable 4.2")

    # calibration data
    # TODO: there are some negative values, which they took as non-saturated, should we do the same? here we are doing the same as them
    cal_saturated = cal_saturated > 0
    print(f"    Saturated Sample Count: {cal_saturated.sum()}")
    cal_sw = cal_data_quality[:, 3] != 0
    print(f"    Microprocessor Status Word: {(cal_sw).sum()}")
    cal_glitch_rate = cal_data_quality[:, 7] > 0
    print(f"    Glitch Rate: {cal_glitch_rate.sum()}")
    # TODO: these totals don't correspond to the tables in the explanatory supplement
    print(
        f"    Cal Records Failed by FDQ: {(cal_saturated | cal_sw | cal_glitch_rate).sum()}"
    )
    print(
        f"    Cal records Passed by FDQ: {(~(cal_saturated | cal_sw | cal_glitch_rate)).sum()}"
    )

    # sky data
    sky_saturated = sky_saturated > 0
    print(f"    Saturated Sample Count: {sky_saturated.sum()}")
    sky_sw = sky_data_quality[:, 3] != 0
    print(f"    Microprocessor Status Word: {(sky_sw).sum()}")
    sky_glitch_rate = sky_data_quality[:, 7] > 0
    print(f"    Glitch Rate: {sky_glitch_rate.sum()}")
    limb = sky_data_quality[:, 5] != 0
    print(f"    Sun, Moon, or Earth Limb: {limb.sum()}")
    no_solution = solution == 0
    print(f"    No Attitude Solution: {no_solution.sum()}")
    print(
        f"    Sky Records Failed by FDQ: {(sky_saturated | sky_sw | sky_glitch_rate | limb | no_solution).sum()}"
    )
    print(
        f"    Sky records Passed by FDQ: {(~(sky_saturated | sky_sw | sky_glitch_rate | limb | no_solution)).sum()}"
    )

    return (
        cal_saturated,
        cal_sw,
        cal_glitch_rate,
        sky_saturated,
        sky_sw,
        sky_glitch_rate,
        limb,
        no_solution,
    )


def table4_5(
    earth_limb,
    midpoint_time_cal,
    midpoint_time_sky,
    ical_cal,
    ical_sky,
    sun_angle,
    upmode,
    dihedral_cal,
    dihedral_sky,
):
    print("\nTable 4.5")

    earth_limb = earth_limb < 87
    print(f"    Earth Limb Angle < 87.0: {earth_limb.sum()}")

    start_times = [
        "89-326-1130",
        "89-328-0000",
        "89-343-0152",
        "90-019-0205",
        "90-080-0115",
        "90-129-0000",
        "90-139-1535",
        "90-193-1850",
        "90-207-1104",
        "90-208-1120",
        "90-220-0500",
    ]
    # Vectorized date parsing using list comprehension (faster than loop)
    start_times = [
        data_utils.parse_date_string(t.replace("-", ""), include_s=False)
        for t in start_times
    ]

    end_times = [
        "89-327-2359",
        "89-343-0151",
        "90-019-0204",
        "90-080-0114",
        "90-128-2359",
        "90-139-1534",
        "90-193-1849",
        "90-207-1103",
        "90-208-1119",
        "90-220-0459",
        "90-264-0936",
    ]
    end_times = [
        data_utils.parse_date_string(t.replace("-", ""), include_s=False)
        for t in end_times
    ]

    allowed_ical_temps = [
        [2.789],
        [2.758, 2.763, 2.789],
        [2.759, 2.771],
        [2.758, 2.771],
        [2.758, 2.771],
        [2.758, 2.770],
        [2.7455, 2.755, 2.768],
        [2.746, 2.757, 2.769],
        [2.757, 2.769],
        [2.758, 2.770],
        [2.758, 2.771],
    ]

    wrong_ical_temp_cal = np.zeros(len(midpoint_time_cal), dtype=bool)
    for i, (start, end, temps) in enumerate(
        zip(start_times, end_times, allowed_ical_temps)
    ):
        time_mask = (midpoint_time_cal > start) & (midpoint_time_cal < end)
        temp_mask = np.ones(len(midpoint_time_cal), dtype=bool)
        for temp in temps:
            temp_mask &= (ical_cal > temp + 0.002) | (ical_cal < temp - 0.002)
        wrong_ical_temp_cal |= time_mask & temp_mask

    wrong_ical_temp_sky = np.zeros(len(midpoint_time_sky), dtype=bool)
    for i, (start, end, temps) in enumerate(
        zip(start_times, end_times, allowed_ical_temps)
    ):
        time_mask = (midpoint_time_sky > start) & (midpoint_time_sky < end)
        temp_mask = np.ones(len(midpoint_time_sky), dtype=bool)
        for temp in temps:
            temp_mask &= (ical_sky > temp + 0.002) | (ical_sky < temp - 0.002)
        wrong_ical_temp_sky |= time_mask & temp_mask
    print(f"    Wrong ICAL Temperature: {sum(wrong_ical_temp_sky)}")

    sun_angle = sun_angle < 91.2
    print(f"    Sun Angle < 91.2: {sun_angle.sum()}")

    wrong_sci_mode = upmode != 4
    print(f"    Wrong Science Mode: {wrong_sci_mode.sum()}")

    dihedral_temp_cal = dihedral_cal > 5.5
    dihedral_temp_sky = dihedral_sky > 5.5
    print(f"    Dihedral Temperature > 5.5: {dihedral_temp_sky.sum()}")

    print(
        f"    Sky Records Failed by FSS: {(earth_limb | wrong_ical_temp_sky | sun_angle | wrong_sci_mode | dihedral_temp_sky).sum()}"
    )
    print(
        f"Sky Records Passed by FSS: {len(midpoint_time_sky) - (earth_limb | wrong_ical_temp_sky | sun_angle | wrong_sci_mode | dihedral_temp_sky).sum()}"
    )

    return (
        earth_limb,
        wrong_ical_temp_cal,
        wrong_ical_temp_sky,
        sun_angle,
        wrong_sci_mode,
        dihedral_temp_cal,
        dihedral_temp_sky,
    )
