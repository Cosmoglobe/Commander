import matplotlib.pyplot as plt
import numpy as np
from matplotlib import markers
from scipy.stats import norm

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

    # start_times = [
    #     "89-326-1130",
    #     "89-328-0000",
    #     "89-343-0152",
    #     "90-019-0205",
    #     "90-080-0115",
    #     "90-129-0000",
    #     "90-139-1535",
    #     "90-193-1850",
    #     "90-207-1104",
    #     "90-208-1120",
    #     "90-220-0500",
    # ]
    # # Vectorized date parsing using list comprehension (faster than loop)
    # start_times = [
    #     data_utils.parse_date_string(t.replace("-", ""), include_s=False)
    #     for t in start_times
    # ]

    # end_times = [
    #     "89-327-2359",
    #     "89-343-0151",
    #     "90-019-0204",
    #     "90-080-0114",
    #     "90-128-2359",
    #     "90-139-1534",
    #     "90-193-1849",
    #     "90-207-1103",
    #     "90-208-1119",
    #     "90-220-0459",
    #     "90-264-0936",
    # ]
    # end_times = [
    #     data_utils.parse_date_string(t.replace("-", ""), include_s=False)
    #     for t in end_times
    # ]

    # allowed_ical_temps = [
    #     [2.789],
    #     [2.758, 2.763, 2.789],
    #     [2.759, 2.771],
    #     [2.758, 2.771],
    #     [2.758, 2.771],
    #     [2.758, 2.770],
    #     [2.7455, 2.755, 2.768],
    #     [2.746, 2.757, 2.769],
    #     [2.757, 2.769],
    #     [2.758, 2.770],
    #     [2.758, 2.771],
    # ]

    # wrong_ical_temp_cal = np.zeros(len(midpoint_time_cal), dtype=bool)
    # for i, (start, end, temps) in enumerate(
    #     zip(start_times, end_times, allowed_ical_temps)
    # ):
    #     time_mask = (midpoint_time_cal > start) & (midpoint_time_cal < end)
    #     temp_mask = np.ones(len(midpoint_time_cal), dtype=bool)
    #     for temp in temps:
    #         temp_mask &= (ical_cal > temp + 0.002) | (ical_cal < temp - 0.002)
    #     wrong_ical_temp_cal |= time_mask & temp_mask

    # wrong_ical_temp_sky = np.zeros(len(midpoint_time_sky), dtype=bool)
    # for i, (start, end, temps) in enumerate(
    #     zip(start_times, end_times, allowed_ical_temps)
    # ):
    #     time_mask = (midpoint_time_sky > start) & (midpoint_time_sky < end)
    #     temp_mask = np.ones(len(midpoint_time_sky), dtype=bool)
    #     for temp in temps:
    #         temp_mask &= (ical_sky > temp + 0.002) | (ical_sky < temp - 0.002)
    #     wrong_ical_temp_sky |= time_mask & temp_mask
    # print(f"    Wrong ICAL Temperature: {sum(wrong_ical_temp_sky)}")

    sun_angle = sun_angle < 91.2
    print(f"    Sun Angle < 91.2: {sun_angle.sum()}")

    wrong_sci_mode = upmode != 4
    print(f"    Wrong Science Mode: {wrong_sci_mode.sum()}")

    # dihedral_temp_cal = dihedral_cal > 5.5
    # dihedral_temp_sky = dihedral_sky > 5.5
    # print(f"    Dihedral Temperature > 5.5: {dihedral_temp_sky.sum()}")

    # print(
    #     f"    Sky Records Failed by FSS: {(earth_limb | wrong_ical_temp_sky | sun_angle | wrong_sci_mode | dihedral_temp_sky).sum()}"
    # )
    print(f"    Sky Records Failed by FSS: {(earth_limb | sun_angle | wrong_sci_mode).sum()}")
    # print(
    #     f"Sky Records Passed by FSS: {len(midpoint_time_sky) - (earth_limb | wrong_ical_temp_sky | sun_angle | wrong_sci_mode | dihedral_temp_sky).sum()}"
    # )
    print(
        f"    Sky Records Passed by FSS: {len(midpoint_time_sky) - (earth_limb | sun_angle | wrong_sci_mode).sum()}"
    )
    
    return (
        earth_limb,
        # wrong_ical_temp_cal,
        # wrong_ical_temp_sky,
        sun_angle,
        wrong_sci_mode,
        # dihedral_temp_cal,
        # dihedral_temp_sky,
    )

def hilo_stats(hi, lo, channel, element, side):
    diff = hi - lo
    plt.hist(diff, bins=np.linspace(-0.05, 0.05, 1000))
    plt.xlabel("High - Low Detector Temperature (K)")
    plt.ylabel("Counts")
    plt.title(f"High - Low Detector Temperature Difference ({channel}) for " + element + " (" + side.upper() + " side)")
    plt.yscale('log')
    plt.savefig(f"data/output/hilo_temp_difference_histogram/{element}_{side}_{channel}.png")
    plt.close()

def fit_gaussian(hi, lo, channel, element, side, ngaussians=2, calsky="sky"):
    """Fit a Gaussian to the difference between hi and lo temperatures."""

    if element == "ical":
        if side == "a":
            lims1 = (0.0105, 0.0155)
            lims2 = (0.031, 0.033)

        elif side == "b":
            lims1 = (0.0009, 0.0070)
            lims2 = (0.008, 0.011)

    elif element == "dihedral":
        return
        if side == "a":
            lims1 = (0.0050, 0.0100)
            lims2 = (0.020, 0.025)

        elif side == "b":
            lims1 = (0.011, 0.016)
            lims2 = (0.018, 0.023)
    else:
        return
        quit()

    # fit two specific gaussians
    # set up plot
    plt.yscale('log')
    plt.title(f"{channel.upper()} for " + element + " (" + side.upper() + " side)")
    plt.xlabel("High - Low Detector Temperature (K)")
    plt.ylabel("Density")

    diff = hi - lo

    # plot the full distribution
    plt.hist(diff, bins=np.linspace(-0.05, 0.05, 1000), density=True)

    # cut and plot the biggest gaussian
    if ngaussians > 0:
        diff_g1 = diff[(diff > lims1[0]) & (diff < lims1[1])]
        plt.hist(diff_g1, bins=np.linspace(-0.05, 0.05, 1000), density=True)
    
    if ngaussians > 1:
        diff_g2 = diff[(diff > lims2[0]) & (diff < lims2[1])]
        plt.hist(diff_g2, bins=np.linspace(-0.05, 0.05, 1000), density=True)
    
    # set axis from the full distribution
    xmin, xmax = plt.xlim()
    # xmin, xmax = 0.01, 0.015
    x = np.linspace(xmin, xmax, 100)
    plt.ylim(1, None)

    # fit the normal distribution to the data
    if ngaussians > 0:
        mu, std = norm.fit(diff_g1, loc=0.015, scale=0.005)

        # plot the histogram and the fitted Gaussian
        p = norm.pdf(x, mu, std)
        plt.plot(x, p, ls='dashed', label='Gaussian 1')
        print(f"Fitted Gaussian 1 for {element} {side} {channel}: mu={mu:.04f}, std={std:.04f}")
        explained = len(diff[(diff > mu - 5*std) & (diff < mu + 5*std)]) / len(diff)
        print(f"Explained fraction by Gaussian 1 (5 std): {explained*100:.02f}%")

    if ngaussians > 1:
        mu2, std2 = norm.fit(diff_g2, loc=0.02, scale=0.005)
        p2 = norm.pdf(x, mu2, std2)
        plt.plot(x, p2, ls='dashed', label='Gaussian 2')
        print(f"Fitted Gaussian 2 for {element} {side} {channel}: mu={mu2:.04f}, std={std2:.04f}")
        explained2 = len(diff[(diff > mu2 - 5*std2) & (diff < mu2 + 5*std2)]) / len(diff)
        print(f"Explained fraction by Gaussian 2 (5 std): {explained2*100:.02f}%")
    
    plt.legend()
    plt.savefig(f"data/output/hilo_temp_difference_gaussian_fit/{element}_{side}_{channel}_{calsky}.png")
    plt.close()

    # plot all temperatures over time for visual inspection
    plt.figure(figsize=(12, 6))
    plt.xlabel("Record Index")
    plt.ylabel("Temperature (K)")
    plt.title(f"Detector Temperatures Over Time for {element} ({side.upper()}) - {channel.upper()}")

    x = np.arange(len(hi))
    
    plt.plot(x, hi, label='High Temp')
    plt.plot(x, lo, label='Low Temp')
    plt.savefig(f"data/output/hilo_temp_timeseries/{element}_{side}_{channel}_{calsky}_g0.png")
    if ngaussians > 0:
        plt.plot(x[(diff > mu - 5*std) & (diff < mu + 5*std)],
                 hi[(diff > mu - 5*std) & (diff < mu + 5*std)], label='Explained by Gaussian 1',
                    linestyle='None', marker='.')
        plt.plot(x[(diff > mu - 5*std) & (diff < mu + 5*std)],
                 lo[(diff > mu - 5*std) & (diff < mu + 5*std)], label='Explained by Gaussian 1',
                    linestyle='None', marker='.')
        plt.legend()
        plt.savefig(f"data/output/hilo_temp_timeseries/{element}_{side}_{channel}_{calsky}_g1.png")
    if ngaussians > 1:
        plt.plot(x[(diff > mu2 - 5*std2) & (diff < mu2 + 5*std2)],
                 hi[(diff > mu2 - 5*std2) & (diff < mu2 + 5*std2)], label='Explained by Gaussian 2',
                    linestyle='None', marker='.')
        plt.plot(x[(diff > mu2 - 5*std2) & (diff < mu2 + 5*std2)],
                 lo[(diff > mu2 - 5*std2) & (diff < mu2 + 5*std2)], label='Explained by Gaussian 2',
                    linestyle='None', marker='.')
        plt.legend()
        plt.savefig(f"data/output/hilo_temp_timeseries/{element}_{side}_{channel}_{calsky}_g2.png")
    
    # plt.savefig(f"data/output/hilo_temp_timeseries/{element}_{side}_{channel}.png")
    plt.close()
    if ngaussians == 1:
        return mu, std
    elif ngaussians == 2:
        return mu, std, mu2, std2
    
def debiase_hi(mu, std, mu2, std2, hi, lo, element, side, channel, calsky):
    """Debiase the high temperature using the fitted Gaussian parameters."""
    diff = hi - lo
    debiased_hi = np.copy(hi)

    # Debias using the first Gaussian
    mask1 = (diff > mu - 5*std) & (diff < mu + 5*std)
    debiased_hi[mask1] = hi[mask1] - mu

    # Debias using the second Gaussian
    mask2 = (diff > mu2 - 5*std2) & (diff < mu2 + 5*std2)
    debiased_hi[mask2] = hi[mask2] - mu2

    # plot before and after
    plt.plot(lo, label='Low Current')
    plt.plot(hi, label='High Current', ls='None', marker='.')
    plt.savefig(f"data/output/hilo_temp_debiasing/{element}_{side}_{channel}_{calsky}_0.png")
    plt.close()

    x = np.arange(len(hi))
    plt.plot(x, lo, label='Low Current')
    plt.plot(x[mask1 | mask2], debiased_hi[mask1 | mask2], label='Debiased High Current', ls='None',
             marker='.', ms=4)
    plt.savefig(f"data/output/hilo_temp_debiasing/{element}_{side}_{channel}_{calsky}_1.png")
    plt.close()

    # plot before and after zoom-in
    plt.plot(lo, label='Low Current')
    plt.plot(hi, label='High Current', ls='None', marker='.', ms=4)
    plt.ylim(2.5, 3)
    plt.savefig(f"data/output/hilo_temp_debiasing/{element}_{side}_{channel}_{calsky}_zoom_0.png")
    plt.close()

    x = np.arange(len(hi))
    plt.plot(x, lo, label='Low Current')
    plt.plot(x[mask1 | mask2], debiased_hi[mask1 | mask2], label='Debiased High Current', ls='None',
             marker='.', ms=4)
    plt.ylim(2.5, 3)
    plt.savefig(f"data/output/hilo_temp_debiasing/{element}_{side}_{channel}_{calsky}_zoom_1.png")
    plt.close()
    return debiased_hi, mask1 | mask2