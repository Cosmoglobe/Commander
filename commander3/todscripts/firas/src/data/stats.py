import numpy as np


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
