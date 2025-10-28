def filter_junk(
    stat_word_1,
    stat_word_5,
    stat_word_9,
    stat_word_12,
    stat_word_13,
    stat_word_16,
    lvdt_stat_a,
    lvdt_stat_b,
):
    """
    Sets the filters needed to remove junk data according to flags.

    ---
    Parameters

    stat_word_1 : np.ndarray

    stat_word_5 : np.ndarray

    stat_word_9 : np.ndarray

    stat_word_12 : np.ndarray

    stat_word_13 : np.ndarray

    stat_word_16 : np.ndarray

    lvdt_stat_a : np.ndarray

    lvdt_stat_b : np.ndarray
    ---
    Returns
    full_filter : np.ndarray
    """

    filter0 = stat_word_1 != 46  # couldn't find in cal ifgs
    filter1 = (
        (stat_word_5 != 16641) & (stat_word_5 != 17921) & (stat_word_5 != 17217)
    )  # took out because looking at cal ifgs it looks fine
    # & (stat_word_5 != 19457) # makes the hole!
    filter2 = (stat_word_9 != 15414) & (stat_word_9 != 45110)
    filter25 = (
        (
            stat_word_12 != 6714
        )  # hot horn season - ifgs look fine, problem is calibration
        & (
            stat_word_12 != 7866
        )  # hot horn season - ifgs look fine, problem is calibration
        & (stat_word_12 != 18536)
        & (stat_word_12 != 19121)
        & (stat_word_12 != 54906)
        & (stat_word_12 != 63675)
    )
    filter3 = (
        (stat_word_13 != 17281)
        & (stat_word_13 != 17345)
        & (stat_word_13 != 17393)
        & (stat_word_13 != 19585)
        & (stat_word_13 != 23681)  # hot horn season
        & (stat_word_13 != 25153)
        & (stat_word_13 != 25201)
        & (stat_word_13 != 25585)
        & (stat_word_13 != 26945)
        & (stat_word_13 != 27009)
        & (stat_word_13 != 27073)
        & (stat_word_13 != 27649)
        & (stat_word_13 != 27697)
        # & (stat_word_13 != 27777) # makes the hole!
        # & (stat_word_13 != 27825) # takes away a lot of data!
        & (stat_word_13 != 60465)
        & (stat_word_13 != 60545)
        & (stat_word_13 != 60593)
    )
    filter4 = (
        (stat_word_16 != 0)
        & (stat_word_16 != 14372)
        & (stat_word_16 != 35392)  # hot horn season
        & (stat_word_16 != 35584)
        & (stat_word_16 != 35642)  # hot horn season
        & (stat_word_16 != 36032)
        & (stat_word_16 != 52992)
        & (stat_word_16 != 53056)
    )
    filter5 = (
        (lvdt_stat_a != 1)
        & (lvdt_stat_a != 2)
        & (lvdt_stat_a != 4)
        & (lvdt_stat_a != 6)
        & (lvdt_stat_a != 13)
        & (lvdt_stat_a != 17)
        & (lvdt_stat_a != 21)
        & (lvdt_stat_a != 24)
        & (lvdt_stat_a != 26)
        & (lvdt_stat_a != 31)
        & (lvdt_stat_a != 32)
        & (lvdt_stat_a != 35)
        & (lvdt_stat_a != 47)
        & (lvdt_stat_a != 49)
    )
    filter6 = (
        (lvdt_stat_b != -127)
        & (lvdt_stat_b != -124)
        & (lvdt_stat_b != -108)
        & (lvdt_stat_b != -83)
        & (lvdt_stat_b != -79)
        & (lvdt_stat_b != -78)
        & (lvdt_stat_b != -71)
        & (lvdt_stat_b != -70)
        & (lvdt_stat_b != 75)
        & (lvdt_stat_b != 79)
        & (lvdt_stat_b != 82)
        & (lvdt_stat_b != 83)
        & (lvdt_stat_b != 85)
        & (lvdt_stat_b != 86)
        & (lvdt_stat_b != 87)
        & (lvdt_stat_b != 91)
        & (lvdt_stat_b != 93)
        & (lvdt_stat_b != 95)
        & (lvdt_stat_b != 96)
        & (lvdt_stat_b != 99)
        & (lvdt_stat_b != 100)
        & (lvdt_stat_b != 103)
        & (lvdt_stat_b != 107)
        & (lvdt_stat_b != 111)
        & (lvdt_stat_b != 121)
    )

    full_filter = (
        filter0 & filter1 & filter2 & filter25 & filter3 & filter4 & filter5 & filter6
    )

    return full_filter


def flag(stat_word_9, lvdt_stat_b):
    """
    This function returns a mask with the bad flags based on looking at the calibration interferograms.
    """

    sw9_filter = stat_word_9 != 15417
    print(
        f"sw9_filter: Removing {(sw9_filter.size- sw9_filter.sum())/sw9_filter.size * 100} % of the data"
    )
    lsb_filter = lvdt_stat_b != 96
    print(
        f"lsb_filter: Removing {(lsb_filter.size - lsb_filter.sum()) / lsb_filter.size * 100} % of the data"
    )

    return sw9_filter & lsb_filter


def filter_bol(
    a_bol_assem_rh,
    a_bol_assem_rl,
    a_bol_assem_lh,
    a_bol_assem_ll,
    b_bol_assem_rh,
    b_bol_assem_rl,
    b_bol_assem_lh,
    b_bol_assem_ll,
    bol_cmd_bias_rh,
    bol_cmd_bias_rl,
    bol_cmd_bias_lh,
    bol_cmd_bias_ll,
):
    filtergrt = (
        (a_bol_assem_rh > 0)
        & (a_bol_assem_rl > 0)
        & (a_bol_assem_lh > 0)
        & (a_bol_assem_ll > 0)
        & (b_bol_assem_rh > 0)
        & (b_bol_assem_rl > 0)
        & (b_bol_assem_lh > 0)
        & (b_bol_assem_ll > 0)
    )
    print(
        f"filtergrt: Removing {(filtergrt.size - filtergrt.sum()) / filtergrt.size * 100} % of the data"
    )
    filtercmd = (
        (bol_cmd_bias_rh > 0)
        & (bol_cmd_bias_rl > 0)
        & (bol_cmd_bias_lh > 0)
        & (bol_cmd_bias_ll > 0)
    )
    print(
        f"filtercmd: Removing {(filtercmd.size - filtercmd.sum()) / filtercmd.size * 100} % of the data"
    )

    return filtergrt & filtercmd


def flag_bad_ifgs():
    """
    These are the IDs given by the pre-processing script so it might change depending on if we change the way of associating IFGs in the preprocessing.
    """
    # idx_cal = [3490]
    # return idx_cal
    pass
