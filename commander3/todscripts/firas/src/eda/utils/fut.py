# FUT is FIRAS Utilities which contains all routines common between pipeline facilities

from __future__ import print_function

import warnings

import numpy as np

import utils.data_types as dt

# This ignores some warnings in the temperature calculation code
# Could probably figure out some indexing to remove the warnings
warnings.simplefilter("ignore")


def temp_list2(rawtemps, rawsigs, grtwt, grttrans, combswitch, singlifg):
    """Calculates temperatures useful in linearization and calibration given the Engineering
    Analog record
    """

    # Reduce the list of 64 GRT readings (A and B sides and high and low current readings for
    # 12 GRTs + 4 cal resistors) down to a short list of 10 "definitive" temperatures. For each
    # object (e.g., XCAL) form a weighted average of the high and low current readings for each
    # GRT to get a reliable reading for that GRT, then form a weighted average of all the GRTs
    # on each object to get a definitive temperature for the object

    # The routine uses the coaddition of sigmas a part of the averaging weights. It handles individual readings
    # as well as multiple ones from coadds. It is also capable of combining the high/low currents alone
    # without averaging the GRTs themselves

    # combswitch - 0=combine A&B, 1=output A only, 2=output B only

    # outtemps - xcal, ical, skyhorn, refhorn, dihedral, structure temp (mirror + collimator),
    #           bol rh, bol rl, bol lh, bol ll

    inout = [0, 3, 1, 2, 4, 9, 5, 6, 7, 8]

    coeff = 1.070e-6

    if rawtemps.ndim == 1:
        nrecords = len(rawtemps)
    else:
        nrecords = 1

    # if combswitch < 0 or combswitch > 2:
    if not 0 <= combswitch <= 2:
        raise ValueError("combswitch must be 0, 1, 2:", combswitch)

    # Extract GRT weights from the input GRTWT record
    # Extract GRT transition and half-width temperatures from the input GRTTRANS record
    gwta = np.tile(np.copy(grtwt["grt_a_weight"]), [nrecords, 1])
    gwtb = np.tile(np.copy(grtwt["grt_b_weight"]), [nrecords, 1])
    trantempa = grttrans["grt_a_trans_temp"]
    trantempb = grttrans["grt_b_trans_temp"]
    tranhwa = grttrans["grt_a_trans_hwid"]
    tranhwb = grttrans["grt_b_trans_hwid"]

    outtemps = np.empty([nrecords, 10])
    outsigs = np.empty([nrecords, 10])

    offgrts = -1 * np.ones(nrecords, dtype=int)

    scratsigs = np.copy(rawsigs)

    # Set the 8 hi-curr bols to infinity (hi_bol_assem) (idx = 5:9 of hi_grt)
    scratsigs["sig_a_hi_bol_assem"] = np.inf
    scratsigs["sig_b_hi_bol_assem"] = np.inf

    # Set any positive sigmas that are too small to a minimum value
    idx = scratsigs["sig_grt"] > 0
    scratsigs["sig_grt"][idx] = np.max(
        [scratsigs["sig_grt"][idx], coeff * rawtemps["grt"][idx]], axis=0
    )

    # A-side temps unless "B side only" was requested, combine Hi and Lo readings of each A-side GRT
    if combswitch != 2:
        tlo = rawtemps["a_lo_grt"]
        thi = rawtemps["a_hi_grt"]
        slo = scratsigs["sig_a_lo_grt"]
        shi = scratsigs["sig_a_hi_grt"]
        tlo[np.isnan(tlo)] = -1
        thi[np.isnan(thi)] = -1
        wt = combine_hilo(thi, tlo, trantempa, tranhwa)
        tempa = (1.0 - wt) * tlo + wt * thi
        siga = np.sqrt(((1.0 - wt) * slo) ** 2 + (wt * shi) ** 2)
        if nrecords == 1:
            tempa.shape = (1, 16)
            siga.shape = (1, 16)
            wt.shape = (1, 16)
        idx = wt == -9999
        tempa[idx] = -9999
        siga[idx] = -9999
        gwta[idx] = 0

    # B-side temps unless "A side only" was requested, combine Hi and Lo readings of each B-side GRT
    if combswitch != 1:
        tlo = rawtemps["b_lo_grt"]
        thi = rawtemps["b_hi_grt"]
        slo = scratsigs["sig_b_lo_grt"]
        shi = scratsigs["sig_b_hi_grt"]
        tlo[np.isnan(tlo)] = -1
        thi[np.isnan(thi)] = -1
        wt = combine_hilo(thi, tlo, trantempb, tranhwb)
        tempb = (1.0 - wt) * tlo + wt * thi
        sigb = np.sqrt(((1.0 - wt) * slo) ** 2 + (wt * shi) ** 2)
        if nrecords == 1:
            tempb.shape = (1, 16)
            sigb.shape = (1, 16)
            wt.shape = (1, 16)
        idx = wt == -9999
        tempb[idx] = -9999
        sigb[idx] = -9999
        gwtb[idx] = 0

    # If "combine A&B sides" was requested, take the merged current readings and combine A and B
    # sides to get a single definitive temperature for each object. Otherwise, supply just the side
    # requested, A or B

    # Fortran code has the same code replicated 3 times, just with certain values emitted based
    # on the combswitch value. It looks like I can just set a variable to 1 or 0 and then just
    # use the code when combining A&B sides to reduce the number of lines
    if combswitch == 0:
        fa = 1
        fb = 1
    elif combswitch == 1:
        fa = 1
        fb = 0
        tempb = np.zeros_like(tempa)
    elif combswitch == 2:
        fa = 0
        fb = 1
        tempa = np.zeros_like(tempb)

    # xcal temp: out#1 = in#1 + in#15 (except I start at 0)
    sumwts = fa * gwta[:, 0] + fb * gwtb[:, 0] + fa * gwta[:, 14] + fb * gwtb[:, 14]

    idx = sumwts > 0
    outtemps[idx, 0] = (
        fa * gwta[idx, 0] * tempa[idx, 0]
        + fb * gwtb[idx, 0] * tempb[idx, 0]
        + fa * gwta[idx, 14] * tempa[idx, 14]
        + fb * gwtb[idx, 14] * tempb[idx, 14]
    ) / sumwts
    if not singlifg:
        outsigs[idx, 0] = (
            np.sqrt(
                (fa * gwta[idx, 0] * siga[idx, 0]) ** 2
                + (fb * gwtb[idx, 0] * sigb[idx, 0]) ** 2
                + (fa * gwta[idx, 14] * siga[idx, 14]) ** 2
                + (fb * gwtb[idx, 14] * sigb[idx, 14]) ** 2
            )
            / sumwts
        )
    else:
        outsigs[idx, 0] = -9999

    outtemps[~idx, 0] = -9999
    outsigs[~idx, 0] = -9999
    offgrts[~idx] = 0

    # structure temp (out#6) = mirror mounts (in#10, aka folding flats) + collimator (in#16)
    # (B-side collimator is unavailable)
    sumwts = fa * gwta[:, 9] + fb * gwtb[:, 9] + fa * gwta[:, 15]

    idx = sumwts > 0
    outtemps[idx, 5] = (
        fa * gwta[idx, 9] * tempa[idx, 9]
        + fb * gwtb[idx, 9] * tempb[idx, 9]
        + fa * gwta[idx, 15] * tempa[idx, 15]
    ) / sumwts
    if not singlifg:
        outsigs[idx, 5] = (
            np.sqrt(
                (fa * gwta[idx, 9] * siga[idx, 9]) ** 2
                + (fb * gwtb[idx, 9] * sigb[idx, 9]) ** 2
                + (fa * gwta[idx, 15] * siga[idx, 15]) ** 2
            )
            / sumwts
        )
    else:
        outsigs[idx, 5] = -9999
    outtemps[~idx, 5] = -9999
    outsigs[~idx, 5] = -9999
    idx = np.logical_and(offgrts == -1, ~idx)
    offgrts[idx] = 5

    # ical, sky horn, ref horn, dihedral, bolometer RH, RL, LH, LL temps
    for k in range(1, 10):
        j = inout[k]
        if j != 0 and k != 5:
            sumwts = fa * gwta[:, j] + fb * gwtb[:, j]

            idx = sumwts > 0
            outtemps[idx, k] = (
                fa * gwta[idx, j] * tempa[idx, j] + fb * gwtb[idx, j] * tempb[idx, j]
            ) / sumwts

            if not singlifg:
                outsigs[idx, k] = (
                    np.sqrt(
                        (fa * gwta[idx, j] * siga[idx, j]) ** 2
                        + (fb * gwtb[idx, j] * sigb[idx, j]) ** 2
                    )
                    / sumwts
                )
            else:
                outsigs[idx, k] = -9999
            outtemps[~idx, k] = -9999
            outsigs[~idx, k] = -9999

            idx1 = np.logical_or(k < offgrts, offgrts == -1)
            idx = np.logical_and(~idx, idx1)
            offgrts[idx] = k

    # Along with the temperatures and sigs, return the lowest index number in the out-array where
    # a 'flag' value was returned because all the GRT switches in its make-up were off; if none, -1

    if nrecords == 1:
        outtemps.shape = (10,)
        outsigs.shape = (10,)
        offgrts = offgrts[0]

    return outtemps, outsigs, offgrts


def temp_list(rawtemps, rawsigs, grtwt, grttrans, combswitch, singlifg):
    """Calculates temperatures useful in linearization and calibration given the Engineering
    Analog record
    """

    # Reduce the list of 64 GRT readings (A and B sides and high and low current readings for
    # 12 GRTs + 4 cal resistors) down to a short list of 10 "definitive" temperatures. For each
    # object (e.g., XCAL) form a weighted average of the high and low current readings for each
    # GRT to get a reliable reading for that GRT, then form a weighted average of all the GRTs
    # on each object to get a definitive temperature for the object

    # The routine uses the coaddition of sigmas a part of the averaging weights. It handles individual readings
    # as well as multiple ones from coadds. It is also capable of combining the high/low currents alone
    # without averaging the GRTs themselves

    # combswitch - 0=combine A&B, 1=output A only, 2=output B only

    # outtemps - xcal, ical, skyhorn, refhorn, dihedral, structure temp (mirror + collimator),
    #           bol rh, bol rl, bol lh, bol ll

    # inout = [0, 4, 2, 3, 5, 10, 6, 7, 8, 9]
    inout = [0, 3, 1, 2, 4, 9, 5, 6, 7, 8]

    coeff = 1.070e-6

    if combswitch < 0 or combswitch > 2:
        raise ValueError("combswitch must be 0, 1, or 2")

    # Extract GRT weights from the input GRTWT record
    # Extract GRT transition and half-width temperatures from the input GRTTRANS record
    # gwta = np.copy(grtwt.grt_a_weight)
    # grtwtb = np.copy(grtwt.grt_b_weight)
    gwta = np.copy(grtwt["grt_a_weight"])
    gwtb = np.copy(grtwt["grt_b_weight"])
    # trantempa = grttrans.grt_a_trans_temp
    # trantempb = grttrans.grt_b_trans_temp
    # tranhwa = grttrans.grt_a_trans_hwid
    # tranhwb = grttrans.grt_b_trans_hwid
    trantempa = grttrans["grt_a_trans_temp"]
    trantempb = grttrans["grt_b_trans_temp"]
    tranhwa = grttrans["grt_a_trans_hwid"]
    tranhwb = grttrans["grt_b_trans_hwid"]
    # tempb = np.empty(16)
    # sigb = np.empty(16)

    outtemps = np.empty(10)
    outsigs = np.empty(10)

    offgrts = -1

    scratsigs = np.copy(rawsigs)

    # Set the 8 hi-curr bols to infinity (hi_bol_assem) (idx = 5:9 of hi_grt)
    # scratsigs['a_hi_bol_assem'] = np.inf
    # scratsigs['b_hi_bol_assem'] = np.inf
    scratsigs["sig_a_hi_bol_assem"] = np.inf
    scratsigs["sig_b_hi_bol_assem"] = np.inf

    # Set any positive sigmas that are too small to a minimum value
    idx = scratsigs["sig_grt"] > 0
    scratsigs["sig_grt"][idx] = np.max(
        [scratsigs["sig_grt"][idx], coeff * rawtemps["grt"][idx]], axis=0
    )

    # A-side temps unless "B side only" was requested, combine Hi and Lo readings of each A-side GRT
    if combswitch != 2:
        tlo = rawtemps["a_lo_grt"]
        thi = rawtemps["a_hi_grt"]
        slo = scratsigs["sig_a_lo_grt"]
        shi = scratsigs["sig_a_hi_grt"]
        tlo[np.isnan(tlo)] = -1
        thi[np.isnan(thi)] = -1
        wt = combine_hilo(thi, tlo, trantempa, tranhwa)
        tempa = (1.0 - wt) * tlo + wt * thi
        siga = np.sqrt(((1.0 - wt) * slo) ** 2 + (wt * shi) ** 2)
        idx = wt == -9999
        tempa[idx] = -9999
        siga[idx] = -9999
        gwta[idx] = 0

    # B-side temps unless "A side only" was requested, combine Hi and Lo readings of each B-side GRT
    if combswitch != 1:
        tlo = rawtemps["b_lo_grt"]
        thi = rawtemps["b_hi_grt"]
        slo = scratsigs["sig_b_lo_grt"]
        shi = scratsigs["sig_b_hi_grt"]
        tlo[np.isnan(tlo)] = -1
        thi[np.isnan(thi)] = -1
        wt = combine_hilo(thi, tlo, trantempb, tranhwb)
        tempb = (1.0 - wt) * tlo + wt * thi
        sigb = np.sqrt(((1.0 - wt) * slo) ** 2 + (wt * shi) ** 2)
        idx = wt == -9999
        tempb[idx] = -9999
        sigb[idx] = -9999
        gwtb[idx] = 0

    # If "combine A&B sides" was requested, take the merged current readings and combine A and B
    # sides to get a single definitive temperature for each object. Otherwise, supply just the side
    # requested, A or B

    # Fortran code has the same code replicated 3 times, just with certain values emitted based
    # on the combswitch value. It looks like I can just set a variable to 1 or 0 and then just
    # use the code when combining A&B sides to reduce the number of lines
    if combswitch == 0:
        fa = 1
        fb = 1
    elif combswitch == 1:
        fa = 1
        fb = 0
        tempb = np.zeros_like(tempa)
    elif combswitch == 2:
        fa = 0
        fb = 1
        tempa = np.zeros_like(tempb)

    # xcal temp: out#1 = in#1 + in#15 (except I start at 0)
    sumwts = fa * gwta[0] + fb * gwtb[0] + fa * gwta[14] + fb * gwtb[14]

    if sumwts > 0:
        outtemps[0] = (
            fa * gwta[0] * tempa[0]
            + fb * gwtb[0] * tempb[0]
            + fa * gwta[14] * tempa[14]
            + fb * gwtb[14] * tempb[14]
        ) / sumwts
        if not singlifg:
            outsigs[0] = (
                np.sqrt(
                    (fa * gwta[0] * siga[0]) ** 2
                    + (fb * gwtb[0] * sigb[0]) ** 2
                    + (fa * gwta[14] * siga[14]) ** 2
                    + (fb * gwtb[14] * sigb[14]) ** 2
                )
                / sumwts
            )

        else:
            outsigs[0] = -9999
    else:
        outtemps[0] = -9999
        outsigs[0] = -9999
        offgrts = 0

    # structure temp (out#6) = mirror mounts (in#10, aka folding flats) + collimator (in#16)
    # (B-side collimator is unavailable)
    sumwts = fa * gwta[9] + fb * gwtb[0] + fa * gwta[15]

    if sumwts > 0:
        outtemps[5] = (
            fa * gwta[9] * tempa[9]
            + fb * gwtb[9] * tempb[9]
            + fa * gwta[15] * tempa[15]
        ) / sumwts
        if not singlifg:
            outsigs[5] = (
                np.sqrt(
                    (fa * gwta[9] * siga[9]) ** 2
                    + (fb * gwtb[9] * sigb[9]) ** 2
                    + (fa * gwta[15] * siga[15]) ** 2
                )
                / sumwts
            )
        else:
            outsigs[5] = -9999
    else:
        outtemps[5] = -9999
        outsigs[5] = -9999
        if offgrts == 0:
            offgrts = 5

    # ical, sky horn, ref horn, dihedral, bolometer RH, RL, LH, LL temps
    for k in range(1, 10):
        j = inout[k]
        if j != 0 and k != 5:
            sumwts = fa * gwta[j] + fb * gwtb[j]

            if sumwts > 0:
                outtemps[k] = (
                    fa * gwta[j] * tempa[j] + fb * gwtb[j] * tempb[j]
                ) / sumwts

                if not singlifg:
                    outsigs[k] = (
                        np.sqrt(
                            (fa * gwta[j] * siga[j]) ** 2
                            + (fb * gwtb[j] * sigb[j]) ** 2
                        )
                        / sumwts
                    )
                else:
                    outsigs[k] = -9999
            else:
                outtemps[k] = -9999
                outsigs[k] = -9999
                if k < offgrts or offgrts == -1:
                    offgrts = k

    # Along with the temperatures and sigs, return the lowest index number in the out-array where
    # a 'flag' value was returned because all the GRT switches in its make-up were off; if none, -1

    return outtemps, outsigs, offgrts


def combine_hilo(rdghi, rdglo, trantemp, tranhwid, report=False):
    """Function to combine high- and low- current readings of an individual GRT temperature.
    Return weights needed for combining temperatures
    """

    # Test of vectorized version
    is_vect = type(rdghi) == np.ndarray

    if not is_vect:
        rdghi, rdglo = np.array([rdghi]), np.array([rdglo])
        trantemp, tranhwid = np.array([trantemp]), np.array([tranhwid])
    elif rdghi.ndim == 2:
        nrecords = len(rdghi)
        trantemp = np.tile(trantemp, [nrecords, 1])
        tranhwid = np.tile(tranhwid, [nrecords, 1])

    tranlo = trantemp - tranhwid
    tranhi = trantemp + tranhwid

    """
    if rdghi < -999:
        if rdglo < -999:
            wt = -9999
        else:
            wt = 0
    elif rdglo < -999:
        wt = 1.0
    elif rdglo < tranlo :
        wt = 0.0
    elif rdglo > tranhi:
        wt = 1.0
    else:
        wt = 0.5*((rdglo - trantemp)/tranhwid + 1.0)
    """

    wt = np.zeros_like(rdghi)
    idx1 = rdghi < -999
    idx2 = rdglo < -999

    idx = np.logical_and(idx1, idx2)
    wt[idx] = -9999.0
    idx = np.logical_and(idx1, ~idx2)
    wt[idx] = 0.0
    idx = np.logical_and(~idx1, idx2)
    wt[idx] = 1.0
    idxA = np.logical_and(~idx1, ~idx2)
    idx = np.logical_and(idxA, rdglo < tranlo)
    wt[idx] = 0.0
    idxA = np.logical_and(idxA, rdglo >= tranlo)
    idx = np.logical_and(idxA, rdglo > tranhi)
    wt[idx] = 1.0
    idx = np.logical_and(idxA, rdglo <= tranhi)
    wt[idx] = 0.5 * ((rdglo[idx] - trantemp[idx]) / tranhwid[idx] + 1.0)

    return wt


def default_peak(mtmspeed, scanlength, chan, ngroup, sci_mode, linearized):
    """Given the MTM speed and data compression, ths function returns the index number
    of the position of the zero path difference in the interferogram

    Notes
    -----
    A lot more information about the calculation is in fut_default_peak.for
    """

    dx = [1.152e-3, 1.728e-3]
    addon = [18.0, 12.0, 20.0, 14.0, 16.0, 6.0, 8.0, 0.0]
    offset = 1.2

    # Different values for ch based on whether we have a high or low channel
    # Set ch to 1 for high channels, 0 for low channels (RH, RL, LH, LL)
    if chan == 0 or chan == 2:
        ch = 1
    else:
        ch = 0

    # Set science mode/linearization indicator to 1 if digital filters are on OR if linearization has
    # occured, 0 if not. Since the peak position after linearization is not affected by whether the digital
    # filters are on or off, the shift factor is added in regardless when the interferogram is already
    # linearized
    if (sci_mode == 1 or sci_mode == 3) and not linearized:
        sc_m_or_l = 0
    else:
        sc_m_or_l = 1

    # Index of fudge factor is determed by binary arithmetic
    data_indx = ch + 2 * mtmspeed + 4 * scanlength

    # Initial peak location is 2 + distance to peak (1.2 cm) divided by adds per groups times distance per sample
    loc_peak = int(offset / (ngroup * dx[mtmspeed]) + 1)

    # If digital filters are on or if the interferogram has been linearized, shift is the "fudge factor"
    # divided by adds per group times 2/3 for slow speed, 1 for fast speed, plus the intefer c....quotient
    # of 4 divided by adds per group. Shift is 0 if digital filters are off and interferogram is not linearized
    shift = sc_m_or_l * (
        3 * addon[data_indx] / (ngroup * (2.0 + mtmspeed))
    ) + sc_m_or_l * (4 // ngroup)

    # Buffer factor is 3 for short slow scan mode, 2 for short fast scan mode, and 8 for long scan modes
    buffactor = int((scanlength * (5 + mtmspeed) - mtmspeed + 3) // ngroup)

    # Amount to take off if linearized is 2 for long scan length or short slow, high channel. It is
    # 3 for short fast, high channel, 5 for short slow, low channel, and 6 for short fast, low channel
    takeoff = linearized * (2 * scanlength + (1 - scanlength) * (mtmspeed + 5 - 3 * ch))

    # Peak location is the sum of the original plus the shift and buffer factor, less the amount
    # subtracted after linearization
    loc_peak += int(shift + buffactor - takeoff)

    # Peak is good if between 0 and 511, otherwie not good
    if loc_peak >= 0 and loc_peak < 512:
        status = True
    else:
        status = False

    return loc_peak, status


def ave_bin_times(bin_time_array, weights):
    """To add and average quadword binary times (the adt)

    Notes
    -----
    I am following the algorithm used in the FIRAS code where to avoid issues
    with large integers to calculate the average on the differences from the
    first time.
    """

    diffs = [0] + list(bin_time_array[1:] - bin_time_array[0])
    meandiffs = int(round(np.average(diffs, weights=weights)))
    bin_time_avg = bin_time_array[0] + meandiffs

    return bin_time_avg


def apod_rotl(aifg, apod, fft_len, peak):
    """This function apodizes and rotates the coadded interferograms using apodization functions
    read from the binary reference dataset FEX_APODL

    Created for the long spectrum pipeline routines (FIL, FFL, FSL, FCL)
    """

    ndim = aifg.ndim

    if ndim == 1:
        aifg.shape = (1, 512)

    nrecords = len(aifg)

    # Initialize the apodized interfogram and the apodized, rotation interferogram
    if fft_len > 512:
        apifg = np.zeros([nrecords, fft_len])

        # Apodize the interferogram, zero-padding from 513 to 720
        apifg[:, :512] = aifg * apod[None, :]
    else:
        # In this case, the fft_len is smaller than the input 512 array
        apifg = aifg * apod[None, :]
        apifg = apifg[:, :fft_len]

    # Move the peak position to the first point
    rifg = np.roll(apifg, -peak, axis=1)

    if ndim == 1:
        rifg.shape = (fft_len,)
        # rifg.shape = (720, )
        aifg.shape = (512,)

    return rifg


def apod_recnuml(channel, scan_mode, fake_it, sci_mode, adds_per_group, linearized):
    """Return the record number for the appropriate apodization function from the
    direct-access apodization file for revised pipeline facilities (FIL, FSL, FFL, FCL)

    Parameters
    ----------
    channel: int
        Value of channel, 0-3 = RH-LL
    scan_mode: int
        Value of scan mode, 0-5 = SS-FL
    fake_it: bool
        Value of fake-it bit
    sci_mode: int
        Value of science mode, 0-4
    adds_per_group: int
        Value of adds per group, 1-12
    linearized: bool
        Whether apodization function should be for linearized interferograms

    Returns
    -------
    recnum: int
        Record number of appropriate apodization function

    """

    # Determine whether channel is high or low
    if channel == 0 or channel == 2:
        ch = 1
    else:
        ch = 0

    # Determine whether digital filters are on or off
    if sci_mode == 1 or sci_mode == 3:
        sc = 0  # digital filters off
    else:
        sc = 1  # digital filters on

    # Apodization function for fake-it data is the last one (432); functions for FS and FL data are after 383
    if fake_it:
        recnum = 432
    elif scan_mode >= 4:
        recnum = 384 + 4 * (adds_per_group - 1) + 2 * sc + linearized
    else:
        recnum = (
            32 * (adds_per_group - 1)
            + 16 * (scan_mode % 2)
            + 8 * int(scan_mode / 2)
            + 4 * ch
            + 2 * sc
            + linearized
        )

    return recnum


# TODO: Write this like the apodization recnum routine just so the code looks consistent across functions
def get_recnum(fakeit, mtmspeed, ichan, micro_mode, iadds_grp):
    """This routine takes the numerical values associated with the appropriate intstrument
    status bits and returns the corresponding records number in the electronics transfer
    function file. That record holds the transfer function associated with the given
    instrument configuration.

    Parameters
    ----------
    fakeit: bool
        Fakeit pulse status
    mtmspeed: int
        0=slow, 1=fast
    ichan: int
        channel number (0-3 = RH-LL)
    micro_mode: int
        uProc science mode ("block type" 0-4)
    iadds_grp: int
        uP data compression (adds per group)

    Returns
    -------
    irecno: int
        desired record number xfcn file
    """

    # ETF is the same for both mtm speeds
    irecno = 96 * (1 - mtmspeed)
    # irecno = 96 * mtmspeed
    if fakeit:
        irecno = 192
    irecno += 24 * ichan
    if micro_mode == 2 or micro_mode == 4:
        idig_fltr = 0  # digital filter on
    else:
        idig_fltr = 1  # digital filter off

    irecno += 12 * idig_fltr
    irecno += iadds_grp - 1

    if irecno <= 0:
        raise ValueError("Electronics transfer function record number < 0")
    elif irecno >= 288:
        raise ValueError("Electronics transfer functions record number >= 288")

    return irecno


def average_angles(angles):
    """Calculates the average angle from a list of input angles. All
    calculations are performed in radians.
    """

    # Subtract the first value. This will put all the values between +/- 2pi. If they are not
    # between +- pi, get them there by adding or subtracting 2pi.

    corr_angle = angles - angles[0]

    idx = corr_angle > np.pi
    corr_angle[idx] -= 2 * np.pi
    idx = corr_angle < -np.pi
    corr_angle[idx] += 2 * np.pi

    avg_angle = np.mean(corr_angle) + angles[0]

    # If the average is not between +/- pi, add or subtract 2pi
    if avg_angle > np.pi:
        avg_angle -= 2 * np.pi
    if avg_angle < -np.pi:
        avg_angle += 2 * np.pi

    return avg_angle
