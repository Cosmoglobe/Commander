import h5py
import numpy as np
import matplotlib.pyplot as plt
import os

import models.frd as frd
import models.fut as fut

from astropy.modeling.models import BlackBody
from astropy import constants as c
from astropy import units as u

fac_icm_ghz = 29.9792458  # convert cm-1 to GHz
fac_watt_to_mjy = 1.0e15 / fac_icm_ghz
fac_erg_to_mjy = 1.0e8 / fac_icm_ghz
fac_fft_length = [640, 640, 640, 640, 160, 160]  # never used?
fac_spec_length = [321, 321, 321, 321, 81, 81]  # what is this?
fac_nyq_correct = 1.00159


# xcalpos == 1 # external calibrator is in
# xcalpos == 2 # external calibrator is out
# else, it's not legal


# channel numbers:
# ll, lh, rl, rh
# 3   2   1   0


def compute_constants(
    chan: int, scan_mode: int, config, fft_len: int = None, nyquist="default"
):

    # This line is useless?
    fcc_spec_length = fac_spec_length[scan_mode]
    fcc_spec_length = 512

    # Get the nyquist frequencies and frequency intervals in icm and Hz.
    frec = 4 * (chan % 2) + scan_mode
    fnyq_icm = config["nyquistl"]["icm"][frec]
    spec_len = fcc_spec_length - 1.0

    # NOTE: Applying fac_nyq_correct is not done in Fortran code, but I think
    # it is correct based on published data
    if nyquist is None or nyquist == "None":
        print("NO NYQUIST CORRECTION FACTOR")
        nyq_fact = 1
    elif nyquist == "fish":
        print("USING 1.00275 AS NYQUIST CORRECTION FACTOR")
        nyq_fact = 1.00275
    elif nyquist == "default":
        print("USING", fac_nyq_correct, "AS NYQUIST CORRECTION FACTOR")
        nyq_fact = fac_nyq_correct

    df = fnyq_icm / spec_len * nyq_fact
    fnyq_hz = config["nyquistl"]["hz"][frec]
    dw = 2.0 * np.pi * fnyq_hz / spec_len

    # Fill in the frequency arrays
    freq = np.arange(fcc_spec_length) * df
    afreq = np.arange(fcc_spec_length) * dw

    # Calculate the phase shift required to make the high frequency short fast
    # spectra compatable with the high frequency long fast spectra
    if (chan == 0 or chan == 2) and scan_mode == 1:
        phase_shift = np.pi / fft_len * 1j
    else:
        phase_shift = 0.0

    consts = {}
    # Not entirely sure what these constants are?
    consts["df"] = df
    consts["dw"] = dw
    consts["freq"] = freq
    consts["afreq"] = afreq
    consts["phase_shift"] = phase_shift
    # consts['spec_norm'] = spec_norm

    return consts


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


def read_samprate(samprate_type="mission"):
    fn = "/mn/stornext/d5/data/duncanwa/FIRAS/firas_analysis_original/firas_pipeline/reference/fex_samprate.txt"
    data = np.loadtxt(fn, comments="!")

    if samprate_type == "int":
        return data[0]
    else:
        return data[1]


def read_cmdgain():
    data = np.loadtxt(
        "/mn/stornext/d5/data/duncanwa/FIRAS/firas_analysis_original/firas_pipeline/reference/fex_cmdgain.txt",
        comments="C",
    )
    output = np.transpose(data)
    return output


def gen_nyquistl(samprate_type="mission"):
    """This function reads the MTM sampling rates from the samprate text file and the optical
    Nyquist frequency correction from the nyquist fuke and then compute the Nyquist frequencies
    in icm and Hz for all channels and scan modes for either I&T data or flight data
    """
    path = "/mn/stornext/d5/data/duncanwa/FIRAS/firas_analysis_original/firas_pipeline/reference"
    fn_samp = os.path.join(path, "fex_samprate.txt")
    fn_nyquist = os.path.join(path, "fex_nyquist.txt")

    mtm_speed = [0, 1, 0, 1, 0, 1, 0, 1, 1, 1]
    fakeit = [1.0, 2.0, 3.0, 8.0, 12.0]
    multiplier = [6.0, 4.0]
    ngroup = [3.0, 2.0, 3.0, 2.0, 3.0, 2.0, 12.0, 8.0, 8.0, 8.0]

    fringes = 20.00e-4
    optical_path = 4 * np.cos(np.pi / 6)

    fac_fft_length = [640, 640, 640, 640, 160, 160]

    data_samp = np.loadtxt(fn_samp, comments="!")

    if samprate_type == "int":
        sampling_rate = data_samp[0]
    else:
        sampling_rate = data_samp[1]

    freq_shift = np.loadtxt(fn_nyquist, comments="!")

    icm = np.zeros(10)
    hz = np.zeros(15)

    # Calculate the Nyquist frequencies

    # Calculate the optical Nyquist frequencies
    norm = freq_shift / fringes / optical_path / 2.0
    for j in range(10):
        icm[j] = norm * multiplier[mtm_speed[j]] / ngroup[j]

    # Calculate the scan speeds
    speed = sampling_rate / icm[0] / multiplier

    # Calculate the electronic Nyquist frequencies
    for j in range(10):
        hz[j] = icm[j] * speed[mtm_speed[j]]

    # Calculate the fakeit Nyquist frequencies
    for j in range(10, 15):
        hz[j] = fac_fft_length[0] / fakeit[j - 10]

    output = {}
    output["icm"] = icm
    output["hz"] = hz

    return output


if __name__ == "__main__":
    # Gets the MTM sampling rate.
    samprate = read_samprate()

    # Generates apodization functions. Shape is 433x512. 433 is the number of records (what does that mean?) and 512 is the number of data points in one IFG.
    apodl_all = frd.apodl()
    # Creates array with post-detector (?) electronics transfer function, derived from the analytical expressions of the analog and digital filters. Shape is 288x512. 288 is the number of records (what does that mean?) and 512 is the number of data points in one IFG. Found in page 40 of Explanatory Suplement.
    # Where did Duncan get the actual values for the functions?
    etfl_all = frd.elex_transfcnl(samprate, 512)
    print(f"shape of etfl_all: {etfl_all.shape}")
    plt.plot(etfl_all)
    plt.savefig("../output/etfl_all.png")

    # Gain values for the four bolometer preamps. Used to normalise the IFG data.
    gains = read_cmdgain()

    config = {}
    config["samprate"] = samprate
    config["cmdgain"] = gains

    data_orig = h5py.File(
        "/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5"
    )
    data_orig_eng = h5py.File(
        "/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5"
    )

    # Data investigation
    # import pandas as pd

    # for i in range(50):
    #     print(np.sort(pd.unique(data_orig_eng["en_xcal"]["xcal_spares"][:, i])))
    # # print(data_orig_eng["en_xcal"]["pos"].keys())
    # print((data_orig_eng["en_xcal"]["xcal_spares"]).shape)
    # # get how many values are one in fake-it data
    # # print(
    # #     len(data_orig_eng["chan"]["fakeit"][data_orig_eng["chan"]["fakeit"][:, 3] == 1])
    # # )
    # # print(np.sort(pd.unique(data_orig_eng["en_tempdiff"]["bol_assem"])))

    # Get science data from the left long channel
    data_ll = data_orig["fdq_sdf_ll"]
    # Setting the channel number to be 3 (LL)
    channel = 3
    # print(f"Keys of data_ll: {data_ll.keys()}")

    # Picking an IFG to investigate
    ifg_ind = 55_000  # this is a pretty good one

    print(f"xcal_pos of {ifg_ind}: {data_ll["dq_data/xcal_pos"][ifg_ind]}")
    # Result is 1, so theoretically the XCAL is in?

    # 50_000 has some prominent glitches

    # for ifg_ind in [5_000, 55_000, 100_000, 550_000]:
    for ifg_ind in [550_000]:
        print(f"xcal_pos of {ifg_ind}: {data_ll["dq_data/xcal_pos"][ifg_ind]}")
        # Result is 2, so theoretically the XCAL is out?

        t_ifg = data_ll["collect_time/midpoint_time"][ifg_ind]
        print(f"time in engineering data: {data_orig_eng["ct_head/time"][()]}")
        print(np.count_nonzero(np.abs(t_ifg - data_orig_eng["ct_head/time"][()]) == 0))
        # Engineering index is the time of the recorded IFG - all of the times in the engineering data.
        # One value will be zero and we want to retrieve that index.
        # I don't understand how this works because there are less engineering data records than IFG records?
        eng_ind = np.argmin(np.abs(t_ifg - data_orig_eng["ct_head/time"][()]))
        print(f"engineering index: {eng_ind}")

        # Getting the MTM length. There are only two values, 0 and 1. Which one corresponds to which? In this IFG, the value is 0.
        mtm_length = data_ll["sci_head/mtm_length"][ifg_ind]
        print(f"mtm_length: {mtm_length}")
        # Same for here.
        mtm_speed = data_ll["sci_head/mtm_speed"][ifg_ind]
        print(f"mtm_speed: {mtm_speed}")
        # Here I am not sure what this gain refers to or what it affects? In this iFG, it is 4
        gain = data_ll["sci_head/gain"][ifg_ind]
        print(f"gain: {gain}")
        # Not sure what this means? In this IFG, it is 4, out of the possibilities {0, 2, 4}.
        sci_mode = data_ll["sci_head/sc_head1a"][ifg_ind]
        print(f"sci_mode: {sci_mode}")
        # Not sure what this means? In this iFG, it is 3, out of the possibilities {-1, 0, 1, 2, 3}
        adds_per_group = data_ll["sci_head/sc_head9"][ifg_ind]
        print(f"adds_per_group: {adds_per_group}")

        # Not used anywhere else.
        # pts_per_sweep = data_ll["sci_head/sc_head10"][ifg_ind]
        # Not sure what this means? In this IFG it is 16, out of the possibilities: {-21845, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 2052}. Used to get counts to voltage in the IFG.
        sweeps = data_ll["sci_head/sc_head11"][ifg_ind]
        print(f"sweeps: {sweeps}")
        # Actual IFG data
        ifg = data_ll["ifg_data/ifg"][ifg_ind]
        # Check if it is in fake-it mode. Here, it is 0, so it should be in normal mode, assuming 0 is off.
        fake_it = data_ll["dq_data/fake"][ifg_ind]
        print(f"fake_it: {fake_it}")

        xcal_pos = data_ll["dq_data/xcal_pos"][ifg_ind]

        # temperature readings
        # Temperature of the ICAL.
        # Not sure what this temperature refers to, because they mention in the Explanatory Supplement that they use 90% cone and 10% tip temperature (page 124)?
        T_ical = data_orig_eng["en_analog/grt/a_lo_ical"][ifg_ind, 0]
        print(f"T_ical: {T_ical}")
        # Temperature of the XCAL tip. In the Explanatory Supplement they mention that they don't take this into account.
        # For this sample the XCAL is out, so this temperature is not used in this specific sample.
        T_xcal_tip = data_orig_eng["en_analog/grt/a_lo_xcal_tip"][ifg_ind, 0]
        print(f"T_xcal_tip: {T_xcal_tip}")
        # Temeprature of the XCAL cone. In the Explanatory Supplement they mention they take this to be the temperature of the XCAL.
        # For this sample the XCAL is out, so this temperature is not used in this specific sample.
        T_xcal_cone = data_orig_eng["en_analog/grt/a_lo_xcal_cone"][ifg_ind, 0]
        print(f"T_xcal_cone: {T_xcal_cone}")
        # Temperature of the mirror. Not used anywhere else.
    #     T_mirror = data_orig_eng["en_analog/grt/a_lo_mirror"][ifg_ind, 0]
        # Temperature of the collimator. Not used anywhere else.
    #     T_collimator = data_orig_eng["en_analog/grt/a_lo_collimator"][ifg_ind, 0]
        # Temperature of the bolometers. Not used anywhere else.
    #     T_bolometers = data_orig_eng["en_analog/grt/a_lo_bol_assem"][
    #         ifg_ind, :
    #     ]  # 4 detectors
        # Temperature of the sky horn. Used as observation temperature if the XCAL is out. Used to get SI units of the intensity.
        T_sky = data_orig_eng["en_analog/grt/a_lo_skyhorn"][ifg_ind, 0]
        print(f"T_sky: {T_sky}")
        # Temperature of the reference horn. Not used anywhere else.
    #     T_refhorn = data_orig_eng["en_analog/grt/a_lo_refhorn"][ifg_ind, 0]
        # Temperature of the dihedral. Not used anywhere else.
    #     T_dihedral = data_orig_eng["en_analog/grt/a_lo_dihedral"][ifg_ind, 0]

    #     print("xcal_pos")
    #     print(xcal_pos)

    #     print("T_ical, T_xcal_tip, T_xcal_cone, T_sky, T_refhorn")
    #     print(T_ical, T_xcal_tip, T_xcal_cone, T_sky, T_refhorn)

        if xcal_pos == 1:
            T_obs = T_xcal_tip
        elif xcal_pos == 2:
            T_obs = T_sky
        else:
            print("Something weird happening.")
            # What is this?
            # asdf

        # Gain of what? Used to get the counts to voltage in the IFG.
        gain = gains[channel, gain]
        print (f"gain: {gain}")

        # Mode used to get the record number of the ETF and the apodization function. This line is not actually needed.
        upmode = sci_mode
        # Not sure what this is? This line is not actually needed.
        ngpetf = adds_per_group

        # Generates the Nyquist frequencies for all channels and scan modes in cm^-1 and Hz. 10 values for cm^-1 and 15 for hz.
        nyquistl = gen_nyquistl()
        print(f"shape of nyquistl in icm: {nyquistl['icm'].shape}")
        print(f"shape of nyquistl in hz: {nyquistl['hz'].shape}")
        # Not sure what this is? Used to get the index for the correct Nyquist frequency. Value is 0.
        sm = 2 * mtm_length + mtm_speed
        print(f"sm: {sm}")
        # Index for the Nyquist frequency. Value is 4.
        index = 4 * (channel % 2) + sm
        print(f"index for nyquist frequency: {index}")
    #     # print(nyquistl['icm'][index])
        # Nyquist frequency in cm^-1.
        fnyq_icm = nyquistl["icm"][index]  # cm-1
        print(f"fnyq_icm: {fnyq_icm}")
        # Etendu factor.
        fac_etendu = 1.5  # cm2 sr
        # Analog-to-digital conversion scale factor.
        fac_adc_scale = 204.75
        # Normalization factor for the spectra.
        spec_norm = fnyq_icm * fac_etendu * fac_adc_scale
        print(f"spec_norm: {spec_norm}")
        if channel < 2:
            spec_norm *= -1
        # The signs of the spectra for the right
        # side channels is flipped so the right side transfer functions will be positive

        # Copying the nyquistl dictionary to the config dictionary.
        config["nyquistl"] = nyquistl

        consts = compute_constants(channel, sm, config)

        print(consts.keys())

        # ETF record number
        erecno = get_recnum(fake_it, mtm_speed, channel, upmode, ngpetf)
        print(f"erecno: {erecno}")
        # Apodization record number
        arecno = apod_recnuml(channel, sm, fake_it, upmode, adds_per_group, 0)
        print(f"arecno: {arecno}")
        etf = etfl_all[erecno, :]
        apod = apodl_all[arecno, :]
        plt.figure()
        plt.plot(etfl_all[erecno,:].real, label="real")
        plt.plot(etfl_all[erecno,:].imag, label="imag")
        plt.title("ETF")
        plt.ylabel("ETF")
        plt.legend(loc="best")
        plt.savefig("../output/etf.png")
        plt.clf()

    #     # print(mtm_length)
    #     # print(mtm_speed)
    #     # print(sweeps)
    #     # print(pts_per_sweep)
    #     # print(adds_per_group)
    #     # print(gain)

        # Speed of what? Not used anywhere else.
        # speed = 0.8  # cm/s
        # Length of what? Total scan length (page 14 of the Explanatory Supplement)
        length = 1.76  # cm
        # Distance between data points. Not used anywhere else.
        # f_samp = 681.43  # Hz
        # Length differencial? 512 for the 512 points in each interferogram.
        dx = length / 512 * u.cm
        # Frequency differencial?
        dGHz = 17 * u.GHz

        x = np.arange(512) * dx
        # I am guessing that 0 corresponds to short/slow and 1 corresponds to long/fast.
        # The status is not used. This function gets the index of the zero path difference in the IFG. (I thought this was always in the same place? Page 37 of the Explanatory Supplement has table with the peak positions for each of the modes.)
        # peak_pos, status = fut.default_peak(
        peak_pos, _ = fut.default_peak(
            mtm_speed, mtm_length, channel, adds_per_group, sci_mode, 0
        )
        print(f"peak_pos: {peak_pos}")
        # print(peak_pos)

        # Dither subtraction. See Section 4.2 of Fixsen 1994
        # Median of the IFG values. Why are the first 20 values not used?
        med = np.median(ifg[20:]).astype("int16")
        print(f"median: {med}")
        ifg -= med
        plt.figure(1)
        plt.clf()
        plt.plot(x, ifg)
        plt.xlabel("Distance [cm]")
        plt.ylabel("Dither-subtracted IFG [counts]")
        plt.title("Dither-subtracted IFG")
        plt.savefig("../output/ifg.png")
        plt.clf()

        plt.figure(2)
        plt.plot(x, ifg / sweeps / gain, label="Non-apodized")
        plt.plot(x, ifg / sweeps / gain * np.roll(apod, peak_pos), label="Apodized")
        plt.xlabel("Distance [cm]")
        plt.ylabel("IFG [volts]")
        plt.title("IFG, pre and post-apodized and rolled")
        plt.legend(loc="best")
        plt.savefig("../output/ifg_apod.png")
        plt.clf()
        
        # The apodization function is designed to be one at the IFG peak and to fall to zero at the ends (page 36 of the Explanatory Supplement).
        plt.figure(3)
        # rfft = real fast fourier transform
        plt.plot(np.fft.rfft(ifg / sweeps / gain), label="Non-apodized")
        plt.plot(np.fft.rfft(ifg / sweeps / gain * np.roll(apod, peak_pos)), label="Apodized")
        plt.xlabel(r"Frequency [$\mathrm{cm^{-1}}$]")
        plt.ylabel("Spectra [V cm]")
        plt.title("Real FFT")
        plt.legend(loc="best")
        plt.savefig("../output/fft.png")
        plt.clf()

    ifg = np.roll(ifg, -peak_pos)
    spec_nonapod = np.fft.rfft(ifg / sweeps / gain)
    spec_apod = np.fft.rfft(ifg / sweeps / gain * apod)
    plt.figure()
    plt.plot(spec_nonapod.real / spec_norm, label="Non-apodized")
    plt.plot(spec_apod.real / spec_norm, label="Apodized")
    plt.legend(loc="best")
    plt.xlabel(r"Frequency [$\mathrm{cm^{-1}}$]")
    plt.ylabel("Spectra")
    plt.title("Spectra")
    plt.savefig("../output/spectra.png")
    plt.clf()

    plt.figure()
    plt.plot(
        spec_apod.real / spec_norm / etf[: len(spec_apod)] * fac_erg_to_mjy,
        label="ETF-removed",
    )
    plt.legend(loc="best")
    plt.xlabel(r"Frequency [$\mathrm{cm^{-1}}$]")
    plt.ylabel("Spectra")
    plt.title("Spectra with ETF removed")
    plt.savefig("../output/spectra_etf_removed.png")

    nu = np.arange(1, 513) * dGHz

    bb_sky = BlackBody(temperature=T_obs * u.K)
    bb_ical = BlackBody(temperature=T_ical * u.K)
    I_sky = bb_sky(nu)  # .to('erg s-1 cm-2 sr-1 Hz-1') 
    I_sky = (I_sky * c.c).to("erg s-1 cm-2 sr-1 cm")
    I_ical = (bb_ical(nu) * c.c).to("erg s-1 cm-2 sr-1 cm")

    plt.figure()
    plt.plot(nu, I_sky, label="Sky")
    plt.plot(nu, I_ical, label="ICAL")
    plt.legend(loc="best")
    plt.xlabel("Frequency [GHz]")
    plt.ylabel("Intensity [erg s-1 cm-2 sr-1 cm]")
    plt.title("Blackbody intensities")
    plt.tight_layout()
    plt.savefig("../output/blackbody.png")
    plt.clf()


    # # emiss = model['emissivity']

    # #  bol_cmd_bias = cspec_recs['coad_spec_data']['bol_cmd_bias']
    # #  bol_cmd_bias[bol_cmd_bias < 0] += 256
    # #  cmd_bias = np.double(bol_cmd_bias) / 25.5
    # #  tdets = temps[:, chan+6]

    # #  bols, qrad, S0, tau = calc_responsivity(bol_volts, cmd_bias, tdets, model['bolparm'])

    # # B = 1.0 + consts['afreq'][None, :]*tau[:, None]*1j
    # # B[:, 0] = 0.0

    # #  emiss = model['emissivity']

    # #  denom = np.outer(S0, emiss[0, idxs])

    # #  S0 is of order -1.5 to -2

    # #  emissivities are something like 0.01 to 0.07, but also complex.

    S0 = -2  # volt/erg
    emiss = 0.05
    tau = 0.04 * u.s
    denom = S0 * emiss

    # What is B?
    # # B = 1 + 2*np.pi*(nu.value*tau.value)*1j
    B = 1 + consts["afreq"] * tau.value * 1j
    plt.plot(np.arange(1, 513), B.imag)
    plt.ylabel("Imaginary part of B")
    plt.title("Imaginary part of B")
    plt.savefig("../output/B_imag.png")
    plt.clf()
    # # B[0] = 0

    # # print(B)

    nu_in_icm = (nu / c.c).to("cm-1")

    b_x = np.zeros(512, dtype="complex")
    # Not used anywhere else.
    # b_x_ZB = np.zeros(512, dtype="complex")
    x0 = np.diff(x)[0].value * peak_pos
    # Why is the ICAL intensity multiplied by 0.95?
    I_tot = I_sky - I_ical * 0.95
    # Why do we have to cut the ETF to the length of the spectra?
    I_tot = I_tot * etf[:512] * denom / B
    for i in range(len(x)):
        b_x[i] = np.sum(
            I_tot * np.cos(2 * np.pi * (x[i].value - x0) * nu_in_icm.value)
        ).value
    # # b_x = np.fft.irfft(

    plt.plot(x, b_x.real)
    plt.xlabel("Distance [cm]")
    plt.ylabel("Real part of b_x")
    plt.title("Real part of b_x")
    plt.savefig("../output/real_b_x.png")
    plt.clf()

    plt.figure()
    # Gets the angle of the complex part.
    phase = np.angle(b_x)
    # # plt.figure()
    plt.plot(x, phase)
    plt.xlabel("Distance [cm]")
    plt.ylabel("Phase")
    plt.title("Phase")
    plt.savefig("../output/phase.png")
    plt.clf()
    phi = -1.26
    plt.figure()
    plt.plot(x, b_x.real * sweeps * gain * spec_norm, label="Real")
    plt.plot(x, b_x.imag * sweeps * gain * spec_norm, label="Imaginary")
    plt.legend(loc="best")
    plt.xlabel("Distance [cm]")
    plt.ylabel("?")
    plt.title("?")
    plt.savefig("../output/real_imag.png")
    plt.clf()
    plt.figure()
    plt.plot(x, (np.exp(-1j*phi)*b_x).real * gain * spec_norm)
    plt.xlabel("Distance [cm]")
    plt.ylabel("?")
    plt.title("?")
    plt.savefig("../output/real.png")
    plt.clf()

    # plt.show()
