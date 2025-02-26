# ================================================================================
#
# Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
#
# This file is part of Commander3.
#
# Commander3 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Commander3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Commander3. If not, see <https://www.gnu.org/licenses/>.
#
# ================================================================================

import numpy as np
from astropy.io import fits
from astropy.time import Time
from scipy.io import readsav
from scipy.interpolate import BSpline, interp1d

import matplotlib.pyplot as plt

from glob import glob

prefix = "/mn/stornext/d16/cmbco/ola/wmap/"
from tqdm import tqdm

"""
The gain model as described by the explanatory supplement requires the total
radiometer bias power, T_FPA, and T_RXB as inputs.

These quantities are indexed by mnemonics, enumerated in table C.2.2, C.2.4, and
C.2.5, respectively. I have listed them in the files t_rxb_mnem.txt,
t_fpa_mnem.txt, and rfbias.txt

AIHK_mnemonic.pro lists the conversion between mnemonic, and the location in the
fits table ("index" and "arr").


AIHK_Arch2Mnemonic.pro returns the physical value associated with an analog
    instrument housekeeping mnemonic, extracting the data from an AEU sweep.
AIHK_GetMnemonic.pro returns the physical value associated with a mnemonic,
    extracting the data out of a sweep of analog instrument housekeeping
    telemetry data.
AIHK_Pckt2Mnemonic.pro returns the physical value associated with an analog
    housekeeping mnemonic, extracting the data from a DEU telemetry sweep.


I think AIHK_Arch2Mnemonic is the thing closest to what I want, since I want to
read in a TOD and just have that data available to me.


in the wmap pro file, execute
file = '/mn/stornext/d16/cmbco/ola/wmap/tods/uncalibrated/wmap_tod_20052022356_20052032356_uncalibrated_v5.fits'
fits_read_tod, file, arch
.r aihk_arch2mnemonic
rf_bias = AIHK_Arch2Mnemonic(arch, 'DRV113RFBI32', 0)
T_RXB = AIHK_Arch2Mnemonic(arch, 'DRV111RXBAMPT', 0)
T_FPA = AIHK_Arch2Mnemonic(arch, 'DFV11FPATEET', 0)



Check
aihk_mnem2serial
get_prt_temp


Comparing against the ExSupp values, it seems like decent approximations to the
reported temperature FPA are
DFQ1AFEEDT
DFW3BFEEDT
DFV22FPATEET
DFW3AOMTT
(although they are all scattering around the full range...)

while the RXB good approximations are...
DRW221RXBAMPT
DRQ2RXBRIBT


I suspect that I have some residual issues here, since when I run the IDL
scripts, all of the values are within roughly 1 degree of each other.
"""


mnem_list = np.loadtxt("gain_params/mnem_indices.txt", delimiter=",", dtype=str)

# Temp_FPA = np.loadtxt('gain_params/fpa_temps.txt')[:,8]
# Temp_RXB = np.loadtxt('gain_params/rxb_temps.txt')[:,1] # Compare to the ExSupp
Temp_FPAs = np.loadtxt("gain_params/fpa_temps.txt")
Temp_RXBs = np.loadtxt("gain_params/rxb_temps.txt")
rng = np.random.default_rng()

Temp_FPA = Temp_FPAs.mean(axis=1)
Temp_RXB = Temp_RXBs.mean(axis=1)
# Temp_FPA = Temp_FPAs[:,6]
# Temp_RXB = Temp_RXBs[:,0]
# Temp_FPA = rng.choice(Temp_FPAs, axis=1)
# Temp_RXB = rng.choice(Temp_RXBs, axis=1)

fnames = glob("/mn/stornext/d16/cmbco/ola/wmap/tods/uncalibrated/*.fits")
fnames.sort()
t1s = []
t2s = []
for f in fnames:
    t1 = f.split("_")[2]
    t2 = f.split("_")[3]
    t1 = t1[:4] + ":" + t1[4:7] + ":" + t1[7:9] + ":" + t1[9:11]
    t2 = t2[:4] + ":" + t2[4:7] + ":" + t2[7:9] + ":" + t2[9:11]
    t1 = Time(t1, format="yday")
    t2 = Time(t2, format="yday")
    t1s.append(t1.jd)
    t2s.append(t2.jd)

t1s = np.array(t1s)
t2s = np.array(t2s)
times = (t1s + t2s) / 2

inds = Temp_FPA > 50

fpa_func = interp1d(times[inds], Temp_FPA[inds], fill_value="extrapolate")
rxb_func = interp1d(times[inds], Temp_RXB[inds], fill_value="extrapolate")


def aihk_mnemonic(mnem):
    for i in range(len(mnem_list)):
        if mnem in mnem_list[i, 0]:
            ind = i
    index, arr = mnem_list[ind, 1:]
    # index - 1 because the indexing is from Fortran
    return int(index) - 1, int(arr)


def get_resist(data, index, arr):
    Win = 0
    Wflg = 0
    if arr == 3:
        Tmp = data[3].data["AEU2"][:, index]
        if (index >= 1) and (index <= 31):
            Wflg = 1
            Win = data[3].data["AWIN2"][:, index // 2]
            # This is some kind of hexadecimal code, but it seems like it's a
            # correction on the order of 30 K. For the sake of my sanity, I'm going
            # to set terms involving Win to zero.
            if index % 2 == 0:
                Win = np.array([w & 0x00FF for w in Win])
            else:
                Win = np.array([(w & 0xFF00) * 2**-8 for w in Win])
            Slp = 254.968244e-6
            YInt = 319.5004
            WInt = 129
            WSlp = 256.0
            Rmax = 650.25838
            Rmax = 0
    elif arr == 2:
        Tmp = data[3].data["AEU1"][:, index]
        # For AEU1:
        if (index >= 1) and (index <= 31):
            Wflg = 2
            Win = data[3].data["AWIN1"][:, index // 2]
            if index % 2 == 0:
                Win = np.array([w & 0x00FF for w in Win])
            else:
                Win = np.array([(w & 0xFF00) * 2**-8 for w in Win])
            Slp = 255.381467e-6
            YInt = 319.5197
            WInt = 129
            WSlp = 256.0
            Rmax = 650.58226
    else:
        Tmp = data[3].data["PDU"][index]
    if Wflg == 0:
        return Tmp, Wflg
    else:
        Res = Tmp * Slp + YInt + Rmax * (Win - WInt) / WSlp
        return Res, Wflg


def mnem2serial(mnem):
    # aihk_mnem2serial_list, info
    mnem_list = np.loadtxt("gain_params/mnem_serial2list.txt", dtype="str")
    ind = np.where(mnem_list[:, 0] == mnem)[0]
    return mnem_list[ind][0][1]


def get_temp(Res, mnem, SerialNum="UG89", SensorID="IHK1_AC_17"):
    # from aihk_mnem2serial_list.pro
    SerialNum = mnem2serial(mnem)
    # from get_prt_temp.pro, converts resistance in ohms to temperature in Kelvin.
    # This requires using the prt_splinedata.xdr file in the ref directory.A

    Temp = get_prt_temp(Res, SerialNum)
    return Temp


def get_prt_temp(resist, SerialNum):
    table = readsav(f"{prefix}/software/ref/prt_data/prt_splinedata.xdr")
    splinedata = table["splinedata"]
    serial_nos = []
    for b in splinedata["serial"]:
        serial_nos.append(bytearray(b).decode("utf8").strip())
    element = np.where(np.array(serial_nos) == SerialNum)[0]
    if len(element) == 0:
        print(f"Serial number {SerialNum} not found, recheck")
        return

    # spline structure for this PRT
    a = splinedata[element][0]
    # using notation of
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.BSpline.html
    t = a["xknot"]
    c = a["bscoef"]
    k = a["korder"]
    maxind = min(sum(c != 0), sum(t != 0))
    t = t[:maxind]
    c = c[:maxind]
    spl = BSpline(t, c, k - 1)
    temperature = spl(resist)

    return temperature


def G(V, T_RXB, T_FPA, t, pars):
    T0, V0, beta, alpha, m, c = pars
    G = alpha * (V - V0 - beta * (T_RXB - 290)) / (T_FPA - T0) + m * t + c
    return G


def get_val_from_mnem(data, mnem):
    # there are several steps to get from the raw TOD to physical units.

    # First read in the data themselves, which are indexed according to the
    # mnemonics mentioned before.

    # Second, convert the raw value into physical units using the polynomial
    # conversion in line 122 of aihk_get_mnemonic.pro.

    # Finally, you need to get the serial number of the thermistor used to measure
    # the temperature, and get that thermistor's conversion to temperature.

    # I'll do this manually as an example, for K1 A side feed.

    # Thinking of doing something like K1A and then.
    # mnem = get_mnem(band, variable)

    # index = 18 in fortran indexing, so 17 in idl. Also, AEU1.

    index, arr = aihk_mnemonic(mnem)
    # index, arr = 17, 2

    # Get temperature from index
    res, wflg = get_resist(data, index, arr)
    if wflg == 0:
        temp = res
    else:
        temp = get_temp(res, mnem)

    coefs_array = np.loadtxt("gain_params/mnem_coefs.txt", dtype=str, delimiter=",")
    for i in range(len(coefs_array)):
        if mnem in coefs_array[i, 0]:
            ind = i
    coefs = coefs_array[ind, 1:].astype("float")
    T_pow = np.array([temp**i for i in range(len(coefs))])
    Val = coefs.dot(T_pow)

    return Val


def gain_tests():
    # to convert to JD, add 2.4500e06
    files = glob(prefix + "tods/uncalibrated/*.fits")
    files.sort()
    files = np.array(files)
    data = fits.open(
        "/mn/stornext/d16/cmbco/ola/wmap/tods/uncalibrated/wmap_tod_20052022356_20052032356_uncalibrated_v5.fits"
    )

    # The only published gain model I've seen in the WMAP paper is Figure 2 of
    # Jarosik et al. 2007, which plots the V223 detector.
    par_array = np.loadtxt("gain_params/T0_sols.txt", dtype=str)
    ind = np.where(par_array[:, 0] == "V2")[0][0]
    pars = par_array[ind, 1 + 2 :: 4].astype("float")

    mnem_rfb = "DRV223RFBI14"
    mnem_tfpa = "DFV22FPATEET"
    mnem_trxb = "DRV222RXBAMPT"

    T_RXB = get_val_from_mnem(data, mnem_trxb)
    T_FPA = get_val_from_mnem(data, mnem_tfpa)

    # This is confirmed correct by comparing to the IDL code.
    RF_bias = get_val_from_mnem(data, mnem_rfb)
    t_JD = data[3].data["TIME"] + 2.45e6
    # 2452131.5 is 0:00 of day 222 of 2001 (August 10)
    t = t_JD - 2452131.50000

    G_V223 = G(RF_bias, T_RXB, T_FPA, t, pars)

    fig, axes = plt.subplots(sharex=True, nrows=2, ncols=2)
    axs = axes.flatten()
    axs[0].plot(t, T_RXB, ".")
    axs[0].set_title(r"$T_\mathrm{RXB}$")
    axs[1].plot(t, T_FPA, ".")
    axs[1].set_title(r"$T_\mathrm{FPA}$")
    axs[2].plot(t, RF_bias, ".")
    axs[2].set_title("RF Bias")
    axs[3].plot(t, G_V223, ".")
    axs[3].set_title("Gain")
    # axs[3].axhline(band_dict['V223'], color='r')
    fig.suptitle("V223")

    # Figure 4 in Hinshaw et al has a couple of figures. K113 and V113;
    # K113
    ind = np.where(par_array[:, 0] == "K1")[0][0]
    pars = par_array[ind, 1::4].astype("float")
    mnem_trxb = "DRK12RXBRIBT"
    mnem_tfpa = "DFK1AFEEDT"
    mnem_rfb = "DRK113RFBI0"

    T_RXB = get_val_from_mnem(data, mnem_trxb)
    T_FPA = get_val_from_mnem(data, mnem_tfpa)
    RF_bias = get_val_from_mnem(data, mnem_rfb)
    G_K113 = G(RF_bias, T_RXB, T_FPA, t, pars)

    fig, axes = plt.subplots(sharex=True, nrows=2, ncols=2)
    axs = axes.flatten()
    axs[0].plot(t, T_RXB, ".")
    axs[0].set_title(r"$T_\mathrm{RXB}$")
    axs[1].plot(t, T_FPA, ".")
    axs[1].set_title(r"$T_\mathrm{FPA}$")
    axs[2].plot(t, RF_bias, ".")
    axs[2].set_title("RF Bias")
    axs[3].plot(t, G_K113, ".")
    # axs[3].axhline(band_dict['K113'], color='r')
    axs[3].set_title("Gain")
    fig.suptitle("K113")

    # V113
    ind = np.where(par_array[:, 0] == "V1")[0][0]
    pars = par_array[ind, 1::4].astype("float")
    mnem_trxb = "DRV111RXBAMPT"
    mnem_tfpa = "DFV11FPATEET"
    mnem_rfb = "DRV113RFBI32"

    T_RXB = get_val_from_mnem(data, mnem_trxb)
    T_FPA = get_val_from_mnem(data, mnem_tfpa)
    RF_bias = get_val_from_mnem(data, mnem_rfb)
    G_V113 = G(RF_bias, T_RXB, T_FPA, t, pars)

    fig, axes = plt.subplots(sharex=True, nrows=2, ncols=2)
    axs = axes.flatten()
    axs[0].plot(t, T_RXB, ".")
    # axs[0].axhline(290, color='r')
    axs[0].set_title(r"$T_\mathrm{RXB}$")
    axs[1].plot(t, T_FPA, ".")
    # axs[0].axhline(90, color='r')
    axs[1].set_title(r"$T_\mathrm{FPA}$")
    axs[2].plot(t, RF_bias, ".")
    axs[2].set_title("RF Bias")
    axs[3].plot(t, G_V113, ".")
    # axs[3].axhline(band_dict['V113'], color='r')
    axs[3].set_title("Gain")
    fig.suptitle("V113")

    # Some good consistency checks are in Figure 1.6 of the supplemental
    # materials. But for now, maybe let's just say that these temperatures are
    # essentially constant.
    T_FPA = 90
    T_RXB = 290
    G_V113 = G(RF_bias, T_RXB, T_FPA, t, pars)
    plt.figure()
    plt.plot(t, G_V113, ".")

    plt.show()


def rfb_tests():
    files = glob(prefix + "tods/uncalibrated/*.fits")
    files.sort()
    files = np.array(files)
    data = fits.open(files[0])
    bands = ["K1", "Ka1", "Q1", "Q2", "V1", "V2", "W1", "W2", "W3", "W4"]
    bands = [[b + "13", b + "14", b + "23", b + "24"] for b in bands]
    bands = np.array(bands).flatten()

    rf_mnems = np.loadtxt("gain_params/rfbias_mnem.txt", dtype=str)
    i = 0
    fig, axes = plt.subplots(nrows=8, ncols=5, sharex=True)
    axs = axes.flatten()
    for b in bands:
        ind = b == rf_mnems[:, 0]
        mnem = rf_mnems[ind][0][1]
        rf_bias = get_val_from_mnem(data, mnem)
        axs[i].plot(rf_bias, "k.", ms=0.5, alpha=0.5)
        axs[i].set_title(b)
        i += 1

    plt.suptitle("RF Bias")
    plt.show()


def fpa_tests():
    files = glob(prefix + "tods/uncalibrated/*.fits")
    files.sort()
    files = np.array(files)
    data = fits.open(
        "/mn/stornext/d16/cmbco/ola/wmap/tods/uncalibrated/wmap_tod_20052022356_20052032356_uncalibrated_v5.fits"
    )

    fpa_mnems = np.loadtxt("gain_params/t_fpa_mnem.txt", dtype=str, delimiter=",")
    fig, axes = plt.subplots(sharex=True, sharey=False, nrows=4, ncols=5)
    fig2, ax2 = plt.subplots()
    axs = axes.flatten()
    for i, mnem in enumerate(fpa_mnems[:, 0]):
        T_FPA = get_val_from_mnem(data, mnem)
        axs[i].plot(T_FPA, ".", ms=1)
        axs[i].set_title(mnem)
        # axs[i].set_ylim([88,92])
        ax2.plot(T_FPA, ".", ms=1)
    plt.show()


def rxb_tests():
    files = glob(prefix + "tods/uncalibrated/*.fits")
    files.sort()
    files = np.array(files)
    data = fits.open(files[0])

    rxb_mnems = np.loadtxt("gain_params/t_rxb.txt", dtype=str, delimiter=",")
    fig, axes = plt.subplots(sharex=True, sharey=False, nrows=3, ncols=5)
    fig2, ax2 = plt.subplots()
    axs = axes.flatten()
    for i, mnem in enumerate(rxb_mnems[:, 0]):
        T_RXB = get_val_from_mnem(data, mnem)
        axs[i].plot(T_RXB, ".", ms=1)
        axs[i].set_title(mnem)
        # axs[i].set_ylim([286, 290])
        ax2.plot(T_RXB, ".", ms=1)
    plt.show()


def get_mnems(band):
    """
    RFB is simple.
    FPA and RXB need to be the temperature "as measured from the nearest
    thermistor". So I need to get a schematic of the instruments and the various
    detectors.
    """
    rfb_mnems = np.char.upper(np.loadtxt("gain_params/rfbias_mnem.txt", dtype=str))
    rfb_mnem = rfb_mnems[band == rfb_mnems[:, 0]][0][1]
    tfpa_mnem = ""
    trxb_mnem = ""
    return rfb_mnem, tfpa_mnem, trxb_mnem


def get_gain(data, band):
    """
    My goal is to have a function that will take a TOD label, like K113, and
    return gain.

    Part of this requires getting the mnemonics. Some dream functions will be
    included here.
    """
    t_JD = data[3].data["TIME"] + 2.45e6
    # 2452131.5 is 0:00 of day 222 of 2001 (August 10)
    t = t_JD - 2452131.50000
    mnems = get_mnems(band)

    RF_bias = get_val_from_mnem(data, mnems[0])
    # T_FPA = get_val_from_mnem(data, 'DFV11FPATEET')
    # T_FPA = get_val_from_mnem(data, 'DFK1BOMTT')
    T_FPA = get_val_from_mnem(data, "DFQ1AFEEDT")
    # T_FPA = get_val_from_mnem(data, 'DFQ1BOMTT')
    # T_RXB = get_val_from_mnem(data, 'DRK12RXBRIBT')
    T_RXB = get_val_from_mnem(data, "DRQ1RXBRIBT")

    T_FPA = fpa_func(t_JD)
    T_RXB = rxb_func(t_JD)

    par_array = np.char.upper(np.loadtxt("gain_params/T0_sols.txt", dtype=str))
    horn_inds = np.array(["13", "14", "23", "24"])

    ind = np.where(par_array[:, 0] == band[:-2])[0][0]
    da_ind = np.where(band[-2:] == horn_inds)[0][0]
    pars = par_array[ind, 1 + da_ind :: 4].astype("float") * 1.01

    G_band = G(RF_bias, T_RXB, T_FPA, t, pars)

    return t_JD, G_band


def fullgain_tests():
    files = glob(prefix + "tods/uncalibrated/*.fits")
    files.sort()
    files = np.array(files)
    data = fits.open(files[100])
    bands = ["K1", "Ka1", "Q1", "Q2", "V1", "V2", "W1", "W2", "W3", "W4"]
    bands = [[b + "13", b + "14", b + "23", b + "24"] for b in bands]
    bands = np.array(bands).flatten()

    bands = data[2].columns.names[1:-6]
    fig, axes = plt.subplots(nrows=10, ncols=4)
    axs = axes.flatten()
    for i, b in enumerate(bands):
        t, G_b = get_gain(data, b)
        axs[i].plot(t, G_b)
        axs[i].set_title(b)

    plt.show()

    return


def publication_plots(nfiles=100):
    """
    Trying to reproduce the figures from Hinshaw 2003 and Jarosik 2007
    """
    files = glob(prefix + "tods/uncalibrated/*.fits")
    files.sort()
    files = np.array(files)
    files = files[:nfiles]

    fpa_mnems = np.loadtxt("gain_params/t_fpa_mnem.txt", dtype=str, delimiter=",")
    rxb_mnems = np.loadtxt("gain_params/t_rxb.txt", dtype=str, delimiter=",")

    # K113
    # from 221st day of 2001 to 256th day of 2002
    GsK = []
    GsV = []
    GsV223 = []
    ts = []
    T_FPAs = [[] for i in range(len(fpa_mnems))]
    T_RXBs = [[] for i in range(len(rxb_mnems))]
    for j in tqdm(range(len(files))):
        data = fits.open(files[j])
        ti, Gi = get_gain(data, "K113")
        GsK += Gi.tolist()
        ti, Gi = get_gain(data, "V113")
        GsV += Gi.tolist()
        ti, Gi = get_gain(data, "V223")
        GsV223 += Gi.tolist()
        ts += ti.tolist()
        # for i in range(len(fpa_mnems)):
        #    T_FPAs[i] += get_val_from_mnem(data, fpa_mnems[i][0]).tolist()
        # for i in range(len(rxb_mnems)):
        #    T_RXBs[i] += get_val_from_mnem(data, rxb_mnems[i][0]).tolist()
        data.close()
    ts = Time(ts, format="jd")

    # T_FPAs = np.array(T_FPAs)
    # T_RXBs = np.array(T_RXBs)
    # np.save('t_fpa', T_FPAs)
    # np.save('t_rxb', T_RXBs)

    # plt.figure()
    # for i in range(len(T_FPAs)):
    #    plt.plot(ts, T_FPAs[i])

    # plt.figure()
    # for i in range(len(T_RXBs)):
    #    plt.plot(ts, T_RXBs[i])

    plt.figure()
    plt.plot(ts.datetime, GsK, ".", ms=1)
    # plt.axhline(-0.97)
    plt.gcf().autofmt_xdate()
    plt.title("K113")
    plt.ylim([-1.04, -0.9])

    # V113
    plt.figure()
    plt.plot(ts.datetime, GsV, ".", ms=1)
    # plt.axhline(0.4986)
    plt.gcf().autofmt_xdate()
    plt.title("V113")
    plt.ylim([0.4, 0.5])

    # V223
    plt.figure()
    plt.plot(ts.datetime, GsV223, ".", ms=1)
    # plt.axhline(0.4096)
    plt.gcf().autofmt_xdate()
    plt.title("V223")
    plt.ylim([0.365, 0.395])

    plt.show()
    #'''
    return


def temp_tests(nfiles=100, band="V113"):
    par_array = np.char.upper(np.loadtxt("gain_params/T0_sols.txt", dtype=str))
    horn_inds = np.array(["13", "14", "23", "24"])

    ind = np.where(par_array[:, 0] == band[:-2])[0][0]
    da_ind = np.where(band[-2:] == horn_inds)[0][0]
    pars = par_array[ind, 1 + da_ind :: 4].astype("float")

    files = glob(prefix + "tods/uncalibrated/*.fits")
    files.sort()
    files = np.array(files)
    files = files[:nfiles]

    Temp_FPAs = np.loadtxt("gain_params/fpa_temps.txt")
    Temp_RXBs = np.loadtxt("gain_params/rxb_temps.txt")

    rng = np.random.default_rng()

    Temp_FPA = rng.choice(Temp_FPAs, axis=1)
    Temp_RXB = rng.choice(Temp_RXBs, axis=1)

    fnames = glob("/mn/stornext/d16/cmbco/ola/wmap/tods/uncalibrated/*.fits")
    fnames.sort()
    t1s = []
    t2s = []
    for f in fnames:
        t1 = f.split("_")[2]
        t2 = f.split("_")[3]
        t1 = t1[:4] + ":" + t1[4:7] + ":" + t1[7:9] + ":" + t1[9:11]
        t2 = t2[:4] + ":" + t2[4:7] + ":" + t2[7:9] + ":" + t2[9:11]
        t1 = Time(t1, format="yday")
        t2 = Time(t2, format="yday")
        t1s.append(t1.jd)
        t2s.append(t2.jd)

    t1s = np.array(t1s)
    t2s = np.array(t2s)
    times = (t1s + t2s) / 2

    inds = Temp_FPA > 50

    fpa_func = interp1d(times[inds], Temp_FPA[inds], fill_value="extrapolate")
    rxb_func = interp1d(times[inds], Temp_RXB[inds], fill_value="extrapolate")

    mnems = get_mnems(band)
    print(mnems)
    for j in range(len(files)):
        data = fits.open(files[j])
        RF_bias = get_val_from_mnem(data, mnems[0])
        t_JD = data[3].data["TIME"] + 2.45e6
        # 2452131.5 is 0:00 of day 222 of 2001 (August 10)
        t = t_JD - 2452131.50000

        G_V223 = G(RF_bias, rxb_func(t), fpa_func(t), t, pars)
        plt.plot(t, G_V223, "k.", ms=1)

    return


def plot_gain_history(band="Q113"):
    fnames = glob("/mn/stornext/d16/cmbco/ola/wmap/tods/uncalibrated/*.fits")
    fnames.sort()

    fnames = fnames[::7]

    ts = np.zeros(len(fnames))
    Gs = np.zeros(len(fnames))
    from tqdm import tqdm

    for i in tqdm(range(len(fnames))):
        data = fits.open(fnames[i])
        t_JD, G = get_gain(data, band)
        ts[i] = t_JD[0]
        Gs[i] = G[0]
    plt.plot(ts, Gs)
    np.save(f"gain_{band}", np.array([ts, Gs]))
    plt.show()


if __name__ == "__main__":
    # publication_plots(nfiles=400)
    # gain_tests()
    # rfb_tests()
    # fpa_tests()
    # rxb_tests()

    fullgain_tests()
    #for i in range(10):
    #    temp_tests(band="K113")
    #plt.show()
    plot_gain_history(band="K113")
    # plot_gain_history(band="Q114")
    # plot_gain_history(band="Q123")
    # plot_gain_history(band="Q124")
