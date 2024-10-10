from __future__ import print_function

import os
import glob
import datetime

import numpy as np
from scipy import stats

from utils.params import fac_apco_date, fac_vax_year_len, fac_apco_gmt


def pointing_conversion(lon, lat):
    """Convert poninting values in deg to values stored in the data files"""

    lon = np.copy(lon)

    lon[lon > 180.0] -= 360.0

    lat_out = np.array(np.radians(lat) * 1e4, dtype="<i2")
    lon_out = np.array(np.radians(lon) * 1e4, dtype="<i2")

    return lon_out, lat_out


def expand_glitch(gltch):

    if gltch.dtype != np.ushort or len(gltch) != 32:
        raise ValueError("Input glitch array is not valid")

    output = np.zeros(512, dtype=np.byte)

    for i in range(512):

        # Find the correct entry of the array
        idx = i / 16
        gltch_word = gltch[idx]

        # Find the correct byte and set the output to one if it is set
        offset = i % 16
        if gltch_word & (1 << offset):
            output[i] = 1

    return output


def find_engineering(fn_fdq_sdf, root=None):
    """Find the associated engineering data. Assumption is that the engineering data is
    stored in ../fdq_eng if no root path is specified
    """

    fn_split = os.path.split(fn_fdq_sdf)

    if root is None:
        root = fn_split[0]
    tmp = fn_split[1][-9:]

    tmp2 = os.path.join(root, "../fdq_eng", "*" + tmp + "*")

    fn = glob.glob(tmp2)

    if len(fn) != 1:
        raise ValueError("Finding multiple engineering files")

    print(fn)

    return os.path.normpath(fn[0])


def find_sdf(fn_fdq_eng, root=None, chan=None):
    """Find the associated science data. Assumption is that the science data is
    stored in ../fdq_sdf if no root path is specified.
    """

    fn_split = os.path.split(fn_fdq_eng)

    if root is None:
        root = fn_split[0]
    tmp = fn_split[1][-9:]

    tmp2 = os.path.join(root, "../fdq_sdf", "*" + tmp + "*")

    fn = glob.glob(tmp2)

    if len(fn) != 4:
        raise ValueError("Not finding 4 science data files")

    fn_out = []
    for i in range(4):
        if chan is not None:
            if chan in fn[i]:
                fn_out = os.path.normpath(fn[i])
        else:
            fn_out.append(os.path.normpath(fn[i]))

    return fn_out


def compare_bits(var1, var2, position, length):

    bitset = (2**length - 1) << position

    var1_b = var1 & bitset
    var2_b = var2 & bitset

    return var1_b == var2_b


def determine_chan_string(chan):

    chan_str = ["rh", "rl", "lh", "ll"]
    # return chan_str[chan]

    if chan == 3:
        return "ll"
    elif chan == 2:
        return "lh"
    elif chan == 1:
        return "rl"
    elif chan == 0:
        return "rh"
    else:
        raise ValueError("Input chan value needs to be between 0 and 3")


def determine_chan_num(str_chan):

    if type(str_chan) is int:
        return str_chan

    if "ll" in str_chan or "LL" in str_chan:
        return 3
    elif "lh" in str_chan or "LH" in str_chan:
        return 2
    elif "rl" in str_chan or "RL" in str_chan:
        return 1
    elif "rh" in str_chan or "RH" in str_chan:
        return 0
    else:
        raise ValueError("Valid chan string not found in input", str_chan)


def determine_scan_mode(str_scanmode):

    str_scanmode = str_scanmode.lower()

    if "ss" in str_scanmode:
        return 0
    elif "sf" in str_scanmode:
        return 1
    elif "ls" in str_scanmode:
        return 2
    elif "lf" in str_scanmode:
        return 3
    elif "fs" in str_scanmode:
        return 4
    elif "fl" in str_scanmode:
        return 5
    elif "all" in str_scanmode:
        return -1
    else:
        raise ValueError("Valid scan mode not found in input", str_scanmode)


def determine_sm_string(sm):
    sm_str = ["ss", "sf", "ls", "lf", "fs", "fl"]
    return sm_str[sm]


def ibits(num, pos, length):
    """Code to reproduce the Fortran ibits functions
    This takes the bits starting at position pos and going to the left with length length.
    Value is right-justified
    """

    # Should be the first length values set to 1 and the rest 0
    bitset = (1 << length) - 1

    # Shifts all the bits to the right so that the bit at position pos is not at 0. Bits on the
    # left are set to 0
    num_shift = num >> pos

    return np.bitwise_and(num_shift, bitset)


def coadd_angles(lat, lon, glwt, invert=False):

    if invert:
        f1 = np.cos
        f2 = np.sin
    else:
        f1 = np.sin
        f2 = np.cos

    sum_arr = [
        np.sum(f1(lat) * f1(lon) * glwt),
        np.sum(f1(lat) * f2(lon) * glwt),
        np.sum(f2(lat) * glwt),
    ]
    norm = np.linalg.norm(sum_arr)

    if norm != 0.0:
        avg_arr = sum_arr / norm

    if invert:
        lat_avg = np.arctan2(avg_arr[2], np.sqrt(avg_arr[0] ** 2 + avg_arr[1] ** 2))
    else:
        lat_avg = np.arctan2(avg_arr[0] ** 2 + avg_arr[1] ** 2, avg_arr[2])

    # ??????
    if avg_arr[0] != 0.0:
        lon_avg = np.arctan2(avg_arr[1], avg_arr[0])
    else:
        lon_avg = 0.0

    return lat_avg, lon_avg


def binary_to_gmt(binary):
    """Converts input ADT tme to 14-element string"""

    # ADT is a 64-bit quadword containing the number of 100-nanosecond ticks
    # since November 17, 1858

    time0 = datetime.datetime(1858, 11, 17)

    deltat = datetime.timedelta(microseconds=0.1 * binary)
    time1 = time0 + deltat

    return time1.strftime("%y%j%H%M%S%f")[:14]


def gmt_to_binary(gmt):
    """Convers an input string value to a VAX binary time"""

    time0 = datetime.datetime(1858, 11, 17)
    time1 = datetime.datetime.strptime(gmt, "%y%j%H%M%S%f")

    deltat = time1 - time0

    daystons = np.array(1, dtype=np.int64)
    daystox = 24 * 60 * 60 * 10000000
    secondstox = 10000000
    microsecondstox = 10
    binary = np.array(0, dtype=np.int64)
    binary += deltat.days * daystox
    binary += deltat.seconds * secondstox
    binary += deltat.microseconds * microsecondstox

    return binary


def adt_to_fits(binary):
    """Converts input ADT tme to value stored in the published FITS files"""

    # ADT is a 64-bit quadword containing the number of 100-nanosecond ticks
    # since November 17, 1858
    # Fits file time is number of seconds since Jan 1, 1981
    time0 = datetime.datetime(1858, 11, 17)
    time1 = datetime.datetime(1981, 1, 1)

    deltat2 = time1 - time0

    deltat = datetime.timedelta(microseconds=0.1 * binary)

    deltat3 = deltat - deltat2

    daystox = 24 * 60 * 60
    secondstox = 1
    binary_out = np.array(0, dtype=np.int64)
    binary_out += deltat3.days * daystox
    binary_out += deltat3.seconds * secondstox

    return binary_out


def fits_to_adt(binary):

    time0 = datetime.datetime(1858, 11, 17)
    time1 = datetime.datetime(1981, 1, 1)

    deltat2 = time1 - time0

    deltat = datetime.timedelta(seconds=1.0 * binary)

    deltat3 = deltat + deltat2

    daystox = 24 * 60 * 60 * 10000000
    secondstox = 10000000
    microsecondstox = 10
    binary = np.array(0, dtype=np.int64)
    binary += deltat3.days * daystox
    binary += deltat3.seconds * secondstox
    binary += deltat3.microseconds * microsecondstox

    return binary


def encode_label(first_gmt, last_gmt, num_ifgs, scan_mode, ical, xcal_pos, pixel_no):

    label = ""
    label += str(first_gmt)
    label += "_"
    label += str(last_gmt)
    label += "_"
    label += str(int(num_ifgs)).zfill(3)
    label += "_"
    label += str(scan_mode)
    label += "_"
    label += str(int(ical)).zfill(5)
    label += str(xcal_pos)
    if pixel_no >= 0:
        label += "_"
        label += str(int(pixel_no)).zfill(4)

    return label


def planck_dist(freq, temp, sigma_t):
    """This routine returns the brightness at one frequency from a distribution of
    blackbodies with a known average temperature and standard deviation. It approximates
    a continuous distribution using five Planck curves, evenly spaced in temperature such
    that the average and rms equals that of the desired distribution. For five curves, this
    turns out to be delta_T = sigma/sqrt(2). The returned value has the units ergs/sec/cm**2/sr/icm.
    """

    dt = sigma_t / np.sqrt(2)
    total_flux = np.zeros_like(freq)
    for j in range(-2, 3):
        total_flux += planck(freq, temp + dt * j)
    dplanck_dist = total_flux / 5.0

    return dplanck_dist


def planck_dist2(freq, temp, sigma_t):
    """planck_dist except temperature can be a vector"""

    dt = sigma_t / np.sqrt(2)

    total_flux = 0.0

    for j in range(-2, 3):
        total_flux += planck2(freq, temp + dt * j)

    dplanck_dist = total_flux / 5.0

    return dplanck_dist


def planck(anu, T):
    """This routine calculates the Planck function for a specified frequency and
    temperature. The returned value has units of ergs/sec/cm**2/sr/icm.
    """

    hck = 1.438769
    thcc = 1.1910439e-5

    # Check for temperature less than or equal to 0
    if T <= 0:
        return 0.0

    # Calculate the exponential term
    x = hck * anu / T

    x = np.atleast_1d(x)
    anu = np.atleast_1d(anu)

    # Calculate the exponential denominator
    d = np.zeros_like(x)
    idx = np.logical_and(0.01 < x, x <= 70.0)
    d[idx] = np.exp(x[idx]) - 1.0
    idx = x <= 0.01
    d[idx] = x[idx] * (1.0 + x[idx] * (0.5 + x[idx] / 6.0))

    # Calculate the value of the Planck function
    planckval = np.zeros_like(d)
    idx = d != 0
    planckval[idx] = thcc * (anu[idx] ** 3 / d[idx])

    return np.squeeze(planckval)


def planck2(anu, T):

    hck = 1.438769
    thcc = 1.1910439e-5

    # Calculate the exponential term
    x = hck * np.outer(1.0 / T, anu)

    # Calculate the exponential denominator
    d = np.zeros_like(x)
    idx = np.logical_and(0.01 < x, x <= 70.0)
    d[idx] = np.exp(x[idx]) - 1.0
    idx = x <= 0.01
    d[idx] = x[idx] * (1.0 + x[idx] * (0.5 + x[idx] / 6.0))

    # Calculate the value of the Planck function
    planckval = np.zeros_like(d)
    idx = d != 0
    n = np.shape(x)[0]
    anu2 = np.tile(anu, n).reshape(n, -1)
    planckval[idx] = thcc * (anu2[idx] ** 3 / d[idx])

    # Check for temperature less than or equal to 0
    idx = T <= 0
    planckval[idx, :] = 0.0

    return np.squeeze(planckval)


def planck_deriv(anu, T):
    """This routine calculates the derivative of the Planck function with respect
    to temperature for a specified frequency and temperature. The returned value
    has units of ergs/cm**2/sr/icm/K
    """

    hck = 1.438769
    thcc = 1.1910439e-5

    # Check for temperature less than or equal to 0
    if T <= 0:
        return 0.0

    # Calculate the exponential term
    x = hck * anu / T

    x = np.atleast_1d(x)

    # Calculate the exponental denominator
    d = np.zeros_like(x)
    idx = np.logical_and(0.01 < x, x <= 70)
    d[idx] = np.exp(x[idx]) + np.exp(-x[idx]) - 2.0
    idx = x <= 0.01
    d[idx] = x[idx] ** 2 + (1.0 + x[idx] ** 2 / 12.0)

    # Calculate the derivative of the Planck function
    planckderiv_val = np.zeros_like(d)
    idx = d != 0

    planckderiv_val[idx] = thcc * (x[idx] / T) * (anu[idx] ** 3 / d[idx])

    if type(anu) is not np.ndarray:
        return planckderiv_val[0]
    else:
        return planckderiv_val


def get_mode(array, default=-1):
    """Get the mode of the array. If the array has no length, return the default"""

    if len(array) == 0:
        return default
    else:
        vals, counts = np.unique(array, return_counts=True)
        return vals[counts.argmax()]


def get_mode2(array, default=-1, num_cgr_good=None):

    # If all entries have already been marked bad, return default value
    if len(array) == 0:
        return default

    arr_unique = np.unique(array)
    counts = np.zeros_like(arr_unique, dtype=int)

    for i, val in enumerate(arr_unique):
        counts[i] = np.sum(array == val)

    idx = np.where(counts == np.max(counts))[0]

    if len(idx) > 1 and num_cgr_good is not None:
        val_out = get_mode2(array[:num_cgr_good])
    else:
        val_out = arr_unique[idx[0]]

    return val_out


def adt2t68(bin_time):
    """Convert UTC time in ADT format to T68 (seconds from 00:00:00.0 UT on 24 May 1968 = JD 2440000.5)"""

    gmt = binary_to_gmt(bin_time)

    return gmt2t68(gmt)


def gmt2t68(gmt):
    """Convert GMT string to T68 (seconds from 00:00:00.0 UT on 24 May 1968 = JD 2440000.5)

    Notes
    -----
    Output seems to be different from Fortran code by several seconds. For example date corresponding
    to 1-JAN-1972 outputs 113788800.0 whereas the comments in the Fortran code says it should be
    113788803.52, which seems weird since it is a fraction of a second different and both are at
    00:00:00.0 UT
    """

    time0 = datetime.datetime(1968, 5, 24)
    time1 = datetime.datetime.strptime(gmt, "%y%j%H%M%S%f")

    deltat = time1 - time0

    return deltat.total_seconds()


def time_since_apco_NOTUSED(spec_time):
    """This function computes the time in years since the aperture cover ejection."""

    # NOTE: This is not used because the bitshift loses resolution
    spec_time = spec_time >> 32

    # Calculate the time since aperture cover ejection for the spectrum
    time = (spec_time - fac_apco_date) / fac_vax_year_len

    return time


def time_since_apco(spec_time):
    apco_eject_time = gmt_to_binary(fac_apco_gmt + "00000")
    year_len = 31556925.19

    time = (spec_time - apco_eject_time) / year_len / 1e7

    return time


def bintime_to_mytime(bintime):
    """My time is just some modified time (probably in terms of hours since near the
    start of FIRAS observations
    """

    mytime = bintime / 1e10 - 4136834

    return mytime


def gmt_to_mytime(gmt):

    if type(gmt) is str:
        bintime = gmt_to_binary(gmt)
    else:
        bintime = np.empty(len(gmt))
        for i, time in enumerate(gmt):
            bintime[i] = gmt_to_binary(time)

    mytime = bintime_to_mytime(bintime)

    return mytime


def fitstime_to_mytime(fitstime):
    """Input is time stored in the published FIRAS fits files. Time is
    number of seconds since 1 January 1981
    """

    bintime = fitstime_to_bintime(fitstime)

    mytime = bintime_to_mytime(bintime)

    return mytime


def fitstime_to_bintime(fitstime):

    # ADT time is 100 ns tics since Nov 17, 1858
    time0 = datetime.datetime(1858, 11, 17)

    time1 = datetime.datetime(1981, 1, 1)

    dt = time1 - time0

    bintime = (fitstime + dt.total_seconds()) * 1e7

    return bintime


def fitstime_to_mytime2(fitstime):
    """My time is in hours since start of FIRAS observations"""

    # bintime of first entry in fdq_eng
    bintime0 = 41342940987990000

    # ADT time is 100 ns tics since Nov 17, 1858
    time0 = datetime.datetime(1858, 11, 17)
    time1 = datetime.datetime(1981, 1, 1)

    dt = time1 - time0

    # To seconds
    bintime0 /= 1e7

    # Get first time in FITS time
    mytime0 = bintime0 - dt.total_seconds()

    mytime = (fitstime - mytime0) / 3600

    return mytime


def extrap1d(interpolator):
    """Code to apply extrapolation to data points outside the interpolator range when using
    scipy.interpolate.interp1d. Needed for VPM template stuff"""

    xs = interpolator.x
    ys = interpolator.y

    # Quadratic for x < xs
    a = np.array(
        [[xs[0] ** 2, xs[0], 1], [xs[1] ** 2, xs[1], 1], [xs[2] ** 2, xs[2], 1]]
    )
    b = ys[:3]
    p1 = np.linalg.solve(a, b)

    # Quadratic for x > xs
    a = np.array(
        [[xs[-3] ** 2, xs[-3], 1], [xs[-2] ** 2, xs[-2], 1], [xs[-1] ** 2, xs[-1], 1]]
    )
    b = ys[-3:]
    p2 = np.linalg.solve(a, b)

    def pointwise_linear(x):
        if x < xs[0]:
            return ys[0] + (x - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0])
        elif x > xs[-1]:
            return ys[-1] + (x - xs[-1]) * (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
        else:
            return interpolator(x)

    def pointwise_quadratic(x):
        if x < xs[0]:
            return p1[0] * x**2 + p1[1] * x + p1[2]
        elif x > xs[-1]:
            return p2[0] * x**2 + p2[1] * x + p2[2]
        else:
            return interpolator(x)

    def ufunclike_linear(x):
        data = interpolator(x)
        idx1 = x < xs[0]
        idx2 = x > xs[-1]

        data[idx1] = ys[0] + (x[idx1] - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0])
        data[idx2] = ys[-1] + (x[idx2] - xs[-1]) * (ys[-1] - ys[-1]) / (xs[-1] - xs[-2])

        return data

    def ufunclike_quadratic(x):
        data = interpolator(x)
        idx1 = x < xs[0]
        idx2 = x > xs[-1]

        data[idx1] = p1[0] * x[idx1] ** 2 + p1[1] * x[idx1] + p1[2]
        data[idx2] = p2[0] * x[idx2] ** 2 + p2[1] * x[idx2] + p2[2]

        return data

    def ufunclike(xs):
        return np.array(map(pointwise_quadratic, np.array(xs)))

    return ufunclike_quadratic
