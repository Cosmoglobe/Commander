from __future__ import print_function

import os
import struct

import numpy as np
import h5py

import data_types as dt
import vax_conversion as vax_conversion
from params import *

fac_xcalin = 1
fac_xcalout = 2
fac_both = 3

fac_att_conv = 1.0e-4 * 180.0 / 3.1415926  # convert data file value to degrees


def read_fdq_sdf(
    fn,
    galaxy=0,
    domask=True,
    minsun=91.2,
    maxsun=180.0,
    minmoon=22,
    maxmoon=180,
    minearth=87,
    maxearth=180.0,
    check_pixels=False,
    report=False,
    fixfloats=False,
    chan=None,
):

    if fn[-2:] == "h5":
        h5file = h5py.File(fn, "r")
        idx = "fdq_sdf_" + chan.lower()
        data_out = h5file[idx][:]
        h5file.close()
        fixfloats = False
        nrecords = len(data_out)
    else:
        statinfo = os.stat(fn)

        nbytes = statinfo.st_size
        bytes_per_record = 1536

        nrecords = nbytes / bytes_per_record

        if nbytes % bytes_per_record != 0:
            raise ValueError("Issue with size of file", fn, nbytes, bytes_per_record)

        data_out = np.fromfile(fn, dtype=dt.dt_fdq_sdf)

    if fixfloats:
        fix_floats(data_out)

    if not domask:
        return data_out

    # Tests to check if data should be used or whether we should mask it
    # This set of tests is taken from fss_read_science
    latitude = data_out["attitude"]["galactic_latitude"]
    sun = data_out["attitude"]["sun_angle"] * fac_att_conv
    moon = data_out["attitude"]["moon_angle"] * fac_att_conv
    earth = data_out["attitude"]["earth_limb"] * fac_att_conv

    sample_mask = 4095  # 0FFF
    start_sample = data_out["sci_head"]["sc_head16"]
    start_sample = np.bitwise_and(data_out["sci_head"]["sc_head16"], sample_mask)

    if check_pixels and pix_type <= 2:
        pixel_in_list = np.zeros(nrecords, dtype=np.bool)
        if pix_type == 1:
            pix_temp = data_out["attitude"]["pixel_no"] / (4 ** (6 - big_pixel))
        elif pix_type == 2:
            pix_temp = data_out["attitude"]["terr_pixel_no"] / (4 ** (6 - big_pixel))
        for i in range(num_pixels):
            pixel_in_list[pix_temp == plist[i]] = True
    else:
        pixel_in_list = True

    tmp = sun < minsun
    tmp2 = sun >= maxsun
    incsun = ~np.logical_or(tmp, tmp2)
    out_minsun = np.sum(tmp, dtype=np.int)
    out_maxsun = np.sum(tmp2, dtype=np.int)

    tmp = moon < minmoon
    tmp2 = moon >= maxmoon
    incmoon = ~np.logical_or(tmp, tmp2)
    out_minmoon = np.sum(tmp, dtype=np.int)
    out_maxmoon = np.sum(tmp2, dtype=np.int)

    tmp = earth < minearth
    tmp2 = earth >= maxearth
    incearth = ~np.logical_or(tmp, tmp2)
    out_minearth = np.sum(tmp, dtype=np.int)
    out_maxearth = np.sum(tmp2, dtype=np.int)

    good_records = np.abs(latitude) >= galaxy
    good_records = np.logical_and(good_records, pixel_in_list)
    good_records = np.logical_and(good_records, start_sample == 1)
    good_records = np.logical_and(good_records, incsun)
    good_records = np.logical_and(good_records, incmoon)
    good_records = np.logical_and(good_records, incearth)

    nrecords_good = np.sum(good_records)

    if report:
        print("Channel = ", chan, "Total science records = ", nrecords)
        print("Records accepted = ", nrecords_good)
        print("Records rejected for (sets may overlap):")
        # print("Low Gal:", out_gal)
        # print("Not in pixel:", out_pix)
        # print("Start sample:", out_samp)
        print("Sun min/max", out_minsun, out_maxsun)
        print("Moon min/max", out_minmoon, out_maxmoon)
        print("Earth min/max", out_minearth, out_maxearth)

    return data_out[good_records]


def read_fdq_eng(fn, xcal_pos=fac_xcalin, fakeit=False, fixfloats=False, domask=True):

    if fn[-2:] == "h5":
        h5file = h5py.File(fn, "r")
        data_tmp = h5file["fdq_eng"][:]
        data_out = data_tmp.view(dtype=dt.dt_fdq_eng)
        h5file.close()
        fixfloats = False
    else:
        statinfo = os.stat(fn)

        nbytes = statinfo.st_size
        bytes_per_record = 1024

        nrecords = nbytes / bytes_per_record

        if nbytes % bytes_per_record != 0:
            raise ValueError("Issue with size of file", fn, nbytes, bytes_per_record)

        data_out = np.fromfile(fn, dtype=dt.dt_fdq_eng)

    if fixfloats:
        fix_floats(data_out)

    if not domask:
        return data_out

    # Lots of tests to check if data should be used or whether we should mask it
    # This set of masking was taken from the fec_load code and might be different
    # for the other scripts
    fakeitvals = data_out["chan"]["fakeit"]
    if fakeit:
        mask_fakeit = np.any(fakeitvals != 1, axis=1)
    else:
        mask_fakeit = np.any(fakeitvals != 0, axis=1)

    pos = data_out["en_xcal"]["pos"]
    if xcal_pos == fac_xcalin:
        mask_xcalpos = np.any(pos != fac_xcalin, axis=1)
    elif xcal_pos == fac_xcalout:
        mask_xcalpos = np.any(pos != fac_xcalout, axis=1)
    elif xcal_pos == fac_both:
        mask_1 = np.any(pos != fac_xcalin, axis=1)
        mask_2 = np.any(pos != fac_xcalout, axis=1)
        mask_xcalpos = np.logical_and(mask_1, mask_2)

    dq_summary_flag = data_out["en_head"]["dq_summary_flag"]
    mask_dqsummary = np.any(dq_summary_flag == 127, axis=1)

    mask = np.logical_or(mask_fakeit, mask_xcalpos)
    mask = np.logical_or(mask, mask_dqsummary)

    # Because if mask is true we want to ignore the data
    mask = ~mask

    return data_out[mask]


def read_short_science(fn, fixfloats=False):

    statinfo = os.stat(fn)

    nbytes = statinfo.st_size
    bytes_per_record = 64

    nrecords = nbytes / bytes_per_record

    if nbytes % bytes_per_record != 0:
        raise ValueError("Issue with size of file", fn, nbytes, bytes_per_record)

    data_out = np.fromfile(fn, dtype=dt.dt_short_science)

    if fixfloats:
        fix_floats(data_out)

    return data_out


def read_fil_data(fn, quality=None):

    statinfo = os.stat(fn)

    nbytes = statinfo.st_size
    bytes_per_record = 11246

    nrecords = nbytes / bytes_per_record

    if nbytes % bytes_per_record != 0:
        raise ValueError("Issue with size of file", fn, nbytes, bytes_per_record)

    data_out = np.fromfile(fn, dtype=dt.dt_fil_sky)

    instr_qual = None
    attit_qual = None
    if quality is not None:
        if type(quality) is int:
            instr_qual = quality
            attit_qual = quality
        else:
            instr_qual = quality[0]
            attit_qual = quality[1]

    if instr_qual is not None:
        idx = data_out["coad_spec_data"]["dq_summary_flag"] <= instr_qual
        data_out = data_out[idx]

    if attit_qual is not None:
        idx = data_out["coad_spec_data"]["att_summary_flag"] <= attit_qual
        data_out = data_out[idx]

    return data_out


def read_fsl_data(fn):

    statinfo = os.stat(fn)

    nbytes = statinfo.st_size
    # bytes_per_record = 11246
    bytes_per_record = dt.dt_fsl_sky.itemsize

    nrecords = nbytes / bytes_per_record

    if nbytes % bytes_per_record != 0:
        raise ValueError("Issue with size of file", fn, nbytes, bytes_per_record)

    data_out = np.fromfile(fn, dtype=dt.dt_fsl_sky)

    return data_out


def save_fsl_data(fn, fsl_recs):

    fsl_recs.tofile(fn)


def fix_floats(data):
    """Fix VAX floats by running binary string through vax_conversion.

    Notes
    -----
    The way VAX float data is stored is different than the IEEE standard so anything
    we read in as a float in our code will be nonsense. This code converts every float
    in the numpy array back into a 4/8-bit binary string and then decodes the float
    value correctly
    """

    names = data.dtype.names

    for name in names:
        dtype = data[name].dtype
        if dtype == np.float32:
            for x in np.nditer(data[name], op_flags=["readwrite"]):
                stream = struct.pack("<f", x)
                x[...] = vax_conversion.vr4tof(stream)
            # ndata = data[name].size
            # tmp = data[name].reshape(ndata,)
            # for i in range(ndata):
            #    stream = struct.pack('<f', tmp[i])
            #    tmp[i] = vax_conversion.vr4tof(stream)
        elif dtype == np.float64:
            for x in np.nditer(data[name], op_flags=["readwrite"]):
                stream = struct.pack("<d", x)
                x[...] = vax_conversion.vd8tod(stream)
            # ndata = data[name].size
            # tmp = data[name].reshape(ndata,)
            # for i in range(ndata):
            #    stream = struct.pack('<d', tmp[i])
            #    tmp[i] = vax_conversion.vd8tod(stream)
        elif (
            dtype == np.byte
            or dtype == np.int16
            or dtype == np.int32
            or dtype == np.int64
            or dtype == np.uint16
        ):
            pass
        elif dtype.type is np.string_:
            pass
        elif dtype.type == np.complex64:
            for x in np.nditer(data[name], op_flags=["readwrite"]):
                xr = x.real
                xi = x.imag
                streamr = struct.pack("<f", xr)
                streami = struct.pack("<f", xi)
                tmpr = vax_conversion.vr4tof(steamr)
                tmpi = vax_conversion.vr4tof(steami)
                x[...] = tmpr + 1j * tmpi
        elif dtype.type == np.complex128:
            for x in np.nditer(data[name], op_flags=["readwrite"]):
                xr = x.real
                xi = x.imag
                streamr = struct.pack("<d", xr)
                streami = struct.pack("<d", xi)
                tmpr = vax_conversion.vd8tod(streamr)
                tmpi = vax_conversion.vd8tod(streami)
                x[...] = tmpr + 1j * tmpi
        else:
            # A dtype that we have defined and need to use recursion
            fix_floats(data[name])


def test_eng():
    import glob
    import datetime
    import matplotlib.pyplot as plt

    fns = glob.glob("/Users/njmille2/COBE/FIRAS/data/fdq_eng/*")

    time0 = datetime.datetime.now()
    data = read_fdq_eng(fns[0])
    time1 = datetime.datetime.now()

    print("Time to run:", time1 - time0)
    print("TEST1:", data["ct_head"]["gmt"][0:10])

    temps_init = data["en_analog"]["a_lo_xcal_tip"]
    temps_fin = []
    import struct

    print(len(temps_init))
    for i in range(len(temps_init)):
        stream = struct.pack("<f", temps_init[i])
        temps_fin.append(vax_conversion.vr4tof(stream))

    fix_floats(data)

    plt.figure()
    plt.plot(temps_init)
    plt.figure()
    plt.plot(temps_fin)
    plt.show()


def test_sdf():
    import glob
    import datetime

    fns = glob.glob("/Users/njmille2/COBE/FIRAS/data/fdq_sdf/*")

    time0 = datetime.datetime.now()
    data = read_fdq_sdf(fns[700])
    time1 = datetime.datetime.now()

    print("Time to run:", time1 - time0)
    print("TEST1:", data["ct_head"]["gmt"][0:10])


def test_ssi():
    import glob
    import datetime
    import matplotlib.pyplot as plt

    fns = glob.glob("/Users/njmille2/COBE/FIRAS/data/fss_sssky/FSS*")

    data = read_short_science(fns[15])

    print(np.shape(data))

    plt.figure()
    plt.plot(data["mtm_scan_speed"], "b")
    plt.plot(data["mtm_scan_length"] + 2, "r")
    plt.show()


if __name__ == "__main__":

    test_eng()
    # test_sdf()
    # test_ssi()
