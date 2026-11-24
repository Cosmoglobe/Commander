import h5py
from cosmoglobe.tod_tools import commander_tod
import random
import glob
import numpy as np

BANDS = ['065', '090', '140', '160']
NEW2OLD_MAP = {'065': 'N60', '090': 'WIDE-S', '140': 'WIDE-L', '160': 'N160'}
OLD2NEW_MAP = {'N60': '065', 'WIDE-S': '090', 'WIDE-L': '140', 'N160': '160'}

NUM_DETECTORS = {
    '065': 40, #N60
    '090': 60, #WIDE-S
    '140': 45, #WIDE-L 
    '160': 30  #N160
}

OLD_START_DET_INDS = { 
    '065': 0,
    '090': 40,
    '140': 0,
    '160': 45
}


def select_random_old_chunk(h5file):
    selection = list(h5file)[:1]
    return random.choice(selection)

def get_random_band():
    return random.choice(BANDS)


def get_random_detector(band):

    detind = random.randrange(0, NUM_DETECTORS[band])
    return detind


def transform_detector_to_old_format(band, detector):
    if band in ('140', '160'):
        return NEW2OLD_MAP[band], detector

    if band == '065':
        old_bandname = 'WIDE-S'
        detector = detector + 20
    elif band == '090':
        if detector <= 40:
            old_bandname = 'N60'
        else:
            old_bandname = NEW2OLD_MAP[band]
            detector = detector - 40

    return old_bandname, detector


def get_random_chunk(h5file):
    selection = list(h5file)[:-1]

    return random.choice(selection)


def resolve_hdffile_band(fname):
    band = None
    for testband in BANDS:
        if f'_{testband}_' in fname:
            band = testband
            break

    if band is None:
        for key, value in NEW2OLD_MAP.items():
            if f'_{value}_' in fname:
                band = value
                break
    
    return band


def load_flags(hdf_fname, flagnums=None, detector='random', chunk='random'):
    h5file = h5py.File(hdf_fname)

    tod_reader = commander_tod.TODLoader(hdf_fname, 'dummy')

    tod_reader.outFile = h5file

    band = resolve_hdffile_band(hdf_fname)

    if detector == 'random':
        detector = str(get_random_detector(band))

    if chunk == 'random':
        chunk = get_random_chunk(h5file)

    fieldname = '/' + chunk + '/' + detector + '/flag'
    print(hdf_fname)
    print(fieldname)

    flags = tod_reader.load_field(fieldname)

    outflags = []
    for i, flagnum in enumerate(flagnums):
        outflags.append([])
        for flag in flags:
            boolstr = bin(flag)[2:]
            if len(boolstr) <= flagnum:
                outflags[i].append(False)
                continue
            outflags[i].append(bool(int(boolstr[::-1][flagnum])))
    return outflags


def load_tods(hdf_fname, band, detector='random', chunk='random'):
    h5file = h5py.File(hdf_fname)

    tod_reader = commander_tod.TODLoader(hdf_fname, 'dummy')

    tod_reader.outFile = h5file

    if detector == 'random':
        detector = str(get_random_detector(band))

    if chunk == 'random':
        chunk = get_random_chunk(h5file)


    print(detector)
    print(chunk)
    fieldname = '/' + chunk + '/' + detector + '/tod'
    print(hdf_fname)
    print(fieldname)

    tod = tod_reader.load_field(fieldname)
    return tod

#    outflags = []
#    for i, flagnum in enumerate(flagnums):
#        outflags.append([])
#        for flag in flags:
#            boolstr = bin(flag)[2:]
#            if len(boolstr) <= flagnum:
#                outflags[i].append(False)
#                continue
#            outflags[i].append(bool(int(boolstr[::-1][flagnum])))
#    return outflags



def get_random_file(path):
    files = glob.glob(path + '*.h5')
    return random.choice(files)


def get_old_file_properties(fname):
#    print(fname)
    h5file = h5py.File(fname)
    chunk = get_random_chunk(h5file)
    file_lastname = fname.split('/')[-1]
    orig_oldband = resolve_hdffile_band(file_lastname)
    h5file.close()
    newdet, olddet, oldband = get_random_detector(orig_oldband)
    if orig_oldband != oldband:
        new_fname = fname.replace(orig_oldband, oldband)

    return chunk, OLD2NEW_MAP[oldband], newdet, olddet, oldband, new_fname


def get_file_given_chunk(fdir, chunk, band, mode='new'):
    if mode == 'new':
        chunknum = int(chunk)
        filenum = int((chunknum - 1) / 200) + 1
        filename = fdir + f'AKARI_{band}_{filenum:06d}.h5'
    return filename
