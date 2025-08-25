from __future__ import annotations
from dataclasses import dataclass
import itertools
from pathlib import Path
from datetime import timedelta
import multiprocessing
import time
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import healpy as hp
from matplotlib import pyplot as plt
import quadcube
import numpy as np
import akari_utils
from scipy.interpolate import interp1d
from astropy.time import Time, TimeDelta
from cosmoglobe.tod_tools import TODLoader
import pickle
import os
from multiprocessing import Pool
from glob import glob
import multiprocessing as mp

import multi_ori_3 as multi
import make_hdf5_ori_3 as make_hdf5

BAND = 'WIDE-S'
AKARI_DATA_PATH = Path(f"/mn/stornext/d23/cmbco/cg/AKARI/tamakim3/data/data_{BAND}")
CIO_PATH = Path('/mn/stornext/d23/cmbco/cg/AKARI/akari_TSD_pkl')
src_dir = f"/mn/stornext/d23/cmbco/cg/AKARI/tamakim3/data"

multi.multi(BAND, AKARI_DATA_PATH, CIO_PATH) # make h5 files for each 1h pkl
make_hdf5.make_hdf5(BAND, src_dir) # combine h5 files to a file for each band
