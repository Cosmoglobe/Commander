import os

import healpy as hp

# analysis parameters
OFFSET = 0

user = os.environ["USER"]

# data paths
ORIGINAL_DATA = "/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5"
ORIGINAL_DATA_ENG = "/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5"
PREPROCESSED_DATA_PATH_SKY = f"/mn/stornext/u3/{user}/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v4.4.h5"
PREPROCESSED_DATA_PATH_CAL = f"/mn/stornext/u3/{user}/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/cal_v4.4.h5"
PROCESSED_DATA_PATH = f"/mn/stornext/u3/{user}/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/data/processed_sky_offset_{OFFSET}.npz"
PROCESSED_DATA_PATH_CAL = f"/mn/stornext/u3/{user}/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/data/processed_cal.npz"

SAVE_PATH = f"/mn/stornext/d16/www_cmb/{user}/firas/"

# original pipeline parameters
PUB_MODEL = "/mn/stornext/d16/cmbco/ola/firas/pub_calibration_model/"

# constants
T_CMB = 2.72548  # Fixsen 2009

# map parameters
COORDINATES = "G"
NSIDE = 32
NPIX = hp.nside2npix(NSIDE)

# save parameters
PNG = True
FITS = False

# plotting parameters
CHANNELS = ["rh"]
MODES = ["ss"]
JOINT = False
SCANUPDOWN = False