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

FAC_ICM_GHZ = 29.9792458
FAC_ERG_TO_MJY = 1.0e8 / FAC_ICM_GHZ

FAC_ETENDU = 1.5  # nathan's pipeline
FAC_ADC_SCALE = 204.75  # nathan's pipeline

# map parameters
COORDINATES = "G"
NSIDE = 32
NPIX = hp.nside2npix(NSIDE)

# save parameters
PNG = False
FITS = True

# plotting parameters
CHANNELS_PLOT = ["ll"]
MODES_PLOT = ["ss"]
JOINT = False
SCANUPDOWN = False

PEAK_POSITIONS = {
    "lh_ss": 357,
    "rh_ss": 357,
    "ll_ss": 360,
    "rl_ss": 360,
    "ll_lf": 90,
    "rl_lf": 90,
}

IFG_SIZE = 512
SPEC_SIZE = 257

CHANNELS = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}
