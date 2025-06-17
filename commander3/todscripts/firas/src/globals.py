# analysis parameters
OFFSET = 0

# data paths
PROCESSED_DATA_PATH = f"/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/data/processed_sky_offset_{OFFSET}.npz"
PREPROCESSED_DATA_PATH_SKY = "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v4.4.h5"
PREPROCESSED_DATA_PATH_CAL = "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/cal_v4.4.h5"

# constants
T_CMB = 2.72548  # Fixsen 2009

# map parameters
COORDINATES = "G"
NSIDE = 32

# save parameters
PNG = True
FITS = False

# plotting parameters
CHANNELS = ["ll"]
MODES = ["ss"]
JOINT = False
SCANUPDOWN = False