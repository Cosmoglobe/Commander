# constants
T_CMB = 2.72548  # Fixsen 2009

# analysis parameters
OFFSET = 0.5

# map parameters
COORDINATES = "G"
NSIDE = 32

# data paths
if OFFSET == 0:
    PROCESSED_DATA_PATH = "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/data/processed_sky.npz"
else:
    PROCESSED_DATA_PATH = f"/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/data/processed_sky_offset_{OFFSET}.npz"

# save parameters
PNG = True
FITS = False