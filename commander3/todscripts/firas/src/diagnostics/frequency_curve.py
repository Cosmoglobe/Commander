"""
This script is meant to plot a curve of the intensity of the galaxy and outside the galaxy over the frequencies.
"""

import globals as g
import numpy as np

data = np.load(g.PROCESSED_DATA_PATH)
