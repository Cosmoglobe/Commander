import detection
import globals as g
import numpy as np

res = np.load(f"{g.DATA_PATH}143-2a_simulations.npy")
detection.matched_filter(res)