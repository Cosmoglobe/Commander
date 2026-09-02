import classification
import detection
import globals as g
import numpy as np

res = np.load(f"{g.DATA_PATH}143-2a_simulations.npy")
seconds = np.linspace(0, len(res) / g.SAMPRATE, len(res))

glitch_idx = detection.matched_filter(res)
# quit()

classification.classify_glitches(glitch_idx[glitch_idx < 1000], res, seconds)