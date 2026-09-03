import classification
import detection
import globals as g
import matplotlib.pyplot as plt
import numpy as np
import subtraction

res = np.load(f"{g.DATA_PATH}143-2a_simulations.npy")
seconds = np.linspace(0, len(res) / g.SAMPRATE, len(res))

glitch_idx = detection.matched_filter(res, seconds)
# quit()

glitch_labels, glitch_amps = classification.classify_glitches(glitch_idx[glitch_idx < 1000], res,
                                                              seconds)

T = subtraction.build_template_matrix(glitch_idx, seconds, glitch_labels, glitch_amps)
plt.imshow(T, aspect='auto', extent=[0, len(glitch_idx), seconds[-1], seconds[0]])
plt.show()