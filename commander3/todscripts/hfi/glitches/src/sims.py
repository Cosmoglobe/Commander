"""
This script generates glitch simulations that last 1 hour (1 scan for Planck).
"""

import Commander.commander3.todscripts.hfi.glitches.globals as g
import matplotlib.pyplot as plt
import numpy as np
import Commander.commander3.todscripts.hfi.glitches.templates as templates

seconds = np.arange(0, 3600, 1 / g.SAMPRATE)  # 1 hour of data at 180.3737 Hz

res = np.random.normal(0, g.SIGMA, len(seconds))  # white noise

# there is about 1 glitch per second, so per hour we have about
nglitches = 3600
# get that number of random points to inser the glitches

glitch_time = np.linspace(0, 2, 2 * 180, endpoint=False)  # glitch model lasts 2 seconds
glitch_indices = np.random.choice(len(seconds), nglitches, replace=False)
glitch_amps = np.random.uniform(0.1, 1.0, nglitches)  # random amplitudes for the glitches
glitch_types = np.random.choice(g.GLITCH_TYPES, nglitches, replace=True)  # random glitch types

for glitch_i, glitch_idx in enumerate(glitch_indices):
    # get the glitch model for this glitch
    glitch_model = templates.glitch_model_func(glitch_time, glitch_amps[glitch_i], "143-2a",
                                               glitch_types[glitch_i])
    # add the glitch to the data
    if glitch_idx + 2 * 180 < len(res):
        res[glitch_idx:glitch_idx + 2 * 180] += glitch_model
    else:
        res[glitch_idx:] += glitch_model[:len(res) - glitch_idx]

# save simulations
np.save(f"{g.DATA_PATH}143-2a_simulations.npy", res)

plot_window = 1000  # number of samples to plot
plt.figure(figsize=(10, 5))
for i in range(0, len(res), plot_window):
    glt_idx = glitch_indices[(glitch_indices >= i) & (glitch_indices < i + plot_window)]
    print(f"Glitch indices in window {i}-{i + plot_window}: {glt_idx}")
    glt_lbl = glitch_types[glt_idx]

    plt.plot(seconds[i:i + plot_window], res[i:i + plot_window])
    plt.scatter(seconds[i:i + plot_window][glt_idx],
                res[i:i + plot_window][glt_idx],
                color='red', label='Glitches')
    for j, glt in enumerate(glt_idx):
        plt.text(seconds[i:i + plot_window][glt],
                 res[i:i + plot_window][glt],
                 glt_lbl[j], fontsize=8, color='red', rotation=45)
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.title("Glitch Simulation")
    plt.savefig(f"{g.FIGURES_PATH}sims/{i}_{i+plot_window}.png")
    plt.close()
