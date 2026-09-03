"""
This script calculates a chi2 from the glitch templates from Guillaume and chooses which type of
glitch it is based on the lowest chi2. It also fits an overall amplitude to each event.
"""
from turtle import color

import globals as g
import matplotlib.pyplot as plt
import numpy as np
import templates
from scipy.optimize import curve_fit


def chi2(data, model):
    return np.sum((data - model) ** 2) / g.SIGMA**2

def classify_glitches(glitch_idx, res, seconds):
    first_samples = int(g.FAST_PART * g.SAMPRATE)
    total_samples = int(g.NSECS * g.SAMPRATE)
    slowpart = np.arange(first_samples, total_samples) / g.SAMPRATE

    _, ax = plt.subplots(2, 1, figsize=(10, 15), sharex=True)
    ax[0].plot(seconds[:1000], res[:1000], label='Original')
    ax[0].scatter(seconds[glitch_idx], res[glitch_idx], color='red', label='Glitches')
    ax[0].set_title("Original TOD with detected Glitches")
    ax[1].set_xlabel("Time (s)")
    ax[1].set_ylim(-0.1, 1.1)
    residual = res.copy()

    glitch_labels = []
    glitch_amps = []

    for glitch_i in glitch_idx:
        data = res[glitch_i + first_samples:glitch_i + total_samples]
        if len(data) != len(slowpart):
            continue
        # seconds_cut = seconds[glitch_i + first_samples:glitch_i + 2 * 180]
        popt = {}
        popt["short"], _ = curve_fit(templates.short_glitch, slowpart, data,
                                  p0=[1], bounds=(0, np.inf))
        popt["long"], _ = curve_fit(templates.long_glitch, slowpart, data, p0=[1],
                                 bounds=(0, np.inf))
        popt["slow"], _ = curve_fit(templates.slow_glitch, slowpart, data, p0=[1],
                                 bounds=(0, np.inf))

        print(f"Glitch at index {glitch_i}: popt_short={popt['short']}, popt_long={popt['long']}, popt_slow={popt['slow']}")

        chi2val = {}
        chi2val["short"] = chi2(data, templates.short_glitch(slowpart, *popt["short"]))
        chi2val["long"] = chi2(data, templates.long_glitch(slowpart, *popt["long"]))
        chi2val["slow"] = chi2(data, templates.slow_glitch(slowpart, *popt["slow"]))

        print(f"Glitch at index {glitch_i}: chi2_short={chi2val['short']}, chi2_long={chi2val['long']}, chi2_slow={chi2val['slow']}")

        glitch_label = min(chi2val, key=chi2val.get)
        glitch_labels.append(glitch_label)
        glitch_amps.append(popt[glitch_label][0])
        print(f"Glitch at index {glitch_i} classified as: {glitch_label}")

        ax[0].plot(slowpart + seconds[glitch_i],
                 templates.glitch_model_func(slowpart, *popt[glitch_label],
                                             glitch_type=glitch_label), color='green', alpha=0.5)
        ax[0].text(seconds[glitch_i], res[glitch_i], glitch_label, fontsize=8, color='red',
                 rotation=45)
        

        
        residual[glitch_i:glitch_i + total_samples] -= templates.glitch_model_func(slowpart, *popt[glitch_label],
                                                     glitch_type=glitch_label)

    ax[1].plot(seconds, residual)
    ax[1].set_title("Residual after subtracting the best-fit glitch model")
    ax[1].set_ylim(-0.1, 1.1)

    plt.xlim(seconds[0], seconds[1000])
    plt.savefig(f"{g.FIGURES_PATH}classification/classified_glitches.png")
    plt.close()

    return glitch_labels, glitch_amps
        
        
        

        

    

        
