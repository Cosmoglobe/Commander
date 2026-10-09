"""
This script calculates a chi2 from the glitch templates from Guillaume and chooses which type of
glitch it is based on the lowest chi2. It also fits an overall amplitude to each event.
"""
from functools import partial

import globals as g
import matplotlib.pyplot as plt
import numpy as np
import templates
from scipy.optimize import curve_fit


def _template(t, A, glitch_type, glitch_params=None):
    # Fixed template (from glitch_params if given) so curve_fit only fits the amplitude A
    return templates.glitch_model_func(t, A, glitch_type=glitch_type,
                                       glitch_params=glitch_params)

def chi2(data, model):
    return np.sum((data - model) ** 2) / g.SIGMA**2

def classify_glitches(glitch_idx, res, seconds, glitch_params=None, iter=0, prev_amps=None,
                      prev_labels=None):
    """glitch_params: optional {glitch_type: template params} from a previous iteration.
    If None, the default templates are used."""
    first_samples = int(g.FAST_PART * g.SAMPRATE)

    # if iter == 0 , we use only 1 second after the glitch is detected to fit
    if glitch_params is None and iter == 0:
        total_samples = int(0.05 * g.SAMPRATE)
    elif glitch_params is not None and iter > 0:
        total_samples = int(g.NSECS * g.SAMPRATE)
    else:
        raise ValueError("Invalid combination of glitch_params and iter")
    slowpart = np.arange(first_samples, total_samples) / g.SAMPRATE
    oneminute = np.arange(0, int(g.NSECS * g.SAMPRATE)) / g.SAMPRATE

    # Drop glitches that are too close to the end of the TOD to provide a full
    # window of data. This must happen before any fitting so that glitch_idx
    # stays aligned, index-for-index, with glitch_labels/glitch_amps below.
    valid_mask = (glitch_idx + total_samples) <= len(res)
    n_dropped = int(np.size(glitch_idx) - np.count_nonzero(valid_mask))
    if n_dropped:
        print(f"Dropping {n_dropped} glitch(es) too close to the end of the TOD "
              "to fit a full template window.")
    glitch_idx = np.asarray(glitch_idx)[valid_mask]

    if g.PLOTS:
        _, ax = plt.subplots(2, 1, figsize=(10, 15), sharex=True)
        ax[0].plot(seconds[:1000], res[:1000], label='Original')
        ax[0].scatter(seconds[glitch_idx], res[glitch_idx], color='red', label='Glitches')
        ax[0].set_title("Original TOD with detected Glitches")
        ax[1].set_xlabel("Time (s)")
        ax[1].set_ylim(-0.1, 1.1)
    residual = res.copy()

    glitch_labels = []
    glitch_amps = []

    glitch_types = ("short", "long", "slow")
    models = {gtype: partial(_template, glitch_type=gtype, glitch_params=glitch_params[gtype] if glitch_params else None)
              for gtype in glitch_types}

    if iter > 0 and prev_amps is not None and prev_labels is not None:
        # subtreact all of the detected glitches from the data before fitting the next iteration
        for glitch_i, glitch_label, glitch_amp in zip(glitch_idx, prev_labels, prev_amps):
            res[glitch_i:glitch_i + len(oneminute)] -= models[glitch_label](oneminute, glitch_amp)

        if g.PLOTS:
            fig_sub, ax_sub = plt.subplots(1, 1)
            ax_sub.plot(seconds, residual, label="Residual after subtracting all detected glitches")
            ax_sub.legend()
            ax_sub.set_ylabel("Residual")
            ax_sub.set_xlabel("Time (s)")
            plt.xlim(seconds[0], seconds[1000])
            fig_sub.savefig(f"{g.FIGURES_PATH}classification/residual_after_subtraction_{iter}.png")
            plt.close(fig_sub)
            quit()

    for glitch_i in glitch_idx:
        if iter > 0 and prev_amps is not None and prev_labels is not None:
            # add back the current glitch so that we can look at the residual timestream with only the current glitch
            res[glitch_i:glitch_i + len(oneminute)] += models[prev_labels[glitch_idx.index(glitch_i)]](oneminute, prev_amps[glitch_idx.index(glitch_i)])

        data = res[glitch_i + first_samples:glitch_i + total_samples]
        popt = {}
        for gtype in glitch_types:
            popt[gtype], _ = curve_fit(models[gtype], slowpart, data, p0=[1],
                                       bounds=(0, np.inf))

        chi2val = {gtype: chi2(data, models[gtype](slowpart, *popt[gtype]))
                   for gtype in glitch_types}

        glitch_label = min(chi2val, key=chi2val.get)
        glitch_labels.append(glitch_label)
        glitch_amps.append(popt[glitch_label][0])

        if g.PLOTS:
            ax[0].plot(oneminute + seconds[glitch_i],
                    models[glitch_label](oneminute, *popt[glitch_label]), color='green', alpha=0.5)
            ax[0].text(seconds[glitch_i], res[glitch_i], glitch_label, fontsize=8, color='red',
                    rotation=45)
        
        residual[glitch_i:glitch_i + len(oneminute)] -= models[glitch_label](oneminute, *popt[glitch_label])

    if g.PLOTS:
        ax[1].plot(seconds, residual)
        ax[1].set_title("Residual after subtracting the best-fit glitch model")
        ax[1].set_ylim(-0.1, 1.1)
        ax[1].set_xlabel("Time (s)")

        plt.xlim(seconds[0], seconds[1000])
        plt.savefig(f"{g.FIGURES_PATH}classification/classified_glitches_{iter}.png")
        plt.close()

    glitch_labels = np.asarray(glitch_labels)
    glitch_idx = np.asarray(glitch_idx)

    return glitch_idx, glitch_labels, glitch_amps