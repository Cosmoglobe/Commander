from functools import cache

import globals as g
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


@cache
def _load_glitch_params():
    return np.load("/mn/stornext/d23/cmbco/globe/planck/glitch/glitch_params.npy",
                   allow_pickle=True).item()


@cache
def glitch_model_samples(sample_count, sample_rate, band, glitch_type):
    """Return a unit-amplitude glitch template on the sample grid."""
    t = np.arange(sample_count) / sample_rate
    glitch_params = _load_glitch_params()[band][glitch_type]
    model = np.zeros(sample_count)
    for amp_i in range(1, 9):
        amplitude = glitch_params[f"Amplitude{amp_i}"]
        if amplitude == 0:
            break
        model += amplitude * np.exp(-t / glitch_params[f"Tau{amp_i}"])
    return model / np.max(model)


def short_glitch(t, A):
    return glitch_model_func(t, A, "143-2a", "short")

def long_glitch(t, A):
    return glitch_model_func(t, A, "143-2a", "long")

def slow_glitch(t, A):
    return glitch_model_func(t, A, "143-2a", "slow")

def glitch_model_func(t, A=1, band="143-2a", glitch_type = "short"):
    t = np.asarray(t)
    glitch_params = _load_glitch_params()[band][glitch_type]

    amp = []
    tau = []
    glitch_model = np.zeros_like(t, dtype=float)
    for amp_i in range(1, 9):
        amp.append(glitch_params[f"Amplitude{amp_i}"])

        if amp[-1] == 0:
            break

        tau.append(glitch_params[f"Tau{amp_i}"])

        # print(f"Band: {band}, Type: {glitch_type}, Amplitude{amp_i}: {amp[-1]}, "
        #       f"Tau{amp_i}: {tau[-1]}")

        glitch_model += amp[-1] * np.exp(-t / tau[-1])

    glitch_model = glitch_model / np.max(glitch_model)  # normalize the glitch model to have a max of 1

    return glitch_model * A

def fit_glitch(res, glitch, secs):
    adjusted_res = res - np.median(res)

    popt, _, _, mesg, _ = curve_fit(glitch_model_func, secs[glitch: glitch + 2 * 180],
                           adjusted_res[glitch: glitch + 2 * 180], p0=[1.0], bounds=(0, np.inf),
                           full_output=True, nan_policy = "raise")

    print(f"Fitted parameters: {popt}")
    # print(f"Infodict: {infodict}")
    print(f"Message: {mesg}")

    # if g.PLOTS:
    #     plt.scatter(secs[glitch], adjusted_res[glitch], color='red', label='Glitch Sample')
    #     plt.plot(secs[glitch - 1 * 180:glitch + 3 * 180],
    #             adjusted_res[glitch - 1 * 180:glitch + 3 * 180],
    #             label="Data")
    #     plt.plot(secs[glitch: glitch + 2 * 180],
    #             glitch_model_func(secs[glitch: glitch + 2 * 180] - secs[glitch], 1),
    #             label="Glitch model")
    #     plt.title("Glitch Model")
    #     plt.legend()
    #     plt.xlabel("Time (s)")
    #     plt.ylabel("Amplitude")
    #     plt.legend()
    #     plt.savefig(f"{g.SAVE_PATH}glitch_model.png")
    # plt.close()

    # subtract the glitch model from the data
    nr_flag = 5
    adjusted_res[glitch + nr_flag: glitch + 2 * 180] -= glitch_model_func(secs[glitch + nr_flag: glitch + 2 * 180] - secs[glitch], 1)
    adjusted_res[glitch: glitch + nr_flag] = 0

    # if g.PLOTS:
    #     plt.vlines(secs[glitch:glitch+nr_flag], ymin=np.min(adjusted_res[glitch - 1 * 180:glitch + 3 * 180]),
    #                 ymax=np.max(adjusted_res[glitch - 1 * 180:glitch + 3 * 180]),
    #                 color='grey', alpha=0.1, label=f'First {nr_flag} samples are flagged')
    #     plt.plot(secs[glitch - 1 * 180:glitch + 3 * 180],
    #             adjusted_res[glitch - 1 * 180:glitch + 3 * 180],
    #             label="Data - Glitch Model")
        
    #     plt.title("Data - Glitch Model")
    #     plt.xlabel("Time (s)")
    #     plt.ylabel("Amplitude")
    #     plt.ylim(-0.1, 0.2)
    #     plt.legend()
    #     plt.savefig(f"{g.SAVE_PATH}data_minus_glitch_model.png")
    #     plt.close()


if __name__ == "__main__":
    twosecs = np.linspace(0, 2, 2 * 180, endpoint=False)  # glitch model lasts 2 seconds
    tensecs = np.linspace(0, 10, 10 * 180, endpoint=False)  # glitch model lasts 10 seconds

    short = glitch_model_func(tensecs, glitch_type = "short")
    long = glitch_model_func(tensecs, glitch_type = "long")
    slow = glitch_model_func(tensecs, glitch_type = "slow")

    # if g.PLOTS:
    #     plt.plot(tensecs, short, label="Short Glitch")
    #     plt.plot(tensecs, long, label="Long Glitch")
    #     plt.plot(tensecs, slow, label="Slow Glitch")
    #     plt.title("Glitch Models")
    #     plt.xlabel("Time (s)")
    #     plt.ylabel("Amplitude")
    #     plt.legend()
    #     plt.savefig(f"{g.FIGURES_PATH}debug/glitch_models.png")

def stacking(result, glitch_idx, glitch_labels, glitch_amps, seconds):
    window = int(g.NSECS * g.SAMPRATE)
    stack = {'short': np.zeros(window), 'long': np.zeros(window), 'slow': np.zeros(window)}
    count = {'short': 0, 'long': 0, 'slow': 0}

    for i, idx in enumerate(glitch_idx):
        if idx + window > len(result):
            continue
        model = glitch_model_func(seconds[:window], glitch_amps[i], glitch_type=glitch_labels[i])

        stack[glitch_labels[i]] += result[idx:idx + window] + model
        count[glitch_labels[i]] += 1

        if i == 0:
            plot_res = result.copy()
            plot_res[idx:idx + window] += model
            plt.plot(seconds[:window], plot_res[:window], label="Data + Glitch Model")
            plt.scatter(seconds[idx], plot_res[idx], color='red',
                        label='First Glitch Sample')
            plt.title("First Glitch and Model")
            plt.xlabel("Time (s)")
            plt.ylabel("Amplitude")
            plt.legend()
            plt.xlim(0, seconds[1000])
            plt.savefig(f"{g.FIGURES_PATH}templates/first_glitch_and_model.png")
            plt.close()

    for glitch_type in ['short', 'long', 'slow']:
        if count[glitch_type] == 0:
            continue

        stack[glitch_type] /= count[glitch_type]
        plt.plot(seconds[:window], stack[glitch_type], label=glitch_type)
        plt.title("Stacked Average")
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude")
        plt.legend()
        plt.xlim(-0.1, 2)
    plt.savefig(f"{g.FIGURES_PATH}templates/stacked_average_{glitch_type}.png")
    plt.close()