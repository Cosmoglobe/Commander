import globals as g
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def short_glitch(t, A, a0):
    return glitch_model_func(t, A, a0, "143-2a", "short")

def long_glitch(t, A, a0):
    return glitch_model_func(t, A, a0, "143-2a", "long")

def slow_glitch(t, A, a0):
    return glitch_model_func(t, A, a0, "143-2a", "slow")

def glitch_model_func(t, A=1, a0=0, band="143-2a", glitch_type = "short"):
    glitches = np.load("/mn/stornext/d23/cmbco/globe/planck/glitch/glitch_params.npy",
                        allow_pickle=True).item()
    
    glitch_params = glitches[band]
    # print(f"Band: {band}, Glitch Params: {glitch_params}")

    amp = []
    tau = []
    glitch_model = 0
    for amp_i in range(1, 9):
        amp.append(glitch_params[glitch_type][f"Amplitude{amp_i}"])

        if amp[-1] == 0:
            break

        tau.append(glitch_params[glitch_type][f"Tau{amp_i}"])

        # print(f"Band: {band}, Type: {glitch_type}, Amplitude{amp_i}: {amp[-1]}, "
        #       f"Tau{amp_i}: {tau[-1]}")

        glitch_model += amp[-1] * np.exp(-t / tau[-1])

    glitch_model = glitch_model / np.max(glitch_model)  # normalize the glitch model to have a max of 1

    return glitch_model * A + a0

def fit_glitch(res, glitch, secs):
    adjusted_res = res - np.median(res)

    popt, _, _, mesg, _ = curve_fit(glitch_model_func, secs[glitch: glitch + 2 * 180],
                           adjusted_res[glitch: glitch + 2 * 180], p0=[1.0], bounds=(0, np.inf),
                           full_output=True, nan_policy = "raise")

    print(f"Fitted parameters: {popt}")
    # print(f"Infodict: {infodict}")
    print(f"Message: {mesg}")

    plt.scatter(secs[glitch], adjusted_res[glitch], color='red', label='Glitch Sample')
    plt.plot(secs[glitch - 1 * 180:glitch + 3 * 180],
             adjusted_res[glitch - 1 * 180:glitch + 3 * 180],
             label="Data")
    plt.plot(secs[glitch: glitch + 2 * 180],
             glitch_model_func(secs[glitch: glitch + 2 * 180] - secs[glitch], 1),
             label="Glitch model")
    plt.title("Glitch Model")
    plt.legend()
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.savefig(f"{g.SAVE_PATH}glitch_model.png")
    # plt.close()

    # subtract the glitch model from the data
    nr_flag = 5
    adjusted_res[glitch + nr_flag: glitch + 2 * 180] -= glitch_model_func(secs[glitch + nr_flag: glitch + 2 * 180] - secs[glitch], 1)
    adjusted_res[glitch: glitch + nr_flag] = 0

    plt.vlines(secs[glitch:glitch+nr_flag], ymin=np.min(adjusted_res[glitch - 1 * 180:glitch + 3 * 180]),
                   ymax=np.max(adjusted_res[glitch - 1 * 180:glitch + 3 * 180]),
                   color='grey', alpha=0.1, label=f'First {nr_flag} samples are flagged')
    plt.plot(secs[glitch - 1 * 180:glitch + 3 * 180],
             adjusted_res[glitch - 1 * 180:glitch + 3 * 180],
             label="Data - Glitch Model")
    
    plt.title("Data - Glitch Model")
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.ylim(-0.1, 0.2)
    plt.legend()
    plt.savefig(f"{g.SAVE_PATH}data_minus_glitch_model.png")
    plt.close()


if __name__ == "__main__":
    twosecs = np.linspace(0, 2, 2 * 180, endpoint=False)  # glitch model lasts 2 seconds
    tensecs = np.linspace(0, 10, 10 * 180, endpoint=False)  # glitch model lasts 10 seconds

    short = glitch_model_func(tensecs, glitch_type = "short")
    long = glitch_model_func(tensecs, glitch_type = "long")
    slow = glitch_model_func(tensecs, glitch_type = "slow")

    plt.plot(tensecs, short, label="Short Glitch")
    plt.plot(tensecs, long, label="Long Glitch")
    plt.plot(tensecs, slow, label="Slow Glitch")
    plt.title("Glitch Models")
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.savefig(f"{g.FIGURES_PATH}debug/glitch_models.png")