"""
In this script, I want to try out the way that the Planck team used to remove the glitches.
"""

import globals as g
import matplotlib.pyplot as plt
import numpy as np
import templates as templates
from scipy.signal import find_peaks


def threepointfilter(data, kernel=(0.25, 0.5, 0.25)):
    filtered = np.convolve(data, kernel, mode='same')

    # plt.plot(data[:windowsize], label='Original')
    # plt.plot(filtered[:windowsize], label='Filtered')

    # plt.legend()
    # plt.savefig(save_path + "filtered_tod.png")
    # plt.clf()
    return filtered

def search(filtered, x_axis, flags=None, windowsize=1000, threshold=3.2):
    window = filtered[:windowsize]
    rms = np.std(window)
    if np.any(np.abs(window) > threshold * rms):
        print("Glitch detected in the first window!")

        # save which samples are above the threshold
        glitch_samples = np.where(np.abs(window) > threshold * rms)[0]
        for sample in glitch_samples:
            if window[sample] <= 0:
                glitch_samples = glitch_samples[glitch_samples != sample]
        print("Glitch samples:", glitch_samples)

        if flags is not None:
            plt.vlines(flags[:windowsize], ymin=np.min(window), ymax=np.max(window),
                       color='grey', alpha=0.5, label='Flags')
        plt.scatter(x_axis[glitch_samples], window[glitch_samples], color='red', label='Glitch Samples')
        plt.plot(x_axis[:windowsize], window, label='Filtered')
        plt.xlabel('Time (s)')
        plt.ylabel('Signal')
        plt.title('Glitch Detection')
        plt.legend()
        plt.savefig(g.SAVE_PATH + "glitch_samples.png")
        plt.clf()

        biggest_glitch = np.argmax(window[glitch_samples])
        biggest_glitch = glitch_samples[biggest_glitch]
        print("Biggest glitch sample:", biggest_glitch)
        # first_glitch = glitch_samples[biggest_glitch]
        # plot first glitch and previous 1 sec and next 3 seconds
        plt.scatter(x_axis[biggest_glitch],
                    window[biggest_glitch], color='red',
                    label='Glitch Samples')
        plt.plot(x_axis[biggest_glitch - 1 * 180:biggest_glitch + 3 * 180],
                 window[biggest_glitch - 1 * 180:biggest_glitch + 3 * 180], label = 'Filtered')
        plt.xlabel('Time (s)')
        plt.ylabel('Signal')
        plt.legend()
        plt.title('First Glitch and previous 1 Second and Next 3 Seconds')
        plt.savefig(g.SAVE_PATH + "biggest_glitch.png")
        plt.close()

        return biggest_glitch


def matched_filter(res, seconds):
    res = np.asarray(res, dtype=float)
    glitch_model = templates.glitch_model_func(np.arange(0, g.NSECS, 1 / g.SAMPRATE), 1)
    template_norm = np.linalg.norm(glitch_model)
    if len(res) < len(glitch_model):
        raise ValueError("The TOD must be longer than the glitch template")

    # Correlate each possible onset with the causal template.  Reversing the
    # template converts correlation into convolution; the valid slice has
    # one score for each possible template start.
    fft_length = len(res) + len(glitch_model) - 1
    correlation = np.fft.irfft(np.fft.rfft(res, fft_length) * np.fft.rfft(glitch_model[::-1],
                                                                          fft_length), fft_length)
    score = correlation[len(glitch_model) - 1:len(res)] / (g.SIGMA * template_norm)

    peak_indices, _ = find_peaks(score, height=g.CUT_OFF, distance=max(1, int(0.05 * g.SAMPRATE)))
    glitch_idx = peak_indices
    print(f"Glitch indices: {(glitch_idx[glitch_idx < 1000])}")

    if g.PLOTS:
        plt.plot(seconds[:1000], res[:1000], label='Original')
        visible = glitch_idx[glitch_idx < 1000]
        plt.scatter(seconds[visible], res[visible], color='red', label='Glitches')
        plt.title("Original TOD with detected Glitches")
        plt.savefig(g.FIGURES_PATH + "debug/original_tod.png")
        plt.close()

        plt.plot(seconds[:1000], score[:1000])
        visible = glitch_idx[glitch_idx < 1000]
        plt.scatter(seconds[visible], score[visible], color='red', label='Glitches')
        plt.legend()
        plt.title("Matched Filter Result")
        plt.ylabel("Amplitude")
        plt.xlabel("Time (s)")
        plt.savefig(g.FIGURES_PATH + "debug/matched_filter_result.png")
        plt.close()

    # Keep the historical np.where-style return shape used by main.py.
    return glitch_idx

if __name__ == "__main__":
    res = np.load(f"{g.DATA_PATH}143-2a_simulations.npy")
    matched_filter(res)
    