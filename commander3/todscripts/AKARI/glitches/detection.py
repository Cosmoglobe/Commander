"""
In this script, I want to try out the way that the Planck team used to remove the glitches.
"""

import matplotlib.pyplot as plt
import numpy as np

save_path = "/mn/stornext/d23/cmbco/akari/aimartin/figures/"

def threepointfilter(data, kernel=(0.25, 0.5, 0.25)):
    filtered = np.convolve(data, kernel, mode='same')

    # plt.plot(data[:windowsize], label='Original')
    # plt.plot(filtered[:windowsize], label='Filtered')

    # plt.legend()
    # plt.savefig(save_path + "filtered_tod.png")
    # plt.clf()
    return filtered

def search(filtered, flags, windowsize=1000, threshold=3.2):
    window = filtered[:windowsize]
    rms = np.std(window)
    if np.any(np.abs(window) > threshold * rms):
        print("Glitch detected in the first window!")

        # save which samples are above the threshold
        glitch_samples = np.where(np.abs(window) > threshold * rms)[0]
        print("Glitch samples:", glitch_samples)

        plt.vlines(flags[:windowsize], ymin=np.min(window), ymax=np.max(window),
                   color='grey', alpha=0.5, label='Flags')
        plt.scatter(glitch_samples, window[glitch_samples], color='red', label='Glitch Samples')
        plt.plot(window, label='Filtered')
        plt.legend()
        plt.savefig(save_path + "glitch_samples.png")
        plt.clf()