"""
In this script, I want to try out the way that the Planck team used to remove the glitches.
"""

import Commander.commander3.todscripts.hfi.glitches.globals as g
import matplotlib.pyplot as plt
import numpy as np
import Commander.commander3.todscripts.hfi.glitches.templates as templates


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


def matched_filter(res):
    fourier = np.fft.fft(res)

    glitch_model_short = templates.glitch_model_func(np.arange(0, 2, 1 / g.SAMPRATE), 1)
    glitch_fourier = np.fft.fft(glitch_model_short)

    matched_filter_fourier = np.convolve(fourier, glitch_fourier)
    result = np.fft.ifft(matched_filter_fourier)

    glitch_idx = np.where(result[:1000] > 5)
    print(f"Glitch indices: {glitch_idx}")

    plt.plot(res[:1000], label='Original')
    plt.scatter(glitch_idx, res[glitch_idx], color='red', label='Glitches')
    plt.title("Original TOD with detected Glitches")
    plt.savefig(g.FIGURES_PATH + "debug/original_tod.png")
    plt.close()

    plt.plot(result[:1000])
    plt.scatter(glitch_idx, result[glitch_idx], color='red', label='Glitches')
    plt.legend()
    plt.title("Matched Filter Result")
    plt.ylabel("Amplitude")
    plt.savefig(g.FIGURES_PATH + "debug/matched_filter_result.png")
    plt.close()

    

    return result

if __name__ == "__main__":
    res = np.load(f"{g.DATA_PATH}143-2a_simulations.npy")
    matched_filter(res)
    