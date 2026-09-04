import cProfile

import classification
import detection
import globals as g
import matplotlib.pyplot as plt
import numpy as np
import subtraction
import templates


def main():
    res = np.load(f"{g.DATA_PATH}143-2a_simulations.npy")
    seconds = np.linspace(0, len(res) / g.SAMPRATE, len(res))

    print("Iteration 0")

    print("Starting glitch detection...")
    glitch_idx = detection.matched_filter(res, seconds)
    print(f"Detected {len(glitch_idx)} glitches at indices: {glitch_idx}")

    print("Starting glitch classification...")
    glitch_labels, glitch_amps = classification.classify_glitches(glitch_idx, res, seconds)

    print("Starting glitch subtraction...")
    result, fitted_amps = subtraction.subtract_glitches_from_data(
        glitch_idx, seconds, glitch_labels, glitch_amps, res)

    if g.PLOTS:
        plt.plot(seconds[:1000], res[:1000], label='Original Data')
        plt.plot(seconds[:1000], result[:1000], label='Data after Glitch Subtraction')
        plt.scatter(seconds[glitch_idx], res[glitch_idx], color='red', label='Detected Glitches')
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.title('Glitch Subtraction Result')
        plt.legend()
        plt.xlim(0, seconds[1000])
        plt.savefig(f"{g.FIGURES_PATH}debug/glitch_subtraction_result_0.png")
        plt.close()

    final_stack = templates.stacking(result, glitch_idx, glitch_labels, fitted_amps, seconds)

    for glitch_type in ['short', 'long', 'slow']:
        amps, taus = templates.glitch_estimation(seconds[:int(g.NSECS * g.SAMPRATE)],
                                           final_stack[glitch_type])
        # print(f"Estimated parameters for {glitch_type} glitch: amplitudes={amps}, taus={taus}")
        if g.PLOTS:
            plt.plot(seconds[:int(2*g.SAMPRATE)], final_stack[glitch_type][:int(2*g.SAMPRATE)],
                    label="Stacked median")
            plt.plot(seconds[:int(2*g.SAMPRATE)],
                    templates.glitch_model(seconds[:int(g.NSECS * g.SAMPRATE)], *amps, *taus)[:int(2*g.SAMPRATE)],
                    label="Fitted Model")

            plt.title(f"Glitch Stacking and Fitting for {glitch_type} Glitches")
            plt.xlabel("Time (s)")
            plt.ylabel("Amplitude")
            plt.legend()
            plt.savefig(f"{g.FIGURES_PATH}templates/glitch_fitting_{glitch_type}.png")
            plt.close()

    maxiter = 1
    for i in range(maxiter):
        print(f"Iteration {i+1}/{maxiter}")
        print("Starting glitch detection...")
        glitch_idx = detection.matched_filter(res, seconds)

if __name__ == "__main__":
    # cProfile.run('main()', sort='time')
    main()