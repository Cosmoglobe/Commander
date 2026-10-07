import classification
import detection
import globals as g
import matplotlib.pyplot as plt
import numpy as np
import subtraction
import templates
import utils


def main():
    res = np.load(f"{g.DATA_PATH}143-2a_simulations.npy")
    sim_types = np.load(f"{g.DATA_PATH}143-2a_simulations_types.npy")
    sim_indices = np.load(f"{g.DATA_PATH}143-2a_simulations_indices.npy")
    seconds = np.linspace(0, len(res) / g.SAMPRATE, len(res))

    # for now, we are excluding the classification etc of the last min of data, so let's exclude those from the sims
    sim_types = sim_types[sim_indices < len(res) - g.SAMPRATE * 60]
    sim_indices = sim_indices[sim_indices < len(res) - g.SAMPRATE * 60]

    print("Iteration 0")

    # print("Starting glitch detection...")
    glitch_idx, _ = detection.matched_filter(res)
    # print(f"Detected {len(glitch_idx)} glitches at indices: {glitch_idx}")
    # match the detected glitches with the simulated ones
    matched_indices = np.intersect1d(glitch_idx, sim_indices)
    print(f"[Iteration 0] Detection accuracy: {(len(matched_indices)) * 100 / len(sim_indices):.2f}%")

    # print("Starting glitch classification...")
    glitch_idx, glitch_labels, glitch_amps = classification.classify_glitches(glitch_idx, res, seconds)

    # count how many wrong classifications
    wrong_classifications = np.sum(glitch_labels[np.isin(glitch_idx, matched_indices)] != sim_types[np.isin(sim_indices, matched_indices)])
    print(f"[Iteration 0] Classification accuracy: {100 * (1 - wrong_classifications / len(glitch_idx)):.2f}%")

    # print("Starting glitch subtraction...")
    result, fitted_amps = subtraction.subtract_glitches_from_data(glitch_idx, seconds,
                                                                  glitch_labels, glitch_amps, res)

    # calculate normalized chi2
    chi2_value = utils.chi2(result)
    print(f"[Iteration 0] Chi2: {int(chi2_value)}")

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

    glitch_params = {}
    for glitch_type in ['short', 'long', 'slow']:
        glitch_params[glitch_type] = templates.glitch_estimation(seconds[:int(g.NSECS * g.SAMPRATE)],
                                           final_stack[glitch_type])
        # print(f"Estimated parameters for {glitch_type} glitch: amplitudes={amps}, taus={taus}")
        if g.PLOTS:
            plt.plot(seconds[:int(2*g.SAMPRATE)], final_stack[glitch_type][:int(2*g.SAMPRATE)],
                    label="Stacked median")
            plt.plot(seconds[:int(2*g.SAMPRATE)],
                    templates.glitch_model(seconds[:int(g.NSECS * g.SAMPRATE)],
                                           *glitch_params[glitch_type]['Amplitude1': 'Amplitude8'],
                                           *glitch_params[glitch_type]['Tau1': 'Tau8'])[:int(2*g.SAMPRATE)],
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
        # print("Starting glitch detection...")
        glitch_idx, score = detection.matched_filter(res, glitch_params['short'])

        matched_indices = np.intersect1d(glitch_idx, sim_indices)
        print(f"[Iteration {i+1}] Detection accuracy: {(len(matched_indices)) * 100 / len(sim_indices):.2f}%")

        if g.PLOTS:
            plt.plot(seconds[:1000], score[:1000])
            visible = glitch_idx[glitch_idx < 1000]
            plt.scatter(seconds[visible], score[visible], color='red', label='Glitches')
            plt.legend()
            plt.title("Matched Filter Result")
            plt.ylabel("Amplitude")
            plt.xlabel("Time (s)")
            plt.savefig(g.FIGURES_PATH + "detection/matched_filter_result_" + str(i) + ".png")
            plt.close()

        glitch_idx, glitch_labels, glitch_amps = classification.classify_glitches(glitch_idx, res,
                                                                                  seconds,
                                                                                  glitch_params,
                                                                                  iter=i+1)

        # print(f"DEBUG: matched_indices: {matched_indices}, glitch_idx: {glitch_idx}, sim_types: {sim_types}, sim_indices: {sim_indices}")

        wrong_classifications = np.sum(glitch_labels[np.isin(glitch_idx, matched_indices)] != sim_types[np.isin(sim_indices, matched_indices)])
        print(f"[Iteration {i+1}] Classification accuracy: {100 * (1 - wrong_classifications / len(glitch_idx)):.2f}%")

if __name__ == "__main__":
    main()