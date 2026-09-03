import classification
import detection
import globals as g
import matplotlib.pyplot as plt
import numpy as np
import subtraction
import templates

res = np.load(f"{g.DATA_PATH}143-2a_simulations.npy")
seconds = np.linspace(0, len(res) / g.SAMPRATE, len(res))

print("Starting glitch detection...")
glitch_idx = detection.matched_filter(res, seconds)
print(f"Detected {len(glitch_idx)} glitches at indices: {glitch_idx}")

print("Starting glitch classification...")
glitch_labels, glitch_amps = classification.classify_glitches(glitch_idx, res, seconds)

# if g.PLOTS:
    # T = subtraction.build_template_matrix(glitch_idx, seconds, glitch_labels, glitch_amps)

    # im = plt.imshow(T, aspect='auto', norm='log')
    # plt.colorbar(im, label='Amplitude')
    # plt.savefig(f"{g.FIGURES_PATH}debug/template_matrix.png")
    # plt.close()

    # b = subtraction.calculate_b(res, T)
    # A = subtraction.calculate_A(T)

    # plt.imshow(A, aspect='auto', norm='log')
    # plt.colorbar(label='Amplitude')
    # plt.title("Matrix A = T^T N^-1 T")
    # plt.savefig(f"{g.FIGURES_PATH}debug/matrix_A.png")
    # plt.close()

print("Starting glitch subtraction...")
result = subtraction.subtract_glitches_from_data(glitch_idx, seconds, glitch_labels, glitch_amps,
                                                 res)

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

# maxiter = 1
# for i in range(maxiter):

templates.stacking(result, glitch_idx, glitch_labels, glitch_amps, seconds)