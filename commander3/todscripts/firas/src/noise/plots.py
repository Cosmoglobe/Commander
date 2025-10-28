import matplotlib.pyplot as plt
import numpy as np


def plot_cov_matrix(cov, channel):
    plt.imshow(cov, cmap="RdBu_r", vmax=1, vmin=-1)
    plt.colorbar()
    plt.title(f"Correlation Coefficient Matrix of IFG {channel.upper()}")
    plt.xlabel("Sample Index")
    plt.ylabel("Sample Index")
    plt.savefig(f"./noise/output/cov_{channel}.png")
    plt.close()


def plot_psd(freqs, psd, channel):
    plt.plot(freqs, psd.T)
    plt.xscale("log")
    plt.yscale("log")
    plt.title(f"PSD of IFG {channel.upper()}")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power Spectral Density")
    plt.savefig(f"./noise/output/psd_{channel}.png")
    plt.close()


def plot_mean_psd(psd, freqs, channel):
    sd_mean = psd.mean(axis=0)
    psd_std = psd.std(axis=0)
    plt.fill_between(
        freqs,
        sd_mean - psd_std / (2 * np.sqrt(psd.shape[0])),
        sd_mean + psd_std / (2 * np.sqrt(psd.shape[0])),
        alpha=0.2,
    )
    plt.plot(freqs, sd_mean, color="orange", label="Mean PSD with Std Dev")
    plt.xscale("log")
    plt.yscale("log")
    plt.title(f"Mean PSD of IFG {channel.upper()} with Std Dev")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power Spectral Density")
    plt.legend()
    plt.savefig(f"./noise/output/psd_{channel}_mean_std.png")
    plt.close()


def plot_psd_waterfall(psd, freqs, channel):
    plt.imshow(
        psd,
        aspect="auto",
        extent=[freqs[0], freqs[-1], 0, psd.shape[0]],
        origin="lower",
        cmap="viridis",
    )
    plt.colorbar(label="Power Spectral Density")
    plt.xscale("log")
    plt.yscale("linear")
    plt.title(f"Waterfall Plot of PSD of IFG {channel.upper()}")
    plt.xlabel("Frequency (Hz)")
    plt.show()
    plt.savefig(f"./noise/output/psd_{channel}_waterfall.png")
    plt.close()
