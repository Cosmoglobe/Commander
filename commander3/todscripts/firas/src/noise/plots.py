import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d


def plot_cov_matrix(cov, channel):
    plt.imshow(cov, cmap="RdBu_r", vmax=1, vmin=-1)
    plt.colorbar()
    plt.title(f"Correlation Coefficient Matrix of IFG {channel.upper()}")
    plt.xlabel("Sample Index")
    plt.ylabel("Sample Index")
    plt.savefig(f"./noise/output/cov_{channel}.png")
    plt.close()


def plot_psd(freqs_hz, freqs_icm, psd, channel, mode, path):
    _, ax = plt.subplots(layout="constrained")
    ax.plot(freqs_hz, psd.T)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(f"PSD of IFG {channel.upper()}{mode.upper()}")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power Spectral Density")

    # Create interpolation for conversion between Hz and icm
    hz_to_icm_func = interp1d(freqs_hz, freqs_icm, fill_value="extrapolate")
    icm_to_hz_func = interp1d(freqs_icm, freqs_hz, fill_value="extrapolate")
    secax = ax.secondary_xaxis("top", functions=(hz_to_icm_func, icm_to_hz_func))
    secax.set_xlabel("Wavenumber (1/cm)")

    # plt.show()
    plt.savefig(f"{path}psd_{channel}{mode}.png")
    plt.close()


def plot_mean_psd(psd, freqs_hz, freqs_icm, channel, mode, path):
    sd_mean = psd.mean(axis=0)
    psd_std = psd.std(axis=0)

    _, ax = plt.subplots(layout="constrained")
    ax.fill_between(
        freqs_hz,
        sd_mean - psd_std / (2 * np.sqrt(psd.shape[0])),
        sd_mean + psd_std / (2 * np.sqrt(psd.shape[0])),
        alpha=0.2,
    )
    ax.plot(freqs_hz, sd_mean, color="orange", label="Mean PSD with Std Dev")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(f"Mean PSD of IFG {channel.upper()}{mode.upper()} with Std Dev")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power Spectral Density")
    ax.legend()

    # Create interpolation for conversion between Hz and icm
    hz_to_icm_func = interp1d(freqs_hz, freqs_icm, fill_value="extrapolate")
    icm_to_hz_func = interp1d(freqs_icm, freqs_hz, fill_value="extrapolate")
    secax = ax.secondary_xaxis("top", functions=(hz_to_icm_func, icm_to_hz_func))
    secax.set_xlabel("Wavenumber (1/cm)")
    plt.savefig(f"{path}psd_{channel}{mode}_mean_std.png")
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
    # plt.show()
    plt.savefig(f"./noise/output/psd_{channel}_waterfall.png")
    plt.close()
