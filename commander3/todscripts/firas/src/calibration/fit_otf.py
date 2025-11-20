"""
Script to fit an optical transfer function to the calibration data. For now, it follows the model:
S_sky/XCAL = 1/OTF * 1/ETF * 1/Bol * Y - R,
where Y is the Fourier transformed interferogram and R is defined by
R = 1/OTF * sum over all emitters i (excluding XCAL) of E_i * P(T_i).
"""

import argparse
import time
from functools import partial
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fmin_l_bfgs_b, minimize

import globals as g
from pipeline import ifg_spec
from utils import my_utils as utils
from utils.config import gen_nyquistl


def D(
    Ei,
    nui,
    ifg,
    channel,
    mode,
    gain,
    sweeps,
    bol_cmd_bias,
    bol_volt,
    temps,
    adds_per_group,
    fnyq_icm,
    apod=False,
):
    # if apod:
    #     fits_data = fits.open(
    #         f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
    #     )
    #     apod_func = fits_data[1].data["APODIZAT"][0]

    # # subtract dither
    # if ifg.ndim == 1:
    #     ifg = ifg - np.median(ifg)

    #     ifg = ifg / gain / sweeps

    #     if apod:
    #         ifg = ifg * apod_func

    #     ifg = np.roll(ifg, -g.PEAK_POSITIONS[f"{channel}_{mode}"])
    #     Y = np.fft.rfft(ifg)
    #     Y = Y[nui]
    # else:
    #     ifg = ifg - np.median(ifg, axis=1)[:, None]

    #     # "normalize" ifg by gain and sweeps
    #     ifg = ifg / gain[:, np.newaxis] / sweeps[:, np.newaxis]

    #     if apod:
    #         ifg = ifg * apod_func

    #     ifg = np.roll(ifg, -g.PEAK_POSITIONS[f"{channel}_{mode}"], axis=1)
    #     Y = np.fft.rfft(ifg, axis=1)
    #     Y = Y[:, nui]

    # B = bolometer.get_bolometer_response_function(
    #     channel,
    #     mode,
    #     bol_cmd_bias / 25.5,
    #     bol_volt,
    #     Tbol=temps[6],  # TODO: generalize for other channels
    # )
    # B = B[:, nui]

    # # electronics transfer function
    # samprate = 681.43  # from fex_samprate
    # Z = electronics.etfunction(channel, adds_per_group, samprate, nui)
    # # plt.plot(Z.real)
    # # plt.plot(Z.imag)
    # # plt.show()

    # H = Ei[0]

    # spec_norm = fnyq_icm * g.FAC_ETENDU * g.FAC_ADC_SCALE

    # print(f"Y: {Y}")
    # print(f"H: {H}")
    # print(f"Z: {Z[Z == 0]}")
    # print(f"Z: {Z.shape}")
    # print(f"B: {B}")
    # return Y / H / Z / B / spec_norm * g.FAC_ERG_TO_MJY

    spec = ifg_spec.ifg_to_spec(
        ifg,
        channel,
        mode,
        adds_per_group,
        bol_cmd_bias,
        bol_volt,
        fnyq_icm,
        otf=Ei[0],
        Tbol=temps[6],  # TODO: generalize for other channels
        apod=apod,
        gain=gain,
        sweeps=sweeps,
        nui=nui,
    )

    return spec


def R(Ei, nui, temps, frequencies):
    """
    R is the function that weights each of the emitters.

    Parameters
    ----------
    H : array_like
        Optical transfer function a.k.a. emissivity of the XCAL.
    Ei : array_like
        2D array of the emissivities over frequency for each of the other nine emitters.
        1: ICAL, 2: dihedral, 3: refhorn, 4: skyhorn, 5: collimator, 6: bolometer_rh, 7: bolometer_rl, 8: bolometer_lh, 9: bolometer_ll
    """
    sum = np.zeros_like(temps[0], dtype=complex)
    for i in range(1, Ei.shape[0]):
        sum += Ei[i] * utils.planck(frequencies[nui], temps[i])

    H = Ei[0]

    return sum / H


def S(nui, frequencies, temps):
    return utils.planck(frequencies[nui], temps[0])


def full_function(
    Ei,
    nui,
    ifg,
    channel,
    mode,
    gain,
    sweeps,
    bol_cmd_bias,
    bol_volt,
    temps,
    adds_per_group,
    fnyq_icm,
    frequencies,
):
    result = (
        D(
            Ei,
            nui,
            ifg,
            channel,
            mode,
            gain,
            sweeps,
            bol_cmd_bias,
            bol_volt,
            temps,
            adds_per_group,
            fnyq_icm,
        )
        - R(Ei, nui, temps, frequencies)
        - S(nui, frequencies, temps)
    )

    return np.sum(np.abs(result) ** 2)  # Use abs to handle complex numbers properly


def real_to_complex(z):
    """
    Real vector of length 2n -> complex of length n.
    Taken from https://stackoverflow.com/questions/51211055/can-scipy-optimize-minimize-functions-of-complex-variables-at-all-and-how.
    """
    return z[: len(z) // 2] + 1j * z[len(z) // 2 :]


def complex_to_real(z):
    """
    Complex vector of length n -> real of length 2n
    Taken from https://stackoverflow.com/questions/51211055/can-scipy-optimize-minimize-functions-of-complex-variables-at-all-and-how.
    """
    return np.concatenate((np.real(z), np.imag(z)))


def full_function_real(
    Ei_real,
    nui,
    ifg,
    channel,
    mode,
    gain,
    sweeps,
    bol_cmd_bias,
    bol_volt,
    temps,
    adds_per_group,
    fnyq_icm,
    frequencies,
):
    """Wrapper that converts real parameters to complex for optimization."""
    Ei = real_to_complex(Ei_real)
    return full_function(
        Ei,
        nui,
        ifg,
        channel,
        mode,
        gain,
        sweeps,
        bol_cmd_bias,
        bol_volt,
        temps,
        adds_per_group,
        fnyq_icm,
        frequencies,
    )


def fit_single_frequency(
    nui,
    ifg,
    channel,
    mode,
    gain,
    sweeps,
    bol_cmd_bias,
    bol_volt,
    temps,
    adds_per_group,
    fnyq_icm,
    frequencies,
):
    """
    Fit emissivities for a single frequency index.
    This function is designed to be called in parallel.
    """
    # Convert complex initial guess to real representation
    Ei_complex = np.ones(temps.shape[0], dtype=complex)
    Ei_real = complex_to_real(Ei_complex)

    # First try L-BFGS-B with better epsilon for gradient approximation
    x, f, d = fmin_l_bfgs_b(
        full_function_real,
        Ei_real,
        args=(
            nui,
            ifg,
            channel,
            mode,
            gain,
            sweeps,
            bol_cmd_bias,
            bol_volt,
            temps,
            adds_per_group,
            fnyq_icm,
            frequencies,
        ),
        approx_grad=True,
        epsilon=1e-8,  # Larger epsilon for more stable gradient approximation
        maxiter=500,
        factr=1e7,  # Moderate tolerance
        pgtol=1e-5,
    )

    # If L-BFGS-B fails with line search error, try Powell method (derivative-free)
    if d.get("warnflag") == 2:
        result = minimize(
            full_function_real,
            Ei_real,
            args=(
                nui,
                ifg,
                channel,
                mode,
                gain,
                sweeps,
                bol_cmd_bias,
                bol_volt,
                temps,
                adds_per_group,
                fnyq_icm,
                frequencies,
            ),
            method="Powell",
            options={"maxiter": 1000, "ftol": 1e-6},
        )
        x = result.x
        success = result.success
        if not success and nui % 10 == 0:
            print(f"Warning: Frequency {nui} - Powell method failed to converge")
    elif d.get("warnflag") != 0 and nui % 10 == 0:
        print(
            f"Warning: Frequency {nui} with code {d.get('warnflag')} - {d.get('task', '')}"
        )

    # Convert result back to complex
    Ei_solution = real_to_complex(x)

    return nui, Ei_solution


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fit optical transfer function.")
    parser.add_argument(
        "--n",
        type=int,
        default=None,
        help="IFG index to fit (if not provided, fit all frequencies).",
    )
    args = parser.parse_args()
    n = args.n

    for channel in g.CHANNELS.keys():
        data = np.load(f"{g.PREPROCESSED_DATA_PATH}cal_{channel}.npz", "r")

        print(f"Data loaded from {g.PREPROCESSED_DATA_PATH}cal_{channel}.npz")
        print("Preparing data for fitting...")

        mtm_length = data["mtm_length"]
        mtm_speed = data["mtm_speed"]

        for mode in g.MODES.keys():
            if mode == "lf" and (channel[1] == "h"):
                continue
            if mode[0] == "s":
                length_filter = mtm_length == 0
            else:
                length_filter = mtm_length == 1
            if mode[1] == "s":
                speed_filter = mtm_speed == 0
            else:
                speed_filter = mtm_speed == 1

            mode_filter = length_filter & speed_filter

            # if n is chosen, load only the needed data
            if n is None:
                ifg = data[f"ifg"][mode_filter]
                xcal = data["xcal_cone"][mode_filter]  # TODO: update this
                ical = data["ical"][mode_filter]
                dihedral = data["dihedral"][mode_filter]
                refhorn = data["refhorn"][mode_filter]
                skyhorn = data["skyhorn"][mode_filter]
                collimator = data["collimator"][mode_filter]
                # TODO: fit for all bolometers
                bol = data["bolometer"][mode_filter]

                adds_per_group = data["adds_per_group"][mode_filter]
                sweeps = data["sweeps"][mode_filter]

                bol_cmd_bias = data[f"bol_cmd_bias"][mode_filter]
                bol_volt = data[f"bol_volt"][mode_filter]
                gain = data[f"gain"][mode_filter]
            else:
                ifg = data[f"ifg"][mode_filter][n]
                xcal = data["xcal_cone"][mode_filter][n]  # TODO: update this
                ical = data["ical"][mode_filter][n]
                dihedral = data["dihedral"][mode_filter][n]
                refhorn = data["refhorn"][mode_filter][n]
                skyhorn = data["skyhorn"][mode_filter][n]
                collimator = data["collimator"][mode_filter][n]
                bol = data["bolometer"][mode_filter][n]

            # plot a random ifg to check
            n = np.random.randint(ifg.shape[0]) if ifg.ndim == 2 else n

            plt.plot(ifg[n] if ifg.ndim == 2 else ifg)
            plt.title(f"Random IFG at index {n}")
            plt.savefig(f"calibration/output/random_ifgs/{n}.png")

            # print(f"data keys: {list(data.keys())}")

            # TODO: fit for temp weight coefficients too
            # xcal = (data["a_xcal"][:] + data["b_xcal"][:]) / 2
            # ical = (data["a_ical"][:] + data["b_ical"][:]) / 2
            # dihedral = (data["a_dihedral"][:] + data["b_dihedral"][:]) / 2
            # refhorn = (data["a_refhorn"][:] + data["b_refhorn"][:]) / 2
            # skyhorn = (data["a_skyhorn"][:] + data["b_skyhorn"][:]) / 2
            # collimator = (data["a_collimator"][:] + data["b_collimator"][:]) / 2
            # bolometer_ll = (data["a_bol_assem_ll"][:] + data["b_bol_assem_ll"][:]) / 2
            # bolometer_lh = (data["a_bol_assem_lh"][:] + data["b_bol_assem_lh"][:]) / 2
            # bolometer_rl = (data["a_bol_assem_rl"][:] + data["b_bol_assem_rl"][:]) / 2
            # bolometer_rh = (data["a_bol_assem_rh"][:] + data["b_bol_assem_rh"][:]) / 2

            temps = np.array([xcal, ical, dihedral, refhorn, skyhorn, collimator, bol])

            frequencies = utils.generate_frequencies(channel, mode, 257)

            fnyq = gen_nyquistl(
                "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
            )
            frec = 4 * (g.CHANNELS[channel] % 2) + g.MODES[mode]
            fnyq_icm = fnyq["icm"][frec]

            solution = np.zeros(
                (g.SPEC_SIZE, temps.shape[0]), dtype=complex
            )  # 10 emissivities to fit

            # Determine number of processes to use
            n_processes = cpu_count()
            print(f"Using {n_processes} CPU cores for parallel processing")

            # Create partial function with fixed arguments
            fit_func = partial(
                fit_single_frequency,
                ifg=ifg,
                channel=channel,
                mode=mode,
                gain=gain,
                sweeps=sweeps,
                bol_cmd_bias=bol_cmd_bias,
                bol_volt=bol_volt,
                temps=temps,
                adds_per_group=adds_per_group,
                fnyq_icm=fnyq_icm,
                frequencies=frequencies,
            )

            start0 = time.time()
            print(f"Starting parallel optimization for {g.SPEC_SIZE} frequencies...")

            # Parallel processing
            with Pool(processes=n_processes) as pool:
                results = pool.map(fit_func, range(g.SPEC_SIZE))

            # Collect results
            failed_fits = 0
            for nui, emissivities in results:
                solution[nui] = emissivities

            end_time = time.time()
            total_seconds = end_time - start0

            print(f"\n{'='*60}")
            print(f"Optimization Complete!")
            print(f"{'='*60}")
            print(
                f"Total time taken: {total_seconds:.2f} seconds ({total_seconds/60:.2f} minutes)"
            )
            print(
                f"Average time per frequency: {total_seconds/g.SPEC_SIZE:.2f} seconds"
            )
            print(
                f"Speed improvement: {149.78/(total_seconds/g.SPEC_SIZE):.1f}x faster than before"
            )

            # save solution
            np.save(
                f"calibration/output/fitted_emissivities_{channel}_{mode}.npy", solution
            )
            n = None
