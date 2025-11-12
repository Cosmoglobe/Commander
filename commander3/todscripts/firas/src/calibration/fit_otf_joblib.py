"""
Advanced parallelized script to fit optical transfer function with joblib.
This version provides better memory management and more control over parallelization.

Usage:
    python fit_otf_joblib.py [--n-jobs N] [--batch-size B] [--backend BACKEND]

Arguments:
    --n-jobs N: Number of parallel jobs (default: all available cores, -1)
    --batch-size B: Batch size for processing (default: 1)
    --backend BACKEND: Joblib backend ('loky', 'multiprocessing', 'threading')
"""

import argparse
import time
from functools import partial

import h5py
import numpy as np
from astropy.io import fits
from joblib import Parallel, delayed
from scipy.optimize import minimize

import globals as g
from calibration import bolometer, electronics
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
    if apod:
        fits_data = fits.open(
            f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
        )
        apod_func = fits_data[1].data["APODIZAT"][0]

    # subtract dither
    if ifg.ndim == 1:
        ifg = ifg - np.median(ifg)

        ifg = ifg / gain / sweeps

        if apod:
            ifg = ifg * apod_func

        ifg = np.roll(ifg, -g.PEAK_POSITIONS[f"{channel}_{mode}"])
        Y = np.fft.rfft(ifg)
        Y = Y[nui]
    else:
        ifg = ifg - np.median(ifg, axis=1)[:, None]

        # "normalize" ifg by gain and sweeps
        ifg = ifg / gain[:, np.newaxis] / sweeps[:, np.newaxis]

        if apod:
            ifg = ifg * apod_func

        ifg = np.roll(ifg, -g.PEAK_POSITIONS[f"{channel}_{mode}"], axis=1)
        Y = np.fft.rfft(ifg, axis=1)
        Y = Y[:, nui]

    B = bolometer.get_bolometer_response_function(
        channel,
        mode,
        bol_cmd_bias / 25.5,
        bol_volt,
        Tbol=temps[9],  # TODO: generalize for other channels
    )
    B = B[:, nui]

    # electronics transfer function
    samprate = 681.43  # from fex_samprate
    Z = electronics.etfunction(channel, adds_per_group, samprate, nui)

    H = Ei[0]

    spec_norm = fnyq_icm * g.FAC_ETENDU * g.FAC_ADC_SCALE

    return Y / H / Z / B / spec_norm * g.FAC_ERG_TO_MJY


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
    return np.sum(result**2)


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
    verbose=False,
):
    """
    Fit emissivities for a single frequency index.
    This function is designed to be called in parallel.
    """
    Ei = np.ones(10, dtype=complex)  # First guess for the emissivities of all emitters

    result = minimize(
        full_function,
        Ei,
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
    )

    if verbose and nui % 10 == 0:
        print(f"Completed frequency {nui}")

    return nui, result.x


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fit optical transfer function in parallel"
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="Number of parallel jobs (-1 for all cores)",
    )
    parser.add_argument(
        "--batch-size", type=int, default=1, help="Batch size for processing"
    )
    parser.add_argument(
        "--backend",
        type=str,
        default="loky",
        choices=["loky", "multiprocessing", "threading"],
        help="Joblib backend to use",
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Print progress messages"
    )
    args = parser.parse_args()

    print("=" * 60)
    print("FIRAS Optical Transfer Function Fitting (Parallel)")
    print("=" * 60)

    data = h5py.File(g.PREPROCESSED_DATA_PATH_CAL, "r")["df_data"]

    # TODO: generalize
    channel = "ll"
    mode = "ss"

    print(f"\nLoading data for channel {channel}, mode {mode}...")
    ifg = data[f"ifg_{channel}"][:]

    # TODO: fit for temp weight coefficients too
    xcal = (data["a_xcal"][:] + data["b_xcal"][:]) / 2
    ical = (data["a_ical"][:] + data["b_ical"][:]) / 2
    dihedral = (data["a_dihedral"][:] + data["b_dihedral"][:]) / 2
    refhorn = (data["a_refhorn"][:] + data["b_refhorn"][:]) / 2
    skyhorn = (data["a_skyhorn"][:] + data["b_skyhorn"][:]) / 2
    collimator = (data["a_collimator"][:] + data["b_collimator"][:]) / 2
    bolometer_ll = (data["a_bol_assem_ll"][:] + data["b_bol_assem_ll"][:]) / 2
    bolometer_lh = (data["a_bol_assem_lh"][:] + data["b_bol_assem_lh"][:]) / 2
    bolometer_rl = (data["a_bol_assem_rl"][:] + data["b_bol_assem_rl"][:]) / 2
    bolometer_rh = (data["a_bol_assem_rh"][:] + data["b_bol_assem_rh"][:]) / 2

    adds_per_group = data["adds_per_group"][:]
    sweeps = data["sweeps"][:]

    bol_cmd_bias = data[f"bol_cmd_bias_{channel}"][:]
    bol_volt = data[f"bol_volt_{channel}"][:]
    gain = data[f"gain_{channel}"][:]

    # we don't want the records where the bol_cmd_bias is lower than 0
    print("Filtering valid records...")
    valid = bol_cmd_bias > 0
    ifg = ifg[valid]
    xcal = xcal[valid]
    ical = ical[valid]
    dihedral = dihedral[valid]
    refhorn = refhorn[valid]
    skyhorn = skyhorn[valid]
    collimator = collimator[valid]
    bolometer_rh = bolometer_rh[valid]
    bolometer_rl = bolometer_rl[valid]
    bolometer_lh = bolometer_lh[valid]
    bolometer_ll = bolometer_ll[valid]
    adds_per_group = adds_per_group[valid]
    sweeps = sweeps[valid]
    bol_cmd_bias = bol_cmd_bias[valid]
    bol_volt = bol_volt[valid]
    gain = gain[valid]

    temps = [
        xcal,
        ical,
        dihedral,
        refhorn,
        skyhorn,
        collimator,
        bolometer_rh,
        bolometer_rl,
        bolometer_lh,
        bolometer_ll,
    ]

    print("Generating frequencies...")
    frequencies = utils.generate_frequencies(channel, mode, 257)

    fnyq = gen_nyquistl(
        "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
    )
    frec = 4 * (g.CHANNELS[channel] % 2) + g.MODES[mode]
    fnyq_icm = fnyq["icm"][frec]

    solution = np.zeros((g.SPEC_SIZE, 10), dtype=complex)  # 10 emissivities to fit

    # Print configuration
    print(f"\nConfiguration:")
    print(f"  Total frequencies: {g.SPEC_SIZE}")
    print(
        f"  Number of jobs: {args.n_jobs} ({'all cores' if args.n_jobs == -1 else f'{args.n_jobs} cores'})"
    )
    print(f"  Batch size: {args.batch_size}")
    print(f"  Backend: {args.backend}")
    print(f"  Verbose: {args.verbose}")
    print()

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
        verbose=args.verbose,
    )

    start0 = time.time()

    # Parallel processing with joblib
    print("Starting parallel fitting process...")
    results = Parallel(
        n_jobs=args.n_jobs, backend=args.backend, batch_size=args.batch_size, verbose=10
    )(delayed(fit_func)(nui) for nui in range(g.SPEC_SIZE))

    # Collect results
    print("\nCollecting results...")
    for nui, emissivities in results:
        solution[nui] = emissivities

    total_time = time.time() - start0
    print(f"\n{'='*60}")
    print(f"Total time taken: {total_time:.2f} seconds ({total_time/60:.2f} minutes)")
    print(f"Average time per frequency: {total_time/g.SPEC_SIZE:.2f} seconds")
    print(
        f"Speedup estimate: {g.SPEC_SIZE * (total_time/g.SPEC_SIZE) / total_time:.1f}x"
    )
    print(f"{'='*60}\n")

    # save solution
    output_file = "fitted_emissivities.npy"
    np.save(output_file, solution)
    print(f"Solution saved to {output_file}")
