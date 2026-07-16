from __future__ import annotations
import sys
import os
import numpy as np
from functools import cache

# sys.path.insert(0, "/mn/stornext/d16/cmbco/bp/metins/Commander/commander3/python")
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "python"))
from commander_tools.tod_tools.commander_instrument import commander_instrument # type: ignore
from numpy.typing import NDArray
import astropy.units as u

from iras_native_tod_reader import BANDS, BAND_DETS
import matplotlib.pyplot as plt

"""
Write IRAS instrument file 
Essentially modifications of Eirik Gjerløw's instrument script for AKARI
[akari branch]: Commander/commander3/todscripts/AKARI/write_instrument.py
"""

"""
13.07.26 UPDATE
 - Needs to be revisited. Not sure if everything here is used by commander atm
 - nside 512 files written april 29 2026 use lmax=0->n_alms=1. 

"""

OUTPUT_PATH = "/mn/stornext/d23/cmbco/globe/iras/bp"
NSIDE = 512

# temporary values that needs to be updated
TEMP_LMAX       = NSIDE * 3
TEMP_MMAX       = 100
# TEMP_LMAX       = 0
TEMP_MMAX       = 0
TEMP_ELIP       = 1
TEMP_PSI_ELL    = 0
TEMP_MBEAM_EFF  = 0

BEAM_SIZES = {
    # Values taken from IRIS: Table 1, row 2, in Miville-Deschenes, 2004 (https://arxiv.org/abs/astro-ph/0412216)
    "012": 3.8 * u.arcmin,
    "025": 3.8 * u.arcmin,
    "060": 4.0 * u.arcmin,
    "100": 4.3 * u.arcmin,
}

WAVELENGTHS     = {}
CENTER_FREQS    = {}
for band in BANDS:
    WAVELENGTHS[band] = float(band)
    CENTER_FREQS[band] = (float(band)*u.micron).to(u.GHz, equivalencies=u.spectral()).value


def get_iras_fwhm():
    fwhms = {}
    for band, band_det in BAND_DETS.items():
        for det in band_det:
            fwhms[det] = BEAM_SIZES[band].value
    return fwhms

def get_iras_beams(nside: int, lmax: int) -> dict[str, NDArray[np.floating]]:
    """
    Returns a dictionary mapping the DIRBE bands to beams.
    NOTE: Currently only returns a sequence of 0's.
    """
    # NSIDE = 128
    # LMAX = 3 * nside
    N_ALMS = lmax**2 + 2 * lmax + 1
    DEFAULT_BEAM = np.zeros((3, N_ALMS))  # Update this with actual beams
    beams: dict[str, NDArray[np.floating]] = {}
    for band_det_list in BAND_DETS.values():
        for det in band_det_list:
            beams[det] = DEFAULT_BEAM
    return beams

def get_iras_sidelobes(nside: int, lmax: int) -> dict[str, NDArray[np.floating]]:
    """
    Returns a dictionary mapping the DIRBE bands to sidelobes.
    NOTE: We dont have dirbe sidelobes so we just returns a sequence of 0's.
    """
    # NSIDE = 128
    # LMAX = 3 * NSIDE
    N_ALMS = lmax**2 + 2 * lmax + 1
    DEFAULT_SIDELOBE = np.zeros((3, N_ALMS))  # Update this with actual sidelobes
    sidelobes: dict[str, NDArray[np.floating]] = {}
    for band_det_list in BAND_DETS.values():
        for det in band_det_list:
            sidelobes[det] = DEFAULT_SIDELOBE
    return sidelobes


@cache
def get_bandpass_dicts():
    """
    Get freq and weight dicts for each band
    Currently using a delta function bandpass
    with all zeros except for 1 at center wl/freq
    """
    N = 100
    k = N // 2
    def get_bp_wl(center_wl, wl_range):
        start, end = wl_range
        step = (end - start) / (N - 1)
        arr_start = center_wl - k * step
        return arr_start + step * np.arange(N)
    bandpass_lims = {
        # Eyeball estimates of bandpass range for each band
        # from Figure II.C.9 in IRAS expl.suppl.
        "012": [5, 12],
        "025": [12, 40],
        "060": [25, 100],
        "100": [50, 150],
    }

    freqs = {}
    weights = {}
    for band in BANDS:
        center_wl       = WAVELENGTHS[band]
        bp_wavelengths  = get_bp_wl(center_wl, bandpass_lims[band])
        freqs[band]     = (bp_wavelengths * u.micron).to(u.GHz, equivalencies=u.spectral()).value
        weights[band]   = np.where(bp_wavelengths == center_wl, 1, 0)
        # Flip direction of bandpass arrays
        # due to going from wl -> freq.
        freqs[band]     = freqs[band][::-1]
        weights[band]   = weights[band][::-1]

    return freqs, weights

@cache
def get_bandpass(band):
    bp_freqs, bp_weights = get_bandpass_dicts()
    return bp_freqs[band], bp_weights[band]


def write_iras_instrument_file(output_path: str, version: int, nside: int = NSIDE) -> None:
    """Writes the iras filelists for Commander3 using Mathew's script."""

    filename = f"iras_instrument_v{version:02}.h5"

    instrument_file = commander_instrument(output_path, filename, version, "w")

    lmax = nside * 3

    # dict keys for fwhms, beams, sidelobes: "IRAS_{band}-{det}"
    fwhms       = get_iras_fwhm() # fwhm in arcmin. Same for every band 
    beams       = get_iras_beams(nside, lmax) # zeros(3,1) 
    sidelobes   = get_iras_sidelobes(nside, lmax) # zeros(3,1)
    for band, det_list in BAND_DETS.items():
        center_frequency        = CENTER_FREQS[band]
        bp_freqs, bp_weights    = get_bandpass(band)
        instrument_file.add_bandpass(f"IRAS_{band}", bp_freqs, bp_weights)

        for det_label in det_list:
            band_group_name = det_label
            instrument_file.add_bandpass(det_label, bp_freqs, bp_weights)
            _add_fields(
                instrument_file = instrument_file,
                band_label      = det_label,
                beam            = beams[band_group_name],
                sidelobe        = sidelobes[band_group_name],
                fwhm            = fwhms[band_group_name],
                elip            = TEMP_ELIP,
                psi_ell         = TEMP_PSI_ELL,
                mbeam_Eff       = TEMP_MBEAM_EFF,
                central_freq    = center_frequency,
            )

    instrument_file.finalize()


def _add_fields(
    instrument_file: commander_instrument,
    band_label: str,
    beam: NDArray[np.floating],
    sidelobe: NDArray[np.floating],
    fwhm: float,
    elip: float,
    psi_ell: float,
    mbeam_Eff: float,
    central_freq: float,
) -> None:
    """Adds various required fields to the instrument file for a band."""

    # Add beam information
    instrument_file.add_alms(band_label, "beam", TEMP_LMAX, TEMP_MMAX, *beam)

    # Add sidelobe information
    instrument_file.add_alms(band_label, "sl", TEMP_LMAX, TEMP_MMAX, *sidelobe)

    # Add beam parameters
    instrument_file.add_field(band_label + "/fwhm", data=[fwhm])
    instrument_file.add_field(band_label + "/elip", data=[elip])
    instrument_file.add_field(band_label + "/psi_ell", data=[psi_ell])
    instrument_file.add_field(band_label + "/mbeam_eff", data=[mbeam_Eff])

    # Add central wavelength
    instrument_file.add_field(band_label + "/centFreq", data=[central_freq])


def main() -> None:

    version = 0
    # Version 1: First attempt.
    write_iras_instrument_file(output_path=OUTPUT_PATH, version=version)


if __name__ == "__main__":
    main()
