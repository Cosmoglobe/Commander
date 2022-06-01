"""
Script for generating the DIRBE instrument file.

NOTE: Currently, most values related to the beams are set to 0 or 1. We also dont 
know what to do about the sidelobes, so these are all zeroes aswell. 
"""

from __future__ import annotations

import h5py
import numpy as np
from numpy.typing import NDArray

from dirbe_utils import get_dirbe_bandpass, get_dirbe_fwhm, get_dirbe_beams, get_dirbe_sidelobes


DEFAULT_OUTPUT_FILE = "DIRBE_instrument_v1.h5py"
TEMP_OUTPUT_PATH = "/mn/stornext/d16/cmbco/bp/metins/dirbe/data/"

NSIDE = 128
DETECTORS = range(1, 11)

DEFAULT_LMAX = int(2 * NSIDE)
DEFAULT_MMAX = 100
DEFAULT_MBEAM_EFF = 1
DEFAULT_PSI_ELL = 0
DEFAULT_ELIP = 0
DEFAULT_SLLMAX = 0
DEFAULT_SLMMAX = 0


def write_dirbe_instrumentfile(output_path: str = TEMP_OUTPUT_PATH) -> None:
    """Generates the DIRBE instrumentfile."""

    wavelengths, bandpasses = get_dirbe_bandpass()
    fwhms = get_dirbe_fwhm()
    beams = get_dirbe_beams()
    sidelobes = get_dirbe_sidelobes()

    with h5py.File(output_path + "DIRBE_instrument.h5", "w") as file:
        for idx, detector in enumerate(DETECTORS):
            bandpass = bandpasses[idx]
            central_wavelength = np.trapz(wavelengths*bandpass, wavelengths)/np.trapz(bandpass, wavelengths)

            # Create detector groups
            detector_group = file.create_group(f"{detector:02}")
            detector_group.create_dataset("bandpass", data=bandpass)
            detector_group.create_dataset("bandpassx", data=wavelengths)

            # Create band groups
            if detector <= 3:
                for idx, band in enumerate(["A", "B", "C"]):
                    group = file.create_group(f"{detector:02}{band}")
                    _write_common_groups(
                        group=group, 
                        wavelengths=wavelengths, 
                        bandpass=bandpass, 
                        central_wavelength=central_wavelength,
                        fwhms=fwhms,
                        beams=beams,
                        sidelobes=sidelobes,
                    )
            else:
                group = file.create_group(f"{detector:02}A")
                _write_common_groups(
                    group=group, 
                    wavelengths=wavelengths, 
                    bandpass=bandpass, 
                    central_wavelength=central_wavelength,
                    fwhms=fwhms,
                    beams=beams,
                    sidelobes=sidelobes,
                )


def _write_common_groups(
    group: h5py.Group, 
    wavelengths: NDArray[np.floating], 
    bandpass: NDArray[np.floating],
    central_wavelength: float,
    fwhms: dict[str, float],
    beams: dict[str, NDArray[np.floating]],
    sidelobes: dict[str, NDArray[np.floating]],
    beam_lmax: int = DEFAULT_LMAX,
    beam_mmax: int = DEFAULT_MMAX,
    mbeam_eff: int = DEFAULT_MBEAM_EFF,
    psi_ell: int = DEFAULT_PSI_ELL,
    elip: int = DEFAULT_ELIP,
    sllmax: int = DEFAULT_SLLMAX,
    slmmax: int = DEFAULT_SLMMAX,

) -> None:
    """Writes the common datasets to a sub-band group."""

    detector = group.name.lstrip("/")
    group.create_dataset("bandpass", data=bandpass)
    group.create_dataset("bandpassx", data=wavelengths)
    group.create_dataset("centFreq", data=[central_wavelength])
    group.create_dataset("fwhm", data=[fwhms[detector]])
    group.create_dataset("beamlmax", data=[beam_lmax])
    group.create_dataset("beammmax", data=[beam_mmax])
    group.create_dataset("mbeam_eff", data=[mbeam_eff])
    group.create_dataset("psi_ell", data=[psi_ell])
    group.create_dataset("elip", data=[elip])
    group.create_dataset("sllmax", data=[sllmax])
    group.create_dataset("slmmax", data=[slmmax])

    beam = beams[detector]
    sl = sidelobes[detector]
    beam_group = group.create_group("beam")
    sl_group = group.create_group("sl")
    for idx, signal in enumerate(["T", "E", "B"]):
        beam_group.create_dataset(signal, data=beam[idx])
        sl_group.create_dataset(signal, data=sl[idx])

if __name__ == "__main__":
    write_dirbe_instrumentfile()
