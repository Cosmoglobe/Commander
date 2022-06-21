from __future__ import annotations

import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ...python.commander_tools.tod_tools.commander_instrument import \
        commander_instrument

sys.path.insert(0, "/mn/stornext/d16/cmbco/bp/metins/Commander/commander3/python")
import numpy as np
from commander_tools.tod_tools.commander_instrument import commander_instrument
from numpy.typing import NDArray

import dirbe_utils

TEMP_OUTPUT_PATH = "/mn/stornext/d16/cmbco/bp/metins/dirbe/data/"
NSIDE = 128

# temporary values that needs to be updated
TEMP_LMAX = NSIDE * 3
TEMP_MMAX = 100
TEMP_ELIP = 1
TEMP_PSI_ELL = 0
TEMP_MBEAM_EFF = 0


def write_dirbe_instrument_file(output_path: str, version: int) -> None:
    """Writes the DIRBE filelists for Commander3 using Mathew's script."""

    filename = f"DIRBE_instrument_{version:02}.h5"

    instrument_file = commander_instrument(output_path, filename, version, "w")

    wavelengths, bandpasses = dirbe_utils.get_dirbe_bandpass()
    fwhms = dirbe_utils.get_dirbe_fwhm()
    beams = dirbe_utils.get_dirbe_beams()
    sidelobes = dirbe_utils.get_dirbe_sidelobes()
    for idx, detector in enumerate(dirbe_utils.DETECTORS):
        wavelength = dirbe_utils.WAVELENGHTS[idx]
        detector_group_name = f"{detector:02}_{wavelength}um"
        instrument_file.add_bandpass(
            detector_group_name, wavelengths, bandpasses[idx]
        )

        for band in dirbe_utils.BANDS_LABELS:
            if idx > 2 and band in dirbe_utils.BANDS_LABELS[1:]:
                break

            band_group_name = f"{detector:02}_{band}"
            instrument_file.add_bandpass(
                band_group_name, wavelengths, bandpasses[idx]
            )
            _add_fields(
                instrument_file=instrument_file,
                band_label=band_group_name,
                beam=beams[band_group_name],
                sidelobe=sidelobes[band_group_name],
                fwhm=fwhms[band_group_name],
                elip=TEMP_ELIP,
                psi_ell=TEMP_PSI_ELL,
                mbeam_Eff=TEMP_MBEAM_EFF,
                central_wavelength=wavelength,
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
    central_wavelength: float,
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
    instrument_file.add_field(band_label + "/centFreq", data=[central_wavelength])


def main() -> None:

    version = 1
    write_dirbe_instrument_file(output_path=TEMP_OUTPUT_PATH, version=version)


if __name__ == "__main__":
    main()
