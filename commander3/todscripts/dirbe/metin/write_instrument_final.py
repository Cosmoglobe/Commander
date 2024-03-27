from __future__ import annotations

import sys
import numpy as np
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from ....python.commander_tools.tod_tools.commander_instrument import commander_instrument

sys.path.insert(0, "/mn/stornext/d16/cmbco/bp/metins/Commander/commander3/python")
from commander_tools.tod_tools.commander_instrument import commander_instrument
from numpy.typing import NDArray
import astropy.units as u
import dirbe_utils

TEMP_OUTPUT_PATH = "/mn/stornext/d5/data/metins/dirbe/data"
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

    fwhms = dirbe_utils.get_dirbe_fwhm()
    beams = dirbe_utils.get_dirbe_beams()
    sidelobes = dirbe_utils.get_dirbe_sidelobes()
    for band, detector in enumerate(dirbe_utils.BANDS, start=1):
        center_frequency = dirbe_utils.WAVELENGHTS[band - 1]
        center_frequency = (center_frequency*u.micron).to(u.GHz, equivalencies=u.spectral()).value
        wavelengths, weights = dirbe_utils.get_bandpass(band)
        frequencies = wavelengths.to(u.GHz, equivalencies=u.spectral())[::-1].value
        weights = weights[::-1]
        detector_group_name = f"{detector:02}"
        instrument_file.add_bandpass(
            detector_group_name, frequencies, weights
        )

     
        band_group_name = f"{detector:02}_A"
        instrument_file.add_bandpass(
            band_group_name, frequencies, weights
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
            central_wavelength=center_frequency,
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

    version = 2
    write_dirbe_instrument_file(output_path=TEMP_OUTPUT_PATH, version=version)


if __name__ == "__main__":
    main()
