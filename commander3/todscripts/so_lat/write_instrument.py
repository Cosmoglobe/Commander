from __future__ import annotations

import sys
import os
import numpy as np

# from typing import TYPE_CHECKING
# if TYPE_CHECKING:
#    from ....python.commander_tools.tod_tools.commander_instrument import commander_instrument

# sys.path.insert(0, "/mn/stornext/d16/cmbco/bp/metins/Commander/commander3/python")
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "python"))
from commander_tools.tod_tools.commander_instrument import commander_instrument
from numpy.typing import NDArray
import astropy.units as u

TEMP_OUTPUT_PATH = "/mn/stornext/d5/data/duncanwa/SO/"
NSIDE = 8192

# temporary values that needs to be updated
TEMP_LMAX = NSIDE * 3
TEMP_MMAX = 100
TEMP_LMAX = 0
TEMP_MMAX = 0
TEMP_ELIP = 1
TEMP_PSI_ELL = 0
TEMP_MBEAM_EFF = 0

BANDS = [
        '090'
        ]


# in arcmin
FWHMS = {
        '030' : 7.4,
        '040' : 5.1,
        '090' : 2.2,
        '145' : 1.4,
        '225' : 1.0,
        '280' : 0.9
        }

FREQS = {
        '090': 93.
        }



prefix = 'Mv28_f090'
det_suffix = [
'Ar00c03A', 'Ar00c04B', 'Ar00c05A', 'Ar00c05B', 'Ar00c06A', 'Ar00c06B', 'Ar00c07A', 'Ar00c07B', 'Ar00c08A', 'Ar00c08B', 'Ar00c09A', 'Ar00c09B', 'Ar00c10B', 'Ar00c11A', 'Ar00c11B', 'Ar01c01A', 'Ar01c03B', 'Ar01c04A', 'Ar01c05B', 'Ar01c06A', 'Ar01c07A', 'Ar01c07B', 'Ar01c08A', 'Ar01c08B', 'Ar01c09A', 'Ar01c09B', 'Ar01c10A', 'Ar01c10B', 'Ar01c11A', 'Ar01c11B', 'Ar02c04A', 'Ar02c05A', 'Ar02c06A', 'Ar02c06B', 'Ar02c07A', 'Ar02c08A', 'Ar02c08B', 'Ar02c09B', 'Ar02c10B', 'Ar02c11A', 'Ar02c11B', 'Ar03c03B', 'Ar03c05A', 'Ar03c05B', 'Ar03c06A', 'Ar03c06B', 'Ar03c07A', 'Ar03c07B', 'Ar03c08A', 'Ar03c08B', 'Ar03c09A', 'Ar03c09B', 'Ar03c10A', 'Ar03c10B', 'Ar03c11B', 'Ar04c00A', 'Ar04c04A', 'Ar04c05B', 'Ar04c06A', 'Ar04c06B', 'Ar04c07A', 'Ar04c07B', 'Ar04c08A', 'Ar04c08B', 'Ar04c09A', 'Ar04c09B', 'Ar04c10A', 'Ar04c10B', 'Ar04c11A', 'Ar04c11B', 'Ar05c00A', 'Ar05c01A', 'Ar05c01B', 'Ar05c03A', 'Ar05c03B', 'Ar05c04A', 'Ar05c05A', 'Ar05c05B', 'Ar05c06A', 'Ar05c06B', 'Ar05c07A', 'Ar05c08A', 'Ar05c08B',
'Ar05c09A', 'Ar05c10A', 'Ar05c11A', 'Ar05c11B', 'Ar06c00A', 'Ar06c00B', 'Ar06c01A', 'Ar06c02A', 'Ar06c02B', 'Ar06c03A', 'Ar06c03B', 'Ar06c04A', 'Ar06c04B', 'Ar06c05A', 'Ar06c05B', 'Ar06c06A', 'Ar06c06B', 'Ar06c07A', 'Ar06c07B', 'Ar06c08A', 'Ar06c08B', 'Ar06c09A', 'Ar06c09B', 'Ar06c10A', 'Ar06c10B', 'Ar06c11A', 'Ar06c11B', 'Ar07c00A', 'Ar07c00B', 'Ar07c01A', 'Ar07c01B', 'Ar07c02A', 'Ar07c02B', 'Ar07c03B', 'Ar07c04A', 'Ar07c05A', 'Ar07c05B', 'Ar07c06B', 'Ar07c07A', 'Ar07c08A', 'Ar07c09B', 'Ar07c10A', 'Ar07c10B', 'Ar07c11A', 'Ar07c11B', 'Ar08c00A', 'Ar08c00B', 'Ar08c01A', 'Ar08c01B', 'Ar08c02A', 'Ar08c02B', 'Ar08c03A', 'Ar08c03B', 'Ar08c04A', 'Ar08c04B', 'Ar08c05A', 'Ar08c06A', 'Ar08c06B', 'Ar08c07A', 'Ar08c07B', 'Ar08c08A', 'Ar08c08B', 'Ar08c09A', 'Ar08c09B', 'Ar08c10A', 'Ar08c10B', 'Ar08c11A', 'Ar08c11B', 'Ar09c00A', 'Ar09c00B', 'Ar09c01A', 'Ar09c01B', 'Ar09c02A', 'Ar09c02B', 'Ar09c03B', 'Ar09c04A', 'Ar09c05A', 'Ar09c05B', 'Ar09c06A', 'Ar09c07A', 'Ar09c07B', 'Ar09c08A', 'Ar09c08B',
'Ar09c09A', 'Ar09c09B', 'Ar09c10A', 'Ar09c10B', 'Ar09c11A', 'Ar09c11B', 'Ar10c00A', 'Ar10c00B', 'Ar10c01A', 'Ar10c02A', 'Ar10c02B', 'Ar10c03B', 'Ar10c04A', 'Ar10c05A', 'Ar10c05B', 'Ar10c06A', 'Ar10c06B', 'Ar10c07A', 'Ar10c07B', 'Ar10c08A', 'Ar10c09A', 'Ar10c09B', 'Ar10c11B', 'Ar11c00A', 'Ar11c00B', 'Ar11c01A', 'Ar11c01B', 'Ar11c02B', 'Ar11c03A', 'Ar11c04A', 'Ar11c04B', 'Ar11c05B', 'Ar11c06A', 'Ar11c07B', 'Ar11c08A', 'Ar11c08B', 'Ar11c09A', 'Ar11c09B', 'Ar11c10A', 'Ar11c11B', 'Br00c02A', 'Br00c04B', 'Br00c05A', 'Br00c05B', 'Br00c06A', 'Br00c06B', 'Br00c07A', 'Br00c07B', 'Br00c08A', 'Br00c08B', 'Br00c09A', 'Br00c09B', 'Br00c10A', 'Br00c10B', 'Br00c11A', 'Br00c11B', 'Br01c01A', 'Br01c02B', 'Br01c03A', 'Br01c04A', 'Br01c04B', 'Br01c05A', 'Br01c05B', 'Br01c06B', 'Br01c07A', 'Br01c07B', 'Br01c08A', 'Br01c08B', 'Br01c09A', 'Br01c10A', 'Br01c10B', 'Br01c11A', 'Br01c11B', 'Br02c00A', 'Br02c02B', 'Br02c03A', 'Br02c04A', 'Br02c04B', 'Br02c05A', 'Br02c05B', 'Br02c06A', 'Br02c06B', 'Br02c07A',
'Br02c07B', 'Br02c08A', 'Br02c08B', 'Br02c09A', 'Br02c09B', 'Br02c10A', 'Br02c10B', 'Br02c11A', 'Br02c11B', 'Br03c00A', 'Br03c02A', 'Br03c03A', 'Br03c03B', 'Br03c04A', 'Br03c04B', 'Br03c05A', 'Br03c05B', 'Br03c06A', 'Br03c06B', 'Br03c07A', 'Br03c07B', 'Br03c08A', 'Br03c08B', 'Br03c09A', 'Br03c10A', 'Br03c11A', 'Br03c11B', 'Br04c00A', 'Br04c01A', 'Br04c02B', 'Br04c04A', 'Br04c04B', 'Br04c05A', 'Br04c05B', 'Br04c06A', 'Br04c06B', 'Br04c07A', 'Br04c07B', 'Br04c08A', 'Br04c08B', 'Br04c09A', 'Br04c09B', 'Br04c10A', 'Br04c10B', 'Br04c11A', 'Br04c11B', 'Br05c00A', 'Br05c01A', 'Br05c01B', 'Br05c02A', 'Br05c02B', 'Br05c03A', 'Br05c04A', 'Br05c04B', 'Br05c05A', 'Br05c05B', 'Br05c06A', 'Br05c06B', 'Br05c07A', 'Br05c07B', 'Br05c08A', 'Br05c08B', 'Br05c09A', 'Br05c09B', 'Br05c10A', 'Br05c10B', 'Br05c11A', 'Br05c11B', 'Br06c00A', 'Br06c00B', 'Br06c01A', 'Br06c02A', 'Br06c02B', 'Br06c03A', 'Br06c03B', 'Br06c04A', 'Br06c04B', 'Br06c05A', 'Br06c05B', 'Br06c07A', 'Br06c07B', 'Br06c08A', 'Br06c09A',
'Br06c09B', 'Br06c11A', 'Br06c11B', 'Br07c00A', 'Br07c00B', 'Br07c01A', 'Br07c01B', 'Br07c02A', 'Br07c03A', 'Br07c04A', 'Br07c04B', 'Br07c05A', 'Br07c06A', 'Br07c06B', 'Br07c07A', 'Br07c07B', 'Br07c08A', 'Br07c08B', 'Br07c09A', 'Br07c09B', 'Br07c10A', 'Br07c10B', 'Br07c11A', 'Br07c11B', 'Br08c00A', 'Br08c00B', 'Br08c01A', 'Br08c01B', 'Br08c02A', 'Br08c03A', 'Br08c04A', 'Br08c05A', 'Br08c05B', 'Br08c06A', 'Br08c07A', 'Br08c07B', 'Br08c08A', 'Br08c08B', 'Br08c09A', 'Br08c09B', 'Br08c10A', 'Br08c10B', 'Br08c11A', 'Br08c11B', 'Br09c00A', 'Br09c00B', 'Br09c01A', 'Br09c02A', 'Br09c02B', 'Br09c03A', 'Br09c03B', 'Br09c04A', 'Br09c04B', 'Br09c05A', 'Br09c05B', 'Br09c06A', 'Br09c06B', 'Br09c07A', 'Br09c07B', 'Br09c08A', 'Br09c08B', 'Br09c09A', 'Br09c09B', 'Br09c10A', 'Br09c10B', 'Br09c11A', 'Br09c11B', 'Br10c00A', 'Br10c00B', 'Br10c01A', 'Br10c01B', 'Br10c02A', 'Br10c02B', 'Br10c03A', 'Br10c03B', 'Br10c04A', 'Br10c04B', 'Br10c05B', 'Br10c06A', 'Br10c07A', 'Br10c07B', 'Br10c08A', 'Br10c08B',
'Br10c09B', 'Br10c10A', 'Br10c10B', 'Br10c11A', 'Br10c11B', 'Br11c01B', 'Br11c02A', 'Br11c02B', 'Br11c03A', 'Br11c03B', 'Br11c04A', 'Br11c04B', 'Br11c05A', 'Br11c05B', 'Br11c06A', 'Br11c06B', 'Br11c07A', 'Br11c07B', 'Br11c08A', 'Br11c08B', 'Br11c10A', 'Br11c10B', 'Br11c11B', 'Cr00c02A', 'Cr00c02B', 'Cr00c04A', 'Cr00c04B', 'Cr00c05A', 'Cr00c05B', 'Cr00c06A', 'Cr00c06B', 'Cr00c07A', 'Cr00c07B', 'Cr00c08A', 'Cr00c08B', 'Cr00c09A', 'Cr00c09B', 'Cr00c10A', 'Cr00c11A', 'Cr00c11B', 'Cr01c03B', 'Cr01c04B', 'Cr01c05A', 'Cr01c05B', 'Cr01c06A', 'Cr01c06B', 'Cr01c07B', 'Cr01c08A', 'Cr01c09A', 'Cr01c10A', 'Cr01c10B', 'Cr01c11A', 'Cr01c11B', 'Cr02c00A', 'Cr02c02A', 'Cr02c03A', 'Cr02c04A', 'Cr02c04B', 'Cr02c05B', 'Cr02c06A', 'Cr02c06B', 'Cr02c07A', 'Cr02c07B', 'Cr02c08A', 'Cr02c08B', 'Cr02c09A', 'Cr02c09B', 'Cr02c10A', 'Cr02c11A', 'Cr02c11B', 'Cr03c00A', 'Cr03c00B', 'Cr03c01A', 'Cr03c01B', 'Cr03c02B', 'Cr03c03A', 'Cr03c04A', 'Cr03c04B', 'Cr03c05A', 'Cr03c05B', 'Cr03c06A', 'Cr03c06B', 'Cr03c07A',
'Cr03c07B', 'Cr03c08A', 'Cr03c09A', 'Cr03c09B', 'Cr03c10A', 'Cr03c10B', 'Cr03c11A', 'Cr03c11B', 'Cr04c00A', 'Cr04c00B', 'Cr04c01A', 'Cr04c01B', 'Cr04c02A', 'Cr04c02B', 'Cr04c03A', 'Cr04c03B', 'Cr04c04A', 'Cr04c04B', 'Cr04c05A', 'Cr04c05B', 'Cr04c06A', 'Cr04c06B', 'Cr04c07A', 'Cr04c07B', 'Cr04c08A', 'Cr04c08B', 'Cr04c09A', 'Cr04c10A', 'Cr04c11A', 'Cr04c11B', 'Cr05c00A', 'Cr05c01A', 'Cr05c01B', 'Cr05c02A', 'Cr05c02B', 'Cr05c03A', 'Cr05c03B', 'Cr05c04A', 'Cr05c04B', 'Cr05c05A', 'Cr05c05B', 'Cr05c06A', 'Cr05c06B', 'Cr05c07A', 'Cr05c07B', 'Cr05c08A', 'Cr05c08B', 'Cr05c09A', 'Cr05c09B', 'Cr05c10A', 'Cr05c10B', 'Cr05c11A', 'Cr05c11B', 'Cr06c00B', 'Cr06c01A', 'Cr06c01B', 'Cr06c02B', 'Cr06c03A', 'Cr06c03B', 'Cr06c04A', 'Cr06c04B', 'Cr06c05A', 'Cr06c05B', 'Cr06c06A', 'Cr06c06B', 'Cr06c07A', 'Cr06c07B', 'Cr06c08A', 'Cr06c08B', 'Cr06c09A', 'Cr06c09B', 'Cr06c10A', 'Cr06c10B', 'Cr06c11A', 'Cr06c11B', 'Cr07c00A', 'Cr07c01A', 'Cr07c02A', 'Cr07c02B', 'Cr07c03A', 'Cr07c03B', 'Cr07c04A', 'Cr07c04B',
'Cr07c05A', 'Cr07c05B', 'Cr07c06A', 'Cr07c06B', 'Cr07c07A', 'Cr07c08A', 'Cr07c08B', 'Cr07c09A', 'Cr07c09B', 'Cr07c10B', 'Cr07c11A', 'Cr07c11B', 'Cr08c00A', 'Cr08c00B', 'Cr08c01A', 'Cr08c01B', 'Cr08c02A', 'Cr08c02B', 'Cr08c03A', 'Cr08c03B', 'Cr08c04A', 'Cr08c05A', 'Cr08c05B', 'Cr08c06A', 'Cr08c07A', 'Cr08c07B', 'Cr08c08A', 'Cr08c08B', 'Cr08c09B', 'Cr08c10A', 'Cr08c10B', 'Cr08c11A', 'Cr08c11B', 'Cr09c00A', 'Cr09c00B', 'Cr09c01A', 'Cr09c01B', 'Cr09c02A', 'Cr09c02B', 'Cr09c03A', 'Cr09c03B', 'Cr09c04A', 'Cr09c04B', 'Cr09c05A', 'Cr09c05B', 'Cr09c06A', 'Cr09c06B', 'Cr09c07A', 'Cr09c08A', 'Cr09c09A', 'Cr09c09B', 'Cr09c10A', 'Cr09c10B', 'Cr09c11A', 'Cr09c11B', 'Cr10c00A', 'Cr10c00B', 'Cr10c01A', 'Cr10c01B', 'Cr10c02A', 'Cr10c02B', 'Cr10c03A', 'Cr10c03B', 'Cr10c05A', 'Cr10c05B', 'Cr10c06A', 'Cr10c07A', 'Cr10c08A', 'Cr10c08B', 'Cr10c09A', 'Cr10c09B', 'Cr10c10A', 'Cr10c10B', 'Cr10c11A', 'Cr10c11B', 'Cr11c00A', 'Cr11c00B', 'Cr11c01A', 'Cr11c01B', 'Cr11c02A', 'Cr11c02B', 'Cr11c03A', 'Cr11c03B',
'Cr11c04A', 'Cr11c04B', 'Cr11c05A', 'Cr11c05B', 'Cr11c06A', 'Cr11c06B', 'Cr11c07A', 'Cr11c07B', 'Cr11c08A', 'Cr11c08B', 'Cr11c09A', 'Cr11c09B', 'Cr11c10A', 'Cr11c10B'

        ]

DETECTORS = [f'{prefix}_{d}' for d in det_suffix]



def write_solat_instrument_file(output_path: str, version: int) -> None:
    """Writes the SO LAT filelists for Commander3 using Mathew's script."""

    filename = f"SO-LAT_instrument_v{version:02}.h5"

    instrument_file = commander_instrument(output_path, filename, version, "w")

    DEFAULT_BEAM = np.zeros((3, TEMP_LMAX**2 + 2*TEMP_LMAX + 1))

    #beams = akari_utils.get_akari_beams(NSIDE, TEMP_LMAX)
    #sidelobes = akari_utils.get_akari_sidelobes(NSIDE, TEMP_LMAX)
    for band, detector in enumerate(BANDS):
        center_frequency = FREQS[BANDS[band]]
        frequencies = np.array([center_frequency])
        weights = np.array([1.])
        instrument_file.add_bandpass(f"SO-LAT_{detector}", frequencies, weights)

        fwhm = FWHMS[BANDS[band]]

        for i in range(len(DETECTORS)):
            band_group_name = DETECTORS[i]
            instrument_file.add_bandpass(
                f"{band_group_name}", frequencies, weights
            )
            _add_fields(
                instrument_file=instrument_file,
                band_label=DETECTORS[i],
                beam=DEFAULT_BEAM,
                sidelobe=DEFAULT_BEAM,
                fwhm=fwhm,
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

    instrument_file.add_field(band_label + "/polEff", data=[1.0])


def main() -> None:

    version = 1
    # version 1 just takes fiducial values.
    write_solat_instrument_file(output_path=TEMP_OUTPUT_PATH, version=version)


if __name__ == "__main__":
    main()
