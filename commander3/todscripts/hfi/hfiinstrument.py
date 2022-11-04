#================================================================================
#
# Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
#
# This file is part of Commander3.
#
# Commander3 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Commander3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Commander3. If not, see <https://www.gnu.org/licenses/>.
#
#================================================================================

import h5py
import os
import healpy as hp
import numpy as np
import math
import argparse
import sys
from astropy.io import fits
sys.path.insert(0, "/mn/stornext/u3/hke/git/Commander_hfi/commander3/python")
from commander_tools.tod_tools.hfi import hfi
from commander_tools.tod_tools.lfi import lfi
from commander_tools.tod_tools import commander_instrument as inst


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--out-dir', type=str, action='store', help='output directory', default='/mn/stornext/d16/cmbco/bp/mathew/hfi')

    parser.add_argument('--rimo', type=str, action='store', help='path to the RIMO file', default='/mn/stornext/d16/cmbco/bp/HFI/aux/RIMO_npipe2.fits')

    parser.add_argument('--beam-dir', type=str, action='store', help='path to the directory containing the sidelobe alms', default='/mn/stornext/d16/cmbco/bp/HFI/aux/beams')

    args = parser.parse_args()
    outDir = args.out_dir

    version = 1

    rimo = fits.open(args.rimo)
    
    inst_file = inst.commander_instrument(outDir, hfi.instrument_filename(version), version, 'w')

    for freq in hfi.freqs:
        bandNo = rimo.index_of('BANDPASS_F' + str(freq))
        inst_file.add_hfi_bandpass(freq, rimo[bandNo].data.field('WAVENUMBER'), rimo[bandNo].data.field('TRANSMISSION'))

        for det in hfi.dets[freq]:
            prefix = str(freq) + '-' + det
            bandNo = rimo.index_of('bandpass_' + str(freq) + '-' + det)
            inst_file.add_hfi_bandpass(prefix, rimo[bandNo].data.field('wavenumber'), rimo[bandNo].data.field('transmission'))
            
            beamData, mmax_b = hp.read_alm(os.path.join(args.beam_dir, 'blm_' + str(freq) + '-' + det + '.fits'), return_mmax=True)

            #beamData_E = None

            #beamData_B = None

            #These should be in beam_fsl_Pxx_100-1a.fits and so on
            if(freq < 545):
                slData, mmax_s = hp.read_alm(os.path.join(args.beam_dir, 'fsl_alms_' + str(freq) + '-' + det + '.fits'), return_mmax=True)
           
                slData_E, mmax_s = hp.read_alm(os.path.join(args.beam_dir, 'fsl_alms_' + str(freq) + '-' + det + '.fits'), return_mmax=True, hdu=1)

                slData_B, mmax_s = hp.read_alm(os.path.join(args.beam_dir, 'fsl_alms_' + str(freq) + '-' + det + '.fits'), return_mmax=True, hdu=2)

                inst_file.add_alms(prefix, 'sl', lfi.getLmax(len(slData), mmax_s), mmax_s, lfi.complex2realAlms(slData, mmax_s), lfi.complex2realAlms(slData_E, mmax_s), lfi.complex2realAlms(slData_B, mmax_s))



            inst_file.add_alms(prefix, 'beam', lfi.getLmax(len(beamData), mmax_b), mmax_b, lfi.complex2realAlms(beamData, mmax_b), None, None)


            #beam parameters
            detnames = rimo[1].data['detector']

            fwhm = rimo[1].data['fwhm'][detnames == str(freq) +'-'+det]
            inst_file.add_field(prefix + '/fwhm', data=fwhm)

            elip = rimo[1].data['ellipticity'][detnames == str(freq) + '-'+det]
            inst_file.add_field(prefix + '/elip', data=elip)

            psi_ell = rimo[1].data['posang'][detnames == str(freq) + '-' + det]
            inst_file.add_field(prefix + '/psi_ell', data=psi_ell)

            #central frequency
            inst_file.add_field(prefix + '/centFreq', data=[hfi.cent_freqs[str(freq) + '-' + det]])

            print(prefix)

    inst_file.finalize()
    hfi.verify_instrument_file(outDir, version)

if __name__ == '__main__':
    main()
