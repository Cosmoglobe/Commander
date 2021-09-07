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
from astropy.io import fits
from commander_tools.tod_tools.lfi import lfi
from commander_tools.tod_tools import commander_instrument as inst


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--out-dir', type=str, action='store', help='output directory', default='/mn/stornext/d16/cmbco/bp/mathew/test')

    parser.add_argument('--rimo', type=str, action='store', help='path to the RIMO file', default='/mn/stornext/d16/cmbco/bp/data/auxiliary_data/LFI_RIMO_R3.31.fits')

    parser.add_argument('--sl-dir', type=str, action='store', help='path to the directory containing the sidelobe alms', default='/mn/stornext/d16/cmbco/bp/data/beamalms/sl')

    parser.add_argument('--beam-dir', type=str, action='store', help='path to the directory containing the beam alms', default='/mn/stornext/d16/cmbco/bp/data/beamalms/mainbeams')

    args = parser.parse_args()
    outDir = args.out_dir

    version = 5

    rimo = fits.open(args.rimo)

    slDir = args.sl_dir
    beamDir = args.beam_dir
    
    inst_file = inst.commander_instrument(outDir, lfi.instrument_filename(version), version, 'w')

    for freq in lfi.freqs:
        bandNo = rimo.index_of('BANDPASS_0' + str(freq))
        inst_file.add_bandpass(freq, rimo[bandNo].data.field('wavenumber'), rimo[bandNo].data.field('transmission'))

        for horn in lfi.horns[freq]:
            for hornType in ['S', 'M']:
                prefix = str(horn) + hornType
                bandNo = rimo.index_of('BANDPASS_0' + str(freq) + '-' + str(horn) + hornType)
                inst_file.add_bandpass(prefix, rimo[bandNo].data.field('wavenumber'), rimo[bandNo].data.field('transmission'))
                beamData, mmax_b = hp.read_alm(os.path.join(beamDir, 'mbib_DX12_LFI' + str(horn) + hornType + '.fits'), return_mmax=True)

                beamData_E = hp.read_alm(os.path.join(beamDir, 'mbib_DX12_LFI' + str(horn) + hornType + '.fits'), hdu=2)

                beamData_B = hp.read_alm(os.path.join(beamDir, 'mbib_DX12_LFI' + str(horn) + hornType + '.fits'), hdu=3)

                beamType = 'y'
                if hornType is 'S':
                    beamType = 'x'

                slData, mmax_s = hp.read_alm(os.path.join(slDir, 'sl_0' + str(freq) + '_' + str(horn) + '_' + beamType + '_qucs-raa.alm'), return_mmax=True)
           
                slData_E = hp.read_alm(os.path.join(slDir, 'sl_0' + str(freq) + '_' + str(horn) + '_' + beamType + '_qucs-raa.alm'), hdu=2)
   
                slData_B = hp.read_alm(os.path.join(slDir, 'sl_0' + str(freq) + '_' + str(horn) + '_' + beamType + '_qucs-raa.alm'), hdu=3) 

                inst_file.add_alms(prefix, 'beam', lfi.getLmax(len(beamData), mmax_b), mmax_b, lfi.complex2realAlms(beamData, mmax_b), lfi.complex2realAlms(beamData_E, mmax_b), lfi.complex2realAlms(beamData_B, mmax_b))

                inst_file.add_alms(prefix, 'sl', lfi.getLmax(len(slData), mmax_s), mmax_s, lfi.complex2realAlms(slData, mmax_s), lfi.complex2realAlms(slData_E, mmax_s), lfi.complex2realAlms(slData_B, mmax_s))

                #beam parameters
                inst_file.add_field(prefix + '/fwhm', data=[lfi.fwhms[str(horn) + hornType]])
                inst_file.add_field(prefix + '/elip', data=[lfi.elips[str(horn) + hornType]])
                inst_file.add_field(prefix + '/psi_ell', data=[math.radians(lfi.psis[str(horn) + hornType])])
                inst_file.add_field(prefix + '/mbeam_eff', data=[lfi.b_effs[str(horn) + hornType]])

                #central frequency
                inst_file.add_field(prefix + '/centFreq', data=[lfi.cent_freqs[str(horn) + hornType]])

                print(prefix)

    inst_file.finalize()
    lfi.verify_instrument_file(outDir, version)

if __name__ == '__main__':
    main()
