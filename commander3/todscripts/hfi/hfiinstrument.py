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
from astropy.table import Table

## Fill in path to the Commander python below!! ###
sys.path.append('../Commander/commander3/python')
####################
sys.path.append('/mn/stornext/u3/raelynsu/code/Commander/commander3/python')
#sys.path.insert(0, "/mn/stornext/u3/hke/git/Commander_hfi/commander3/python")
from commander_tools.tod_tools.hfi import hfi
from commander_tools.tod_tools.lfi import lfi
from commander_tools.tod_tools import commander_instrument as inst


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--out-dir', type=str, action='store', help='output directory', default='/mn/stornext/d23/cmbco/globe/planck/rimo')

    parser.add_argument('--rimo', type=str, action='store', help='path to the RIMO file', default='/mn/stornext/d23/cmbco/globe/orig/planck/aux/RIMO_npipe2.fits')

    parser.add_argument('--beam-dir', type=str, action='store', help='path to the directory containing the beam alms', default='/mn/stornext/d23/cmbco/globe/orig/planck/aux/beams')

    parser.add_argument('--fsl-dir', type=str, action='store', help='path to the directory containing the sidelobe alms', default='/mn/stornext/d23/cmbco/globe/planck/beam/fsl')

    tbol_file="/mn/stornext/d23/cmbco/globe/orig/planck/aux/Tbol_SI.txt"

    args = parser.parse_args()
    outDir = args.out_dir

    version = 6

    rimo = fits.open(args.rimo)
    
    inst_file = inst.commander_instrument(outDir, hfi.instrument_filename(version), version, 'w')

    ##from https://www.aanda.org/articles/aa/full_html/2016/10/aa25844-15/T9.html
    tbol = Table.read(tbol_file,format='ascii')

    for freq in hfi.freqs:
        bandNo = rimo.index_of('BANDPASS_F' + str(freq))
        inst_file.add_hfi_bandpass(freq, rimo[bandNo].data.field('WAVENUMBER'), rimo[bandNo].data.field('TRANSMISSION'))

        for det in hfi.dets[freq]:
            prefix = str(freq) + '-' + det
            bandNo = rimo.index_of('bandpass_' + str(freq) + '-' + det)
            inst_file.add_hfi_bandpass(prefix, rimo[bandNo].data.field('wavenumber'), rimo[bandNo].data.field('transmission'))
            
            # beamData, mmax_b = hp.read_alm(os.path.join(args.beam_dir, 'blm_' + str(freq) + '-' + det + '.fits'), return_mmax=True)

            #beamData_E = None

            #beamData_B = None

            #These should be in beam_fsl_Pxx_100-1a.fits and so on
            # if(freq < 545):
            slData, mmax_s = hp.read_alm(os.path.join(args.fsl_dir, 'fsl_alms_' + str(freq) + '-' + det + '_v2.fits'), return_mmax=True)
            mbData, mmax_mb = hp.read_alm(os.path.join(args.fsl_dir, 'mb_alms_' + str(freq) + '-' + det + '_v2.fits'), return_mmax=True)
            fpibData, mmax_fpib = hp.read_alm(os.path.join(args.fsl_dir, 'fpib_alms_' + str(freq) + '-' + det + '_v2.fits'), return_mmax=True)

                # slData_E, mmax_s = hp.read_alm(os.path.join(args.beam_dir, 'fsl_alms_' + str(freq) + '-' + det + '.fits'), return_mmax=True, hdu=1)

                # slData_B, mmax_s = hp.read_alm(os.path.join(args.beam_dir, 'fsl_alms_' + str(freq) + '-' + det + '.fits'), return_mmax=True, hdu=2)

            # inst_file.add_alms(prefix, 'sl', lfi.getLmax(len(slData), mmax_s), mmax_s, lfi.complex2realAlms(slData, mmax_s), lfi.complex2realAlms(slData_E, mmax_s), lfi.complex2realAlms(slData_B, mmax_s))
            out=lfi.complex2realAlms(mbData, mmax_mb)
            # inst_file.add_alms(prefix, 'sl', lfi.getLmax(len(slData), mmax_s), mmax_s, lfi.complex2realAlms(slData, mmax_s), lfi.complex2realAlms(slData_E, mmax_s), lfi.complex2realAlms(slData_B, mmax_s))
            inst_file.add_alms(prefix, 'mb', lfi.getLmax(len(mbData), mmax_mb), mmax_mb, lfi.complex2realAlms(mbData, mmax_mb), None, None)#, [0], [0])
            inst_file.add_alms(prefix, 'beam', lfi.getLmax(len(fpibData), mmax_fpib), mmax_fpib, lfi.complex2realAlms(fpibData, mmax_fpib), None, None)#, [0], [0])
            inst_file.add_alms(prefix, 'sl', lfi.getLmax(len(slData), mmax_s), mmax_s, lfi.complex2realAlms(slData, mmax_s), None, None)#, [0], [0])
            
            
            # else:
                #we need sidelobe models for 545 and 857
                # inst_file.add_alms(prefix, 'sl', 0, 0, [0], [0], [0])


            # inst_file.add_alms(prefix, 'beam', lfi.getLmax(len(beamData), mmax_b), mmax_b, lfi.complex2realAlms(beamData, mmax_b), None, None)


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

            inst_file.add_field(prefix + '/polEff', data=[hfi.pol_effs[str(freq) + '-' + det]])

            inst_file.add_field(prefix+ '/Tbol', data=np.asarray(tbol[tbol['Bolometer'] == str(freq) + '-' + det][0][1:]))

            print(prefix)

    inst_file.finalize()
    hfi.verify_instrument_file(outDir, version)

if __name__ == '__main__':
    main()
