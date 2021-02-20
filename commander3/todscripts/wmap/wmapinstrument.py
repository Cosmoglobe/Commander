import h5py
import os
import healpy as hp
import numpy as np
import math
import argparse
from astropy.io import fits
import sys
sys.path.insert(0, "/mn/stornext/u3/duncanwa/Commander/commander3/python")
from commander_tools.tod_tools.wmap import wmap
from commander_tools.tod_tools import commander_instrument as inst


def main():

    #parser = argparse.ArgumentParser()

    #parser.add_argument('--out-dir', type=str, action='store', help='output directory', default='/mn/stornext/d16/cmbco/bp/mathew/test')

    #parser.add_argument('--rimo', type=str, action='store', help='path to the RIMO file', default='/mn/stornext/d16/cmbco/bp/data/auxiliary_data/LFI_RIMO_R3.31.fits')

    #parser.add_argument('--sl-dir', type=str, action='store', help='path to the directory containing the sidelobe alms', default='/mn/stornext/d16/cmbco/bp/data/beamalms/sl')

    #parser.add_argument('--beam-dir', type=str, action='store', help='path to the directory containing the beam alms', default='/mn/stornext/d16/cmbco/bp/data/beamalms/totalAlm')

    #args = parser.parse_args()
    #outDir = args.out_dir

    version = 0

    #rimo = fits.open(args.rimo)

    #slDir = args.sl_dir
    #beamDir = args.beam_dir

    outDir = 'data'
    
    inst_file = inst.commander_instrument(outDir, wmap.instrument_filename(version), version, 'w')

    for freq in wmap.freqs:

        for horn in wmap.horns[freq]:
            for hornType in ['1', '2']:
                nu, B3, B4 = np.loadtxt(f'data/wmap_bandpass_{freq}{horn}{hornType}_v5.cbp').T
                for B, pairType in zip([B3, B4], ['3', '4']):
                    print(freq + str(horn) + hornType + pairType)
                    prefix = freq + str(horn) + hornType + pairType
                    inst_file.add_bandpass(prefix, nu, B)

                    #central frequency
                    inst_file.add_field(prefix + '/centFreq',
                            data=[wmap.cent_freqs[freq + str(horn) + hornType + pairType]])
                prefix = freq + str(horn) + hornType

                #beamData_E = hp.read_alm(os.path.join(beamDir, 'totalAlm_0' + str(freq) + '_' + str(horn) + '_' + beamType + '_qucs-raa.alm'), hdu=2)

                #beamData_B = hp.read_alm(os.path.join(beamDir, 'totalAlm_0' + str(freq) + '_' + str(horn) + '_' + beamType + '_qucs-raa.alm'), hdu=3)

                #slData, mmax_s = hp.read_alm(os.path.join(slDir, 'sl_0' + str(freq) + '_' + str(horn) + '_' + beamType + '_qucs-raa.alm'), return_mmax=True)
           
                #slData_E = hp.read_alm(os.path.join(slDir, 'sl_0' + str(freq) + '_' + str(horn) + '_' + beamType + '_qucs-raa.alm'), hdu=2)
   
                #slData_B = hp.read_alm(os.path.join(slDir, 'sl_0' + str(freq) + '_' + str(horn) + '_' + beamType + '_qucs-raa.alm'), hdu=3) 

                #inst_file.add_alms(prefix, 'beam', wmap.getLmax(len(beamData),
                #    mmax_b), mmax_b, wmap.complex2realAlms(beamData, mmax_b),
                #    wmap.complex2realAlms(beamData_E, mmax_b), wmap.complex2realAlms(beamData_B, mmax_b))

                #inst_file.add_alms(prefix, 'sl', wmap.getLmax(len(slData),
                #    mmax_s), mmax_s, wmap.complex2realAlms(slData, mmax_s),
                #    wmap.complex2realAlms(slData_E, mmax_s), wmap.complex2realAlms(slData_B, mmax_s))
    
                #beam parameters
                inst_file.add_field(prefix + '/fwhm', data=[wmap.fwhms[freq + str(horn)]])
                inst_file.add_field(prefix + '/elip', data=[1])
                inst_file.add_field(prefix + '/psi_ell', data=[0])
                inst_file.add_field(prefix + '/mbeam_eff', data=[1])



    inst_file.finalize()
    wmap.verify_instrument_file(outDir, version)

if __name__ == '__main__':
    main()
