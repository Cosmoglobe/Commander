import h5py
import os
from astropy.io import fits

def main():

    outDir = '/mn/stornext/d14/bp/mathew/test'

    rimo = fits.open('/mn/stornext/d14/bp/data/auxiliary_data/LFI_RIMO_R3.31.fits')

    horns = {30:[27, 28], 44:[24, 25, 26], 70:[18, 19, 20, 21, 22, 23]}

    outFile = h5py.File(os.path.join(outDir, 'LFI_bandpasses.h5'), 'w')

    for freq in [30, 44, 70]:
        for horn in horns[freq]:
            for hornType in ['S', 'M']:
                prefix = str(horn) + hornType
                bandNo = rimo.index_of('BANDPASS_0' + str(freq) + '-' + str(horn) + hornType)
                outFile.create_dataset(prefix + '/bandpassx', data=rimo[bandNo].data.field('wavenumber'))
                outFile.create_dataset(prefix + '/bandpass', data=rimo[bandNo].data.field('transmission'))



if __name__ == '__main__':
    main()
