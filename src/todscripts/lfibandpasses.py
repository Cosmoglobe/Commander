import h5py
import os
import healpy as hp
import numpy as np
from astropy.io import fits

def main():

    outDir = '/mn/stornext/d16/cmbco/bp/mathew/test'

    rimo = fits.open('/mn/stornext/d16/cmbco/bp/data/auxiliary_data/LFI_RIMO_R3.31.fits')

    slDir = '/mn/stornext/d16/cmbco/bp/data/beamalms/sl'
    beamDir = '/mn/stornext/d16/cmbco/bp/data/beamalms/totalAlm'

    horns = {30:[27, 28], 44:[24, 25, 26], 70:[18, 19, 20, 21, 22, 23]}

    outFile = h5py.File(os.path.join(outDir, 'LFI_instrument.h5'), 'w')

    for freq in [30, 44, 70]:
        bandNo = rimo.index_of('BANDPASS_0' + str(freq))
        outFile.create_dataset('0' + str(freq) + '/bandpassx', data=rimo[bandNo].data.field('wavenumber'))
        outFile.create_dataset('0' + str(freq) + '/bandpass', data=rimo[bandNo].data.field('transmission'))

        for horn in horns[freq]:
            for hornType in ['S', 'M']:
                prefix = str(horn) + hornType
                bandNo = rimo.index_of('BANDPASS_0' + str(freq) + '-' + str(horn) + hornType)
                outFile.create_dataset(prefix + '/bandpassx', data=rimo[bandNo].data.field('wavenumber'))
                outFile.create_dataset(prefix + '/bandpass', data=rimo[bandNo].data.field('transmission'))
                beamType = 'y'
                if hornType is 'S':
                    beamType = 'x'
                beamData, mmax_b = hp.read_alm(os.path.join(beamDir, 'totalAlm_0' + str(freq) + '_' + str(horn) + '_' + beamType + '_qucs-raa.alm'), return_mmax=True)

                beamData_E = hp.read_alm(os.path.join(beamDir, 'totalAlm_0' + str(freq) + '_' + str(horn) + '_' + beamType + '_qucs-raa.alm'), hdu=2)

                beamData_B = hp.read_alm(os.path.join(beamDir, 'totalAlm_0' + str(freq) + '_' + str(horn) + '_' + beamType + '_qucs-raa.alm'), hdu=3)

                slData, mmax_s = hp.read_alm(os.path.join(slDir, 'sl_0' + str(freq) + '_' + str(horn) + '_' + beamType + '_qucs-raa.alm'), return_mmax=True)
           
                slData_E = hp.read_alm(os.path.join(slDir, 'sl_0' + str(freq) + '_' + str(horn) + '_' + beamType + '_qucs-raa.alm'), hdu=2)
   
                slData_B = hp.read_alm(os.path.join(slDir, 'sl_0' + str(freq) + '_' + str(horn) + '_' + beamType + '_qucs-raa.alm'), hdu=3) 

                outFile.create_dataset(prefix + '/beam/T', data=complex2realAlms(beamData, mmax_b))
                outFile.create_dataset(prefix + '/beam/E', data=complex2realAlms(beamData_E, mmax_b))
                outFile.create_dataset(prefix + '/beam/B', data=complex2realAlms(beamData_B, mmax_b))


                outFile.create_dataset(prefix + '/sl/T', data=complex2realAlms(slData, mmax_s))
                outFile.create_dataset(prefix + '/sl/E', data=complex2realAlms(slData_E, mmax_s))
                outFile.create_dataset(prefix + '/sl/B', data=complex2realAlms(slData_B, mmax_s))

                outFile.create_dataset(prefix + '/beammmax', data=[mmax_b])
                outFile.create_dataset(prefix + '/beamlmax', data=[getLmax(len(beamData), mmax_b)])
                outFile.create_dataset(prefix + '/slmmax', data=[mmax_s])
                outFile.create_dataset(prefix + '/sllmax', data=[getLmax(len(slData), mmax_s)]) 

                print(prefix)
 
def complex2realAlms(data, mmax):
    lmax = getLmax(len(data), mmax)
    outData = np.zeros((lmax+1)**2)
    
    for l in range(0, lmax):
        for m in range(0, mmax):
            if(m > l):
                continue
            #TODO: figure this out
            scaling = np.sqrt(2)
            if(m == 0):
                scaling = 1
            healpixI = hp.sphtfunc.Alm.getidx(lmax, l, m) 
            outI = getOutidx(l, m)
            outJ = getOutidx(l, -1*m)
            outData[outI] = np.real(data[healpixI]) * scaling
            if(m is not 0):
                outData[outJ] = np.imag(data[healpixI]) * scaling

    return outData

def getLmax(N, mmax):

    return int((2.0*(N - 1.0)/mmax + mmax -1)/(2.0*(1.0 + 1.0/mmax)))

def getOutidx(l, m):
    return l**2 + l + m

if __name__ == '__main__':
    main()
