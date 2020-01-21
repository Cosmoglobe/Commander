import h5py
import os
import healpy as hp
import numpy as np
import math
from astropy.io import fits

def main():

    outDir = '/mn/stornext/d16/cmbco/bp/mathew/test'

    rimo = fits.open('/mn/stornext/d16/cmbco/bp/data/auxiliary_data/LFI_RIMO_R3.31.fits')

    slDir = '/mn/stornext/d16/cmbco/bp/data/beamalms/sl'
    beamDir = '/mn/stornext/d16/cmbco/bp/data/beamalms/totalAlm'

    horns = {30:[27, 28], 44:[24, 25, 26], 70:[18, 19, 20, 21, 22, 23]}

    #fwhm, elipticity and psi_ell from https://www.aanda.org/articles/aa/full_html/2016/10/aa25809-15/T6.html
    fwhms = {'18M':13.44, '18S':13.5, '19M':13.14, '19S':13.07, '20M':12.84, '20S':12.84, '21M':12.77, '21S':12.87, '22M':12.92, '22S':12.97, '23M':13.35, '23S':13.36, '24M':23.18, '24S':23.04, '25M':30.23, '25S':30.94, '26M':30.29, '26S':30.64, '27M':32.02, '27S':33.11, '28M':33.1, '28S':33.09}

    elips = {'18M':1.23, '18S':1.27, '19M':1.25, '19S':1.28, '20M':1.27, '20S':1.29, '21M':1.28, '21S':1.29, '22M':1.27, '22S':1.28, '23M':1.23, '23S':1.28, '24M':1.39, '24S':1.34, '25M':1.19, '25S':1.19, '26M':1.19, '26S':1.19, '27M':1.37, '27S':1.38, '28M':1.37, '28S':1.37}

    psis = {'18M':85, '18S':86, '19M':78, '19S':79, '20M':71, '20S':72, '21M':107, '21S':106, '22M':101, '22S':101, '23M':92, '23S':92, '24M':89, '24S':89, '25M':114, '25S':117, '26M':62, '26S':61, '27M':101, '27S':101, '28M':78, '28S':78}

    b_effs = {'18M':0.9921, '18S':0.9887, '19M':0.9883, '19S':0.9898, '20M':0.9885, '20S':0.9881, '21M':0.9894, '21S':0.9882, '22M':0.9916, '22S':0.9915, '23M':0.9926, '23S':0.9919, '24M':0.9972, '24S':0.9973, '25M':0.9975, '25S':0.9976, '26M':0.9974, '26S':0.9977, '27M':0.9904, '27S':0.9889, '28M':0.9907, '28S':0.9879}

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
    
                #beam parameters
                outFile.create_dataset(prefix + '/fwhm', data=[fwhms[str(horn) + hornType]])
                outFile.create_dataset(prefix + '/elip', data=[elips[str(horn) + hornType]])
                outFile.create_dataset(prefix + '/psi_ell', data=[math.radians(psis[str(horn) + hornType])])
                outFile.create_dataset(prefix + '/mbeam_eff', data=[b_effs[str(horn) + hornType]])

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
