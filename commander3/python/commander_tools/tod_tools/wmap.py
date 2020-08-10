import h5py
import commander_tools.tod_tools.commander_instrument as inst
import os
import numpy as np
import healpy as hp
import random
import math

class wmap(object):   
 
    freqs = []
    horns = {}
    hornTypes = []
    npsi = 4096
    ntodsigma = 100
    nsides = {}
    #compression arrays 
    huffman = ['huffman', {'dictNum':1}]
    huffTod = ['huffman', {'dictNum':2}]
    psiDigitize = ['digitize', {'min':0, 'max':2*np.pi,'nbins':npsi}]
    todDytpe = ['dtype', {'dtype':'f4'}]
    todSigma = ['sigma', {'sigma0':None, 'nsigma':ntodsigma}] 
    #fwhm, elipticity and psi_ell from https://www.aanda.org/articles/aa/full_html/2016/10/aa25809-15/T6.html
    fwhms = {}

    elips = {}

    psis = {}

    b_effs = {}

    cent_freqs = {}

    #psi_uv from https://www.aanda.org/articles/aa/full_html/2016/10/aa25818-15/T5.html
    mbangs = {}

    
    def __init__():
        return

    @staticmethod
    def instrument_filename(version):
        return 'WMAP_instrument_v' + str(version) + '.h5'

    @staticmethod
    def ring_outer_product(theta, phi):
        outAng = [0, 0]
        nsamps = min(100, len(theta))
        pair1 = random.sample(range(len(theta)), nsamps)
        pair2 = random.sample(range(len(theta)), nsamps)
        vecs = hp.pixelfunc.ang2vec(theta, phi)
        for a1, a2 in zip(pair1, pair2):
            crossP = np.cross(vecs[a1], vecs[a2])
            if(crossP[0] < 0):
                crossP *= -1
            theta1, phi1 = hp.vec2ang(crossP)
            if(not math.isnan(theta1) and not math.isnan(phi1)):
                outAng[0] += theta1/nsamps
                outAng[1] += phi1/nsamps
            else:
                outAng[0] += outAng[0]/nsamps
                outAng[1] += outAng[1]/nsamps
        return outAng

    def horn2freq(self, horn):
         return

    def freq2horns(self, freq):
        return self.horns[freq]

    @staticmethod
    def verify_instrument_file(outDir, version):
        f = inst.commander_instrument(outDir, lfi.instrument_filename(version), 'r')

        if version <= 3:
            raise ValueError("Versions of LFI instrument files <= 3 are lost to time")
        if version == 4:
            print("Should check the version here")

        if version > 4:
            raise ValueError("Version " + str(version) + " of LFI instrument file has not yet been defined.")

    @staticmethod
    def mbang(horn):
        mbangs = {}
        if len(str(horn)) == 3:
            horn = int(str(horn[:-1]))
        return mbangs[int(horn)]

    @staticmethod
    def complex2realAlms(data, mmax):
        lmax = lfi.getLmax(len(data), mmax)
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
                outI = lfi.getOutidx(l, m)
                outJ = lfi.getOutidx(l, -1*m)
                outData[outI] = np.real(data[healpixI]) * scaling
                if(m is not 0):
                    outData[outJ] = np.imag(data[healpixI]) * scaling

        return outData

    @staticmethod
    def getLmax(N, mmax):
        return int((2.0*(N - 1.0)/mmax + mmax -1)/(2.0*(1.0 + 1.0/mmax)))

    @staticmethod
    def getOutidx(l, m):
        return l**2 + l + m

