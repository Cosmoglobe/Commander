import h5py
import commander_tools.tod_tools.commander_instrument as inst
import os
import numpy as np
import healpy as hp
import random
import math

class wmap(object):   
 
    freqs = ['K', 'Ka', 'Q', 'V', 'W']
    nhorns = [1, 1, 2, 2, 4]
    # These are really DAs...
    horns = {f:np.arange(nhorn)+1 for nhorn, f in zip(nhorns, freqs)}
    hornTypes = [1, 2] # each DA has two pairs
    pairType = [3,4] # Horns A/B, or 3/4.
    npsi = 4096
    ntodsigma = 100
    nsides = {f:512 for f in freqs}
    #compression arrays 
    huffman = ['huffman', {'dictNum':1}]
    huffTod = ['huffman', {'dictNum':2}]
    psiDigitize = ['digitize', {'min':0, 'max':2*np.pi,'nbins':npsi}]
    todDytpe = ['dtype', {'dtype':'f4'}]
    todSigma = ['sigma', {'sigma0':None, 'nsigma':ntodsigma}] 
    # fwhm for WMAP is surprisingly difficult to find, I am just using the
    # values in table 1 of Planck Collaboration X 2015
    fwhms = {'K1': 51.5, 'Ka1': 44.4, 'Q1': 33.9, 'Q2':34.4, 
            'V1':23.3, 'V2':23.3, 
            'W1': 15.0, 'W2': 13.8, 'W3': 14.4, 'W4': 15.0}

    # Again, hard to find a specific number for this.
    elips = {}

    # I think this doesn't matter for WMAP
    psis = {}

    # Not sure what this is
    b_effs = {}

    cent_freqs = {'K113':22.0, 'K114': 22.6, 'K123':23.2, 'K124': 23.1,
                  'Ka113': 32.7, 'Ka114': 32.9, 'Ka123':33.1, 'Ka124':33.1,
                  'Q113': 40.7, 'Q114': 40.8, 'Q123': 40.9, 'Q124': 40.8,
                  'Q213': 40.2, 'Q214': 40.2, 'Q223': 40.9, 'Q224': 41.0, 
                  'V113': 59.3, 'V114': 59.1, 'V123': 61.1, 'V124': 61.1,
                  'V213': 61.6, 'V214': 61.6, 'V223': 60.7, 'V224': 60.6,
                  'W113': 93.9, 'W114': 93.3, 'W123': 93.1, 'W124': 93.3,
                  'W213': 93.6, 'W214': 93.4, 'W223': 94.0, 'W224': 94.6,
                  'W313': 92.4, 'W314': 92.3, 'W323': 92.9, 'W324': 93.7,
                  'W413': 94.2, 'W414': 94.4, 'W423': 93.2, 'W424': 93.0,
                  }

    # I also think this doesn't matter for WMAP
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
        f = inst.commander_instrument(outDir, wmap.instrument_filename(version), 'r')

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

