import h5py
import os

class commander_instrument:

    def __init__(self, path, fileName, version, mode):
        self.version = version
        self.outPath = path
        self.h5file = h5py.File(os.path.join(path, fileName), mode)
        
    #Raw field writing
    def add_field(self, fieldName, data, compression=None):
        if(compression is None):
            self.h5file.create_dataset(fieldName, data=data)
        else:
            raise ValueError('Compression type ' + compression + ' is not supported')
    
    #add bandpass data for detector det
    def add_bandpass(self, det, freqs, response):
        if len(freqs) is not len(response):
            raise ValueError('Bandpass frequency x and f(x) must be the same length')
        det = self.parse_det(det)

        self.add_field(det + '/bandpassx', freqs)
        self.add_field(det + '/bandpass', response)
       
    def add_alms(self, det, almType, lmax, mmax, T, E, B):
       
        if(mmax > lmax):
            mmax = lmax
 
        nalm = self.nalm(lmax, mmax)
        if (len(T) != nalm) or (len(E) != nalm) or (len(B) != nalm):
            raise ValueError('Length of alms is wrong for lmax=' + str(lmax) + ', mmax=' +str(mmax) + '. Lengths: ' + str(len(T)) + ', ' + str(len(E)) + ', ' + str(len(B)) + ' nalm= ' + str(nalm))
               
        det = self.parse_det(det) 

        self.add_field(det + '/' + almType + '/T', T)
        #polarized experiment
        if(E is not None and B is not None):
            self.add_field(det + '/' + almType + '/E', E)
            self.add_field(det + '/' + almType + '/B', B)

        self.add_field(det + '/'+almType + 'mmax', [mmax])
        self.add_field(det + '/'+almType + 'lmax', [lmax])

    #currently does nothing but maybe it will in the future, who knowssss
    def finalize(self):
        self.add_field('common/version', self.version)

    def parse_det(self, det):
        if type(det) is not str:
            det = str(det).zfill(3)
        return det

    def nalm(self, lmax, mmax):
        return lmax**2 + 2*lmax +1
