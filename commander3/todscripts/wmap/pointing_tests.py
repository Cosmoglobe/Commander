import numpy as np
import matplotlib.pyplot as plt

import healpy as hp
from astropy.io import fits
from glob import glob

from get_gain_model import get_gain


file_num = 30

prefix = '/mn/stornext/d16/cmbco/bp/wmap/'
files = glob(prefix + 'tod/new/*.fits')
files.sort()
data = fits.open(files[file_num])

version=13
nside = 256

#version = 14
#nside = 512

allbands = ['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']
band = 'K1'

time = data[2].data['time']

import h5py
from glob import glob
import huffman
from cg_solver import make_dipole
fnames = glob(f'/mn/stornext/d16/cmbco/bp/wmap/data/wmap_{band}_*v{version}.h5')
fnames.sort()
fname = fnames[file_num]
print(fname)
f= h5py.File(fname, 'r')
obsid = str(list(f.keys())[0])
labels = [f'{band}13', f'{band}14',f'{band}23',f'{band}24']

huffTree = f[obsid+'/common/hufftree']
huffSymb = f[obsid+'/common/huffsymb']
h = huffman.Huffman(tree=huffTree, symb=huffSymb)




DAs = [[], [], [], []]
sigmas = []
npsi = 2048
psiBins = np.linspace(0, 2*np.pi, npsi)
gains = np.zeros(len(labels))
for num, label in enumerate(labels):
    TODs = np.array(f[obsid + '/' + label + '/tod'])
    scalars = f[obsid + '/' + label + '/scalars']
    gains[num] = scalars[0]
    TODs = TODs - np.median(TODs)
    DAs[num] = DAs[num] + TODs.tolist()
    sigmas.append(TODs.std())
    if label == f'{band}13':
        pixA = h.Decoder(np.array(f[obsid + '/' + label + \
            '/pixA'])).astype('int')
        pixB = h.Decoder(np.array(f[obsid + '/' + label + \
            '/pixB'])).astype('int')
        psiA = psiBins[h.Decoder(np.array(f[obsid + '/' + label + \
            '/psiA'])).astype('int')]
        psiB = psiBins[h.Decoder(np.array(f[obsid + '/' + label + \
            '/psiB'])).astype('int')]


hp.mollview(np.zeros(12)+hp.UNSEEN)
theta, phi = hp.pix2ang(nside, pixA[:100000])
hp.projplot(theta, phi, color='C1')
theta, phi = hp.pix2ang(nside, pixB[:100000])
hp.projplot(theta, phi, color='C2')
plt.savefig(f'scan_v{version}.png')

plt.figure()
plt.plot(psiA[:10000])
plt.plot(psiB[:10000])
plt.savefig(f'psi_v{version}.png')


plt.show()
