import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import huffman
import h5py


from glob import glob

version = 24
band = 'K1'

fnames = glob(f'/mn/stornext/d16/cmbco/bp/wmap/data/wmap_{band}_*v{version}.h5')
fnames.sort()
fname = fnames[0]
fname = '/mn/stornext/d16/cmbco/bp/wmap/data/wmap_K1_001588_v24.h5'

labels = [f'{band}13', f'{band}14',f'{band}23',f'{band}24']
f= h5py.File(fname, 'r')
obsid = str(list(f.keys())[0])
obsid = '038101'
huffTree = f[obsid+'/common/hufftree']
huffSymb = f[obsid+'/common/huffsymb']
h = huffman.Huffman(tree=huffTree, symb=huffSymb)


DAs = [[], [], [], []]
flags = [[], [], [], []]
sigmas = []
gains = np.zeros(len(labels))
npsi = 2048
psiBins = np.linspace(0, 2*np.pi, npsi)
for num, label in enumerate(labels):
    TODs = np.array(f[obsid + '/' + label + '/tod'])
    scalars = f[obsid + '/' + label + '/scalars']
    gains[num] = scalars[0]
    flag = h.Decoder(np.array(f[obsid + '/' + label + '/flag']))
    flags[num] = flags[num] + flag.tolist()
    DAs[num] = DAs[num] + TODs.tolist()
    sigmas.append(TODs.std())
    if label == f'{band}13':
        pixA = h.Decoder(np.array(f[obsid + '/' + label + \
            '/pixA'])).astype('int')
        pixB = h.Decoder(np.array(f[obsid + '/' + label + \
            '/pixB'])).astype('int')



nside = 512
thetaA, phiA = hp.pix2ang(nside, pixA)
thetaB, phiB = hp.pix2ang(nside, pixB)
# loop over:
for i in range(7):
    hp.mollview(hp.UNSEEN*np.ones(12), cbar=False, title='', sub=(3,3,i+1))
    ax = plt.gca()
    ax.projscatter(thetaA[4000*i:4000*(i+1)], phiA[4000*i:4000*(i+1)],
            color='r', s=0.1)
    ax.projscatter(thetaB[4000*i:4000*(i+1)], phiB[4000*i:4000*(i+1)],
            color='b', s=0.1)
    plt.title(f'{4000*i}--{4000*(i+1)}')
hp.mollview(hp.UNSEEN*np.ones(12), cbar=False, title='')
hp.projscatter(thetaA[:4000], phiA[:4000], color='r', s=0.5)
hp.projscatter(thetaB[:4000], phiB[:4000], color='b', s=0.5)

hp.mollview(hp.UNSEEN*np.ones(12), cbar=False, title='')
hp.projscatter(thetaA[4000:8000], phiA[4000:8000], color='r', s=0.5)
hp.projscatter(thetaB[4000:8000], phiB[4000:8000], color='b', s=0.5)
hp.mollview(hp.UNSEEN*np.ones(12), cbar=False, title='')
hp.projscatter(thetaA[9000:13000], phiA[9000:13000], color='r', s=0.5)
hp.projscatter(thetaB[9000:13000], phiB[9000:13000], color='b', s=0.5)
plt.title('Scans 9000--13000')
plt.show()
