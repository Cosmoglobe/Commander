import numpy as np
import matplotlib.pyplot as plt

import h5py
import healpy as hp
f= h5py.File('h5_wmap_test.h5', 'r')
obsid = str(list(f.keys())[0])

DAs = []
labels = ['K113', 'K114', 'K123', 'K124']
for label in labels:
    TODs = np.array(f[obsid + '/' + label + '/TOD'])
    DAs.append(TODs)
DAs = np.array(DAs)


d1 = 0.5*(DAs[0] + DAs[1])
d2 = 0.5*(DAs[2] + DAs[3])

d = 0.5*(d1 + d2) # = i_A - i_B
p = 0.5*(d1 - d2) # = q_A*cos(2*g_A) + u_A*sin(2*g_A) - q_B*cos(2*g_B) - u_B*sin(2*g_B)

gal  = np.array(f[obsid + '/common/gal'])
time = np.array(f[obsid + '/common/time'])
vel  = np.array(f[obsid + '/common/vel'])
sin2g  = np.array(f[obsid + '/common/sin_2_g'])
cos2g  = np.array(f[obsid + '/common/cos_2_g'])

print(sin2g.shape)
print(gal.shape)

band_labels = [
'K1A',
'K1B',
'KA1A',
'KA1B',
'Q1A',
'Q1B',
'Q2A',
'Q2B',
'V1A',
'V1B',
'V2A',
'V2B',
'W1A',
'W1B',
'W2A',
'W2B',
'W3A',
'W3B',
'W4A',
'W4B']


for band in range(len(gal)):
    print(band_labels[band])
    hp.mollview(np.zeros(12)+hp.UNSEEN, coord='G', min=time.min(), max=time.max(), unit='Time', title=band_labels[band])
    for i in range(len(gal[band])):
        hp.projscatter(gal[band,i,0], gal[band,i,1], color=plt.cm.viridis(i/len(time)), lonlat=True, s=1)
    
    plt.savefig(f'time_ordered_pointing_{band_labels[band]}.png', bbox_inches='tight')
    
    hp.mollview(np.zeros(12)+hp.UNSEEN, coord='G', min=-1, max=1, unit=r'$\sin2\gamma$', title=band_labels[band])
    for i in range(len(gal[band])):
        hp.projscatter(gal[band,i,0], gal[band,i,1], color=plt.cm.viridis((sin2g[band][i]+1)/2), lonlat=True, s=1)
    
    plt.savefig(f'sin2gamma_stream_{band_labels[band]}.png', bbox_inches='tight')
    
    hp.mollview(np.zeros(12)+hp.UNSEEN, coord='G', min=-1, max=1, unit=r'$\cos2\gamma$', title=band_labels[band])
    for i in range(len(gal[band])):
        hp.projscatter(gal[band,i,0], gal[band,i,1], color=plt.cm.viridis((cos2g[band][i]+1)/2), lonlat=True, s=1)
    
    plt.savefig(f'cos2gamma_stream_{band_labels[band]}.png', bbox_inches='tight')
    
    plt.figure()
    plt.title(band_labels[band])
    plt.scatter(cos2g[band], sin2g[band], c=time, s=1)
    plt.xlabel(r'$\cos2\gamma$')
    plt.ylabel(r'$\sin2\gamma$')
    plt.colorbar(label='Time')
    plt.savefig(f'cos2_v_sin2_{band_labels[band]}.png', bbox_inches='tight')
    plt.close('all')
