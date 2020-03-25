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
print(DAs.shape)


d1 = 0.5*(DAs[0] + DAs[1])
d2 = 0.5*(DAs[2] + DAs[3])

d = 0.5*(d1 + d2) # = i_A - i_B
p = 0.5*(d1 - d2) # = q_A*cos(2*g_A) + u_A*sin(2*g_A) - q_B*cos(2*g_B) - u_B*sin(2*g_B)

print(d.shape)

gal  = np.array(f[obsid + '/common/gal'])
time = np.array(f[obsid + '/common/time'])
vel  = np.array(f[obsid + '/common/vel'])
sin2g  = np.array(f[obsid + '/common/sin_2_g'])
cos2g  = np.array(f[obsid + '/common/cos_2_g'])




hp.mollview(np.zeros(12)+hp.UNSEEN, coord='G', min=time.min(), max=time.max(), unit='Time')
for i in range(len(gal)):
    hp.projscatter(gal[i,0], gal[i,1], color=plt.cm.viridis(i/len(time)), lonlat=True, s=1)

plt.savefig('time_ordered_pointing.png', bbox_inches='tight')

hp.mollview(np.zeros(12)+hp.UNSEEN, coord='G', min=-1, max=1, unit=r'$\sin2\gamma$')
for i in range(len(gal)):
    hp.projscatter(gal[i,0], gal[i,1], color=plt.cm.viridis((sin2g[i][0][0]+1)/2), lonlat=True, s=1)

plt.savefig('sin2gamma_stream.png', bbox_inches='tight')

hp.mollview(np.zeros(12)+hp.UNSEEN, coord='G', min=-1, max=1, unit=r'$\cos2\gamma$')
for i in range(len(gal)):
    hp.projscatter(gal[i,0], gal[i,1], color=plt.cm.viridis((cos2g[i][0][0]+1)/2), lonlat=True, s=1)

plt.savefig('cos2gamma_stream.png', bbox_inches='tight')
