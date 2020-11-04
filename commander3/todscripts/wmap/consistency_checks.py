#================================================================================
#
# Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
#
# This file is part of Commander3.
#
# Commander3 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Commander3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Commander3. If not, see <https://www.gnu.org/licenses/>.
#
#================================================================================

import numpy as np
import matplotlib.pyplot as plt

import h5py
import healpy as hp
f= h5py.File('wmap_tods.h5', 'r')
obsid = str(list(f.keys())[0])

print(obsid)

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

plot_map = True



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


if plot_map:
    for band in range(len(gal)):
        print(band_labels[band])
        hp.mollview(np.zeros(12)+hp.UNSEEN, coord='G', min=time.min(), max=time.max(), unit='Time', title=band_labels[band])
        for i in range(len(gal[band])):
            hp.projscatter(gal[band,i,0], gal[band,i,1], color=plt.cm.viridis(i/len(time)), lonlat=True, s=1)
        
        plt.savefig(f'plots/time_ordered_pointing_{band_labels[band]}.png', bbox_inches='tight')
        
        hp.mollview(np.zeros(12)+hp.UNSEEN, coord='G', min=-1, max=1, unit=r'$\sin2\gamma$', title=band_labels[band])
        for i in range(len(gal[band])):
            hp.projscatter(gal[band,i,0], gal[band,i,1], color=plt.cm.viridis((sin2g[band][i]+1)/2), lonlat=True, s=1)
        
        plt.savefig(f'plots/sin2gamma_stream_{band_labels[band]}.png', bbox_inches='tight')
        
        hp.mollview(np.zeros(12)+hp.UNSEEN, coord='G', min=-1, max=1, unit=r'$\cos2\gamma$', title=band_labels[band])
        for i in range(len(gal[band])):
            hp.projscatter(gal[band,i,0], gal[band,i,1], color=plt.cm.viridis((cos2g[band][i]+1)/2), lonlat=True, s=1)
        
        plt.savefig(f'plots/cos2gamma_stream_{band_labels[band]}.png', bbox_inches='tight')
        
        plt.figure()
        plt.title(band_labels[band])
        plt.scatter(cos2g[band], sin2g[band], c=time, s=1)
        plt.xlabel(r'$\cos2\gamma$')
        plt.ylabel(r'$\sin2\gamma$')
        plt.colorbar(label='Time')
        plt.savefig(f'plots/cos2_v_sin2_{band_labels[band]}.png', bbox_inches='tight')
        plt.close('all')
