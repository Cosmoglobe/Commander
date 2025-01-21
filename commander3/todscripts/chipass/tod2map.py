import numpy as np 
import matplotlib.pyplot as plt
import h5py
import healpy as hp

filename = 'chipass_full.h5'

f1 = h5py.File(filename, 'r')

j = 0
nside = 256                                                                                                                                
npix = hp.nside2npix(nside)                                                                                                                
                                                                                                                                           
# print(indices)                                                                                                                           
# Initate the map and fill it with the values                                                                                              
nhit = np.zeros(npix, dtype=np.float)                                                                                                      
hpxmap = np.zeros(npix, dtype=np.float)  

for key, data in f1.items():
    # print(key)

    point = data['point_gal'][:, :, :2] * np.pi / 180
    # plt.scatter(point[:,:,1].flatten(), point[:,:,0].flatten())
    # plt.show()
    tod = data['tod'][()]
    # print(np.mean((point[:,:,1].flatten() + np.pi / 2.0) * 180 / np.pi), np.mean(point[:,:,0].flatten() * 180 / np.pi))
    # phi in healpix is 90 - phi from the file
    tod = np.sqrt(tod[:, :, 0] ** 2 + tod[:, :, 1] ** 2).mean(2).flatten()
    indices = hp.ang2pix(nside, np.pi / 2.0 - point[:,:,1].flatten(), point[:,:,0].flatten())
    for i in range(len(indices)):
        nhit[indices[i]] += 1
        hpxmap[indices[i]] += tod[i] #

    j += 1    
    # if j > 10000:
    #     break
    if j % 100 == 0:
        print(j)

hpxmap = hpxmap / nhit  
maxval = 20
np.save('map_%i' % nside, hpxmap) 
np.save('nhit_%i' % nside, nhit)                                                                                                                        
hp.mollview(np.log10(hpxmap)) #, min=-maxval, max=maxval) 
plt.show()