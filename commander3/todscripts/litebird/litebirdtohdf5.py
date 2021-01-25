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

import h5py
import numpy as np
import healpy as hp
import huffman
import os

vel_file = "/mn/stornext/u3/ranajoyb/genesys/data/satellite_velocity.fits"
tod_dir = "/mn/stornext/u3/ranajoyb/genesys/output/sim_test_old/scan_test/LFT_40GHz"
det_list = ["0001a", "0001b", "0002a", "0002b", "0003a", "0003b"]
segment = "0000002"
out_dir = "/mn/stornext/u3/ranajoyb/genesys/output/sim_test_old/"
outname = 'LFT_40GHz_' + segment + '.h5'

out_f = h5py.File(os.path.join(out_dir, outname), 'w')

nside = 256
npsi = 4096
fsamp = 31
chunk_size = 3600
nsamp = chunk_size*fsamp
chunk_list = np.arange(24)

prefix = '/common/'
out_f.create_dataset(prefix + 'fsamp', data=fsamp)
out_f.create_dataset(prefix + 'nside', data=nside)
out_f.create_dataset(prefix + 'npsi', data=npsi)
out_f.create_dataset(prefix + 'det', data=np.string_(', '.join(det_list)))

polang = np.radians([0.0, 90.0, 45.0, 135.0, 0.0, 90.0])
mbang = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
out_f.create_dataset(prefix + 'polang', data=polang)
out_f.create_dataset(prefix + 'mbang', data=mbang)
out_f[prefix + 'polang'].attrs['legend'] = ', '.join(det_list)
out_f[prefix + 'mbang'].attrs['legend'] = ', '.join(det_list)


for chunk in chunk_list: 
    prefix = '/' + str(chunk+1).zfill(6) + '/common/'
    out_f.create_dataset(prefix + 'time', data=(chunk-1)*chunk_size)
    out_f[prefix + 'time'].attrs['type'] = 'second'

    out_f.create_dataset(prefix + 'vsun', data=np.random.normal(size=3))
    out_f[prefix + 'vsun'].attrs['info'] = '[x,y,z]'
    out_f[prefix + 'vsun'].attrs['coords'] = 'galactic'
    pix_array = [[],[],[]]
    i_start = chunk*nsamp
    i_stop = (chunk+1)*nsamp
    for det in det_list:
        det_file = h5py.File(os.path.join(tod_dir, det, segment) + '.hdf5', 'r')
        # pixels
        pixels = hp.ang2pix(nside, det_file['theta'][i_start:i_stop], det_file['phi'][i_start:i_stop])
        delta = np.diff(pixels)
        delta = np.insert(delta, 0, pixels[0])
        pix_array[0].append(delta)
        # psi
        psi_bins = np.linspace(0, 2*np.pi, num=npsi)
        psi_index = np.digitize(det_file['psi'][i_start:i_stop], psi_bins)
        delta = np.diff(psi_index)
        delta = np.insert(delta, 0, psi_index[0])
        pix_array[1].append(delta)
        # flag
        flag = np.ones(nsamp)
        delta = np.diff(flag)
        delta = np.insert(delta, 0, flag[0])
        pix_array[2].append(delta)

    h = huffman.Huffman("", nside)
    h.GenerateCode(pix_array)
    huffarray = np.append(np.append(np.array(h.node_max), h.left_nodes), h.right_nodes)
    out_f.create_dataset(prefix + 'hufftree', data=huffarray)
    out_f.create_dataset(prefix + 'huffsymb', data=h.symbols)

    for det in det_list:
        prefix = '/' + str(chunk+1).zfill(6) + '/' + det + '/'
        # signal
        out_f.create_dataset(prefix + 'tod', data=det_file['signal'][i_start:i_stop])
        # flag
        flag = np.ones(nsamp)
        delta = np.diff(flag)
        delta = np.insert(delta, 0, flag[0])
        out_f.create_dataset(prefix + 'flag', data=np.void(bytes(h.byteCode(delta))))
        # pixels
        pixels = hp.ang2pix(nside, det_file['theta'][i_start:i_stop], det_file['phi'][i_start:i_stop])
        delta = np.diff(pixels)
        delta = np.insert(delta, 0, pixels[0])
        out_f.create_dataset(prefix + 'pix', data=np.void(bytes(h.byteCode(delta))))
        # psi
        psi_bins = np.linspace(0, 2*np.pi, num=npsi)
        psi_index = np.digitize(det_file['psi'][i_start:i_stop], psi_bins)
        delta = np.diff(psi_index)
        delta = np.insert(delta, 0, psi_index[0])
        out_f.create_dataset(prefix + 'psi', data=np.void(bytes(h.byteCode(delta))))
        # scalars
        gain = 1.0
        sigma0 = 10.0
        fknee = 0.0
        alpha = 0.0
        out_f.create_dataset(prefix + 'scalars', data=np.array([gain, sigma0, fknee, alpha]).flatten())
        out_f[prefix + 'scalars'].attrs['legend'] = 'gain, sigma0, fknee, alpha'

out_f.close()
