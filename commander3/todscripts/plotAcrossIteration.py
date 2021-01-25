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
import argparse
import h5py
import os

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('data_dir', type=str, action='store', help='path to the data files')

    parser.add_argument('--dets', type=str, nargs='+', action='store', help='detectors to plot')

    parser.add_argument('--fields', type=str, nargs='+', action='store', help='fields to plot')

    parser.add_argument('--ods', type=int, nargs=2, action='store', help='ods to plot from', default=[91, 1604])

    in_args = parser.parse_args()

    freqs = {18:70, 19:70, 20:70, 21:70, 22:70, 23:70, 24:44, 25:44, 26:44, 27:30, 28:30}
    
    used_freqs = []
    for det in in_args.dets:
        freq = freqs[int(det[0:-1])]
        if not freq in used_freqs:
            used_freqs.append(freq)

    fields = {}
    for field in in_args.fields:
        for det in in_args.dets:
            fields[field+det] = []

    print(used_freqs)

    for od in range(in_args.ods[0], in_args.ods[1]):
        for freq in used_freqs:
            f = h5py.File(os.path.join(in_args.data_dir, 'LFI_0' + str(freq) + '_' + str(od).zfill(6) + '.h5'))
            for group in f['/']:
                if 'common' not in group:
                    for field in in_args.fields:
                        newfield = field
                        scalars = ['gain', 'sigma0', 'fknee', 'alpha']
                        if field in scalars:
                            newfield = 'scalars'
                            
                        for det in in_args.dets:
                            #print('/' + group + '/' + det + '/' + newfield)
                            data = f['/' + group + '/' + det + '/' + newfield]
                            if newfield is 'scalars':
                                fields[field+det].append(data[scalars.index(field)])
                            else:
                                fields[field+det].append(data)                

    for field in in_args.fields:
        plt.figure()
        plt.title(field)
        for det in in_args.dets:
            plt.plot(fields[field+det], label=det)
            plt.ylim([0,0.1])
        plt.legend(loc='best')
        plt.savefig(field+'.png')
        plt.close()

if __name__ == '__main__':
    main()
