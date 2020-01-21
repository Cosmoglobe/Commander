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
