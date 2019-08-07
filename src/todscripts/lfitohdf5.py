import h5py
import argparse
import multiprocessing as mp
import os
import numpy as np
import math
from astropy.io import fits
import healpy as hp        
import sys

import huffman

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('planck_dir', type=str, action='store', help='path to the legacy planck data in hdf format')

    parser.add_argument('--gains-dir', type=str, action='store', help='path to a directory with the initial gain estimates', default=None)

    parser.add_argument('--velocity-file', type=str, action='store', help='path to a file with the satelite velocities', default=None)

    parser.add_argument('--rimo', type=str, action='store', help='path to on disk rimo file', default='/mn/stornext/d14/bp/data/auxiliary_data/LFI_RIMO_R3.31.fits')

    parser.add_argument('--out-dir', type=str, action='store', default=os.getcwd(), help='path to output data structure you want to generate')

    parser.add_argument('--num-procs', type=int, action='store', default=1, help='number of processes to use')

    parser.add_argument('--freqs', type=int, nargs='+', default=[30, 44, 70], help='which lfi frequencies to operate on')

    parser.add_argument('--ods', type=int, nargs=2, default=[91, 1540], help='the operational days to operate on')

    in_args = parser.parse_args()

    os.environ['OMP_NUM_THREADS'] = '1'

    pool = mp.Pool(processes=in_args.num_procs)
    manager = mp.Manager()

    ods = range(in_args.ods[0], in_args.ods[1], 1)

    outbufs = {30:manager.dict(), 44:manager.dict(), 70:manager.dict()}

    x = [[pool.apply_async(make_od, args=[freq, od, in_args, outbufs[freq]]) for freq in in_args.freqs] for od in ods]

    for res1 in np.array(x):
        for res in res1:
            res.get()

    pool.close()
    pool.join()

    #write file lists 
    for freq in in_args.freqs:
        outfile = open(os.path.join(in_args.out_dir, 'filelist_' + str(freq) + '.txt'), 'w')
        for buf in outbufs[freq].values():
            outfile.write(buf)

        outfile.close()


def make_od(freq, od, args, outbuf):

    print(freq, od)

    horns = {30:[27, 28], 44:[24, 25, 26], 70:[18, 19, 20, 21, 22, 23]}
    #psi_uv from https://www.aanda.org/articles/aa/full_html/2016/10/aa25818-15/T5.html
    mbangs = {27:-22.46, 28:22.45, 24:0.01, 25:-113.23, 26:113.23, 18:22.15, 19:22.4, 20:22.38, 21:-22.38, 22:-22.34, 23:-22.08}


    nsides = {30:512, 44:512, 70:1024}
    nside = nsides[freq]

    npsi = 4096

    outName = os.path.join(args.out_dir, 'LFI_0' + str(freq) + '_' + str(od).zfill(6) + '.h5')

    outFile = h5py.File(outName, 'w')

    exFile = h5py.File(os.path.join(args.planck_dir, 'LFI_0' + str(freq) + '_' + str(horns[freq][0]) + '_L2_002_OD' + str(od).zfill(4) +'.h5'), 'r')

    rimo = fits.open(args.rimo)

    if args.velocity_file is not None:
        velFile = fits.open(args.velocity_file)

    #make common group for things we only read once
    #polang, mbeamang, nside, fsamp, npsi
    prefix = '/common'
    rimo_i = np.where(rimo[1].data.field('detector').flatten() == 'LFI' + str(horns[freq][0]) + 'M')

    #sampling frequency
    fsamp = rimo[1].data.field('f_samp')[rimo_i]
    outFile.create_dataset(prefix + '/fsamp', data=fsamp)

    #nside
    outFile.create_dataset(prefix + '/nside', data=[nside])

    #psi angle resolution
    outFile.create_dataset(prefix + '/npsi', data=[npsi])    

    detNames = ''
    polangs = []
    mainbeamangs = []
    for horn in horns[freq]:
        for hornType in ['M', 'S']:
            rimo_i = np.where(rimo[1].data.field('detector').flatten() == 'LFI' + str(horn) + hornType)

            detNames += str(horn) + hornType + ', '
            polangs.append(math.radians(rimo[1].data.field('psi_pol')[rimo_i]))
            mainbeamangs.append(math.radians(mbangs[horn]))

    #make detector names lookup
    outFile.create_dataset(prefix + '/det', data=np.string_(detNames[0:-2]))

    #make polarization angle
    outFile.create_dataset(prefix + '/polang', data=polangs)
    outFile[prefix + '/polang'].attrs['legend'] = detNames[0:-2]

    #make main beam angle
    outFile.create_dataset(prefix + '/mbang', data=mainbeamangs)
    outFile[prefix + '/mbang'].attrs['legend'] = detNames[0:-2]


    #huffman coded bits
    for pid, index in zip(exFile['AHF_info/PID'], range(len(exFile['AHF_info/PID']))):
        startIndex = np.where(exFile['Time/OBT'] > exFile['AHF_info/PID_start'][index])
        endIndex = np.where(exFile['Time/OBT'] > exFile['AHF_info/PID_end'][index])
        if len(startIndex[0]) > 0:
            pid_start = startIndex[0][0]
        else:#catch days with no pids
            continue
        if len(endIndex[0]) is not 0:
            pid_end = endIndex[0][0]
        else:#catch final pid per od
            pid_end = len(exFile['Time/OBT'])
        if pid_start == pid_end:#catch chunks with no data like od 1007
            continue

        cut1 = exFile['Time/OBT'][exFile['Time/OBT'] > exFile['AHF_info/PID_start'][index]]


        #common fields
        prefix = str(pid).zfill(6) + '/common'

        #time field
        outFile.create_dataset(prefix + '/time', data=[exFile['Time/MJD'][pid_start]])
        outFile[prefix + '/time'].attrs['type'] = 'MJD'

        #velocity field
        velIndex = np.where(velFile[1].data.scet > exFile['Time/SCET'][pid_start])[0][0]
        #rotate from ecliptic to galactic 
        r = hp.Rotator(coord=['E', 'G'])
        outFile.create_dataset(prefix + '/vsun', data=r([velFile[1].data.xvel[velIndex], velFile[1].data.yvel[velIndex], velFile[1].data.zvel[velIndex]])) 
      
        #add some metadata so someone might be able to figure out what is going on 
        outFile[prefix + '/vsun'].attrs['info'] = '[x, y, z]'
        outFile[prefix + '/vsun'].attrs['coords'] = 'galactic'

        #make huffman code table
        pixArray = [[], [], []]
        for horn in horns[freq]:
            for hornType in ['M', 'S']:
                fileName = h5py.File(os.path.join(args.planck_dir, 'LFI_0' + str(freq) + '_' + str(horn) + '_L2_002_OD' + str(od).zfill(4) +'.h5'), 'r')
                rimo_i = np.where(rimo[1].data.field('detector').flatten() == 'LFI' + str(horn) + hornType)



                #get all pointing data
                newTheta, newPhi = r(fileName[str(horn) + hornType + '/THETA'][pid_start:pid_end], fileName[str(horn) + hornType + '/PHI'][pid_start:pid_end])
                pixels = hp.pixelfunc.ang2pix(nside, newTheta, newPhi)
                if len(pixels > 0):
                    delta = np.diff(pixels)
                    delta = np.insert(delta, 0, pixels[0])
                    pixArray[0].append(delta)

                #get all pol angle data
                psiArray = fileName[str(horn) + hornType + '/PSI'][pid_start:pid_end] + r.angle_ref(fileName[str(horn) + hornType + '/THETA'][pid_start:pid_end], fileName[str(horn) + hornType + '/PHI'][pid_start:pid_end]) + math.radians(rimo[1].data.field('psi_pol')[rimo_i])
                psiArray = np.where(psiArray < 0, 2*np.pi + psiArray, psiArray)
                psiArray = np.where(psiArray >= 2*np.pi, psiArray - 2*np.pi, psiArray)


                psiBins = np.linspace(0, 2*np.pi, num=4096)

                psiIndexes = np.digitize(psiArray, psiBins)
                if(len(psiIndexes) > 0):
                    delta = np.diff(psiIndexes)
                    delta = np.insert(delta, 0, psiIndexes[0])
                    pixArray[1].append(delta)

                #get all flag data
                flagArray = fileName[str(horn) + hornType + '/FLAG'][pid_start:pid_end]           
                if (len(flagArray) > 0):
                    delta = np.diff(flagArray)
                    delta = np.insert(delta, 0, flagArray[0])
                    pixArray[2].append(delta)

        h = huffman.Huffman("", nside)
        h.GenerateCode(pixArray)

        huffarray = np.append(np.append(np.array(h.node_max), h.left_nodes), h.right_nodes)

        outFile.create_dataset(prefix + '/hufftree', data=huffarray)
        outFile.create_dataset(prefix + '/huffsymb', data=h.symbols)

        for horn in horns[freq]:
            fileName = h5py.File(os.path.join(args.planck_dir, 'LFI_0' + str(freq) + '_' + str(horn) + '_L2_002_OD' + str(od).zfill(4) +'.h5'), 'r')
       
            for hornType in ['S', 'M']:
                prefix = str(pid).zfill(6) + '/' + str(horn) + hornType


                #get RIMO index
                #print(rimo[1].data.field('detector').flatten().shape, rimo[1].data.field('detector').flatten(), 'LFI' +str(horn) + hornType)
                rimo_i = np.where(rimo[1].data.field('detector').flatten() == 'LFI' + str(horn) + hornType)

                #make tod data 
                outFile.create_dataset(prefix + '/tod', data=fileName[str(horn) + hornType +'/SIGNAL'][pid_start:pid_end], dtype='f4')

                #undifferenced data? TODO
                #outFile.create_dataset(prefix + '/')
                
                #make flag data
                flagArray = fileName[str(horn) + hornType + '/FLAG'][pid_start:pid_end]
                
                if (len(flagArray) > 0):
                    delta = np.diff(flagArray)
                    delta = np.insert(delta, 0, flagArray[0])
                    outFile.create_dataset(prefix + '/flag', data=np.void(bytes(h.byteCode(delta))))
                
                #outFile.create_dataset(prefix + '/flag', data=flagArray, compression='gzip', shuffle=True)

                #make pixel number
                newTheta, newPhi = r(fileName[str(horn) + hornType + '/THETA'][pid_start:pid_end], fileName[str(horn) + hornType + '/PHI'][pid_start:pid_end])
                pixels = hp.pixelfunc.ang2pix(nside, newTheta, newPhi)
                
                if len(pixels > 0):
                    delta = np.diff(pixels)
                    delta = np.insert(delta, 0, pixels[0])
                    outFile.create_dataset(prefix + '/pix', data=np.void(bytes(h.byteCode(delta))))
                                
                #outFile.create_dataset(prefix + '/pix', data=pixels, compression='gzip', shuffle=True)

                #make pol angle
                psiArray = fileName[str(horn) + hornType + '/PSI'][pid_start:pid_end] + r.angle_ref(fileName[str(horn) + hornType + '/THETA'][pid_start:pid_end], fileName[str(horn) + hornType + '/PHI'][pid_start:pid_end]) + math.radians(rimo[1].data.field('psi_pol')[rimo_i])
                psiArray = np.where(psiArray < 0, 2*np.pi + psiArray, psiArray)
                psiArray = np.where(psiArray >= 2*np.pi, psiArray - 2*np.pi, psiArray)
                psiBins = np.linspace(0, 2*np.pi, num=4096)

                psiIndexes = np.digitize(psiArray, psiBins)
                
                if(pid == 3798 and horn == 28 and hornType == 'M'):
                    print(len(psiIndexes))
                    np.set_printoptions(threshold=sys.maxsize)
                    for i in range(4000):
                        print(i, psiArray[i], psiIndexes[i])

                if(len(psiIndexes) > 0):
                    delta = np.diff(psiIndexes)
                    delta = np.insert(delta, 0, psiIndexes[0]) 
                    outFile.create_dataset(prefix + '/psi', data=np.void(bytes(h.byteCode(delta))))
                
                #outFile.create_dataset(prefix + '/psi', data=psiIndexes, compression='gzip', shuffle=True)

                #scalars
                
                gain = 1
                #make gain
                if args.gains_dir is not None:
                    gainFile = fits.open(os.path.join(args.gains_dir, 'LFI_0' + str(freq) + '_LFI' + str(horn) + hornType + '_001.fits'))
                    gain=1.0/gainFile[1].data.GAIN[np.where(gainFile[1].data.PID == pid)]
                    gainFile.close()
 
                #make white noise
                sigma0 = rimo[1].data.field('net')[rimo_i] * math.sqrt(fsamp)

                #make f_knee
                fknee = rimo[1].data.field('f_knee')[rimo_i]

                #make 1/f noise exponent 
                alpha = rimo[1].data.field('alpha')[rimo_i]

                outFile.create_dataset(prefix + '/scalars', data=np.array([gain, sigma0, fknee, alpha]).flatten())
                outFile[prefix + '/scalars'].attrs['legend'] = 'gain, sigma0, fknee, alpha'

                #make psd noise
                
                #make other
 
                #write to output file
                outbuf['id' + str(pid)] = str(pid) + ' "' + outName + '" ' + '1\n'

if __name__ == '__main__':
    main()
