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

from commander_tools.tod_tools.hfi import hfi
from commander_tools.tod_tools import commander_tod as tod
import argparse
import multiprocessing as mp
import os
import numpy as np
import math
from astropy.io import fits
import healpy as hp        
import sys
import random
import h5py
import sqlite3
import traceback
import cProfile
import glob
from scipy.spatial.transform import Rotation as rot
import pickle
import codecs

import matplotlib.pyplot as plt

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--planck_dir', type=str, action='store', help='path to the legacy planck data in hdf format', default='/mn/stornext/d16/cmbco/bp/HFI/hfi_miss03_adc')

    parser.add_argument('--gains-dir', type=str, action='store', help='path to a directory with the initial gain estimates', default='/mn/stornext/d16/cmbco/bp/HFI/aux/gains')

    parser.add_argument('--rimo', type=str, action='store', help='path to on disk rimo file', default='/mn/stornext/d16/cmbco/bp/HFI/aux/RIMO_npipe2.fits')

    parser.add_argument('--pid-database', type=str, action='store', help='path to the sql database storing the PID info', default='/mn/stornext/d16/cmbco/bp/HFI/aux/hfi_raw_rings_v3.db')

    parser.add_argument('--extra-flags', type=str, action='store', help='path to extra flagging in txt file', default='/mn/stornext/d16/cmbco/bp/HFI/aux/hfi_bad_intervals_15s_elephants.txt')

    #https://github.com/planck-npipe/toast-npipe/blob/master/toast_planck/preproc_modules/transf1_nodemod.py
    parser.add_argument('--calib-params', type=str, action='store', help='hash dump of HFI housekeeping imo from npipe', default='/mn/stornext/d16/cmbco/bp/HFI/aux/hficalibparams.dat')

    parser.add_argument('--out-dir', type=str, action='store', default=os.getcwd(), help='path to output data structure you want to generate')

    parser.add_argument('--num-procs', type=int, action='store', default=1, help='number of processes to use')

    parser.add_argument('--freqs', type=int, nargs='+', default=hfi.freqs, help='which hfi frequencies to operate on')

    parser.add_argument('--ods', type=int, nargs=2, default=[91, 975], help='the operational days to operate on')

    parser.add_argument('--no-compress', action='store_true', default=False, help='Produce uncompressed data output')

    parser.add_argument('--restart', action='store_true', default=False, help="restart from a previous run that didn't finish")

    parser.add_argument('--produce-filelist', action='store_true', default=False, help='force the production of a filelist even if only some files are present')

    in_args = parser.parse_args()

    in_args.version = 1
    
    random.seed()

    os.environ['OMP_NUM_THREADS'] = '1'

    pool = mp.Pool(processes=in_args.num_procs)
    manager = mp.Manager()

    ods = range(in_args.ods[0], in_args.ods[1], 1)

    manager = mp.Manager()
    dicts = {}
    for freq in hfi.freqs:
        dicts[freq] = manager.dict()

    comm_tod = tod.commander_tod(in_args.out_dir, 'HFI', in_args.version, dicts, not in_args.restart)

    x = [[pool.apply_async(make_od, args=[comm_tod, freq, od, in_args]) for freq in in_args.freqs] for od in ods]

    for res1 in np.array(x):
        for res in res1:
            res.get()

    pool.close()
    pool.join()

    if ((in_args.ods[0] == 91 and in_args.ods[1] == 975) or in_args.produce_filelist) :
        comm_tod.make_filelists()
        #write file lists 

def make_od(comm_tod, freq, od, args):

    try:
      
        #load housekeeping imo from file 
        calibfile = open(args.calib_params, mode='r')
        calibparams = calibfile.read()
        calibfile.close()
        calibparams = bytes(calibparams, encoding='utf-8')
        hsk = pickle.loads(codecs.decode(bytes(calibparams), 'base64'), fix_imports=True, encoding='bytes')
 
        nside = hfi.nsides[freq]

        comm_tod.init_file(freq, od, mode='w')

        if(args.restart and comm_tod.exists):
            comm_tod.finalize_file()
            print('Skipping existing file ' + comm_tod.outName)
            return

        rimo = fits.open(args.rimo)

        #make common group for things we only read once
        #polang, mbeamang, nside, fsamp
        prefix = '/common'
        rimo_i = np.where(rimo[1].data.field('detector').flatten() == str(freq) + '-' + hfi.dets[freq][0])

        #sampling frequency
        fsamp = rimo[1].data.field('f_samp')[rimo_i]
        comm_tod.add_field(prefix + '/fsamp', fsamp)

        #nside
        comm_tod.add_field(prefix + '/nside', [nside])

        detNames = ''
        polangs = []
        mainbeamangs = []
        for det in hfi.dets[freq]:
            rimo_i = np.where(rimo[1].data.field('detector').flatten() == str(freq) + '-' + det)

            detNames += str(freq) + '-' + det + ', '
            polangs.append(math.radians(rimo[1].data.field('psi_pol')[rimo_i]))
            mainbeamangs.append(math.radians(rimo[1].data.field('psi_uv')[rimo_i]))

        compArr=[hfi.huffman]
        if args.no_compress:
            compArr = None

        #make detector names lookup
        comm_tod.add_field(prefix + '/det', np.string_(detNames[0:-2]))

        #make polarization angle
        comm_tod.add_field(prefix + '/polang', polangs)
        comm_tod.add_attribute(prefix + '/polang', 'index', detNames[0:-2])

        #make main beam angle
        comm_tod.add_field(prefix + '/mbang', mainbeamangs)
        comm_tod.add_attribute(prefix + '/mbang', 'index', detNames[0:-2])

        dataFileName = glob.glob(os.path.join(args.planck_dir, str(od).zfill(4), 'H' + str(freq) + '_' + str(od).zfill(4) + '_R_201505??.fits'))[0]

        try:
            exFile = fits.open(dataFileName)
        except (OSError):
            print("Failed to open file " + dataFileName)
            return

        try:
            pointingFile = fits.open(os.path.join(args.planck_dir, str(od).zfill(4), 'pointing-' + str(od).zfill(4)  + '.fits'))
        except (OSError):
            print("Failed to open file " + os.path.join(args.planck_dir, str(od).zfill(4), 'pointing-' + str(od).zfill(4)  + '.fits'))
            return

        try:
            extraFlagsFile = h5py.File(os.path.join(args.planck_dir, str(od).zfill(4), 'extra_flags_' + str(od).zfill(4)  + '.h5'), 'r')
        except (OSError):
            print("Failed to open file " + s.path.join(args.planck_dir, str(od).zfill(4), 'extra_flags_' + str(od).zfill(4)  + '.h5'))
            return


        #open database with pid info in it
        conn = sqlite3.connect(args.pid_database) 
        c = conn.cursor()

        starttime = pointingFile[1].data['obt'][0]/1e9
        endtime = pointingFile[1].data['obt'][-1]/1e9

        #per pid
        for dbentry in c.execute("SELECT * FROM ring_times_hfi WHERE stop_time >= '{0}' AND start_time < '{1}'".format(starttime, endtime)):
            pid = dbentry[0]
            start_time = dbentry[2]
            end_time = dbentry[3]        

            startIndex = np.where(exFile[1].data['obt']/1e9 > start_time)
            endIndex = np.where(exFile[1].data['obt']/1e9 > end_time)
            if len(startIndex[0]) > 0:
                pid_start = startIndex[0][0]
            else:#catch days with no pids
                continue
            if len(endIndex[0]) != 0:
                pid_end = endIndex[0][0] 
            else:#catch final pid per od
                pid_end = len(exFile[1].data['obt'])
            if pid_start == pid_end:#catch chunks with no data like od 1007
                continue

            #common fields per pid
            prefix = str(pid).zfill(6) + '/common'

            #time field
            comm_tod.add_field(prefix + '/time', [exFile[1].data['obt'][pid_start]])
            comm_tod.add_attribute(prefix + '/time','index','OBT in s')

            #length of the tod
            comm_tod.add_field(prefix + '/ntod', [pid_end - pid_start])

            #velocity field
            #rotate from ecliptic to galactic
            #might be required, I'm not sure what coordinate system these are in 
            r = hp.Rotator(coord=['E', 'G'])
            comm_tod.add_field(prefix + '/vsun', r([pointingFile[4].data['x_vel'][pid_start], pointingFile[4].data['y_vel'][pid_start], pointingFile[4].data['z_vel'][pid_start]])) 
          
            #add some metadata so someone might be able to figure out what is going on 
            comm_tod.add_attribute(prefix + '/vsun','index', '[x, y, z]')
            comm_tod.add_attribute(prefix + '/vsun','coords','galactic')

            #satelite position
            comm_tod.add_field(prefix + '/satpos', [pointingFile[5].data['x_pos'][pid_start], pointingFile[5].data['y_pos'][pid_start], pointingFile[5].data['z_pos'][pid_start]])
            #add metadata
            comm_tod.add_attribute(prefix + '/satpos','index','X, Y, Z')
            comm_tod.add_attribute(prefix + '/satpos','coords','heliocentric')
            comm_tod.add_attribute(prefix + '/satpos', 'units', 'AU')

            #read boresight pointing for this chunk
            quat_x = pointingFile[3].data['quaternion_x'][pid_start:pid_end]
            quat_y = pointingFile[3].data['quaternion_y'][pid_start:pid_end]
            quat_z = pointingFile[3].data['quaternion_z'][pid_start:pid_end]
            quat_s = pointingFile[3].data['quaternion_s'][pid_start:pid_end]

            quat_arr = np.array([quat_x, quat_y, quat_z, quat_s]).transpose()

            r_boresight = rot.from_quat(quat_arr)

            extra_flags = extraFlagsFile[str(pid).zfill(6) + '/flag_extra']


            #per detector fields
            for det in hfi.dets[freq]:
           
                prefix = str(pid).zfill(6) + '/' + str(freq) + '-' + det

                #get RIMO index
                #print(rimo[1].data.field('detector').flatten().shape, rimo[1].data.field('detector').flatten(), 'LFI' +str(horn) + hornType)
                rimo_i = np.where(rimo[1].data.field('detector').flatten() == str(freq) + '-' + det)
                data_i = exFile.index_of(str(freq) + '-' + det)        
            
                #make flag data
                flagArray = exFile[data_i].data.field('flag')[pid_start:pid_end]

                #add extra flagging from txt files
                ex_flags = extra_flags
                if(det == '353-1'):
                    ex_flags = extraFlagsFile[str(pid).zfill(6) + '/flags_extra_353-1']
                flagArray += ex_flags
                    
                if (len(flagArray) > 0):
                    comm_tod.add_field(prefix + '/flag', flagArray, compArr)

                #make pixel number 
                
                #create detector pointing offset quaternion
                # Follows this npipe function: 
                # https://github.com/planck-npipe/toast-npipe/blob/master/toast_planck/utilities.py#L1503-L1575            

                phi = math.radians(rimo[1].data.field('phi_uv')[rimo_i])
                theta = math.radians(rimo[1].data.field('theta_uv')[rimo_i])
                psi = math.radians(rimo[1].data.field('psi_uv')[rimo_i] + rimo[1].data.field('psi_pol')[rimo_i]) - phi

                det_s = np.cos(0.5 * theta) * np.cos(0.5 * (phi + psi))
                # vector part
                det_x = -np.sin(0.5 * theta) * np.sin(0.5 * (phi - psi))
                det_y = np.sin(0.5 * theta) * np.cos(0.5 * (phi - psi))
                det_z = np.cos(0.5 * theta) * np.sin(0.5 * (phi + psi))

                #convert from boresight pointing to per detector pointing
                r_det = rot.from_quat([det_x, det_y, det_z, det_s])

                r_total = r_boresight * r_det

                #convert to theta, phi, psi
                angs = r_total.as_euler('ZYZ') # could also be 'ZXZ' if we are supposed to be using intrinsic rotations instead of extrinsic
                # idk what the difference is

                phi_array = angs[:,0]
                theta_array = angs[:,1]
                psi_array = angs[:,2]
                
                r = hp.rotator.Rotator(coord=['E', 'G'], deg=False)

                galTheta, galPhi = r(theta_array, phi_array)
                
                pixels = hp.pixelfunc.ang2pix(nside, galTheta, galPhi)

                if(pid == 25083):

                    x = np.arange(0, 20000)

                    theta_, phi_ = hp.pix2ang(nside, pixels)

                    plt.plot(x, theta_[-20000:], label='theta')
                    plt.plot(x, phi_[-20000:], label='phi')
                    plt.plot(x, quat_x[-20000:], label='x')
                    plt.plot(x, quat_y[-20000:], label='y')
                    plt.plot(x, quat_z[-20000:], label='z')
                    plt.plot(x, quat_s[-20000:], label='s')
                    plt.legend(loc='best')

                    plt.savefig('quat_test.pdf')
                    #sys.exit()


                if len(pixels > 0):
                    #compute average outer product
                    outAng = hfi.ring_outer_product(galTheta, galPhi)
                    comm_tod.add_field(prefix + '/outP', data=outAng)
                    if(args.no_compress):
                        comm_tod.add_field(prefix+'/theta', data=galTheta)
                        comm_tod.add_field(prefix+'/phi', data=galPhi) 
                    comm_tod.add_field(prefix + '/pix', pixels, compArr)


                #make pol angle
                if(len(psi_array) > 0):
                    psi_array = np.where(psi_array < 0, 2*np.pi + psi_array, psi_array)
                    psi_array = np.where(psi_array >= 2*np.pi, psi_array - 2*np.pi, psi_array)
                    compArray = None
                    if not args.no_compress:
                        compArray = [hfi.psiDigitize, hfi.huffman]            
                    comm_tod.add_field(prefix + '/psi', psi_array, compArray)
                    
                #scalars
                    
                #make gain

                #TODO: figure this shit out
                # yikes the gain is complicated as hell
                # There is more information needed here from the housekeeping data
                # From Reijo:
                # To convert DSP-valued HFI data into K_CMB, you also need to calibrate twice. First from DSP into pseudo-Volts. Then from Volts to K_CMB.  

                # The first step is performed in this operator: https://github.com/planck-npipe/toast-npipe/blob/master/toast_planck/preproc_modules/transf1_nodemod.py The second step involves multiplying the data with the factors in the gain files you mention. After these, there will be small time-dependent fluctuations coming from residual ADC nonlinearity and (possibly) from real gain fluctuations.

                # Aside from calibration, you will need to demodulate out the square wave modulation between the two calibration steps and correct for bolometric nonlinearity. Please see this part of the preprocessing operator: https://github.com/planck-npipe/toast-npipe/blob/master/toast_planck/preproc.py#L921-L939.

                gainFile = fits.open(os.path.join(args.gains_dir, 'npipe5_gains_' + str(freq) + '.fits'))
                gain = gainFile[gainFile.index_of(str(freq) +'-' + det)].data[0][0]

                gain1, offset = hfi.compute_l1_gain(str(freq) + '-' + det, exFile[1].data['obt'][pid_start:pid_end], hsk)


                #print(gain, gain1, offset, gain/gain1)

                gain = gain/gain1[0]

                #make white noise
                sigma0 = rimo[1].data.field('net')[rimo_i][0] * math.sqrt(fsamp)

                #make f_knee
                #fknee = rimo[1].data.field('f_knee')[rimo_i][0]
                fknee = 0.5

                #make 1/f noise exponent 
                #alpha = rimo[1].data.field('alpha')[rimo_i][0]
                alpha = -1		

                #print(gain, sigma0, fknee, alpha)
                comm_tod.add_field(prefix + '/scalars', np.array([gain, sigma0, fknee, alpha]).flatten())
                comm_tod.add_attribute(prefix + '/scalars','index','gain, sigma0, fknee, alpha')

                #make psd noise
                   
                #make tod data
                tod = exFile[exFile.index_of(str(freq)+'-' + det)].data.field('signal')[pid_start:pid_end] - offset

                #compArray = [hfi.todDtype, hfi.rice]
                compArray = [hfi.todDtype, hfi.huffTod]
                if(args.no_compress):
                    compArray = [hfi.todDtype] 
                comm_tod.add_field(prefix + '/tod', tod, compArray)

            print(freq, od, pid)                 
            comm_tod.finalize_chunk(pid, loadBalance=outAng)
        comm_tod.finalize_file()
    except:
        exc_type, exc_val, exc_tb = sys.exc_info()
        raise Exception("".join(traceback.format_exception(exc_type, exc_val, exc_tb)))


if __name__ == '__main__':
    main()
