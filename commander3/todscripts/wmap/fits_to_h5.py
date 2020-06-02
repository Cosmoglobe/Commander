import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
import h5py
import healpy as hp

from pytest import approx
from glob import glob
import multiprocessing as mp
from multiprocessing import Pool
from joblib import Parallel, delayed

import huffman


from scipy.interpolate import interp1d
from joblib import Parallel, delayed
import os


prefix = '/mn/stornext/d16/cmbco/bp/wmap/'

from time import sleep
from time import time as timer





def write_file_parallel(file_ind, i, obsid, obs_ind, daflags, TODs, gain_guesses,
        band_labels, band, psi_A, psi_B, pix_A, pix_B, fknee, alpha, n_per_day,
        ntodsigma, npsi, psiBins, nside, fsamp, pos, vel, time, compress=False):
    file_out = prefix + f'data/wmap_{band}_{str(file_ind+1).zfill(6)}_v6.h5'
    if os.path.exists(file_out):
        return
    dt0 = np.diff(time).mean()
    det_list = []
    # make huffman code tables
    # Pixel, Psi, Flag
    pixArray_A = [[], [], []]
    pixArray_B = [[], [], []]
    todArray = []
    for j in range(len(band_labels)):
        label = band_labels[j]
        if label[:-2] == band.upper():
            TOD = TODs[j]
            gain = gain_guesses[j]
            sigma_0 = TOD.std()
            scalars = np.array([gain, sigma_0, fknee, alpha])
                    
            tod = np.zeros(TOD.size)
            for n in range(len(TOD[0])):
                tod[n::len(TOD[0])] = TOD[:,n]
            todi = np.array_split(tod, n_per_day)[i]

            todInd = np.int32(ntodsigma*todi/(sigma_0*gain))
            delta = np.diff(todInd)
            delta = np.insert(delta, 0, todInd[0])
            todArray.append(delta)


            pix = np.array_split(pix_A[j//4], n_per_day)[i]
            delta = np.diff(pix)
            delta = np.insert(delta, 0, pix[0])
            pixArray_A[0].append(delta)

            pix = np.array_split(pix_B[j//4], n_per_day)[i]
            delta = np.diff(pix)
            delta = np.insert(delta, 0, pix[0])
            pixArray_B[0].append(delta)


            psi = np.array_split(psi_A[j//4], n_per_day)[i]
            psi = np.where(psi < 0,         2*np.pi+psi,    psi)
            psi = np.where(psi >= 2*np.pi,  psi - 2*np.pi,  psi)
            psiIndexes = np.digitize(psi, psiBins)
            delta = np.diff(psiIndexes)
            delta = np.insert(delta, 0, psiIndexes[0])
            pixArray_A[1].append(delta)

            psi = np.array_split(psi_B[j//4], n_per_day)[i]
            psi = np.where(psi < 0,         2*np.pi+psi,    psi)
            psi = np.where(psi >= 2*np.pi,  psi - 2*np.pi,  psi)
            psiIndexes = np.digitize(psi, psiBins)
            delta = np.diff(psiIndexes)
            delta = np.insert(delta, 0, psiIndexes[0])
            pixArray_B[1].append(delta)

            flags = np.array_split(daflags[:,j//4], n_per_day)[i]
            t0 = np.arange(len(flags))
            t = np.linspace(t0.min(), t0.max(), len(todi))
            func = interp1d(t0, flags, kind='previous')
            flags = func(t)
            delta = np.diff(flags)
            delta = np.insert(delta, 0, flags[0])
            pixArray_A[2].append(delta)
            pixArray_B[2].append(delta)


    h_A = huffman.Huffman("", nside)
    h_A.GenerateCode(pixArray_A)

    h_B = huffman.Huffman("", nside)
    h_B.GenerateCode(pixArray_B)

    h_Tod = huffman.Huffman("", nside)
    h_Tod.GenerateCode(todArray)

    huffarray_A = np.append(np.append(np.array(h_A.node_max), h_A.left_nodes), h_A.right_nodes)
    huffarray_B = np.append(np.append(np.array(h_B.node_max), h_B.left_nodes), h_B.right_nodes)
    huffarray_Tod = np.append(np.append(np.array(h_Tod.node_max), h_Tod.left_nodes), h_Tod.right_nodes)


    #with h5py.File(file_out, 'a') as f:
    #with h5py.File(file_out, 'w') as f:
    f = h5py.File(file_out, 'a')
    for j in range(len(band_labels)):
        label = band_labels[j]
        if label[:-2] == band.upper():
            TOD = TODs[j]
            gain = gain_guesses[j]
            sigma_0 = TOD.std()
            scalars = np.array([gain, sigma_0, fknee, alpha])


            tod = np.zeros(TOD.size)
            for n in range(len(TOD[0])):
                tod[n::len(TOD[0])] = TOD[:,n]
            todi = np.array_split(tod, n_per_day)[i]

            todInd = np.int32(ntodsigma*todi/(sigma_0*gain))
            deltatod = np.diff(todInd)
            deltatod = np.insert(deltatod, 0, todInd[0])


            pixA = np.array_split(pix_A[j//4], n_per_day)[i]
            deltapixA = np.diff(pixA)
            deltapixA = np.insert(deltapixA, 0, pixA[0])


            pixB = np.array_split(pix_B[j//4], n_per_day)[i]
            deltapixB = np.diff(pixB)
            deltapixB = np.insert(deltapixB, 0, pixB[0])


            psiA = np.array_split(psi_A[j//4], n_per_day)[i]
            psiA = np.where(psiA < 0,           2*np.pi+psiA,   psiA)
            psiA = np.where(psiA >= 2*np.pi,    psiA - 2*np.pi, psiA)
            psiIndexesA = np.digitize(psiA, psiBins)
            deltapsiA = np.diff(psiIndexesA)
            deltapsiA = np.insert(deltapsiA, 0, psiIndexesA[0])

            psiB = np.array_split(psi_B[j//4], n_per_day)[i]
            psiB = np.where(psiB < 0,           2*np.pi+psiB,   psiB)
            psiB = np.where(psiB >= 2*np.pi,    psiB - 2*np.pi, psiB)
            psiIndexesB = np.digitize(psiB, psiBins)
            deltapsiB = np.diff(psiIndexesB)
            deltapsiB = np.insert(deltapsiB, 0, psiIndexesB[0])

            flags = np.array_split(daflags[:,j//4], n_per_day)[i]
            t0 = np.arange(len(flags))
            t = np.linspace(t0.min(), t0.max(), len(todi))
            func = interp1d(t0, flags, kind='previous')
            flags = func(t)
            deltaflag = np.diff(flags)
            deltaflag = np.insert(deltaflag, 0, flags[0])


            if compress:
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/tod',
                        data=np.void(bytes(h_Tod.byteCode(deltatod))))

                f.create_dataset(obsid + '/' + label.replace('KA','Ka') + '/flag',
                        data=np.void(bytes(h_A.byteCode(deltaflag))))
                
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/pixA',
                        data=np.void(bytes(h_A.byteCode(deltapixA))))
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/pixB',
                        data=np.void(bytes(h_B.byteCode(deltapixB))))

                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/psiA',
                        data=np.void(bytes(h_A.byteCode(deltapsiA))))
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/psiB',
                        data=np.void(bytes(h_B.byteCode(deltapsiB))))
            else:
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/tod',
                        data=todInd)
                f.create_dataset(obsid + '/' + label.replace('KA','Ka') + '/flag',
                        data=flags)
                
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/pixA',
                        data=pixA)
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/pixB',
                        data=pixB)

                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/psiA',
                        data=psiA)
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/psiB',
                        data=psiB)


            det_list.append(label.replace('KA','Ka'))

            f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/scalars',
                    data=scalars)
            f[obsid + '/' + label.replace('KA','Ka') + '/scalars'].attrs['legend'] = 'gain, sigma0, fknee, alpha'
            # filler 
            f.create_dataset(obsid + '/' + label.replace('KA','Ka') + '/outP',
                    data=np.array([0,0]))



    if compress:
        f.create_dataset(obsid + '/common/todtree', data=huffarray_Tod)
        f.create_dataset(obsid + '/common/todsymb', data=h_Tod.symbols)

        f.create_dataset(obsid + '/common/hufftree_A', data=huffarray_A)
        f.create_dataset(obsid + '/common/huffsymb_A', data=h_A.symbols)

        f.create_dataset(obsid + '/common/hufftree_B', data=huffarray_B)
        f.create_dataset(obsid + '/common/huffsymb_B', data=h_B.symbols)
    
    f.create_dataset(obsid + '/common/satpos',
            data=np.array_split(pos,n_per_day)[i][0])
    f[obsid + '/common/satpos'].attrs['info'] = '[x, y, z]'
    f[obsid + '/common/satpos'].attrs['coords'] = 'galactic'

    f.create_dataset(obsid + '/common/vsun',
            data=np.array_split(vel,n_per_day)[i][0])
    f[obsid + '/common/vsun'].attrs['info'] = '[x, y, z]'
    f[obsid + '/common/vsun'].attrs['coords'] = 'galactic'

    dt = dt0/len(TOD[0])
    time_band = np.arange(time.min(), time.min() + dt*len(tod), dt)
    f.create_dataset(obsid + '/common/time',
            data=[np.array_split(time_band, n_per_day)[i][0],0,0])
    f[obsid + '/common/time'].attrs['type'] = 'MJD, null, null'

    f.create_dataset(obsid + '/common/ntod',
            data=len(np.array_split(tod,n_per_day)[i]))

    if "/common/fsamp" not in f:
        f.create_dataset('/common/fsamp', data=fsamp*len(TOD[0]))
        f.create_dataset('/common/nside', data=nside)
        f.create_dataset('/common/npsi', data=npsi)
        f.create_dataset('/common/det', data=np.string_(', '.join(det_list)))
        f.create_dataset('/common/datatype', data='WMAP')
        # fillers
        #f.create_dataset('/common/mbang', data=0)
        f.create_dataset('/common/ntodsigma', data=100)
        #f.create_dataset('/common/polang', data=0)
    with open(prefix + f'data/filelist_{band}_v6.txt', 'a') as file_list: 
        file_list.write(f'{str(obs_ind).zfill(6)}\t"{file_out}"\t1\t0\t0\n')
    return

def coord_trans(pos_in, coord_in, coord_out, lonlat=False):
    r = hp.rotator.Rotator(coord=[coord_in, coord_out])
    pos_out = r(pos_in.T).T

    if lonlat:
        if pos_out.shape[1] == 2:
            return pos_out
        elif pos_out.shape[1] == 3:
            return hp.vec2dir(pos_out.T, lonlat=True).T
    else:
        return pos_out



def Q2M(Q):
    '''
    PURPOSE:
        Converts quaternions to rotation matrices.
        
    CALLING SEQUENCE:
        M =Q2M(Q)
    
    INPUTS:
        Q - Quaternions.  May be a 2-D array dimensioned 4xN or
            simply a vector dimensioned 4.
    
    OUTPUTS:  
        M - Cube of attitude rotation matrices, 3x3xN (or 3x3
            if only one input quaternion).
    '''
    q1=-Q[0,:]
    q2=-Q[1,:]
    q3=-Q[2,:]
    q4= Q[3,:]
    
    q11=q1**2
    q22=q2**2
    q33=q3**2
    q44=q4**2
    s=q11+q22+q33+q44
    w = (abs(s-1.0) > 1e-5)
    if sum(w) > 0:
        s=np.sqrt(s)
        q1=q1/s
        q2=q2/s
        q3=q3/s
        q4=q4/s
    
    q12=q1*q2
    q13=q1*q3
    q14=q1*q4
    q23=q2*q3
    q24=q2*q4
    q34=q3*q4

    
    M = np.zeros((len(q1), 3,3))
    
    M[:,0,0] = q11 - q22 - q33 + q44
    M[:,0,1] = 2. * ( q12 + q34 )
    M[:,0,2] = 2. * ( q13 - q24 )
    M[:,1,0] = 2. * ( q12 - q34 )
    M[:,1,1] = -q11 + q22 - q33 + q44
    M[:,1,2] = 2. * ( q23 + q14 )
    M[:,2,0] = 2. * ( q13 + q24 )
    M[:,2,1] = 2. * ( q23 - q14 )
    M[:,2,2] = -q11 - q22 + q33 + q44

    M = np.transpose(M, [1,2,0])



    return M

def gamma_from_pol(gal, pol, fixed_basis=False):
    '''
    It should be possible to distribute the inner product among all time
    observations, but the matrix algebra is escaping me right now. In any case,
    for a one time operation this doesn't seem too slow yet.
    '''
    # gal and pol are galactic lonlat vectors
    dir_A_gal = hp.ang2vec(gal[:,0]%np.pi,gal[:,1]%(2*np.pi), lonlat=False)
    dir_A_pol = hp.ang2vec(pol[:,0]%np.pi,pol[:,1]%(2*np.pi), lonlat=False)

    dir_Z = np.array([0,0,1])


    sin_theta_A = np.sqrt(dir_A_gal[:,0]**2 + dir_A_gal[:,1]**2)

    dir_A_west_x = dir_A_gal[:,1]/sin_theta_A
    dir_A_west_y = -dir_A_gal[:,0]/sin_theta_A
    dir_A_west_z = dir_A_gal[:,1]*0
    dir_A_west = np.array([dir_A_west_x, dir_A_west_y, dir_A_west_z]).T
    dir_A_north = (dir_Z - dir_A_gal[2]*dir_A_gal)/sin_theta_A[:,np.newaxis]
    '''
    if sin_theta_A == 0:
        dir_A_west = np.array([1,0,0])
        dir_A_north = np.array([0,1,0])

    assert dir_A_north.dot(dir_A_west) == approx(0), 'Vectors not orthogonal'
    assert dir_A_north.dot(dir_A_north) == approx(1), 'North-vector not normalized'
    assert dir_A_west.dot(dir_A_west) == approx(1), 'North-vector not normalized'
    '''
    sin_gamma_A = dir_A_pol[:,0]*dir_A_west[:,0] + dir_A_pol[:,1]*dir_A_west[:,1] + dir_A_pol[:,2]*dir_A_west[:,2]
    cos_gamma_A = dir_A_pol[:,0]*dir_A_north[:,0] + dir_A_pol[:,1]*dir_A_north[:,1] + dir_A_pol[:,2]*dir_A_north[:,2]

    cos_2_gamma_A = 2*cos_gamma_A**2 - 1
    sin_2_gamma_A = 2*sin_gamma_A*cos_gamma_A

    return sin_2_gamma_A, cos_2_gamma_A

def q_interp(q_arr, t):
    '''
    Copied from interpolate_quaternions.pro

    This is an implementation of Lagrange polynomials.
    ;   input_q  - Set of 4 evenly-spaced quaternions (in a 4x4 array).
    ;          See the COMMENTS section for how this array should
    ;          be arranged.
    ;   offset   - Dimensionless time offset relative to the first quaternion.
    ;   This routine expects a unifomly sampled set of quaternions Q1,Q2,Q3,Q4.
    ;   It interpolate a quaternion for any time between Q1 and Q4, inclusive.
    ;   The output is calculated at a time T_Out, expressed in terms of the
    ;   sampling of the input quaternions:
    ;
    ;                   T_Out - T(Q1)
    ;       Offset = -----------------
    ;                   T(Q2) - T(Q1)
    ;
    ;   where T(Q1) is the time at quaternion Q1, and so forth.  That is,
    ;   the time for the output quaternion (variable OFFSET) should be
    ;   a number in the range -1.000 to 4.000 inclusive.  Input values outside
    ;   that range result in an error.  Input values outside 0.0 to 3.0 result
    ;   in extrapolation instead of interpolation.
    ;
    ;       In other words, Offset is essentially a floating point subscript,
    ;       similar to the those used by the IDL intrinsic routine INTERPOLATE.
    ;
    ;   For optimal results, OFFSET should be in the range [1.0, 2.0] -- that
    ;   is, the input quaternions Q1...Q4 should be arranged such that 2 come
    ;   before the desired output and 2 come after.

    '''
    xp0 = t-1
    xn0 = -xp0
    xp1 = xp0 + 1
    xn1 = xp0 - 1
    xn2 = xp0 - 2
    w = np.array([xn0*xn1*xn2/6, xp1*xn1*xn2/2, xp1*xn0*xn2/2, xp1*xp0*xn1/6])
    Qi = q_arr.dot(w)
    Qi = Qi/np.sum(Qi**2, axis=0)**0.5
    return Qi


def quat_to_sky_coords(quat, center=True):
    Nobs_array = np.array([12, 12, 15, 15, 20, 20, 30, 30, 30, 30])
    '''
    Quaternion is of form (N_frames, 30, 4), with one redundant frame at the
    beginning and two redundant ones at the end, that match the adjacent frames.
    '''
    nt = len(quat)
    Q = np.zeros( (4, 33, nt))
    q0 = quat[:,0::4]
    q1 = quat[:,1::4]
    q2 = quat[:,2::4]
    q3 = quat[:,3::4]
    q0 = np.array([q0[0,0]] + q0[:,1:-2].flatten().tolist() +
            q0[-1,-2:].tolist())
    q1 = np.array([q1[0,0]] + q1[:,1:-2].flatten().tolist() +
            q1[-1,-2:].tolist())
    q2 = np.array([q2[0,0]] + q2[:,1:-2].flatten().tolist() +
            q2[-1,-2:].tolist())
    q3 = np.array([q3[0,0]] + q3[:,1:-2].flatten().tolist() +
            q3[-1,-2:].tolist())
    Q = np.zeros((4, 30*nt + 3))
    Q[0] = q0
    Q[1] = q1
    Q[2] = q2
    Q[3] = q3
    t0 = np.arange(30*nt + 3)


    da_str = ''

    dir_A_los = np.array([
                [ 0.03997405,  0.92447851, -0.37913264],
                [-0.03834152,  0.92543237, -0.37696797],
                [-0.03156996,  0.95219303, -0.30386144],
                [ 0.03194693,  0.95220414, -0.3037872 ],
                [-0.03317037,  0.94156392, -0.33519711],
                [ 0.03336979,  0.94149584, -0.33536851],
                [-0.0091852 ,  0.93943624, -0.34260061],
                [-0.00950387,  0.94586233, -0.32442894],
                [ 0.00980826,  0.9457662 , -0.32470001],
                [ 0.00980739,  0.93934639, -0.34282965]])
    dir_B_los = np.array([
                [ 0.03795967, -0.92391895, -0.38070045],
                [-0.0400215 , -0.92463091, -0.37875581],
                [-0.03340367, -0.95176817, -0.30499432],
                [ 0.03014983, -0.95193039, -0.30482702],
                [-0.03504541, -0.94094355, -0.33674479],
                [ 0.03143652, -0.94113826, -0.33655687],
                [-0.01148033, -0.93883144, -0.3441856 ],
                [-0.01158651, -0.94535168, -0.32584651],
                [ 0.00767888, -0.9454096 , -0.32579398],
                [ 0.00751565, -0.93889159, -0.34413092]])

    dir_A_pol = np.array([  
                [ 0.69487757242271, -0.29835139515692, -0.65431766318192, ],
                [ -0.69545992357813, -0.29560553030986, -0.65494493291187, ],
                [  0.71383872060219, -0.19131247543171, -0.67367189173456, ],
                [ -0.71390969181845, -0.19099503229669, -0.67368675923286, ],
                [ -0.69832280289930, -0.26176968417604, -0.66619959126169, ],
                [  0.69826122350352, -0.26204606404493, -0.66615548040223, ],
                [  0.70944248806767, -0.23532277684296, -0.66431509603747, ],
                [ -0.70476543555624, -0.23649685267332, -0.66886091193973, ],
                [  0.70468980214241, -0.23690904054153, -0.66879472879665, ],
                [ -0.70959923775957, -0.23501806310177, -0.66425554705017]])
    dir_B_pol = np.array([  
                [ 0.69546590081501,  0.29798590641998, -0.65385899120425,],
                [ -0.69486414021667,  0.29814186328140, -0.65442742607568, ],
                [  0.71423586688235,  0.19072845484161, -0.67341650037147, ],
                [ -0.71357469183546,  0.19306390125546, -0.67345192048426, ],
                [ -0.69775710213559,  0.26425762446771, -0.66580998365151, ],
                [  0.69876566230957,  0.26145991550208, -0.66585678772745, ],
                [  0.71002796142313,  0.23471528678222, -0.66390438178103, ],
                [ -0.70422900931886,  0.23906270891214, -0.66851366750529, ],
                [  0.70521159225086,  0.23611413753036, -0.66852578425466, ],
                [ -0.70903152581832,  0.23766935833457, -0.66391834701609]])

    M = Q2M(Q)
    M = np.transpose(M, [2,0,1])

    gal_A = []
    pol_A = []
    gal_B = []
    pol_B = []
    for n, Nobs in enumerate(Nobs_array):
        # for each group from 0--4, the interpolation is valid between 1.5--2.5,
        # which is equivalent to cutting out the first 1.5 time units from the
        # beginning of the total array and the final set of quaternions does not
        # need the last half of the time interval.
        t = np.arange(t0.min() + 1.5, t0.max() - 0.5, 1/Nobs)

        M2 = np.zeros((len(t), 3, 3))
        for i in range(3):
            for j in range(3):
                f = interp1d(t0, M[:,i,j], kind='cubic')
                M2[:,i,j] = f(t)




        Npts = 30*nt*Nobs
        dir_A_los_cel = []
        dir_B_los_cel = []
        dir_A_los_cel = np.sum(M2*np.tile(dir_A_los[n, np.newaxis, np.newaxis,:], (Npts,3,1)),axis=2)
        dir_B_los_cel = np.sum(M2*np.tile(dir_B_los[n, np.newaxis, np.newaxis,:], (Npts,3,1)),axis=2)

        dir_A_los_gal = coord_trans(dir_A_los_cel, 'C', 'G')
        Pll_A = np.array(hp.vec2ang(dir_A_los_gal, lonlat=False))

        dir_B_los_gal = coord_trans(dir_B_los_cel, 'C', 'G')
        Pll_B = np.array(hp.vec2ang(dir_B_los_gal, lonlat=False))
        gal_A.append(Pll_A.T)
        gal_B.append(Pll_B.T)

        dir_A_pol_cel = np.sum(M2*np.tile(dir_A_pol[n, np.newaxis, np.newaxis,:], (Npts,3,1)),axis=2)
        dir_B_pol_cel = np.sum(M2*np.tile(dir_B_pol[n, np.newaxis, np.newaxis,:], (Npts,3,1)),axis=2)

        dir_A_pol_gal = coord_trans(dir_A_pol_cel, 'C', 'G')
        Pll_A = np.array(hp.vec2ang(dir_A_pol_gal, lonlat=False))

        dir_B_pol_gal = coord_trans(dir_B_pol_cel, 'C', 'G')
        Pll_B = np.array(hp.vec2ang(dir_B_pol_gal, lonlat=False))
        pol_A.append(Pll_A.T)
        pol_B.append(Pll_B.T)




    return gal_A, gal_B, pol_A, pol_B


def get_psi(gal, pol, band_labels):
    sin_2_gamma = np.zeros( (len(band_labels), len(gal[0])) )
    cos_2_gamma = np.zeros( (len(band_labels), len(gal[0])) )
    psi = []
    for band in range(len(band_labels)):
        sin2g, cos2g = gamma_from_pol(gal[band], pol[band])
        psi.append(0.5*np.arctan2(sin2g, cos2g))
    return psi

def get_psi_multiprocessing(gal, pol):
    sin_2_gamma = np.zeros(len(gal))
    cos_2_gamma = np.zeros(len(gal))
    for t in range(len(sin_2_gamma)):
        sin_2_gi, cos_2_gi = gamma_from_pol(gal[t], pol[t])
        sin_2_gamma[t] = sin_2_gi
        cos_2_gamma[t] = cos_2_gi
    psi = 0.5*np.arctan2(sin_2_gamma, cos_2_gamma)
    return psi

def get_psi_multiprocessing_2(i):
    gal = gals[i]
    pol = pols[i]
    sin_2_gamma = np.zeros(len(gal))
    cos_2_gamma = np.zeros(len(gal))
    for t in range(len(sin_2_gamma)):
        sin_2_gi, cos_2_gi = gamma_from_pol(gal[t], pol[t])
        sin_2_gamma[t] = sin_2_gi
        cos_2_gamma[t] = cos_2_gi
    psi = 0.5*np.arctan2(sin_2_gamma, cos_2_gamma)
    return psi

def ang2pix_multiprocessing(nside, theta, phi):
    return hp.ang2pix(nside, theta, phi)

def fits_to_h5(file_input, file_ind, plot):
    f_name = file_input.split('/')[-1][:-8]
    # It takes about 30 seconds for the extraction from the fits files, which is
    # very CPU intensive. After that, it maxes out at 1 cpu/process.
    file_out = prefix + f'data/wmap_K1_{str(file_ind+1).zfill(6)}_v6.h5'
    if os.path.exists(file_out):
        return
    t0 = timer()

    # from table 3 of astro-ph/0302222
    gain_guesses = np.array([   -0.974, +0.997,
                                +1.177, -1.122,
                                +0.849, -0.858,
                                -1.071, +0.985,
                                +1.015, -0.948,
                                +0.475, -0.518,
                                -0.958, +0.986,
                                -0.783, +0.760,
                                +0.449, -0.494,
                                -0.532, +0.532,
                                -0.450, +0.443,
                                +0.373, -0.346,
                                +0.311, -0.332,
                                +0.262, -0.239,
                                -0.288, +0.297,
                                +0.293, -0.293,
                                -0.260, +0.281,
                                -0.263, +0.258,
                                +0.226, -0.232,
                                +0.302, -0.286])

    alpha = -1
    fknee = 0.1

    nside = 256
    ntodsigma = 100
    npsi = 2048
    psiBins = np.linspace(0, 2*np.pi, npsi)
    fsamp = 30/1.536 # A single TOD record contains 30 1.536 second major science frames
    chunk_size = 1875
    nsamp = chunk_size*fsamp
    chunk_list = np.arange(25)
    # WMAP data divides evenly into 25 chunks per day...





    bands = ['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']
    #bands = ['K1']

    t2jd = 2.45e6
    jd2mjd = 2400000.5

    data = fits.open(file_input, memmap=False)

    band_labels = data[2].columns.names[1:-6]


    daflags = data[2].data['daflags']

    TODs = []
    for key in data[2].columns.names[1:-6]:
        TODs.append(data[2].data[key])

    
    
   
    # position (and velocity) in km(/s) in Sun-centered coordinates
    pos = data[1].data['POSITION']
    vel = data[1].data['VELOCITY']
    time_aihk = data[1].data['TIME'] + t2jd - jd2mjd
    time = data[2].data['TIME'] + t2jd - jd2mjd

    dt0 = np.diff(time).mean()
    
    if np.any(~np.isfinite(data[1].data['QUATERN'])):
        print(f'{file_input} has NaNs in the quaternion...')
        return
    gal_A, gal_B, pol_A, pol_B = quat_to_sky_coords(data[1].data['QUATERN'])
    # This file has NaNs in the quaternion???
    #/mn/stornext/d16/cmbco/bp/wmap/tod/wmap_tod_20013082358_20013091720_uncalibrated_v5.fits

    data.close()


    if plot:
        # Are gal and pol behaving as expected?
        for i in range(len(gal_A)):
            plt.figure('galtheta')
            plt.plot(np.linspace(0,1,len(gal_A[i][:,0])), gal_A[i][:,0],
                    color=plt.cm.viridis(i/len(gal_A)))
            plt.figure('poltheta')
            plt.plot(np.linspace(0,1,len(gal_A[i][:,0])), pol_A[i][:,0],
                    color=plt.cm.viridis(i/len(gal_A)))
            plt.figure('galphi')
            plt.plot(np.linspace(0,1,len(gal_A[i][:,0])), gal_A[i][:,1],
                    color=plt.cm.viridis(i/len(gal_A)))
            plt.figure('polphi')
            plt.plot(np.linspace(0,1,len(gal_A[i][:,0])), pol_A[i][:,1],
                    color=plt.cm.viridis(i/len(gal_A)))
        plt.figure('galtheta')
        fi = str(file_ind).zfill(5)
        plt.savefig(f'plots/pointing/galtheta_{fi}.png', bbox_inches='tight')
        plt.figure('galphi')
        plt.savefig(f'plots/pointing/galphi_{fi}.png', bbox_inches='tight')
        plt.figure('poltheta')
        plt.savefig(f'plots/pointing/poltheta_{fi}.png', bbox_inches='tight')
        plt.figure('polphi')
        plt.savefig(f'plots/pointing/polphi_{fi}.png', bbox_inches='tight')
        plt.close('all')

    psi_A = get_psi(gal_A, pol_A, band_labels[::4])
    psi_B = get_psi(gal_B, pol_B, band_labels[1::4])


    args_A = [(nside, gal_A[i][:,0], gal_A[i][:,1]) for i in range(len(gal_A))]
    args_B = [(nside, gal_B[i][:,0], gal_B[i][:,1]) for i in range(len(gal_B))]
    pix_A = []
    pix_B = []
    for i in range(len(args_A)):
        pix_A.append(ang2pix_multiprocessing(*args_A[i]))
        pix_B.append(ang2pix_multiprocessing(*args_B[i]))

    n_per_day = 1

    obs_inds = np.arange(n_per_day) + n_per_day*file_ind + 1
    obsids = [str(obs_ind).zfill(6) for obs_ind in obs_inds]
    for band in bands:
        args = [(file_ind, i, obsids[i], obs_inds[i], daflags, TODs, gain_guesses,
                    band_labels, band, psi_A, psi_B, pix_A, pix_B, fknee,
                    alpha, n_per_day, ntodsigma, npsi, psiBins, nside,
                    fsamp, pos, vel, time) for i in range(len(obs_inds))]
        for i in range(n_per_day):
            write_file_parallel(*args[i])


    #print(f'\t{f_name} took {int(timer()-t0)} seconds')

    return

def main(par=True, plot=False):
    '''
    Make 1 hdf5 file for every 10 fits files
    # Actually, just doing 1 hdf5 file for every fits file. Too much clashing is
    # happening right now.
    '''

    files = glob(prefix + 'tod/new/*.fits')
    files.sort()
    files = np.array(files)[:50]

    inds = np.arange(len(files))

    if par:
        nprocs = 90
        os.environ['OMP_NUM_THREADS'] = '1'

        pool = Pool(processes=nprocs)
        x = [pool.apply_async(fits_to_h5, args=[f, i, plot]) for i, f in zip(inds, files)]
        for i, res in enumerate(x):
            res.get()
            #res.wait()
        pool.close()
        pool.join()
    else:
        for i, f in zip(inds, files):
            print(i, f)
            fits_to_h5(f,i,plot)

    bands = ['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']
    for band in bands:
        file_name = prefix + f'data/filelist_{band}.txt'
        data = np.loadtxt(file_name, dtype=str)
        with open(file_name, 'r+') as f:
            content = f.read()
            f.seek(0,0)
            f.write(line.rstrip('\r\n') + '\n' + int(len(data)))
if __name__ == '__main__':
    '''
    If the file exists, skip it!
    '''
    main(par=True, plot=False)
