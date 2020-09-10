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

from tqdm import tqdm


prefix = '/mn/stornext/d16/cmbco/bp/wmap/'

version = 13

from time import sleep
from time import time as timer


from get_gain_model import get_gain



def write_file_parallel(file_ind, i, obsid, obs_ind, daflags, TODs, gain_guesses,
        baseline_guesses,
        band_labels, band, psi_A, psi_B, pix_A, pix_B, fknee, alpha, n_per_day,
        ntodsigma, npsi, psiBins, nside, fsamp, pos, vel, time, compress=False):
    file_out = prefix + f'data/wmap_{band}_{str(file_ind+1).zfill(6)}_v{version}.h5'
    if os.path.exists(file_out):
        return
    dt0 = np.diff(time).mean()
    det_list = []
    # make huffman code tables
    # Pixel, Psi, Flag
    pixArray = [[], [], []]
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
            pixArray[0].append(delta)

            pix = np.array_split(pix_B[j//4], n_per_day)[i]
            delta = np.diff(pix)
            delta = np.insert(delta, 0, pix[0])
            pixArray[0].append(delta)


            psi = np.array_split(psi_A[j//4], n_per_day)[i]
            psi = np.where(psi < 0,         2*np.pi+psi,    psi)
            psi = np.where(psi >= 2*np.pi,  psi - 2*np.pi,  psi)
            psiIndexes = np.digitize(psi, psiBins)
            delta = np.diff(psiIndexes)
            delta = np.insert(delta, 0, psiIndexes[0])
            pixArray[1].append(delta)

            psi = np.array_split(psi_B[j//4], n_per_day)[i]
            psi = np.where(psi < 0,         2*np.pi+psi,    psi)
            psi = np.where(psi >= 2*np.pi,  psi - 2*np.pi,  psi)
            psiIndexes = np.digitize(psi, psiBins)
            delta = np.diff(psiIndexes)
            delta = np.insert(delta, 0, psiIndexes[0])
            pixArray[1].append(delta)

            flags = np.array_split(daflags[:,j//4], n_per_day)[i]
            t0 = np.arange(len(flags))
            t = np.linspace(t0.min(), t0.max(), len(todi))
            func = interp1d(t0, flags, kind='previous')
            flags = func(t)
            delta = np.diff(flags)
            delta = np.insert(delta, 0, flags[0])
            # just necessary to make the array have the correct shape. Redundant
            # info.
            pixArray[2].append(delta)
            pixArray[2].append(delta)


    h = huffman.Huffman("", nside)
    h.GenerateCode(pixArray)


    h_Tod = huffman.Huffman("", nside)
    h_Tod.GenerateCode(todArray)

    huffarray = np.append(np.append(np.array(h.node_max), h.left_nodes), h.right_nodes)
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
            baseline = np.median(todi)
            todi = todi - baseline

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


            f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/tod',
                    data=np.int32(todi))
            if compress:
                f.create_dataset(obsid + '/' + label.replace('KA','Ka') + '/flag',
                        data=np.void(bytes(h.byteCode(deltaflag))))
                
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/pixA',
                        data=np.void(bytes(h.byteCode(deltapixA))))
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/pixB',
                        data=np.void(bytes(h.byteCode(deltapixB))))

                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/psiA',
                        data=np.void(bytes(h.byteCode(deltapsiA))))
                f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/psiB',
                        data=np.void(bytes(h.byteCode(deltapsiB))))
            else:
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
            # Subtracting baseline
            f.create_dataset(obsid + '/' + label.replace('KA','Ka')+ '/baseline',
                    data=np.array([baseline]))
            f[obsid + '/' + label.replace('KA','Ka') + '/baseline'].attrs['legend'] = 'baseline'
            # filler 
            f.create_dataset(obsid + '/' + label.replace('KA','Ka') + '/outP',
                    data=np.array([0,0]))



    if compress:
        f.create_dataset(obsid + '/common/todtree', data=huffarray_Tod)
        f.create_dataset(obsid + '/common/todsymb', data=h_Tod.symbols)

        f.create_dataset(obsid + '/common/hufftree', data=huffarray)
        f.create_dataset(obsid + '/common/huffsymb', data=h.symbols)
    
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
    with open(prefix + f'data/filelist_{band}_v{version}.txt', 'a') as file_list: 
        file_list.write(f'{str(obs_ind).zfill(6)}\t"{file_out}"\t1\t0\t0\n')
    return

def coord_trans(pos_in, coord_in, coord_out, lonlat=False):
    r = hp.rotator.Rotator(coord=[coord_in, coord_out])
    pos_out = r(pos_in.T).T

    if lonlat:
        if pos_out.shape[1] == 2:
            return pos_out
        elif pos_out.shape[1] == 3:
            return hp.vec2dir(pos_out.T, lonlat=lonlat).T
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
    # gal and pol are galactic lonlat vectors
    dir_gal = hp.ang2vec(gal[:,0]%np.pi,gal[:,1]%(2*np.pi), lonlat=False)
    dir_pol = hp.ang2vec(pol[:,0]%np.pi,pol[:,1]%(2*np.pi), lonlat=False)

    dir_Z = np.array([0,0,1])


    sin_theta = np.sqrt(dir_gal[:,0]**2 + dir_gal[:,1]**2)

    dir_west_x = dir_gal[:,1]/sin_theta
    dir_west_y = -dir_gal[:,0]/sin_theta
    dir_west_z = dir_gal[:,1]*0
    dir_west = np.array([dir_west_x, dir_west_y, dir_west_z]).T
    dir_north = (dir_Z - dir_gal[2]*dir_gal)/sin_theta[:,np.newaxis]


    sin_gamma = dir_pol[:,0]*dir_west[:,0] + dir_pol[:,1]*dir_west[:,1] + dir_pol[:,2]*dir_west[:,2]
    cos_gamma = dir_pol[:,0]*dir_north[:,0] + dir_pol[:,1]*dir_north[:,1] + dir_pol[:,2]*dir_north[:,2]

    cos_2_gamma = 2*cos_gamma**2 - 1
    sin_2_gamma = 2*sin_gamma*cos_gamma

    return sin_2_gamma, cos_2_gamma

def q_interp(q_arr, t):
    '''
    Copied from interpolate_quaternions.pro

    This is an implementation of Lagrange polynomials, equation 3.2.1 of
    numerical recipes 3rd edition.


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

    dir_A_los = np.array([
                [  0.03993743194318,  0.92448267167832, -0.37912635267982],
                [ -0.03836350153280,  0.92543717887494, -0.37695393578810],
                [ -0.03157188095163,  0.95219265474988, -0.30386241059657],
                [  0.03193385161530,  0.95220162163922, -0.30379647935526],
                [ -0.03317333754910,  0.94156429439011, -0.33519577742792],
                [  0.03337676771235,  0.94149468374332, -0.33537106592570],
                [ -0.00918939185649,  0.93943847522010, -0.34259437583453],
                [ -0.00950701394255,  0.94586439605663, -0.32442281201900],
                [  0.00980040822398,  0.94576779947882, -0.32469558276581],
                [  0.00980808738477,  0.93934799994236, -0.34282522723123]])
    dir_B_los = np.array([
                [  0.03794083653062, -0.92391755783762, -0.38070571212253],
                [ -0.04002167684949, -0.92463440201100, -0.37874726137612],
                [ -0.03340297596219, -0.95176877819247, -0.30499251475222],
                [  0.03014337784306, -0.95192770480751, -0.30483605690947],
                [ -0.03503633693827, -0.94094544143324, -0.33674045100040],
                [  0.03144454385558, -0.94113854675448, -0.33655530968115],
                [ -0.01147317267740, -0.93883247845653, -0.34418300902847],
                [ -0.01159000320270, -0.94535005109668, -0.32585112047876],
                [  0.00768184749607, -0.94540702221088, -0.32580139897397],
                [  0.00751408106677, -0.93889226303920, -0.34412912836731  ]])


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
        # offset = 1 + (k+0.5)/Nobs
        # or
        # offset = 1 + k/Nobs
        t = np.arange(t0.min() + 1, t0.max()-1, 1/Nobs) + 0.5/Nobs

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

def fits_to_h5(file_input, file_ind, compress, plot):
    f_name = file_input.split('/')[-1][:-8]
    # It takes about 30 seconds for the extraction from the fits files, which is
    # very CPU intensive. After that, it maxes out at 1 cpu/process.
    file_out = prefix + f'data/wmap_K1_{str(file_ind+1).zfill(6)}_v{version}.h5'
    if (os.path.exists(file_out) and file_ind != 1):
        return
    t0 = timer()

    # from programs.pars
    #       Gains and baselines given signs and values from the first
    #       attempt at cal_map for Pass 1.  These are median values
    #       from the hourly calibration files.
    #
    gain_guesses =np.array([ -0.9700,  0.9938,  1.1745, -1.1200, 
                              0.8668, -0.8753, -1.0914,  1.0033, 
                              1.0530, -0.9834,  0.4914, -0.5365, 
                             -0.9882,  1.0173, -0.8135,  0.7896, 
                              0.4896, -0.5380, -0.5840,  0.5840, 
                             -0.4948,  0.4872,  0.4096, -0.3802, 
                              0.3888, -0.4139,  0.3290, -0.3003, 
                             -0.3587,  0.3701,  0.3655, -0.3666, 
                             -0.3255,  0.3517, -0.3291,  0.3225, 
                              0.2841, -0.2918,  0.3796, -0.3591 ])
    baseline = [ 32136.98,  31764.96,  31718.19,  32239.29, 
                 31489.19,  32356.00,  32168.49,  31634.28, 
                 25621.62,  25502.45,  25500.11,  25667.74, 
                 26636.99,  24355.67,  26751.75,  24240.62, 
                 19050.87,  19129.09,  19380.14,  19081.39, 
                 19291.37,  18966.04,  18730.91,  19505.05, 
                 12428.31,  13567.44,  13049.69,  12930.21, 
                 13516.42,  12477.95,  12229.22,  13363.53, 
                 12678.23,  12934.65,  12730.91,  12692.85, 
                 11759.38,  13704.71,  11537.42,  13956.94  ]



    alpha = -1.7
    fknee = 0.1
    fknees = np.array([0.40, 0.51, 0.71, 0.32,
                       1.09, 0.35, 5.76, 8.62,
                       0.09, 1.41, 0.88, 8.35,
                       7.88, 0.66, 9.02, 7.47,
                       0.93, 0.28, 46.5, 26.0]) # mHz
    # From Table 2 of Jarosik et al. 2003, "On-orbit radiometer
    # characterization", alpha=-1.7 from Figure

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

    data = fits.open(file_input, memmap=False)

    band_labels = data[2].columns.names[1:-6]

    # Returns the gain model estimate at the start of each frame.
    gain_guesses = np.array([get_gain(data, b)[1][0] for b in band_labels])


    # If genflags == 1, there is an issue with the spacecraft attitude. Is this
    # the quaternion problem?
    genflags = data[2].data['genflags']
    daflags = data[2].data['daflags']
    print(genflags.shape)
    print(daflags.shape)

    TODs = []
    for key in band_labels:
        TODs.append(data[2].data[key])

    
    
   
    # position (and velocity) in km(/s) in Sun-centered coordinates
    pos = data[1].data['POSITION']
    vel = data[1].data['VELOCITY']
    # time2jd = 2.45e6, comverts table time (modified reduced Julian day) to Julian day, for both...
    time_aihk = data[1].data['TIME'] + t2jd






    time = data[2].data['TIME'] + t2jd

    dt0 = np.median(np.diff(time))

    quat = data[1].data['QUATERN']
    print(quat.shape)
    #if np.any(genflags != 0):
        #return
    if np.any(~np.isfinite(quat)):
        print(f'{file_input} has non-finite quaternions...')
        return
    gal_A, gal_B, pol_A, pol_B = quat_to_sky_coords(quat)

    data.close()


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
        args = [(file_ind, i, obsids[i], obs_inds[i], daflags, TODs, gain_guesses, baseline,
                    band_labels, band, psi_A, psi_B, pix_A, pix_B, fknee,
                    alpha, n_per_day, ntodsigma, npsi, psiBins, nside,
                    fsamp, pos, vel, time, compress) for i in range(len(obs_inds))]
        for i in range(n_per_day):
            write_file_parallel(*args[i])


    #print(f'\t{f_name} took {int(timer()-t0)} seconds')

    return

def main(par=True, plot=False, compress=False, nfiles=-1):
    '''
    Make 1 hdf5 file for every 10 fits files
    # Actually, just doing 1 hdf5 file for every fits file. Too much clashing is
    # happening right now.
    '''

    files = glob(prefix + 'tod/new/*.fits')
    files.sort()
    files = np.array(files)[:nfiles]

    inds = np.arange(len(files))

    if par:
        nprocs = 100
        os.environ['OMP_NUM_THREADS'] = '1'

        pool = Pool(processes=nprocs)
        x = [pool.apply_async(fits_to_h5, args=[f, i, compress, plot]) for i, f in zip(inds, files)]
        for i in tqdm(range(len(x))):
            x[i].get()
            #res.wait()
        pool.close()
        pool.join()
    else:
        for i, f in zip(inds, files):
            #print(i, f)
            fits_to_h5(f,i,compress, plot)

    # I don't know what this line of code was doing, so I will skip it...
    '''
    bands = ['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']
    for band in bands:
        file_name = prefix + f'data/filelist_{band}_v{version}.txt'
        data = np.loadtxt(file_name, dtype=str)
        with open(file_name, 'r+') as f:
            content = f.read()
            f.seek(0,0)
            f.write(line.rstrip('\r\n') + '\n' + int(len(data)))
    '''
if __name__ == '__main__':
    '''
    If the file exists, skip it!
    '''
    main(par=False, plot=False, compress=True)
