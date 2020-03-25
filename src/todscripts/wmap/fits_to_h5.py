import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
import h5py
import healpy as hp

from pytest import approx


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
    
    q11=q1*q1
    q22=q2*q2
    q33=q3*q3
    q44=q4*q4
    s=q11+q22+q33+q44
    w = (abs(s-1.0) > 1.0e-5)
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
    dir_A_gal = hp.ang2vec(gal[0],gal[1], lonlat=True)

    dir_A_pol = np.zeros((len(pol), 3))
    for i in range(len(pol)):
        dir_A_pol[i] = hp.ang2vec(pol[i][0],pol[i][1], lonlat=True)


    dir_Z = np.array([0,0,1])

    sin_theta_A = np.sqrt(dir_A_gal[0]**2 + dir_A_gal[1]**2)

    if sin_theta_A != 0:
        dir_A_west_x = dir_A_gal[1]/sin_theta_A
        dir_A_west_y = -dir_A_gal[0]/sin_theta_A
        dir_A_west_z = dir_A_gal[1]*0
        dir_A_west = np.array([dir_A_west_x, dir_A_west_y, dir_A_west_z])
        dir_A_north = (dir_Z - dir_A_gal[2]*dir_A_gal)/sin_theta_A
    else:
        dir_A_west = np.array([1,0,0])
        dir_A_north = np.array([0,1,0])

    assert dir_A_north.dot(dir_A_west) == approx(0), 'Vectors not orthogonal'
    assert dir_A_north.dot(dir_A_north) == approx(1), 'North-vector not normalized'
    assert dir_A_west.dot(dir_A_west) == approx(1), 'North-vector not normalized'


    sin_gamma_A = dir_A_pol.dot(dir_A_west)
    cos_gamma_A = dir_A_pol.dot(dir_A_north)

    cos_2_gamma_A = 2*cos_gamma_A**2 - 1
    sin_2_gamma_A = 2*sin_gamma_A*cos_gamma_A

    return sin_2_gamma_A, cos_2_gamma_A

def quat_to_sky_coords(quat, center=True, Nobs=12):
    '''
    Quaternion is of form (N_frames, 30, 4)
    '''
    nt = len(quat)
    Q = np.zeros( (4, 30, nt))
    for i in range(nt):
        for j in range(30):
            qt = quat[i,j:j+4]
            qout = qt
            Q[:,j,i] = qt
    Q = np.reshape(Q, (4, 30*nt))

    da_str = ''

    dir_A_LOS = np.array([0, 0.94264149, -0.33380686])#   ; Optical axis
    dir_B_LOS = np.array([0,-0.94264149, -0.33380686])    

    dir_A_POL = np.array([  
                [ 0.69487757242271, -0.29835139515692, -0.65431766318192, ],
                [ -0.69545992357813, -0.29560553030986, -0.65494493291187, ],
                [  0.71383872060219, -0.19131247543171, -0.67367189173456, ],
                [ -0.71390969181845, -0.19099503229669, -0.67368675923286, ],
                [ -0.69832280289930, -0.26176968417604, -0.66619959126169, ],
                [  0.69826122350352, -0.26204606404493, -0.66615548040223, ],
                [  0.70944248806767, -0.23532277684296, -0.66431509603747, ],
                [ -0.70476543555624, -0.23649685267332, -0.66886091193973, ],
                [  0.70468980214241, -0.23690904054153, -0.66879472879665, ],
                [ -0.70959923775957, -0.23501806310177, -0.66425554705017
                    ]])
    dir_B_POL = np.array([  
                [ 0.69546590081501,  0.29798590641998, -0.65385899120425,],
                [ -0.69486414021667,  0.29814186328140, -0.65442742607568, ],
                [  0.71423586688235,  0.19072845484161, -0.67341650037147, ],
                [ -0.71357469183546,  0.19306390125546, -0.67345192048426, ],
                [ -0.69775710213559,  0.26425762446771, -0.66580998365151, ],
                [  0.69876566230957,  0.26145991550208, -0.66585678772745, ],
                [  0.71002796142313,  0.23471528678222, -0.66390438178103, ],
                [ -0.70422900931886,  0.23906270891214, -0.66851366750529, ],
                [  0.70521159225086,  0.23611413753036, -0.66852578425466, ],
                [ -0.70903152581832,  0.23766935833457, -0.66391834701609
                    ]])


    Npts = np.prod(Q.shape)//4
    M = Q2M(Q)
    M = np.transpose(M, [2,0,1])


    dir_A_LOS_cel = np.sum(M*np.tile(dir_A_LOS[np.newaxis, np.newaxis,:], (Npts,3,1)),axis=2)
    dir_B_LOS_cel = np.sum(M*np.tile(dir_B_LOS[np.newaxis, np.newaxis,:], (Npts,3,1)),axis=2)

    dir_A_LOS_gal = coord_trans(dir_A_LOS_cel, 'C', 'G')
    Cll_A = np.array(hp.vec2ang(dir_A_LOS_gal.T, lonlat=True))

    dir_B_LOS_gal = coord_trans(dir_B_LOS_cel, 'C', 'G')
    Cll_B = np.array(hp.vec2ang(dir_B_LOS_gal.T, lonlat=True))

    
    gal = np.vstack((Cll_A, Cll_B)).T


    dir_A_POL_cel = []
    dir_B_POL_cel = []
    for i in range(len(dir_A_POL)):
        dir_A_POL_cel.append(np.sum(M*np.tile(dir_A_POL[i, np.newaxis, np.newaxis,:], (Npts,3,1)),axis=2))
        dir_B_POL_cel.append(np.sum(M*np.tile(dir_B_POL[i, np.newaxis, np.newaxis,:], (Npts,3,1)),axis=2))
    dir_A_POL_cel = np.array(dir_A_POL_cel)
    dir_B_POL_cel = np.array(dir_B_POL_cel)

    dir_A_POL_gal = np.zeros_like(dir_A_POL_cel)
    dir_B_POL_gal = np.zeros_like(dir_A_POL_cel)
    pol = []
    for i in range(len(dir_A_POL)):
        dir_A_POL_gal[i] = coord_trans(dir_A_POL_cel[i], 'C', 'G')
        Pll_A = np.array(hp.vec2ang(dir_A_POL_gal[i].T, lonlat=True))

        dir_B_POL_gal[i] = coord_trans(dir_B_POL_cel[i], 'C', 'G')
        Pll_B = np.array(hp.vec2ang(dir_B_POL_gal[i].T, lonlat=True))
        pol.append(np.vstack((Pll_A, Pll_B)).T)
    pol = np.array(pol)


    return gal, pol


def main():
    file_in = 'wmap.fits'
    file_out = 'h5_wmap_test.h5'

    data = fits.open(file_in)
    f = h5py.File(file_out, 'w')

    obsid = '00001'
    band_labels = ['K113', 'K114', 'K123', 'K124']
    
    
    time = data[2].data['TIME']
    for label in band_labels:
        TOD = data[2].data[label]
        # Write to h5 file.
        dset = f.create_dataset(obsid + '/' + label + '/TOD', data=TOD)
   
    n_records = data[0].header['NUMREC']
    
    # position (and velocity) in km(/s) in Sun-centered coordinates
    pos = data[1].data['POSITION']
    vel = data[1].data['VELOCITY']
    
    
    gal, pol = quat_to_sky_coords(data[1].data['QUATERN'])

    sin_2_gamma = np.zeros((len(gal), 2, 10))
    cos_2_gamma = np.zeros((len(gal), 2, 10))
    for tind in range(len(gal)):
        for hornind in range(2):
            sin_2_gi, cos_2_gi = gamma_from_pol(gal[tind,2*hornind:2*hornind+2], pol[:,tind,2*hornind:2*hornind+2])
            sin_2_gamma[tind, hornind,:] = sin_2_gi
            cos_2_gamma[tind, hornind,:] = cos_2_gi


    plt.figure()
    plt.plot(time, np.radians(gal[:,0]), '.')
    plt.plot(time, sin_2_gamma[:,0,0], '.')
    plt.plot(time, cos_2_gamma[:,0,0], '.')

    plt.figure()
    plt.scatter(gal[:,0], gal[:,1], c=time, s=1)
    plt.colorbar()

    plt.figure()
    plt.scatter(cos_2_gamma[:,0,0], sin_2_gamma[:,0,0], c=time, s=1)
    plt.colorbar()

    plt.show()



    f.create_dataset(obsid + '/common/gal', data=gal)
    f.create_dataset(obsid + '/common/vel', data=vel)
    f.create_dataset(obsid + '/common/time', data=time)
    f.create_dataset(obsid + '/common/sin_2_g', data=sin_2_gamma)
    f.create_dataset(obsid + '/common/cos_2_g', data=cos_2_gamma)


    f.close()

if __name__ == '__main__':
    main()
