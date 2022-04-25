import numpy as np
from tqdm import tqdm
import healpy as hp
import scipy.sparse
import scipy.linalg

import h5py

import matplotlib.pyplot as plt

import sys
sys.path.append('/mn/stornext/u3/duncanwa/Commander/commander3/python/')
from commander_tools.tod_tools import commander_tod as tod
# Takes the directory where .h5 files are located and the dataset in question as
# arguments
comm_tod = tod.commander_tod('/mn/stornext/d16/cmbco/bp/wmap/data_precal', 'wmap')


from scipy.interpolate import interp1d


# kind can be previous, cubic?
def orb_vel(kind='previous'):
    tod_inds = np.arange(1, 468+1)

    v_suns = []

    for tod_ind in tod_inds:
        ind = str(tod_ind).zfill(6)
        comm_tod.init_file('Q1', ind)
        v_sun = comm_tod.load_field(f'{ind}/common/vsun')[()]
        v_suns.append(v_sun)

    v_suns = np.array(v_suns)

    #fig, axes = plt.subplots(nrows=3, sharex=True)
    #axes[0].plot(tod_inds, v_suns[:,0], '.')
    #axes[1].plot(tod_inds, v_suns[:,1], '.')
    #axes[2].plot(tod_inds, v_suns[:,2], '.')



    fx = interp1d(tod_inds, v_suns[:,0], fill_value='extrapolate',
        assume_sorted=True, kind=kind)
    fy = interp1d(tod_inds, v_suns[:,1], fill_value='extrapolate',
        assume_sorted=True, kind=kind)
    fz = interp1d(tod_inds, v_suns[:,2], fill_value='extrapolate',
        assume_sorted=True, kind=kind)

    #t = np.linspace(1, 469, 1000)
    #axes[0].plot(t, fx(t), 'k', lw=1)
    #axes[1].plot(t, fy(t), 'k', lw=1)
    #axes[2].plot(t, fz(t), 'k', lw=1)
    #plt.show()


    return fx, fy, fz


def dot_prod(b_map, pixA, pixB, d, p, psiA, psiB, flags,
    pmask, x1, x2, sigma):
    '''
    Creates synthetic timestream for differential horns,
    d1 = (1+x1) * [T_A + Q_A\cos2\gamma_A + U_A\sin2\gamma_A + S_A]
       - (1-x1) * [T_B + Q_B\cos2\gamma_B + U_B\sin2\gamma_B + S_B]
    d2 = (1+x2) * [T_A - Q_A\cos2\gamma_A - U_A\sin2\gamma_A - S_A]
       - (1-x2) * [T_B - Q_B\cos2\gamma_B - U_B\sin2\gamma_B - S_B]

    and then creates binned map 
    b = P_{am}^T N^-1 d
    for the polarized and unpolarized timesreams
    d = (d1+d2)/2
    and
    p = (d1-d2)/2
    and assuming uniform noise, i.e., N^-1 = diag(1)
    '''

    x = (x1+x2)/2
    dx = (x1-x2)/2




    inds = ((flags % 2) == 0)
    if ~np.any(inds):
        return M

    d = d[inds]
    p = p[inds]
    pixA = pixA[inds]
    pixB = pixB[inds]
    psiA = psiA[inds]
    psiB = psiB[inds]

    f_B = pmask[pixA]
    f_A = pmask[pixB]

    npix = 12*nside_out**2
    
    t_i = np.arange(len(pixA))
    T  = np.ones_like(t_i)
    QA = np.cos(2*psiA)
    UA = np.sin(2*psiA)
    QB = np.cos(2*psiB)
    UB = np.sin(2*psiB)
    SA = T
    SB = T

    t = np.concatenate((t_i,t_i,t_i,t_i))


    data_A = np.concatenate(((1+x)*T*f_A,  dx*QA*f_A,  dx*UA*f_A,  dx*SA*f_A))
    data_B = np.concatenate(((1-x)*T*f_B, -dx*QB*f_B, -dx*UB*f_B, -dx*SB*f_B))
    pixA = np.concatenate((pixA, pixA+npix, pixA+2*npix, pixA+3*npix))
    pixB = np.concatenate((pixB, pixB+npix, pixB+2*npix, pixB+3*npix))


    P_A = scipy.sparse.csr_matrix((data_A, (t, pixA)), shape=(len(t)//4, 4*npix))
    P_B = scipy.sparse.csr_matrix((data_B, (t, pixB)), shape=(len(t)//4, 4*npix))
    P = P_A - P_B

    b_vec = np.zeros(b_map.size)

    b_vec += P.T.dot(d/sigma**2)

    data_A = np.concatenate(( dx*T*f_A, (1+x)*QA*f_A, (1+x)*UA*f_A, (1+x)*SA*f_A))
    data_B = np.concatenate((-dx*T*f_B, (1-x)*QB*f_B, (1-x)*UB*f_B, (1-x)*SB*f_B))

    P_A = scipy.sparse.csr_matrix((data_A, (t, pixA)), shape=(len(t)//4, 4*npix))
    P_B = scipy.sparse.csr_matrix((data_B, (t, pixB)), shape=(len(t)//4, 4*npix))
    P = P_A - P_B
    b_vec += P.T.dot(p)

    b_mapi = np.split(b_vec,4)
    for i in range(4):
      b_map[i] += b_mapi[i]
    return b_map

def make_M(M, Prec, pixA, pixB, psiA, psiB, flags, pmask, x1, x2, sigma):
    '''
    Constructs the asymmetric mapmaking matrix
    M = P_{am}^T N^-1 P
    '''
    inds = ((flags % 2) == 0)
    if ~np.any(inds):
        return M
    pixA = pixA[inds]
    pixB = pixB[inds]
    psiA = psiA[inds]
    psiB = psiB[inds]

    f_B = pmask[pixA]
    f_A = pmask[pixB]
    npix = 12*nside_out**2

    x = (x1+x2)/2
    dx = (x1-x2)/2

    t_i = np.arange(len(pixA))
    T = np.ones_like(t_i)
    QA = np.cos(2*psiA)
    UA = np.sin(2*psiA)
    QB = np.cos(2*psiB)
    UB = np.sin(2*psiB)
    SA = T
    SB = T

    t = np.concatenate((t_i,t_i,t_i,t_i))
    del t_i


    data_A = np.concatenate(((1+x)*T,  dx*QA,  dx*UA,  dx*SA))
    data_B = np.concatenate(((1-x)*T, -dx*QB, -dx*UB, -dx*SB))
    pixA = np.concatenate((pixA, pixA+npix, pixA+2*npix, pixA+3*npix))
    pixB = np.concatenate((pixB, pixB+npix, pixB+2*npix, pixB+3*npix))


    P_A = scipy.sparse.csr_matrix((data_A/sigma**2, (t, pixA)), shape=(len(t)//4, 4*npix))
    P_B = scipy.sparse.csr_matrix((data_B/sigma**2, (t, pixB)), shape=(len(t)//4, 4*npix))
    P = P_A - P_B

    data_A = np.concatenate((f_A*(1+x)*T, f_A*dx*QA,  f_A*dx*UA,  f_A*dx*SA))
    data_B = np.concatenate((f_B*(1-x)*T, f_B*-dx*QB, f_B*-dx*UB, f_B*-dx*SB))
    P_A = scipy.sparse.csr_matrix((data_A, (t, pixA)), shape=(len(t)//4, 4*npix))
    P_B = scipy.sparse.csr_matrix((data_B, (t, pixB)), shape=(len(t)//4, 4*npix))
    P_am = P_A - P_B

    M += P_am.T.dot(P)
    Prec += P_am.T.dot(P_am)


    data_A = np.concatenate(( dx*T, (1+x)*QA, (1+x)*UA, (1+x)*SA))
    data_B = np.concatenate((-dx*T, (1-x)*QB, (1-x)*UB, (1-x)*SB))

    P_A = scipy.sparse.csr_matrix((data_A/sigma**2, (t, pixA)), shape=(len(t)//4, 4*npix))
    P_B = scipy.sparse.csr_matrix((data_B/sigma**2, (t, pixB)), shape=(len(t)//4, 4*npix))
    P = P_A - P_B

    data_A = np.concatenate((f_A*dx*T,  f_A*(1+x)*QA, f_A*(1+x)*UA, f_A*(1+x)*SA))
    data_B = np.concatenate((f_B*-dx*T, f_B*(1-x)*QB, f_B*(1-x)*UB, f_B*(1-x)*SB))
    P_A = scipy.sparse.csr_matrix((data_A, (t, pixA)), shape=(len(t)//4, 4*npix))
    P_B = scipy.sparse.csr_matrix((data_B, (t, pixB)), shape=(len(t)//4, 4*npix))
    P_am = P_A - P_B

    M += P_am.T.dot(P)
    Prec += P_am.T.dot(P_am)


    return M, Prec

def accumulate(tod_ind, x1, x2, x1hat, x2hat, pmask, band):
    '''
    For given week of mission, tod_ind, get timestream given pointing and
    polarization angle for horns A and B, then accumulate this week's
    contributions to M_i = P_{am}^T N^{-1} P and b_i = P_{am}^T N^{-1} d.
    '''
    ind = str(tod_ind).zfill(6)
    comm_tod.init_file('Q1', ind)
    
    pixA_i = comm_tod.load_field(f'{ind}/{band}/pixA').astype('int')
    pixB_i = comm_tod.load_field(f'{ind}/{band}/pixB').astype('int')
    psiA_i = comm_tod.load_field(f'{ind}/{band}/psiA')
    psiB_i = comm_tod.load_field(f'{ind}/{band}/psiB')
    flags = comm_tod.load_field(f'{ind}/{band}/flag')
    d13 = comm_tod.load_field(f'{ind}/{band}13/tod')[()]
    d14 = comm_tod.load_field(f'{ind}/{band}14/tod')[()]
    d23 = comm_tod.load_field(f'{ind}/{band}23/tod')[()]
    d24 = comm_tod.load_field(f'{ind}/{band}24/tod')[()]

    v_sun = comm_tod.load_field(f'{ind}/common/vsun')[()]

    fx, fy, fz = orb_vel()
    t = tod_ind + np.arange(len(pixA_i))/len(pixA_i)

    v_sun = np.array([fx(t), fy(t), fz(t)])

    #print(tod_ind, v_sun)
    c         = 2.99792458e8
    k_B       = 1.3806503e-23
    h         = 1.0545726691251021e-34 * 2*np.pi
    T_CMB     = 2.72548

    beta = v_sun/c
    b = np.sqrt(np.sum(beta**2))
    nu = 61e9
    x = h * nu/(k_B * T_CMB)
    q = (x/2)*(np.exp(x)+1)/(np.exp(x) -1)


    P_A = np.array(hp.pix2vec(512, pixA_i))
    P_B = np.array(hp.pix2vec(512, pixB_i))

    b_dot_A = (beta*P_A).sum(axis=0)
    b_dot_B = (beta*P_B).sum(axis=0)

    s_dip_A = 1e3*T_CMB * (b_dot_A + q*(b_dot_A**2 - b**2/3.))
    s_dip_B = 1e3*T_CMB * (b_dot_B + q*(b_dot_B**2 - b**2/3.))


    dip_d1 = (1 + x1)*s_dip_A - (1-x1)*s_dip_B
    dip_d2 = (1 + x2)*s_dip_A - (1-x2)*s_dip_B

    dip_d = (dip_d1 + dip_d2)/2
    dip_p = (dip_d1 - dip_d2)/2

    #plt.plot(s_dip_A, 'k')
    #plt.plot(s_dip_B, 'r')
    #plt.show()
    #asdfsd


    # b = sqrt(sum(v_ref**2))/c   ! beta for the given scan


    DIR='/mn/stornext/d16/cmbco/bp/wmap/data_sim/'

    with h5py.File(f'{DIR}/wmap_{band}_{ind}.h5', 'r') as data:
        sigma113 = data[f'{ind}/{band}13/scalars'][1]
        sigma114 = data[f'{ind}/{band}14/scalars'][1]
        sigma123 = data[f'{ind}/{band}23/scalars'][1]
        sigma124 = data[f'{ind}/{band}24/scalars'][1]

    sigmas = np.array([sigma113, sigma114, sigma123, sigma124])

    nside_in = 512
    thetaA, phiA = hp.pix2ang(nside_in, pixA_i)
    thetaB, phiB = hp.pix2ang(nside_in, pixB_i)
    pixA = hp.ang2pix(nside_out, thetaA, phiA)
    pixB = hp.ang2pix(nside_out, thetaB, phiB)

   
    '''
    # build pointings
    npnt = len(thetaA)
    #ptg = np.zeros((npnt, 3))
    #ptg[:, 0] = thetaA          # theta
    #ptg[:, 1] = phiA            # phi
    #ptg[:, 2] = psiA_i          # psi
    res_A = I[pixA] + Q[pixA]*np.cos(2*psiA_i) + U[pixA]*np.sin(2*psiA_i)
    
    #ptg[:, 0] = thetaB          # theta
    #ptg[:, 1] = phiB            # phi
    #ptg[:, 2] = psiB_i          # psi
    res_B = I[pixB] + Q[pixB]*np.cos(2*psiB_i) + U[pixB]*np.sin(2*psiB_i)
    d1 = (1+x1)*res_A - (1-x1)*res_B

    res_A = I[pixA] - Q[pixA]*np.cos(2*psiA_i) - U[pixA]*np.sin(2*psiA_i)
    res_B = I[pixB] - Q[pixB]*np.cos(2*psiB_i) - U[pixB]*np.sin(2*psiB_i)
    d2 = (1+x2)*res_A - (1-x2)*res_B
    '''

    d1 = (d13 + d14)/2
    d2 = (d23 + d24)/2
    del d13, d14, d23, d24

    d = 0.5*(d1+d2) - dip_d
    p = 0.5*(d1-d2) - dip_p
    del d1, d2

    d = dip_d
    p = dip_p

    sigma = np.sqrt(sum(sigmas**2)/4)

    M = scipy.sparse.csr_matrix((4*npix, 4*npix))
    Prec = scipy.sparse.csr_matrix((4*npix, 4*npix))
    b_map = np.zeros((4, npix))
    M_i, Prec_i = make_M(M, Prec, pixA, pixB, psiA_i, psiB_i, flags, pmask,
        x1hat, x2hat, sigma)
    b_i = dot_prod(b_map, pixA, pixB, d, p, psiA_i, psiB_i, flags,
        pmask, x1hat, x2hat, sigma)

    return M_i, b_i, Prec_i


def make_dipole_alms(amp=3355, l=263.99, b=48.26, lmax=128):
    ipix = np.arange(12*512**2)
    x, y, z = hp.pix2vec(512, ipix)

    theta, phi = np.pi/180*l, np.pi/2 - np.pi/180*b
    amps = amp*np.array([np.cos(theta)*np.sin(phi), np.sin(theta)*np.sin(phi), np.cos(phi)])

    dipole = x*amps[0] + y*amps[1] + z*amps[2]
    m = hp.read_map('/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_Q1_v5.fits',
        field=(0,1,2))*1e3
    m *= 0
    m[0] = + dipole


    #m = hp.read_map('cmb_plus_dip.fits', field=(0,1,2))

    #m *= 0

    #m[0] += 1000

    slm = hp.map2alm(m, lmax=lmax)

    m = hp.alm2map(slm, 16, lmax=lmax)
    hp.write_map(f'x_input.fits', m, overwrite=True)

    return slm





if __name__ == '__main__':


    # To what extent is this effect due to the mapmaking algorithm, and to what
    # extent is it due to the imbalance parameters in the data model itself?


    MASK_DIR = '/mn/stornext/d16/cmbco/ola/wmap/ancillary_data/masks'
    bands = np.array(['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4'])
    band = bands[2]

    nside_out = 16
    lmax = 128

    # Signal and sidelobe alm model
    slm = make_dipole_alms(lmax=lmax)
    s = hp.alm2map(slm, nside_out)

    # Initializing map structures
    npix = hp.nside2npix(nside_out)
   

    # Low-resolution processing map
    pmask = hp.read_map(f'{MASK_DIR}/wmap_processing_mask_{band[:-1]}_r4_9yr_v5.fits')
    pmask = hp.ud_grade(pmask, nside_out).astype('int')
    
    
   
    # Imbalance parameters from Bennett et al. (2013)
    x_ims = {'K1':  (-0.00067, 0.00536),
             'Ka1': ( 0.00353, 0.00154),
             'Q1':  (-0.00013, 0.00414),
             'Q2':  ( 0.00756, 0.00986),
             'V1':  ( 0.00053, 0.00250),
             'V2':  ( 0.00352, 0.00245),
             'W1':  ( 0.01134, 0.00173),
             'W2':  ( 0.01017, 0.01142),
             'W3':  (-0.00122, 0.00463),
             'W4':  ( 0.02311, 0.02054)}
    x_sigs = {'K1':  (0.00017, 0.00014),
              'Ka1': (0.00014, 0.00008),
              'Q1':  (0.00046, 0.00025),
              'Q2':  (0.00052, 0.00115),
              'V1':  (0.00020, 0.00057),
              'V2':  (0.00033, 0.00098),
              'W1':  (0.00199, 0.00036),
              'W2':  (0.00216, 0.00121),
              'W3':  (0.00062, 0.00041),
              'W4':  (0.00380, 0.00202)}
  
    x1, x2 = x_ims[band]
    sx1, sx2 = x_sigs[band]
    x1hat = x1 + sx1*np.random.randn()
    x2hat = x2 + sx2*np.random.randn()
    tod_inds = np.arange(1, 468+1)
    tod_inds = np.arange(1, 52+1)
    #tod_inds = np.arange(1, 2)

    x1, x2 = 0,0

    x1s = np.arange(x1-sx1, x1+2*sx1, round(sx1/2, 5))[:5]
    x2s = np.arange(x2-sx2, x2+2*sx2, round(sx2/2, 5))[:5]
    x1s = [x1s.mean()]
    x2s = [x2s.mean()]

    x1s = [0]
    x2s = [0.01]

    x1s = [-0.00013]
    x2s = [0.00414]
    # From WMAP 9year Q1 band
    print(x1s)
    print(x2s)
    import multiprocessing
    from functools import partial 
    for i in range(len(x1s)):
      for j in range(len(x2s)):
        M = scipy.sparse.csr_matrix((4*npix, 4*npix))
        Prec = scipy.sparse.csr_matrix((4*npix, 4*npix))
        b_map = np.zeros((4, 12*nside_out**2))
        x1hat = x1s[i]
        x2hat = x2s[j]
        #x1hat = 0.05
        #x2hat = 0.04
        # Maximum number of cpus before my node throws an error
        ncpus = 64
        pool = multiprocessing.Pool(processes=ncpus)
        pool_outputs = list(tqdm(pool.imap(
                           partial(accumulate, x1=x1, x2=x2, x1hat=x1hat,
                             x2hat=x2hat, pmask=pmask, band=band),
                           tod_inds), total=len(tod_inds)))
        pool.close()
        pool.join()
        for k in range(len(pool_outputs)):
            M += pool_outputs[k][0]
            b_map += pool_outputs[k][1]
            Prec += pool_outputs[k][2]


        
        # Useful for visualizing poorly measured modes
        # M = scipy.sparse.load_npz('M.npz')
        
        
        b = np.concatenate((b_map[0], b_map[1], b_map[2], b_map[3]))
        
        print('Solving Mx=b')
        #x = scipy.sparse.linalg.spsolve(M, b)
        #x = np.linalg.inv(np.array(M.todense())).dot(b)
        x = np.linalg.solve(np.array(M.todense()), b)
        print('Solved Mx=b')
        
        I,Q,U,S = np.split(x, 4)

        m = np.array([I,Q,U,S])
        hp.write_map(f'x_pencil_{x1hat}_{x2hat}.fits', m, overwrite=True)
