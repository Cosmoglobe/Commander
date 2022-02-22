import numpy as np
import ducc0.totalconvolve as totalconvolve
from tqdm import tqdm
import healpy as hp
import scipy.sparse
import scipy.linalg

import matplotlib.pyplot as plt

import sys
sys.path.append('/mn/stornext/u3/duncanwa/Commander/commander3/python/')
from commander_tools.tod_tools import commander_tod as tod
# Takes the directory where .h5 files are located and the dataset in question as
# arguments
comm_tod = tod.commander_tod('/mn/stornext/d16/cmbco/bp/wmap/data', 'wmap')


def dot_prod(b_map, pixA, pixB, res_A, res_B, psiA, psiB, flags, npnt,
    pmask, x1, x2):
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

    d1 = (1+x1)*res_A - (1-x1)*res_B
    d2 = (1+x2)*res_A - (1-x2)*res_B
    d = 0.5*(d1+d2)
    p = 0.5*(d1-d2)

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

    b_vec += P.T.dot(d)

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

def make_M(M, Prec, pixA, pixB, psiA, psiB, flags, pmask, x1, x2):
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


    P_A = scipy.sparse.csr_matrix((data_A, (t, pixA)), shape=(len(t)//4, 4*npix))
    P_B = scipy.sparse.csr_matrix((data_B, (t, pixB)), shape=(len(t)//4, 4*npix))
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

    P_A = scipy.sparse.csr_matrix((data_A, (t, pixA)), shape=(len(t)//4, 4*npix))
    P_B = scipy.sparse.csr_matrix((data_B, (t, pixB)), shape=(len(t)//4, 4*npix))
    P = P_A - P_B

    data_A = np.concatenate((f_A*dx*T,  f_A*(1+x)*QA, f_A*(1+x)*UA, f_A*(1+x)*SA))
    data_B = np.concatenate((f_B*-dx*T, f_B*(1-x)*QB, f_B*(1-x)*UB, f_B*(1-x)*SB))
    P_A = scipy.sparse.csr_matrix((data_A, (t, pixA)), shape=(len(t)//4, 4*npix))
    P_B = scipy.sparse.csr_matrix((data_B, (t, pixB)), shape=(len(t)//4, 4*npix))
    P_am = P_A - P_B

    M += P_am.T.dot(P)
    Prec += P_am.T.dot(P_am)


    return M, Prec

def accumulate(tod_ind, x1, x2, pmask):
    '''
    For given week of mission, tod_ind, get timestream given pointing and
    polarization angle for horns A and B, then accumulate this week's
    contributions to M_i = P_{am}^T N^{-1} P and b_i = P_{am}^T N^{-1} d.
    '''
    ind = str(tod_ind).zfill(6)
    comm_tod.init_file('Q1', ind)
    
    pixA_i = comm_tod.load_field(f'{ind}/Q1/pixA').astype('int')
    pixB_i = comm_tod.load_field(f'{ind}/Q1/pixB').astype('int')
    psiA_i = comm_tod.load_field(f'{ind}/Q1/psiA')
    psiB_i = comm_tod.load_field(f'{ind}/Q1/psiB')
    flags = comm_tod.load_field(f'{ind}/Q1/flag')

    nside_in = 512
    thetaA, phiA = hp.pix2ang(nside_in, pixA_i)
    thetaB, phiB = hp.pix2ang(nside_in, pixB_i)
    pixA = hp.ang2pix(nside_out, thetaA, phiA)
    pixB = hp.ang2pix(nside_out, thetaB, phiB)

    
    # build pointings
    npnt = len(thetaA)
    ptg = np.zeros((npnt, 3))
    ptg[:, 0] = thetaA          # theta
    ptg[:, 1] = phiA            # phi
    ptg[:, 2] = psiA_i          # psi
    res_A = interp_A.interpol(ptg)[0]
    
    ptg[:, 0] = thetaB          # theta
    ptg[:, 1] = phiB            # phi
    ptg[:, 2] = psiB_i          # psi
    res_B = interp_B.interpol(ptg)[0]

    M = scipy.sparse.csr_matrix((4*npix, 4*npix))
    Prec = scipy.sparse.csr_matrix((4*npix, 4*npix))
    b_map = np.zeros((4, npix))
    M_i, Prec_i = make_M(M, Prec, pixA, pixB, psiA_i, psiB_i, flags, pmask, x1, x2)
    b_i = dot_prod(b_map, pixA, pixB, res_A, res_B, psiA_i, psiB_i, flags,
        npnt, pmask, x1, x2)

    return M_i, b_i, Prec_i


def make_dipole_alms(amp=3355, l=263.99, b=48.26, lmax=128):
    ipix = np.arange(12*512**2)
    x, y, z = hp.pix2vec(512, ipix)

    theta, phi = np.pi/180*l, np.pi/2 - np.pi/180*b
    amps = amp*np.array([np.cos(theta)*np.sin(phi), np.sin(theta)*np.sin(phi), np.cos(phi)])

    dipole = x*amps[0] + y*amps[1] + z*amps[2]
    slm = hp.map2alm(dipole, lmax=lmax)

    slm = slm[np.newaxis,:].astype('complex128')
    return slm


def get_sidelobe_alms(band='Q1', lmax=128, kmax=100):
    # LOS geometry extracted from program.pars
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

    bands = np.array(['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4'])
    SIDELOBE_DIR = '/mn/stornext/d16/cmbco/ola/wmap/ancillary_data/far_sidelobe_maps'
    # Construct sidelobe model
    sidelobe = hp.read_map(f'{SIDELOBE_DIR}/wmap_sidelobe_map_{band}_9yr_v5.fits')
    
    # Normalized such that \int B_A d\Omega = 1, converting from 
    # sum(abs(sidelobe)) = 2*N_pix normalization
    beam_A = hp.reorder(sidelobe, n2r=True)/(4*np.pi)
    beam_A[beam_A < 0] = 0
    beam_B = hp.reorder(sidelobe, n2r=True)/(4*np.pi)
    beam_B[beam_B > 0] = 0
    beam_B = -beam_B
    
    
    dir_A = dir_A_los[band == bands][0]
    theta = np.arccos(dir_A[2])
    phi = np.arctan2(dir_A[1], dir_A[0])
    
    # Rotate so that main beam A is pointing in the z-direction
    r = hp.rotator.Rotator(rot=(phi, -theta, 0), \
        deg=False, eulertype='ZYX')
    beam_A = r.rotate_map_pixel(beam_A)
    
    dir_B = dir_B_los[band == bands][0]
    theta = np.arccos(dir_B[2])
    phi = np.arctan2(dir_B[1], dir_B[0])
    
    # Rotate so that main beam B is pointing in the z-direction
    r = hp.rotator.Rotator(rot=(phi, -theta, 0), \
        deg=False, eulertype='ZYX')
    beam_B = r.rotate_map_pixel(beam_B)


    blm_A = hp.map2alm(beam_A, lmax=lmax, mmax=kmax)
    blm_B = hp.map2alm(beam_B, lmax=lmax, mmax=kmax)

    blm_A = blm_A[np.newaxis,:].astype('complex128')
    blm_B = blm_B[np.newaxis,:].astype('complex128')

    return blm_A, blm_B



if __name__ == '__main__':


    # To what extent is this effect due to the mapmaking algorithm, and to what
    # extent is it due to the imbalance parameters in the data model itself?

    MASK_DIR = '/mn/stornext/d16/cmbco/ola/wmap/ancillary_data/masks'
    bands = np.array(['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4'])
    bands = [bands[0]]
    for band in bands:
        print(band)
        # Sets maximum lmax, mmax for sidelobe convolution
        lmax = 128
        kmax = 100
   
        # Signal and sidelobe alm model
        slm = make_dipole_alms(lmax=lmax)
        blm_A, blm_B = get_sidelobe_alms(band=band, lmax=lmax, kmax=kmax)
   
        # totalconvolver interpolator, grid in theta,phi,psi
        interp_A = totalconvolve.Interpolator(
            slm, blm_A, separate=False, lmax=lmax, kmax=kmax, epsilon=1e-4, nthreads=0)
        interp_B = totalconvolve.Interpolator(
            slm, blm_B, separate=False, lmax=lmax, kmax=kmax, epsilon=1e-4, nthreads=0)
        
        
   

        # Initializing map structures
        nside_out = 16
        npix = hp.nside2npix(nside_out)
        M = scipy.sparse.csr_matrix((4*npix, 4*npix))
        Prec = scipy.sparse.csr_matrix((4*npix, 4*npix))
        b_map = np.zeros((4, 12*nside_out**2))
   

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
  
        x1, x2 = x_ims[band]
        # Valid weeks from 1--468
        tod_inds = np.arange(1, 468+1)
        
        import multiprocessing
        from functools import partial
        # Maximum number of cpus before my node throws an error
        ncpus = 26
        pool = multiprocessing.Pool(processes=ncpus)
        pool_outputs = list(tqdm(pool.imap(
                           partial(accumulate, x1=x1, x2=x2, pmask=pmask),
                           tod_inds), total=len(tod_inds)))
        pool.close()
        pool.join()
        for i in range(len(pool_outputs)):
            M += pool_outputs[i][0]
            b_map += pool_outputs[i][1]
            Prec += pool_outputs[i][2]


       
        # Useful for visualizing poorly measured modes
        scipy.sparse.save_npz(f'M_{band}.npz', M)
        scipy.sparse.save_npz(f'Precond_{band}.npz', Prec)
        # M = scipy.sparse.load_npz('M.npz')
        
        hp.write_map(f'b_{band}.fits', b_map, overwrite=True)
        
        b = np.concatenate((b_map[0], b_map[1], b_map[2], b_map[3]))
       
        print('Solving Mx=b')
        x = scipy.sparse.linalg.spsolve(M, b)
        print('Solved Mx=b')
        
        I,Q,U,S = np.split(x, 4)
        
        m = np.array([I,Q,U,S])
        hp.write_map(f'x_{band}.fits', m, overwrite=True)
