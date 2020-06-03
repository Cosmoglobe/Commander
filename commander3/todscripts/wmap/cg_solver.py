import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import h5py

import huffman

from multiprocessing import Pool


def cg_solve(A, b, Minv, imax=1000, eps=1e-6):
    x = np.zeros_like(b)
    i = 0
    r = b - A.dot(x)
    d = Minv.dot(r)
    delta_new = r.dot(d)
    delta_0 = r.dot(d)
    while ((i < imax) & (delta_new > eps**2*delta_0)):
        q = A.dot(d)
        alpha = delta_new/d.dot(q)
        x = x + alpha*d
        if (i % 50 == 0):
            r = b - A.dot(x)
        else:
            r = r - alpha*q
        s = Minv.dot(r)
        delta_old = np.copy(delta_new)
        delta_new = r.dot(s)
        beta = delta_new/delta_old
        d = s + beta*d
        i += 1
    return x

def cg_solve_map(A, b, Minv, imax=1000, eps=1e-6):
    '''
    Try to rewrite this so that A is an operator, not a matrix that needs to be
    held in memory.
    '''
    x = np.zeros_like(b)
    i = 0
    r = b - A.dot(x)
    d = Minv.dot(r)
    delta_new = r.dot(d)
    delta_0 = r.dot(d)
    while ((i < imax) & (delta_new > eps**2*delta_0)):
        q = A.dot(d)
        alpha = delta_new/d.dot(q)
        x = x + alpha*d
        if (i % 50 == 0):
            r = b - A.dot(x)
        else:
            r = r - alpha*q
        s = Minv.dot(r)
        delta_old = np.copy(delta_new)
        delta_new = r.dot(s)
        beta = delta_new/delta_old
        d = s + beta*d
        i += 1
    return x


def cg_test():
    A = np.array([[3,2],
                  [2,6]])
    b = np.array([2,-8])
    Minv = np.eye(2)

    x = cg_solve(A, b, Minv)
    assert np.allclose(A.dot(x), b), 'CG solution is not close enough'
    return

def get_Mb(fname, nside=256):
    npix = hp.nside2npix(nside)
    M = np.zeros(npix)
    b = np.zeros(npix)
    labels = ['K113', 'K114', 'K123', 'K124']
    f= h5py.File(fname, 'r')
    obsid = str(list(f.keys())[0])
    TOD0 = np.array(f[obsid + '/' + labels[0] + '/tod'])
    if len(TOD0) != 675000:
        print(f'{fname} has wrong length')
        return M, b
    
    
    DAs = [[], [], [], []]
    pixAs = []
    pixBs = []
    sigmas = []
    for num, label in enumerate(labels):
        TODs = np.array(f[obsid + '/' + label + '/tod'])
        gain, sigma0 = f[obsid + '/' + label + '/scalars'][0:2]
        DAs[num] = DAs[num] + TODs.tolist()
        sigmas.append(sigma0)
        if label == 'K113':
            pixA = np.array(f[obsid + '/' + label + '/pixA']).tolist()
            pixB = np.array(f[obsid + '/' + label + '/pixB']).tolist()

    DAs = np.array(DAs)
    sigma0 = sum(np.array(sigmas)**2)**0.5

    
    
    d1 = 0.5*(DAs[0] + DAs[1])
    d2 = 0.5*(DAs[2] + DAs[3])
    
    d = 0.5*(d1 + d2) # = i_A - i_B
    p = 0.5*(d1 - d2) # = q_A*cos(2*g_A) + u_A*sin(2*g_A) - q_B*cos(2*g_B) - u_B*sin(2*g_B)



    for t in range(len(d)):
        b[pixA[t]] += d[t]/sigma0
        b[pixB[t]] -= d[t]/sigma0


    for t in range(len(d)):
        M[pixA[t]] += sigma0**-2
        M[pixB[t]] += sigma0**-2

    return M, b

def get_cg_K(nside=256):
    npix = hp.nside2npix(nside)
    b = np.zeros(hp.nside2npix(nside))
    M_diag = np.zeros(npix)

    from glob import glob
    fnames = glob('/mn/stornext/d16/cmbco/bp/wmap/data/wmap_K1_*v6.h5')
    fnames.sort()
    #fnames = fnames[:1000]
    pool = Pool(processes=120)
    x = [pool.apply_async(get_Mb, args=[fname]) for fname in fnames]
    for res in x:
        M_i, b_i = res.get()
        M_diag += M_i
        b += b_i
    print('Loaded data')

    # I would like to run an asynchronous pool that will return a sum of all of
    # the arguments.




    i = 0
    x = np.zeros_like(b)
    d = np.zeros_like(b)
    q = np.zeros_like(b)
    s = np.zeros_like(b)
    r = b
    p = (M_diag != 0)
    d[p] = (r/M_diag)[p]
    delta_new = r.dot(d)
    delta_0 = np.copy(delta_new)
    i_max = npix
    eps = 1e-7
    while ((i < i_max) & (delta_new > eps**2*delta_0)):
        q[p] = d[p]/M_diag[p]
        alpha = delta_new/d.dot(q)
        x = x + alpha*d
        if i % 50 == 0:
            r[p] = b[p] - x[p]/M_diag[p]
            #print(f'{str(i).zfill(5)}th iteration, log10(chi2){int(np.log10(delta_new))}')
        else:
            r = r - alpha*q
        if i % 100 == 0:
            # Whichever reaches 1 first wins. The second term seems to go
            # ~cubic.
            print(np.round(i/i_max,3),\
                    np.round(1/np.log10(delta_new/(delta_0*eps**2)),3),\
                    int(delta_new))
        s[p] = r[p]/M_diag[p]
        delta_old = np.copy(delta_new)
        delta_new = r.dot(s)
        beta = delta_new/delta_old
        d = s + beta*d
        i += 1

    hp.mollview(b, cmap='coolwarm', norm='hist')
    hp.mollview(M_diag, norm='hist')
    hp.mollview(x, norm='hist', title='Solution')
    hp.mollview(x, min=-250, max=250, title='Solution')

    hp.write_map('cg_v2.fits', x, overwrite=True)

    plt.show()

    

    return


if __name__ == '__main__':
    #cg_test()
    #get_cg_K()
    check_hdf5()


    # Take official W-band K-band, scan it with the same pointing matrix, divide
    # by gain, subtract from timestream, check the gain and pointing solution
    # directly. Should just be white noise.
