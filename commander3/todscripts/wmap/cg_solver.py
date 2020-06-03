import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import h5py

import huffman


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

def check_hdf5(nside=256):
    npix = hp.nside2npix(nside)
    b = np.zeros(hp.nside2npix(nside))
    M_diag = np.zeros(npix)

    labels = ['K113', 'K114', 'K123', 'K124']
    sigmas = []
    from glob import glob
    fnames = glob('/mn/stornext/d16/cmbco/bp/wmap/data/wmap_K1_*v6.h5')
    fnames.sort()
    for fname in fnames:
        print(fname)
        f= h5py.File(fname, 'r')
        obsid = str(list(f.keys())[0])
        TOD0 = np.array(f[obsid + '/' + labels[0] + '/tod'])
        if len(TOD0) != 675000:
            continue
        
        
        DAs = [[], [], [], []]
        pixAs = []
        pixBs = []
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
            M_diag[pixA[t]] += sigma0**-2
            M_diag[pixB[t]] += sigma0**-2

    #hp.mollview(b, cmap='coolwarm', min=-100, max=100)
    #hp.mollview(M_diag, norm='hist')



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
    i_max = 10000
    eps = 1e-7
    while ((i < i_max) & (delta_new > eps**2*delta_0)):
        print(delta_new)
        q[p] = d[p]/M_diag[p]
        alpha = delta_new/d.dot(q)
        x = x + alpha*d
        if i % 50 == 0:
            r[p] = b[p] - x[p]/M_diag[p]
        else:
            r = r - alpha*q
        s[p] = r[p]/M_diag[p]
        delta_old = np.copy(delta_new)
        delta_new = r.dot(s)
        beta = delta_new/delta_old
        d = s + beta*d
        i += 1

    #hp.mollview(x, norm='hist')

    hp.write_map('cg_v1.fits', x, overwrite=True)

    #plt.show()

    

    return


if __name__ == '__main__':
    cg_test()
    check_hdf5()
