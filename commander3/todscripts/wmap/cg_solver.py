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


def cg_test():
    A = np.array([[3,2],
                  [2,6]])
    b = np.array([2,-8])
    Minv = np.eye(2)

    x = cg_solve(A, b, Minv)
    assert np.allclose(A.dot(x), b), 'CG solution is not close enough'
    return

def check_hdf5():
    fname = '/mn/stornext/d16/cmbco/bp/wmap/data/wmap_K1_000001_v4.h5'

    f= h5py.File(fname, 'r')
    obsid = str(list(f.keys())[0])
    
    print(obsid)
    
    DAs = []
    labels = ['K113', 'K114', 'K123', 'K124']
    for label in labels:
        TODs = np.array(f[obsid + '/' + label + '/tod'])
        gain = f[obsid + '/' + label + '/scalars'][0]
        DAs.append(TODs/gain)
        if label == 'K113':
            pixA = f[obsid + '/' + label + '/pixA']
            pixB = f[obsid + '/' + label + '/pixA']
    DAs = np.array(DAs)



    hufftreeA = f[obsid + '/common/hufftree_A']
    huffsymbA = f[obsid + '/common/huffsymb_A']
    hufftreeB = f[obsid + '/common/hufftree_B']
    huffsymbB = f[obsid + '/common/huffsymb_B']
    
    
    d1 = 0.5*(DAs[0] + DAs[1])
    d2 = 0.5*(DAs[2] + DAs[3])
    
    d = 0.5*(d1 + d2) # = i_A - i_B
    p = 0.5*(d1 - d2) # = q_A*cos(2*g_A) + u_A*sin(2*g_A) - q_B*cos(2*g_B) - u_B*sin(2*g_B)
    
    #plt.plot(d)
    #plt.show()

    return


if __name__ == '__main__':
    cg_test()
    check_hdf5()
