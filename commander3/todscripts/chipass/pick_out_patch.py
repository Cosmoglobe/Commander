import numpy as np 
import matplotlib.pyplot as plt
import h5py
import healpy as hp
import tools


def find_local_scans(lon0, lat0, dlon=10, dlat=10):
    filename = 'chipass_full.h5'
    filenames = []
    j = 0
    f1 = h5py.File(filename, 'r')
    for key, data in f1.items():
        # print(key)

        lon = data['point_gal'][:, :, 0].flatten()
        lat = data['point_gal'][:, :, 1].flatten()

        delta_lon = np.abs(lon - lon0)
        delta_lon[delta_lon > 180] = 360 - delta_lon[delta_lon > 180]
        delta_lat = np.abs(lat - lat0)
        wh = np.where(((delta_lon < dlon) & (delta_lat < dlat)))

        n_hit = len(wh[0])
        if n_hit > 0:
            print(n_hit)
            #filenames.append(key.replace(" ", ""))
            filenames.append(key)

        j = j + 1
        # if j > 10000:
        #     break
    f1.close()
    return filenames

def combine_patch_files(file_list_path, nside=256):
    scans = np.genfromtxt(file_list_path, dtype='str', delimiter=',')
    filename = 'chipass_full.h5'
    f1 = h5py.File(filename, 'r')
    tods = []
    pix = []
    nsamp = []
    for scan in scans:
        data = f1[scan]
        point = data['point_gal'][:, :, :2] * np.pi / 180

        tod = data['tod'][()]
        
        # phi in healpix is 90 - phi from the file
        tod = np.sqrt(tod[:, :, 0] ** 2 + tod[:, :, 1] ** 2).mean(2)
        n = len(tod[:, 0])
        indices = hp.ang2pix(nside, np.pi / 2.0 - point[:,:,1], point[:,:,0])
        # print(tod.shape)
        # print(indices.shape)
        tods.append(tod)
        pix.append(indices)
        nsamp.append(n)
    
    return tods, pix, nsamp


def find_baseline_lenghts(n_t, n_baselines):
    divide = n_t // n_baselines
    lefovers = n_t % n_baselines

    baseline_lengths = np.zeros(n_baselines, dtype=int)
    baseline_lengths[:] = divide 

    if lefovers > 0:
        baseline_lengths[:lefovers] += 1
    
    indices = np.zeros(n_t, dtype=int)
    j = 0
    for i in range(n_baselines):
        indices[j:j + baseline_lengths[i]] = i
        j += baseline_lengths[i]
    return indices

def find_F_indices(nsamp, n_baselines):
    F_indices = np.zeros(n_tot, dtype=int)
    baseline_ind = 0
    samp_ind = 0
    print(n_scans)
    for i in range(n_scans):
        F_indices[samp_ind:samp_ind + nsamp[i]] = find_baseline_lenghts(nsamp[i], n_baselines) + baseline_ind

        samp_ind += nsamp[i]
        baseline_ind += n_baselines
    return F_indices

# def apply_F(x):
#     return x[F_indices, :]

# def apply_FT(z):
#     FTz = np.zeros((n_scans * n_baselines, n_det))
#     for i in range(n_det):
#         FTz[:, i] = np.bincount(F_indices, weights=z[:, i])
#     return FTz

def apply_F(x):
    Fx = np.zeros((n_tot, n_det))
    return Fx + x[:, None]

def apply_FT(z):
    return z.sum(1)

def apply_Z(y):
    my_map = np.bincount(pix.flatten(), weights=y.flatten(), minlength=npix)
    return y - (my_map / nhit)[pix.flatten()].reshape((n_tot, n_det))



def apply_A(x):
    return apply_FT(apply_Z(apply_F(x)))

def run_CG(x, b, max_iter=1, eps=0.0001):
    r = b - apply_A(x)
    p = r
    r2 = np.dot(r.flatten(), r.flatten())
    for k in range(max_iter):
        Ap = apply_A(p)
        alp = r2 / np.dot(p.flatten(), Ap.flatten())
        x = x + alp * p
        r = r - alp * Ap
        r2_new = np.dot(r.flatten(), r.flatten())
        print(np.sqrt(r2_new))
        if np.sqrt(r2_new) < eps:
            return x 
        bet = r2_new / r2
        p = r + bet * p 
        r2 = r2_new
    print("Max iterations reached without convergence")
    return x



if __name__ == '__main__':
    # print(find_baseline_lenghts(13, 3))
    # sys.exit()
    # filenames = find_local_scans(0, 0, dlon=180, dlat=180)
    # np.savetxt('patch_files_full.txt', filenames, fmt="%s")
    nside = 512


    tods, pix, nsamp = combine_patch_files('patch_files_full.txt', nside=nside)

    # Combine all scans (to bet nice numpy arrays)
    tod = np.concatenate(tods, axis=0)
    pix = np.concatenate(pix, axis=0)

    # Global parameters used in the CG method
    n_det = 13
    n_scans = len(nsamp)
    n_tot = np.sum(nsamp)
    n_baselines = 3  # baselines per scan
    F_indices = find_F_indices(nsamp, n_baselines)
    npix = hp.nside2npix(nside)
    nhit = np.bincount(pix.flatten(), minlength=npix)
    # plt.hist(F_indices)
    # plt.show()

    # Initate the map and fill it with the values                                                                                              
    # nhit = np.zeros(npix, dtype=np.float)                                                                                                      
    # hpxmap = np.zeros(npix, dtype=np.float)
    # for i in range(13):
    #     nhit[pix[:, i]] += 1
    #     hpxmap[pix[:, i]] += tod[:, i]

    
    # pix = pix.flatten()
    # tod = tod.flatten()
    # baseline_len = 39
    # a = np.random.randn(len(tod) // baseline_len) * 1e-6

    # a = np.random.randn(n_scans * n_baselines, n_det) * 1e-6


    # a = np.random.randn(n_tot) * 1e-6
    # b = apply_FT(apply_Z(tod))

    # a_out = run_CG(a, b, max_iter=350)

    # nside = 512
    # tods, pix, nsamp = combine_patch_files('patch_files.txt', nside=nside)
    # tod = np.concatenate(tods, axis=0)
    # pix = np.concatenate(pix, axis=0)
    # npix = hp.nside2npix(nside)
    # nhit = np.bincount(pix.flatten(), minlength=npix)
    binmap = np.bincount(pix.flatten(), weights=tod.flatten(), minlength=npix)
    binmap = binmap / nhit

    s_bin = binmap[pix]

    # dstmap = np.bincount(pix.flatten(), weights=(tod.flatten() - apply_F(a_out).flatten()), minlength=npix)
    # dstmap = dstmap / nhit
    
    # med = np.median(tod, axis=1)  # this is the wrong way to do it
    # medmap = np.bincount(pix.flatten(), weights=(tod - med[:, None]).flatten(), minlength=npix)
    # medmap = medmap / nhit
    maxval = 20

    xsize = 600 #330
    ysize = 400 #260
    #hp.gnomview(np.log10(binmap+2), reso=5, xsize=xsize, ysize=ysize, min=0, max=2)
    hp.mollview(np.log10(binmap+2), min=0, max=2)
    plt.savefig('binned.pdf', bbox_inches='tight')
    # hp.gnomview(np.log10(dstmap+2), reso=5, xsize=330, ysize=260, min=0, max=2)
    #hp.mollview(np.log10(hpxmap), xsize=4000) #, min=-maxval, max=maxval) 

    # hp.gnomview(np.log10(medmap+2), reso=5, xsize=330, ysize=260, min=0, max=2)
    print(tod.shape)
    for i in range(3):
        mask = (np.abs(tod - s_bin) < 6.0 + 0.2 * s_bin) * 1.0
        nhit = np.bincount(pix.flatten(), minlength=npix, weights=mask.flatten())
        maskmap = np.bincount(pix.flatten(), weights=tod.flatten() * mask.flatten(), minlength=npix)
        maskmap = maskmap / nhit
        #hp.gnomview(np.log10(maskmap+2), reso=5, xsize=330, ysize=260, min=0, max=2)
        print("Fraction used: ", mask.flatten().sum() / len(mask.flatten()))
        s_bin = maskmap[pix]
    # np.savetxt('maskmap_%i' % nside, maskmap) 
    hf = h5py.File('chipass_maskmap_%i.h5' % nside, 'w')
    hf.create_dataset('map', data=maskmap)
    hf.close()
    hp.mollview(np.log10(maskmap+2), min=0, max=2)
    # hp.gnomview(np.log10(maskmap+2), reso=5, xsize=xsize, ysize=ysize, min=0, max=2)
    plt.savefig('binned_mask.pdf', bbox_inches='tight')
    # hp.gnomview(maskmap - binmap, reso=5, xsize=xsize, ysize=ysize, max=3, min=-3)
    hp.mollview(maskmap - binmap, max=3, min=-3)
    plt.savefig('difference.pdf', bbox_inches='tight')
    nu = 1.3945e9 # Hz
    area = 2.0 * np.pi * (14.0 / 60 / 180 * np.pi / np.sqrt(8 * np.log(2))) ** 2

    dnu = 64e6 / 1024
    jy2k = tools.jy2kcmb(1, nu, area)
    chipass = hp.read_map("lambda_chipass_healpix_r10.fits") * 1e-3
    print("median: ", np.median(chipass.flatten()))
    mono = np.median(chipass.flatten())
    hp.mollview(maskmap - chipass, max=3 - mono, min=-3 - mono)
    plt.savefig('difference_vs_official.pdf', bbox_inches='tight')
    hp.mollview(np.log10(chipass+2), min=0, max=2)
    plt.savefig('official.pdf', bbox_inches='tight')
    # cosmo = hp.read_map("1394MHz.fits")
    # cosmo = hp.smoothing(cosmo, fwhm=np.radians(14.0 / 60))
    # hp.gnomview(np.log10(cosmo+2), reso=5, xsize=xsize, ysize=ysize, min=0, max=2)
    # plt.figure()
    # plt.plot(tod[:1500, 0])
    # plt.plot(s_bin[:1500, 0])
    # plt.plot(mask[:1500, 0]*30)
    plt.show()




# def apply_F(x, baseline_len=39):
#     Fx = np.zeros((len(x), baseline_len))
#     return (Fx + x[:, None]).flatten()

# def apply_Z(y):
#     my_map = np.bincount(pix, weights=y, minlength=npix)
#     return y - (my_map / nhit)[pix]

# def apply_FT(z, baseline_len=39):
#     return z.reshape((len(z) // baseline_len, baseline_len)).sum(1)
