'''
Make a WMAP_instrument.h5 file in /mn/stornext/d16/cmbco/bp/dwatts/wmap/data

compare to the LFI RIMO in
/mn/stornext/d16/cmbco/bp/dwatts/wmap/data/LFI_instrument_v4.h5

Top level includes the averaged channel name, plus the different horns,
[18--28](M/S).

030/040/070 include bandpass and bandpassx, which are B(nu) and nu,
respectively.

Each horn contains:
bandpass                 Dataset {251}
bandpassx                Dataset {251}
beam                     Group
beamlmax                 Dataset {1}
beammmax                 Dataset {1}
centFreq                 Dataset {1}
elip                     Dataset {1}
fwhm                     Dataset {1}
mbeam_eff                Dataset {1}
psi_ell                  Dataset {1}
sl                       Group
sllmax                   Dataset {1}
slmmax                   Dataset {1}


beam and sl are both alm representations of the beam and sidelobes.

mbeam_eff is main beam efficiency, assume it is one.
'''

from scipy import interpolate
from tqdm import tqdm
from scipy.optimize import curve_fit

import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from astropy.io import fits
import h5py

import healpy as hp

from scipy.integrate import trapz
from scipy.optimize import curve_fit, minimize
def getOutidx(l, m):
    return l**2 + l + m

def padAlms(data, lmax, mmax):
    outData = np.zeros((lmax+1)**2, dtype='complex64')

    for l in range(0, lmax):
        for m in range(0, mmax):
            healpixI = hp.sphtfunc.Alm.getidx(lmax, l, m)
            outData[healpixI] = data[healpixI]

    return outData

def complex2realAlms(data, lmax, mmax):
    outData = np.zeros((lmax+1)**2)

    for l in range(0, lmax):
        for m in range(0, mmax):
            if(m > l):
                continue
            scaling = np.sqrt(2)
            if(m == 0):
                scaling = 1
            healpixI = hp.sphtfunc.Alm.getidx(lmax, l, m)
            outI = getOutidx(l, m)
            outJ = getOutidx(l, -1*m)
            outData[outI] = np.real(data[healpixI]) * scaling
            if(m != 0):
                outData[outJ] = np.imag(data[healpixI]) * scaling

    return outData

def cl2realAlms(data, lmax, mmax):
    outData = np.zeros((lmax+1)**2)

    for l in range(0, min(lmax, len(data))):
        for m in range(0, mmax):
            if(m > l):
                continue
            scaling = np.sqrt(2)
            scaling = 1
            if(m == 0):
                scaling = 1
            #healpixI = hp.sphtfunc.Alm.getidx(lmax, l, m)
            outI = getOutidx(l, m)
            outJ = getOutidx(l, -1*m)
            outData[outI] = np.real(data[l]) * scaling
            if(m != 0):
                outData[outJ] = np.imag(data[l]) * scaling

    return outData

def gauss(x, sigma):
    return np.exp(-x**2/(2*sigma**2))


def real2complexAlms(data, lmax, mmax):
    outData = np.zeros((lmax+1)**2, dtype='complex64')
    for l in range(0, lmax):
        for m in range(0, mmax):
            if(m > l):
                continue
            healpixI = hp.sphtfunc.Alm.getidx(lmax, l, m)
            realI = getOutidx(l, m)
            realJ = getOutidx(l, -1*m)
            if (m == 0):
                outData[healpixI] = data[realI]
            else:
                outData[healpixI] = (data[realI] + 1j*data[realJ])/np.sqrt(2)

    return outData


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


rots = np.arange(0, 360, 45)
rots = [0]
for rot in rots:
  #fname_out = f'/mn/stornext/d16/cmbco/bp/dwatts/WMAP/data_WMAP/WMAP_rot{rot}.h5'
  fname_out = '/mn/stornext/d16/cmbco/bp/dwatts/WMAP/data_WMAP/WMAP_instrument_v9.h5'
  #fname_out = 'test.h5'
  #fname_out = '/mn/stornext/d16/cmbco/bp/dwatts/WMAP/data_WMAP/test.h5'
  
  
  with h5py.File(fname_out, 'a') as f:
      labels = ['K', 'Ka', 'Q', 'V', 'W']
      # Bandpasses
      fnames = glob('data/wmap_bandpass_*_v5.cbp')
      fnames.sort()
      
      
      totals = [[],[],[],[],[]]
      for i in range(len(fnames)):
          band = fnames[i].split('_')[2]
          nu, B1, B2 = np.loadtxt(fnames[i]).T
          f.create_dataset(band + '3/bandpassx', data=nu)
          f.create_dataset(band + '4/bandpassx', data=nu)
          f.create_dataset(band + '3/bandpass', data=B1)
          f.create_dataset(band + '4/bandpass', data=B2)
          centFreq1 = trapz(nu*B1, nu)/trapz(B1, nu)
          centFreq2 = trapz(nu*B2, nu)/trapz(B2, nu)
          f.create_dataset(band + '3/centFreq', data=[centFreq1])
          f.create_dataset(band + '4/centFreq', data=[centFreq2])
  
      
          for j in range(len(labels)):
              if labels[j] == band[:-2]:
                  totals[j].append([nu, B1])
                  totals[j].append([nu, B2])
      
      for i in range(len(labels)):
          bp = np.array(totals[i])
          nu, B = np.mean(bp, axis=0)
          f.create_dataset(labels[i] + '/bandpassx', data=nu)
          f.create_dataset(labels[i] + '/bandpass', data=B)
  
      # FWHMs
      ## From radial beam profiles
      DAs = ['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']
      fwhms = [0.93, 0.68, 0.53, 0.53, 0.35, 0.35, 0.23, 0.23, 0.23, 0.23]
      fnames = glob('data/wmap_symm_beam_profile_*_9yr_v5.txt')
      fnames.sort()
      for ind, fname in enumerate(fnames):
          theta, B = np.loadtxt(fname).T
  
          #B = B[1:]
          #theta = theta[1:]
  
          #hwhm_deg = theta[B <= B[5:].max()/2][0]
          #fwhm_deg = 2*hwhm_deg
          fwhm_deg = fwhms[ind]
          fwhm_arcmin= 60*fwhm_deg
          #plt.figure()
          #plt.semilogx(theta, B/B.max())
          #cdf = np.cumsum(B)/sum(B)
          #theta_sigma = cdf[cdf >= 0.34][0]
          #plt.axvline(hwhm_deg)
          #sigma = fwhm_deg/np.sqrt(8*np.log2(2))
          #plt.plot(theta, np.exp(-theta**2/(2*theta_sigma**2)))
      
          DA = fname.split('_')[4]
          #print(DA, fwhm_arcmin)
          #print('\n')
          f.create_dataset(DA + '13/fwhm', data=[fwhm_arcmin])
          f.create_dataset(DA + '14/fwhm', data=[fwhm_arcmin])
          f.create_dataset(DA + '23/fwhm', data=[fwhm_arcmin])
          f.create_dataset(DA + '24/fwhm', data=[fwhm_arcmin])
          """
          Beam ellipticity: WMAP doesn't report the ellipticity directly. They describe
          the beam-symmetrization map-making process in http://arxiv.org/abs/1212.5225.
          The beams for a single pointing are fairly elliptical and non-Gaussian.
          
          For a first pass, I will set these to zero.
          """
          f.create_dataset(DA + '13/elip', data=[0])
          f.create_dataset(DA + '14/elip', data=[0])
          f.create_dataset(DA + '23/elip', data=[0])
          f.create_dataset(DA + '24/elip', data=[0])
  
  
  
  
      fnames = glob('data/wmap_hybrid_beam_maps_*_9yr_v5.fits')
      fnames.sort()
      
      
      hdus = [fits.open(fname)[0] for fname in fnames]
      beamAs = [fits.open(fname)[0].data[0] for fname in fnames]
      beamBs = [fits.open(fname)[0].data[2] for fname in fnames]
      
      
      # pixel size is 0.04 ~ 58.6/nside
      lmax = 1700
      mmax = 100
      
      
      X = np.arange(-11.98, 12, 0.04)*np.pi/180
      Y = np.arange(-11.98, 12, 0.04)*np.pi/180
      X2 = np.linspace(X[0], X[-1], len(X)*5)
      Y2 = np.linspace(Y[0], Y[-1], len(Y)*5)
      xx, yy = np.meshgrid(X,Y)
      theta = 2*np.arcsin(np.sqrt(xx**2+yy**2)/2)
      phi = np.arctan2(yy, xx)



      xx, yy = np.meshgrid(X2,Y2)
      theta = 2*np.arcsin(np.sqrt(xx**2+yy**2)/2)
      phi = np.arctan2(yy, xx)
      
      #data = fits.open(fnames[-1])
      #beamA = data[0].data[0]
      #fig = plt.figure()
      #ax = fig.add_subplot(111)
      #ax.pcolormesh(X,Y,beamA)

      #fig = plt.figure()
      #ax = fig.add_subplot(111)
      #ax.pcolormesh(X2,Y2,beamA_2)

      #plt.show()

     
      #fig, axes = plt.subplots(nrows=5, ncols=2, sharey=True)
      #axs = axes.flatten()
      #lmaxes = [600, 800, 1000, 1000, 1200, 1200, 1500, 1500, 1500, 1500]
      for beam_ind, fname in enumerate(fnames):
          data = fits.open(fname)
          beamA = data[0].data[0]
          beamB = data[0].data[2]
          f = interpolate.interp2d(X, Y, beamA)
          beamA_2 = f(X2, Y2)
          f = interpolate.interp2d(X, Y, beamB)
          beamB_2 = f(X2, Y2)

          nside = 4096//4
          mA = np.zeros(12*nside**2)
          mB = np.zeros(12*nside**2)
          N = np.zeros(12*nside**2)
          #for i in range(len(theta)):
          #    for j in range(len(theta)):
          #        pix = hp.ang2pix(nside, theta[i,j], phi[i,j])
          #        mA[pix] += beamA[i,j]
          #        mB[pix] += beamB[i,j]
          #        N[pix] += 1
          #for i in tqdm(range(len(X2))):
          #    for j in range(len(Y2)):
          pix = hp.ang2pix(nside, theta, phi)
          mA[pix] += beamA_2
          mB[pix] += beamB_2
          N[pix] += 1
          #hp.gnomview(m/N, rot=(-130,86.38,0), reso=1)
          mA = mA/N
          mB = mB/N
          mA[~np.isfinite(mA)] = 0
          mB[~np.isfinite(mB)] = 0

          #hp.write_map(f'freq{beam_ind}_hornA.fits', mA)
          #hp.write_map(f'freq{beam_ind}_hornB.fits', mB)


          ind = np.argmax(mA)
          #fig = plt.figure()
          #ax = fig.add_subplot(121)
          #ax.pcolormesh(X,Y,beamA)
          #ax.set_aspect('equal')
          #ax = fig.add_subplot(122)
          #ax.pcolormesh(X,Y,beamB)
          #ax.set_aspect('equal')


          th, ph = hp.pix2ang(nside, ind)
          r = hp.rotator.Rotator(rot=(ph, -th, 0), \
              deg=False, eulertype='ZYX')
          mA = r.rotate_map_pixel(mA)

          #ind = np.argmax(mA)
          #lon, lat = hp.pix2ang(nside, ind, lonlat=True)
          #print(lon, lat)
          #hp.gnomview(mA, rot=(lon,lat,0), reso=0.5, title='A')

          ind = np.argmax(mB)
          th, ph = hp.pix2ang(nside, ind)
          r = hp.rotator.Rotator(rot=(ph, -th, 0), \
              deg=False, eulertype='ZYX')
          mB = r.rotate_map_pixel(mB)

          #lon, lat = hp.pix2ang(nside, ind, lonlat=True)
          #hp.gnomview(mB, rot=(lon,lat,0), reso=0.5, title='B')
          #plt.show()
          #print('\n')
          
          alm_A = hp.map2alm(mA, lmax=lmax, mmax=mmax)
          b_lm_A = complex2realAlms(alm_A, lmax, mmax)

          # # Test the real2complexAlms function
          # alm_A_test = real2complexAlms(b_lm_A, lmax, mmax)
          # padded_alm = padAlms(alm_A, lmax, mmax)
          # inds = (abs(padded_alm) > 0)
          # print(np.mean(abs(alm_A_test[inds] - padded_alm[inds])))
          # print(np.std(abs(alm_A_test[inds] - padded_alm[inds])))
          # print(np.mean(padded_alm[inds] - alm_A))
          # print(asdf)

          alm_B = hp.map2alm(mB, lmax=lmax, mmax=mmax)
          b_lm_B = complex2realAlms(alm_B, lmax, mmax)

          b_lA = hp.alm2cl(alm_A, lmax=lmax, mmax=mmax)**0.5
          b_lB = hp.alm2cl(alm_B, lmax=lmax, mmax=mmax)**0.5
          b_lm_A /= max(b_lA)
          b_lm_B /= max(b_lB)
          #axs[beam_ind].plot(b_lA/max(b_lA), label='A')
          #axs[beam_ind].plot(b_lB/max(b_lB), ':', label='B')
          #axs[beam_ind].set_ylim([0,1])
          #axs[beam_ind].set_xlim([0, lmaxes[beam_ind]])
          #axs[1].legend()

      fnames = glob('data/wmap_ampl_bl*.txt')
      fnames.sort()
      for beam_ind, fname in enumerate(fnames):

          l, b_lA, sl = np.loadtxt(fname).T

          b_lm_A = cl2realAlms(b_lA, lmax, lmax)

          #b_lm_A = b_lm_A*0 + 1
          #b_lm_B = b_lm_B*0 + 1
      
          DA = fname.split('_')[3]
           
          with h5py.File(fname_out, 'a') as f:
              f.create_dataset(DA + '13/beam/T', data=b_lm_A)
              f.create_dataset(DA + '14/beam/T', data=b_lm_A)
              f.create_dataset(DA + '23/beam/T', data=b_lm_A)
              f.create_dataset(DA + '24/beam/T', data=b_lm_A)
              f.create_dataset(DA + '13/beam/E', data=b_lm_A*0)
              f.create_dataset(DA + '14/beam/E', data=b_lm_A*0)
              f.create_dataset(DA + '23/beam/E', data=b_lm_A*0)
              f.create_dataset(DA + '24/beam/E', data=b_lm_A*0)
              f.create_dataset(DA + '13/beam/B', data=b_lm_A*0)
              f.create_dataset(DA + '14/beam/B', data=b_lm_A*0)
              f.create_dataset(DA + '23/beam/B', data=b_lm_A*0)
              f.create_dataset(DA + '24/beam/B', data=b_lm_A*0)
              f.create_dataset(DA + '13/beamlmax', data=[lmax])
              f.create_dataset(DA + '14/beamlmax', data=[lmax])
              f.create_dataset(DA + '23/beamlmax', data=[lmax])
              f.create_dataset(DA + '24/beamlmax', data=[lmax])
              f.create_dataset(DA + '13/beammmax', data=[mmax])
              f.create_dataset(DA + '14/beammmax', data=[mmax])
              f.create_dataset(DA + '23/beammmax', data=[mmax])
              f.create_dataset(DA + '24/beammmax', data=[mmax])
              
              f.create_dataset(DA + '13/mbeam_eff', data=[1])
              f.create_dataset(DA + '14/mbeam_eff', data=[1])
              f.create_dataset(DA + '23/mbeam_eff', data=[1])
              f.create_dataset(DA + '24/mbeam_eff', data=[1])
              f.create_dataset(DA + '13/psi_ell', data=[0])
              f.create_dataset(DA + '14/psi_ell', data=[0])
              f.create_dataset(DA + '23/psi_ell', data=[0])
              f.create_dataset(DA + '24/psi_ell', data=[0])
      
     
      plt.show()
      nside = 2**7
      sllmax = 512
      slmmax = 100
      labels = ['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']
      fnames = glob('data/wmap_sidelobe*.fits')
      for i in range(len(labels)):
        lab = labels[i]
        for fname in fnames:
          if lab in fname:
            data = hp.read_map(fname, nest=True)
            break
        
        print(lab, fname)
        # Beam is normalized such that sum(slAB) = Npix, or
        #                              \int B(\Omega)\,d\Omega = 4\pi
        # Commander expects \int B\,d\Omega = 1.
        # Not sure if this factor is needed...
        beamtot = hp.reorder(data, n2r=True)
        
        beam_A = hp.reorder(data, n2r=True)/(4*np.pi)
        #beam_A = hp.reorder(data, n2r=True)
        beam_A[beam_A < 0] = 0
        beam_B = hp.reorder(data, n2r=True)/(4*np.pi)
        #beam_B = hp.reorder(data, n2r=True)
        beam_B[beam_B > 0] = 0
        beam_B = -beam_B
  
        dir_A = dir_A_los[i]
        theta = np.arccos(dir_A[2])
        phi = np.arctan2(dir_A[1], dir_A[0])

       
        #hp.mollview(beam_A, min=0, max=0.3)
  
        r = hp.rotator.Rotator(rot=(phi, -theta, 0), \
            deg=False, eulertype='ZYX')
        beam_A_temp = r.rotate_map_pixel(beam_A)

        #hp.mollview(beam_A_temp, min=0, max=0.3)
  
        r = hp.rotator.Rotator(rot=(rot*np.pi/180, 0, 0), \
            deg=False, eulertype='ZYX')
        beam_A = r.rotate_map_pixel(beam_A_temp)
  
        #hp.mollview(beam_A, min=0, max=0.3)
  
        alm_A = hp.map2alm(beam_A, lmax=sllmax, mmax=slmmax)
        s_lm_A = complex2realAlms(alm_A, sllmax, slmmax)
  
        dir_B = dir_B_los[i]
        theta = np.arccos(dir_B[2])
        phi = np.arctan2(dir_B[1], dir_B[0])
        
        r = hp.rotator.Rotator(rot=(phi, -theta, 0), \
            deg=False, eulertype='ZYX')
        beam_B_temp = r.rotate_map_pixel(beam_B)
  
        r = hp.rotator.Rotator(rot=(-rot*np.pi/180, 0, 0), \
            deg=False, eulertype='ZYX')
        beam_B = r.rotate_map_pixel(beam_B_temp)
  
  
        alm_B = hp.map2alm(beam_B, lmax=sllmax, mmax=slmmax)
        s_lm_B = complex2realAlms(alm_B, sllmax, slmmax)
  
        DA = labels[i]
      
        with h5py.File(fname_out, 'a') as f:
            f.create_dataset(DA + '13/sl/T', data=s_lm_A)
            f.create_dataset(DA + '14/sl/T', data=s_lm_A)
            f.create_dataset(DA + '23/sl/T', data=s_lm_B)
            f.create_dataset(DA + '24/sl/T', data=s_lm_B)
            f.create_dataset(DA + '13/sl/E', data=s_lm_A*0)
            f.create_dataset(DA + '14/sl/E', data=s_lm_A*0)
            f.create_dataset(DA + '23/sl/E', data=s_lm_B*0)
            f.create_dataset(DA + '24/sl/E', data=s_lm_B*0)
            f.create_dataset(DA + '13/sl/B', data=s_lm_A*0)
            f.create_dataset(DA + '14/sl/B', data=s_lm_A*0)
            f.create_dataset(DA + '23/sl/B', data=s_lm_B*0)
            f.create_dataset(DA + '24/sl/B', data=s_lm_B*0)
            f.create_dataset(DA + '13/sllmax', data=[sllmax])
            f.create_dataset(DA + '14/sllmax', data=[sllmax])
            f.create_dataset(DA + '23/sllmax', data=[sllmax])
            f.create_dataset(DA + '24/sllmax', data=[sllmax])
            f.create_dataset(DA + '13/slmmax', data=[slmmax])
            f.create_dataset(DA + '14/slmmax', data=[slmmax])
            f.create_dataset(DA + '23/slmmax', data=[slmmax])
            f.create_dataset(DA + '24/slmmax', data=[slmmax])
