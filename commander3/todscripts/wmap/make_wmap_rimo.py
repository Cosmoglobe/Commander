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

def gauss(x, sigma):
    return np.exp(-x**2/(2*sigma**2))

from scipy.optimize import curve_fit

import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import h5py

import healpy as hp

from scipy.integrate import trapz
from scipy.optimize import curve_fit, minimize

import reproject


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
for rot in rots:
  fname_out = f'/mn/stornext/d16/cmbco/bp/dwatts/WMAP/data_WMAP/WMAP_rot{rot}.h5'
  
  
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
      wcs = WCS(hdus[0].header) # hdus are the same
      
      
      # pixel size is 0.04 ~ 58.6/nside
      nside_beam = 512
      lmax = 2*nside_beam
      mmax = 100
      
      
      # I am nearly certain that the projection that they use, 
      # X=2*sin(theta/2)*cos(phi), is zenithal equal area, making the coordinate
      # system centered at the north pole.
      #target_header = fits.Header.fromstring("""
      #NAXIS   =                    2
      #NAXIS1  =                  600
      #NAXIS2  =                  600
      #CTYPE1  = 'GLON-ZEA'
      #CRPIX1  =                0.5
      #CRVAL1  =                -11.98
      #CDELT1  =               0.04
      #CUNIT1  = 'deg     '
      #CTYPE2  = 'GLAT-ZEA'
      #CRPIX2  =                0.5
      #CRVAL2  =                -11.98
      #CDELT2  =                0.04
      #CUNIT2  = 'deg     '
      #COORDSYS= 'icrs    '
      #""", sep='\n')
      
      for fname in fnames:
          data = fits.open(fname)
          
          b_lm_A = np.zeros((lmax+1)**2) + 1
          b_lm_B = np.zeros((lmax+1)**2) + 1
      
          DA = fname.split('_')[4]
      
          with h5py.File(fname_out, 'a') as f:
              f.create_dataset(DA + '13/beam/T', data=b_lm_A)
              f.create_dataset(DA + '14/beam/T', data=b_lm_A)
              f.create_dataset(DA + '23/beam/T', data=b_lm_B)
              f.create_dataset(DA + '24/beam/T', data=b_lm_B)
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
      
      
      nside = 2**7
      sllmax = lmax
      slmmax = 100
      labels = ['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']
      fnames = glob('data/wmap_sidelobe*.fits')
      for i in range(len(labels)):
        lab = labels[i]
        for fname in fnames:
          if lab in fname:
            data = hp.read_map(fname, nest=True)
            break
      
     
        # Beam is normalized such that \int B(\Omega)\,d\Omega = 4\pi, Commander
        # expects \int B\,d\Omega = 1.
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
        
  
        r = hp.rotator.Rotator(rot=(phi, -theta, 0), \
            deg=False, eulertype='ZYX')
        beam_A_temp = r.rotate_map_pixel(beam_A)
  
        r = hp.rotator.Rotator(rot=(rot*np.pi/180, 0, 0), \
            deg=False, eulertype='ZYX')
        beam_A = r.rotate_map_pixel(beam_A_temp)
  
  
        s_lm_A = np.zeros((lmax+1)**2)
        alm_A = hp.map2alm(beam_A)
  
        for j in range(len(alm_A)):
            li, mi = hp.sphtfunc.Alm.getlm(sllmax,i=j)
            #if (li <= lmax) & (mi <= slmmax):
            if (li <= lmax) & (mi <= slmmax) & (mi == 0):
                ind_real = li**2 + li + mi
                ind_imag = li**2 + li - mi
                s_lm_A[ind_real] = alm_A[j].real
                if mi != 0:
                    s_lm_A[ind_imag] = alm_A[j].imag
  
        dir_B = dir_B_los[i]
        theta = np.arccos(dir_B[2])
        phi = np.arctan2(dir_B[1], dir_B[0])
        
        r = hp.rotator.Rotator(rot=(phi, -theta, 0), \
            deg=False, eulertype='ZYX')
        beam_B_temp = r.rotate_map_pixel(beam_B)
  
        r = hp.rotator.Rotator(rot=(-rot*np.pi/180, 0, 0), \
            deg=False, eulertype='ZYX')
        beam_B = r.rotate_map_pixel(beam_B_temp)
  
  
        s_lm_B = np.zeros((lmax+1)**2)
        alm_B = hp.map2alm(beam_B)
        plt.show()
  
        for j in range(len(alm_B)):
            li, mi = hp.sphtfunc.Alm.getlm(sllmax,i=j)
            #if (li <= lmax) & (mi <= slmmax):
            if (li <= lmax) & (mi <= slmmax) & (mi == 0):
                ind_real = li**2 + li + mi
                ind_imag = li**2 + li - mi
                s_lm_B[ind_real] = alm_B[j].real
                if mi != 0:
                    s_lm_B[ind_imag] = alm_B[j].imag
  
  
        DA = labels[i]
      
        with h5py.File(fname_out, 'a') as f:
            f.create_dataset(DA + '13/sl/T', data=s_lm_A)
            f.create_dataset(DA + '14/sl/T', data=s_lm_A)
            f.create_dataset(DA + '23/sl/T', data=s_lm_B)
            f.create_dataset(DA + '24/sl/T', data=s_lm_B)
            f.create_dataset(DA + '13/sllmax', data=[sllmax])
            f.create_dataset(DA + '14/sllmax', data=[sllmax])
            f.create_dataset(DA + '23/sllmax', data=[sllmax])
            f.create_dataset(DA + '24/sllmax', data=[sllmax])
            f.create_dataset(DA + '13/slmmax', data=[slmmax])
            f.create_dataset(DA + '14/slmmax', data=[slmmax])
            f.create_dataset(DA + '23/slmmax', data=[slmmax])
            f.create_dataset(DA + '24/slmmax', data=[slmmax])
