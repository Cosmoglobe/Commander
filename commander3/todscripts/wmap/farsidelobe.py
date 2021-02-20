import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

from glob import glob

from astropy.io import fits
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

# I am nearly certain that the projection that they use, 
# X=2*sin(theta/2)*cos(phi), is zenithal equal area, making the coordinate
# system centered at the north pole.
#The beam coordinates form an equal area rectangular coordinate system centered on the optic axis of the spacecraft. They are related to coordinates theta (elevation from optic axis) and phi (azimuth about optic axis) as follows:
#
#    Xbeam = 2*sin(theta/2) * cos(phi)
#    Ybeam = 2*sin(theta/2) * sin(phi)
#
#The "optic axis" of the spacecraft is elevated by 19.5 degrees from the S/C XY plane and lies within the S/C YZ plane. Although this vector is close to the S/C Y axis (+ or - depending on A or B side), it becomes the Z axis of the focal plane coordinate system. 
target_header = fits.Header.fromstring("""
NAXIS   =                    2
NAXIS1  =                  600
NAXIS2  =                  600
CTYPE1  = 'GLON-ZEA'
CRPIX1  =                0.5
CRVAL1  =                -11.98
CDELT1  =               0.04
CUNIT1  = 'deg     '
CTYPE2  = 'GLAT-ZEA'
CRPIX2  =                0.5
CRVAL2  =                -11.98
CDELT2  =                0.04
CUNIT2  = 'deg     '
COORDSYS= 'icrs    '
""", sep='\n')

fnames = glob('data/wmap_hybrid_beam_maps_*_9yr_v5.fits')
fnames.sort()
fnames = [fnames[0]]
nside_beam = 2**10
for i, fname in enumerate(fnames):
    data = fits.open(fname)
    #plt.figure()
    #plt.imshow(np.log10(data[0].data[0]))
    #plt.figure()
    #plt.imshow(np.log10(data[0].data[2]))
    #plt.show()
    
    m_A, footprint_A = reproject.reproject_to_healpix((data[0].data[0], target_header), 'G', 
                                                   nside=nside_beam,
                                                   order=3)

    m_B, footprint_B = reproject.reproject_to_healpix((data[0].data[2], target_header), 'G', 
                                                   nside=nside_beam,
                                                   order=3)
    m_A[footprint_A==0] = 0
    m_B[footprint_B==0] = 0

    hp.mollview(m_A, min=-0.5, max=0.5, cmap='RdBu_r')
    hp.gnomview(m_A, min=-0.5, max=0.5, cmap='RdBu_r', reso=20)
    hp.gnomview(m_B, min=-0.5, max=0.5, cmap='RdBu_r', reso=20)

    dir_A = dir_A_los[i]
    theta = np.arccos(dir_A[2])
    phi = np.arctan2(dir_A[1], dir_A[0])
    
    hp.orthview(m_A, cmap='Reds', min=0, max=0.5, rot=(0,90,0), title=fname)
    hp.graticule()
plt.show()


labels = ['K1', 'Ka1', 'Q1', 'Q1', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']
fnames = glob('data/wmap_sidelobe*.fits')
fnames.sort()
fnames = np.array(fnames)

labels = [labels[0]]
for i in range(len(labels)):
  lab = labels[i]
  for fname in fnames:
    if lab in fname:
      data = hp.read_map(fname, nest=True)
      break

  print(sum(abs(data)), len(data))

  # The WMAP data are normalized such that
  # \int G\,d\Omega = 4\pi
  # while Commander assumes the full-sky integral is equal to 1.
  # If the data array is G, then we need to divide it by 4\pi to give the
  # quantity that Commander expects.

  beamtot = hp.reorder(data, n2r=True)
  hp.mollview(beamtot, min=-0.5, max=0.5, cmap='RdBu_r', cbar=False,
      title = fname)
  hp.orthview(beamtot, min=-0.5, max=0.5, cmap='RdBu_r', cbar=False,
      title = fname)
  hp.graticule()
  
  beam_A = hp.reorder(data, n2r=True)
  beam_A[beam_A < 0] = 0
  beam_B = hp.reorder(data, n2r=True)
  beam_B[beam_B > 0] = 0
  beam_B = -beam_B
  
  #plt.figure()
  #hp.mollview(beam_A, min=-0.5, max=0.5, cmap='RdBu_r', cbar=False,
  #    title='Beam A', sub=121)
  #hp.mollview(-beam_B, min=-0.5, max=0.5, cmap='RdBu_r', cbar=False,
  #    title='Beam B', sub=122)
  #
  #
  dir_A = dir_A_los[i]
  theta = np.arccos(dir_A[2])
  phi = np.arctan2(dir_A[1], dir_A[0])
  
  r = hp.rotator.Rotator(rot=(0, theta, phi), deg=False, eulertype='X')
  beam_A = r.rotate_map_pixel(beam_A)
  
  
  hp.orthview(beam_A, cmap='Reds', min=0, max=0.5, rot=(0,90,0), title=fname)
  hp.graticule()
  hp.gnomview(beam_A, cmap='Reds', min=0, max=0.5, rot=(0,90,0), title=fname,
      reso=20)
  hp.graticule()
  #
  #dir_B = dir_B_los[i]
  #theta = np.arccos(dir_B[2])
  #phi = np.arctan2(dir_B[1], dir_B[0])
  #
  #r = hp.rotator.Rotator(rot=(0, -theta, -phi), deg=False, eulertype='X')
  #beam_B = r.rotate_map_pixel(beam_B)
  #
  #hp.orthview(beam_B, cmap='Reds', min=0, max=0.5, rot=(0,90,0))
  #hp.graticule()
  #
  #
  #hp.mollview(beam_A, min=0, max=1, cmap='afmhot')
    
print(sum(m_A)*hp.nside2pixarea(nside_beam))
print(sum(data[data>0])*hp.nside2pixarea(hp.npix2nside(len(data))))
print(4*np.pi)
plt.show()
