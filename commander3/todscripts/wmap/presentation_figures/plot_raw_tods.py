import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits


data = fits.open('/mn/stornext/d16/cmbco/ola/wmap/tods/uncalibrated/wmap_tod_20083292351_20083302351_uncalibrated_v5.fits')

data2 = fits.open('/mn/stornext/d16/cmbco/ola/wmap/tods/calibrated/wmap_tod_20083292351_20083302351_calibrated_v5.fits')



d1 = data[2].data['K114']

d1_flat = np.zeros(d1.shape[0]*d1.shape[1])

for i in range(12):
  d1_flat[i::12] = d1[:,i]

plt.figure(figsize=(16, 3))
plt.plot(d1_flat[:5000], 'k.', ms=2)
plt.xlabel(r'Samples')
plt.ylabel('Digital units')
plt.title('Raw Time-ordered Data', size=20)
plt.savefig('raw_tod.png', bbox_inches='tight', transparent=True, dpi=300)


plt.figure(figsize=(16, 3))
d1 = data2[2].data['K114']

d1_flat = np.zeros(d1.shape[0]*d1.shape[1])

for i in range(12):
  d1_flat[i::12] = d1[:,i]

plt.plot(d1_flat[:5000], 'k.', ms=2)
plt.xlabel(r'Samples')
plt.ylabel('Differential Temperature (mK)')
plt.title('Calibrated Time-ordered Data', size=20)
plt.show()
