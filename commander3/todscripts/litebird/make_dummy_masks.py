import fitsio
import numpy as np 




nside = 512
dummy_mask = np.ones((3, 12*nside**2), dtype=np.float32)


filename = "dummy_mask_n512.fits"
fitsio.write(filename, dummy_mask)



