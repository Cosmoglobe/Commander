from astropy.io import fits
import numpy as np

def load_old_fits(filename):
    """
    Load a non-standard (old) FITS-like file, apply BSCALE/BZERO and mask BLANK.

    Parameters
    ----------
    filename : str
        Path to the file.

    Returns
    -------
    data : np.ndarray
        Scaled image data as float array, with blanks set to np.nan.
    header : astropy.io.fits.Header
        FITS header.
    """
    with fits.open(filename,
                   do_not_scale_image_data=True,
                   ignore_missing_end=True) as hdul:
        header = hdul[0].header
        raw = hdul[0].data  # raw array from file

        if raw is None:
            raise ValueError(f"No data found in primary HDU of {filename}")

        # mask blanks before scaling
        blank_val = header.get('BLANK')
        if blank_val is not None:
            mask = (raw == blank_val)
        else:
            mask = None

        # apply scaling
        bscale = header.get('BSCALE', 1.0)
        bzero = header.get('BZERO', 0.0)
        data = raw.astype(float) * bscale + bzero

        if mask is not None:
            data = np.where(mask, np.nan, data)

    return data, header


from glob import glob
import matplotlib.pyplot as plt

DIR = '/mn/stornext/d16/cmbco/ola/IRAS/kester_rawdb'


# spline and lost+found are not used, and/or are useless

band = 4

dt = 1/16 # I think?

# Flash times?
T_cut = 20

unknown_fnames = glob(f'{DIR}/disk*to*/d??/unknown.b{band}/sop.??_/??????')
survey_fnames = glob(f'{DIR}/disk*to*/d??/survey.b{band}/sop.??_/??????')
addobs_fnames = glob(f'{DIR}/disk*to*/d??/AO.b{band}/sop.??_/??????')


addobs_fnames.sort()

#d, h = load_old_fits(addobs_fnames[0])
d, h = load_old_fits(survey_fnames[0])
print(h)

fig, axes = plt.subplots(4, 4, sharex=True, sharey=True)
axs = axes.flatten()

for i in range(len(d)):
    axs[i].plot(d[i], 'k.', ms=1)

plt.show()
