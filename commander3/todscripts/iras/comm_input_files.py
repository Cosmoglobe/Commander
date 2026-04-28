import numpy as np
import healpy
import glob
import astropy.units as u
from pathlib import Path
from astropy.io import fits
import matplotlib.pyplot as plt 

from iras_native_tod_reader import BANDS

OUTDIR = "/mn/stornext/d23/cmbco/globe/iras"
akaridir = "/mn/stornext/d23/cmbco/globe/akari" # Used for testing


def make_beam_files(nside, unit=u.arcsec):
    input("Make beam files? Enter to cnt.")
    # BEAM_SIZES = {
    #     # Values taken from: https://irsa.ipac.caltech.edu/IRASdocs/exp.sup/ch2/C3.html#tabC2
    #     "012": 25  * u.arcsec,
    #     "025": 25  * u.arcsec,
    #     "060": 60  * u.arcsec,
    #     "100": 100 * u.arcsec,
    # }

    BEAM_SIZES = {
        # Values taken from IRIS: Table 1, row 2, in Miville-Deschenes, 2004 (https://arxiv.org/abs/astro-ph/0412216)
        "012": 3.8 * u.arcmin,
        "025": 3.8 * u.arcmin,
        "060": 4.0 * u.arcmin,
        "100": 4.3 * u.arcmin,
    }

    l_max = int(2 * nside)

    for band in BANDS:
        beamsize = BEAM_SIZES[band].to(unit)
        beamsize_value = str(beamsize.value).replace(r".0", "")
        beamsize_unit = str(unit)
        outfile = Path(f"{OUTDIR}/beam/iras_{band}_beam_{beamsize_value}{beamsize_unit}_n{nside}.fits")

        beamsize_rad = beamsize.to(u.rad).value

        beam = healpy.gauss_beam(beamsize_rad, l_max)

        # hdu     = fits.PrimaryHDU()
        # col     = fits.Column(name='TEMPERATURE', array=beam, format='E')
        # tab     = fits.TableHDU.from_columns([col], name="Power spectrum")        
        # hdul    = fits.HDUList([hdu, tab])
        # hdul.writeto(outfile, overwrite=True)

        healpy.write_cl(outfile, 
                        beam, 
                        overwrite=True, 
                        dtype="float32", 
                        extra_header=["PSPEC"],
                        column_names=["TEMPERATURE"])


def make_maps_files(nside=512):
    """
    Use iris maps as the default map
    """
    input("Make maps files? Enter to cnt.")
    iris_path = "/mn/stornext/d16/cmbco/ola/IRAS/iris"
    iris_band_maps = {
        # See /mn/stornext/d16/cmbco/ola/IRAS/iris/iris_info.txt
        "012": f"{iris_path}/IRIS_nohole_1_2048_v2.fits",
        "025": f"{iris_path}/IRIS_nohole_2_2048_v2.fits",
        "060": f"{iris_path}/IRIS_nohole_3_2048_v2.fits",
        "100": f"{iris_path}/IRIS_nohole_4_2048_v2.fits",
    }

    for band in BANDS:
        outfname = f"{OUTDIR}/maps/iras_{band}_1_n{nside}.fits"
        if Path(outfname).exists():
            print(f"File {outfname} exists, skipping...")
            continue
        iris_map = healpy.read_map(iris_band_maps[band])
        # iris_map, header = healpy.read_map(iris_band_maps[band], h=True)
        # print(f"{healpy.npix2nside(len(iris_map))=}")
        # print()
        # print(f"{*header}")
        # print(*header, sep="\n")
        if nside == 2048:
            # This is the nside used in the iris maps.
            healpy.write_map(outfname, iris_map)
        else:
            out_map = healpy.ud_grade(iris_map, nside_out=nside)
            print(f"Storing {outfname}")
            healpy.write_map(outfname, out_map)

def make_mask_files(nside=512):
    input("Make mask files? Enter to cnt.")

    iris_path = "/mn/stornext/d16/cmbco/ola/IRAS/iris"
    iris_band_maps = {
        # See /mn/stornext/d16/cmbco/ola/IRAS/iris/iris_info.txt
        "012": f"{iris_path}/IRIS_1_2048.fits",
        "025": f"{iris_path}/IRIS_2_2048.fits",
        "060": f"{iris_path}/IRIS_3_2048.fits",
        "100": f"{iris_path}/IRIS_4_2048.fits",
    }
    for band in BANDS:
        outfname = f"{OUTDIR}/mask/bitmask_iras_{band}_v1_n{nside}.fits"
        if Path(outfname).exists():
            print(f"File {outfname} exists, skipping...")
            continue
        iris_map, h = healpy.read_map(iris_band_maps[band], h=True)
        # um = healpy.ud_grade(iris_map, nside_out=nside)
        # um[um <= 0] = np.nan
        # um[um > 0] = 1.0

        if nside == 2048:
            # This is the nside used in the iris maps.
            out_map = np.where(iris_map <= 0, 0.0, 1.0)
            healpy.write_map(outfname, out_map)
        else:
            # Setting a high and low value to avoid zeros during ud_grade.
            # Maybe not relevant.
            # iris_map[iris_map > 0] = 100
            # iris_map[iris_map <= 0] = 10
            iris_map = np.where(iris_map <= 0, 10.0, 100.0)
            out_map = healpy.ud_grade(iris_map, nside_out=nside)
            out_map = np.where(out_map < 50, 0.0, 1.0)
            print(f"Storing {outfname}")
            healpy.write_map(outfname, out_map)




def make_rms_files(nside=512):
    input("Make rms files? Enter to cnt.")
    for band in BANDS:
        outfname = f"{OUTDIR}/rms/iras_{band}_temprms_n{nside}.fits"
        rms_map = np.ones(healpy.nside2npix(nside))
        healpy.write_cl(outfname, rms_map)


def make_tod_files(nside=512):
    pass

def make_bp_files(nside=512):
    pass

exit()

def test_map():
    # fname = "/mn/stornext/d16/cmbco/ola/IRAS/iris/IRIS_nohole_2_2048_v2.fits"
    fname = "/mn/stornext/d16/cmbco/ola/IRAS/issa/ISSA_B4H0_healpix_ns1024.fits"
    iris25, h = healpy.read_map(fname, field=None, h=True)

    # print(iris25.shape)
    # nn = len(iris25[1])
    # print(healpy.npix2nside(nn))
    # print(*h, sep='\n')
    healpy.mollview(iris25[0]/iris25[1], min=0, max=50, cmap="jet")#norm="hist")
    plt.show()
    # healpy.mollview(iris25, min=0, max=10, cmap="jet")#norm="hist")
    # plt.show()






# test_map()
