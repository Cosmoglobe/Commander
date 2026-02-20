from astropy.time import Time
from astropy.coordinates import (
    EarthLocation,
    HeliocentricMeanEcliptic,
    Galactic,
    SkyCoord,
    angular_separation,
)
from astropy.coordinates import (
    get_body_barycentric,
    get_body,
    SphericalRepresentation,
    CartesianRepresentation,
)
import astropy.units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy_healpix import HEALPix


from astroquery.jplhorizons import Horizons
from astroplan import Observer
from astroplan.moon import moon_phase_angle

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import h5py
import glob
import os
from typing import Optional

import tqdm

from pixell import enmap, utils
from cosmoglobe.tod_tools import TODLoader


class EphemerisFlagger:
    def __init__(
        self,
        ephemeris_file_dir: Optional[
            str
        ] = "/mn/stornext/d23/cmbco/globe/common/aux/ephemeris",
    ):
        if not os.path.isdir(ephemeris_file_dir):
            raise ValueError(f"Provided path {ephemeris_file_dir} is not a directory.")

        self.ephemeris_file_dir = ephemeris_file_dir
        self.available_ephemeris_files = glob.glob(
            os.path.join(ephemeris_file_dir, "*.txt")
        )

    def load_ephemeris(self, ephemeris_file: str):
        ephemeris_data = np.loadtxt(ephemeris_file, skiprows=2)
        mjd = Time(ephemeris_data[:, 0], format="mjd", scale="tdb")
        pos = SkyCoord(
            x=ephemeris_data[:, 1] * u.au,
            y=ephemeris_data[:, 2] * u.au,
            z=ephemeris_data[:, 3] * u.au,
            representation_type="cartesian",
            frame=HeliocentricMeanEcliptic,
        )
        return mjd, pos

    def linearly_interpolate_ephemeris(
        self, ephemeris_time: Time, ephemeris_position: SkyCoord, target_time: Time
    ):
        if target_time < ephemeris_time[0] or target_time > ephemeris_time[-1]:
            raise ValueError("Target time is out of bounds of the ephemeris data.")

        ephemeris_x = np.interp(
            target_time.tdb.mjd,
            ephemeris_time.tdb.mjd,
            ephemeris_position.x.to(u.km).value,
        )
        ephemeris_y = np.interp(
            target_time.tdb.mjd,
            ephemeris_time.tdb.mjd,
            ephemeris_position.y.to(u.km).value,
        )
        ephemeris_z = np.interp(
            target_time.tdb.mjd,
            ephemeris_time.tdb.mjd,
            ephemeris_position.z.to(u.km).value,
        )

        interpolated_position = SkyCoord(
            x=ephemeris_x * u.km,
            y=ephemeris_y * u.km,
            z=ephemeris_z * u.km,
            representation_type="cartesian",
            frame=HeliocentricMeanEcliptic,
        )

        return interpolated_position
