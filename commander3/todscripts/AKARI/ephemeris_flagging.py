from astropy.time import Time
from astropy.coordinates import (
    EarthLocation,
    HeliocentricMeanEcliptic,
    Galactic,
    SkyCoord,
    angular_separation,
    Angle,
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
            ephemeris_position.cartesian.x.to(u.km).value,
        )
        ephemeris_y = np.interp(
            target_time.tdb.mjd,
            ephemeris_time.tdb.mjd,
            ephemeris_position.cartesian.y.to(u.km).value,
        )
        ephemeris_z = np.interp(
            target_time.tdb.mjd,
            ephemeris_time.tdb.mjd,
            ephemeris_position.cartesian.z.to(u.km).value,
        )

        interpolated_position = SkyCoord(
            x=ephemeris_x * u.km,
            y=ephemeris_y * u.km,
            z=ephemeris_z * u.km,
            representation_type="cartesian",
            frame=HeliocentricMeanEcliptic,
        )

        return interpolated_position

    def generate_tod_flag(
        self,
        celestial_body_names: list[str],
        observatory_pointing: SkyCoord,
        observatory_position: SkyCoord,
        observatory_time: Time,
        proximity_threshold: list[Angle],
    ):

        # Ensure inputs are of the correct types
        if not isinstance(observatory_pointing, SkyCoord):
            raise ValueError("Observatory pointing must be an astropy SkyCoord object.")
        if not isinstance(observatory_position, SkyCoord):
            raise ValueError("Observatory position must be an astropy SkyCoord object.")
        if not isinstance(observatory_time, Time):
            raise ValueError("Observatory time must be an astropy Time object.")

        if observatory_position.name != observatory_pointing.name:
            observatory_position = observatory_position.transform_to(
                observatory_pointing.name
            )

        celestial_body_names = [
            name.casefold() for name in celestial_body_names
        ]  # Normalize names to lowercase for consistent file naming

        proximity_flags = []
        # Looping through each celestial body to generate flags based on proximity to the observatory's pointing
        for body in celestial_body_names:
            if body not in [
                os.path.basename(f).replace("_ephemeris.txt", "")
                for f in self.available_ephemeris_files
            ]:
                raise ValueError(
                    f"""Ephemeris file for {body} not found in directory. 
                    Run ephemeris generator script to generate the required file."""
                )

            ephemeris_file = os.path.join(
                self.ephemeris_file_dir, f"{body}_ephemeris.txt"
            )

            # Interpolating the celestial body's position at the observatory's time using the ephemeris data
            ephemeris_time, ephemeris_position = self.load_ephemeris(ephemeris_file)

            # Ensure all coordinates are in the same frame as the ovservatory pointing for accurate angular separation calculations
            if ephemeris_position.name != observatory_pointing.name:
                ephemeris_position = ephemeris_position.transform_to(
                    observatory_pointing.name
                )

            body_position = self.linearly_interpolate_ephemeris(
                ephemeris_time, ephemeris_position, observatory_time
            )

            # Translating the celestial body's position to the observatory's frame of reference by calculating the relative position
            relative_cartesian_position = CartesianRepresentation(
                body_position.cartesian.xyz.to(u.km)
                - observatory_position.cartesian.xyz.to(u.km)
            )

            # Converting the relative position to spherical coordinates to calculate
            # the apparent longitude and latitude of the celestial body as seen from the observatory
            relative_spherical_position = SphericalRepresentation.from_cartesian(
                relative_cartesian_position
            )

            # Calculating the angular separation between the observatory's pointing
            # and the celestial body's apparent position
            apparent_lon = relative_spherical_position.lon
            apparent_lat = relative_spherical_position.lat

            # Astropy SkyCoord's "longitude" and "latitude" component names depend on the frame,
            # so we need to determine the correct component names for the observatory pointing
            inverse_representation_component_names = {
                v: k
                for k, v in observatory_pointing.representation_component_names.items()
            }

            pointing_lon = getattr(
                observatory_pointing, inverse_representation_component_names["lon"]
            )
            pointing_lat = getattr(
                observatory_pointing, inverse_representation_component_names["lat"]
            )

            separation = angular_separation(
                apparent_lon,
                apparent_lat,
                pointing_lon,
                pointing_lat,
            )

            proximity_flag = (
                separation < proximity_threshold[celestial_body_names.index(body)]
            )

            proximity_flags.append(proximity_flag)

        return proximity_flags


def main():
    # Example usage of the EphemerisFlagger class
    ephemeris_flagger = EphemerisFlagger()

    celestial_bodies = ["Jupiter", "Saturn", "Mars"]
    proximity_thresholds = [Angle(60, unit=u.degree)] * len(celestial_bodies)

    # Example observatory parameters (these would typically come from the TOD data)
    # observatory_pointing = SkyCoord(ra=150 * u.deg, dec=2 * u.deg, frame="icrs")
    observatory_pointing = SkyCoord(
        l=150 * u.deg,
        b=2 * u.deg,
        frame="galactic",
        # lon=150 * u.deg, lat=2 * u.deg, frame="heliocentricmeanecliptic"
    )
    observatory_position = SkyCoord(
        x=1 * u.au,
        y=0 * u.au,
        z=0 * u.au,
        representation_type="cartesian",
        frame="heliocentricmeanecliptic",
    )
    observatory_time = Time("2024-01-01T00:00:00", scale="tdb")

    flags = ephemeris_flagger.generate_tod_flag(
        celestial_body_names=celestial_bodies,
        observatory_pointing=observatory_pointing,
        observatory_position=observatory_position,
        observatory_time=observatory_time,
        proximity_threshold=proximity_thresholds,
    )

    print("Proximity Flags:", flags)


if __name__ == "__main__":
    main()
