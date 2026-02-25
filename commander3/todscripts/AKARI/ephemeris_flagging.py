import glob
import os
from typing import Optional

import astropy.units as u
import numpy as np
import tqdm
from astropy.coordinates import (
    Angle,
    CartesianRepresentation,
    HeliocentricMeanEcliptic,
    SkyCoord,
    SphericalRepresentation,
    angular_separation,
)
from astropy.time import Time


class EphemerisFlagger:
    """Generate flags for celestial bodies based on proximity to observatory pointing.

    This class loads ephemeris data for celestial bodies and generates boolean flags
    indicating when the bodies are within a specified proximity threshold of the
    observatory's pointing direction.
    """

    def __init__(
        self,
        ephemeris_file_dir: Optional[
            str
        ] = "/mn/stornext/d23/cmbco/globe/common/aux/ephemeris",
    ):
        """Initialize EphemerisFlagger with ephemeris data directory.

        Args:
            ephemeris_file_dir: Path to directory containing ephemeris data files.
                Defaults to the common aux ephemeris directory.

        Raises:
            ValueError: If the provided directory path does not exist.
        """
        if not os.path.isdir(ephemeris_file_dir):
            raise ValueError(f"Provided path {ephemeris_file_dir} is not a directory.")

        self.ephemeris_file_dir = ephemeris_file_dir
        self.available_ephemeris_files = glob.glob(
            os.path.join(ephemeris_file_dir, "*.txt")
        )

    def load_ephemeris(self, ephemeris_file: str):
        """Load ephemeris data from file.

        Args:
            ephemeris_file: Path to ephemeris file containing time and position data.

        Returns:
            Tuple[Time, SkyCoord]: Tuple containing:
                - ephemeris time as astropy Time object in TDB scale
                - position as astropy SkyCoord in Heliocentric Mean Ecliptic frame
        """
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
        """Interpolate ephemeris position at a specific time.

        Performs linear interpolation of celestial body position between ephemeris
        data points at the requested time.

        Args:
            ephemeris_time: Array of times from ephemeris data (astropy Time).
            ephemeris_position: Array of positions from ephemeris data (astropy SkyCoord).
            target_time: Time at which to interpolate position (astropy Time).

        Returns:
            SkyCoord: Interpolated position at target time in the same frame as
                ephemeris_position.

        Raises:
            ValueError: If target_time is outside the bounds of ephemeris data.
        """
        if (
            np.min(target_time.tdb.mjd) < ephemeris_time[0].tdb.mjd
            or np.max(target_time.tdb.mjd) > ephemeris_time[-1].tdb.mjd
        ):
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
        """Generate proximity flags for celestial bodies.

        Determines which celestial bodies are within specified proximity thresholds
        of the observatory's pointing direction at a given time.

        Args:
            celestial_body_names: List of celestial body names (e.g., 'Jupiter', 'Saturn').
            observatory_pointing: Observatory's pointing direction (astropy SkyCoord).
            observatory_position: Observatory's position in space (astropy SkyCoord).
            observatory_time: Time of observation (astropy Time).
            proximity_threshold: List of angular separation thresholds for each body
                (astropy Angle, same length as celestial_body_names).

        Returns:
            list: Boolean arrays indicating when each celestial body is within its
                proximity threshold. Each array has same length as input time series.

        Raises:
            ValueError: If input coordinates are not SkyCoord objects or if
                required ephemeris files are not found.
        """
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

        pbar = tqdm.tqdm(
            total=len(celestial_body_names),
            desc="Processing",
            colour="red",
            ncols=80,
            position=0,
            leave=True,
        )

        # Looping through each celestial body to generate flags based on proximity to the observatory's pointing
        for body in celestial_body_names:
            pbar.set_description(f"\033[92mProcessing {body}\033[00m")

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

            # Convert observatory pointing to spherical representation to get lon/lat
            # (needed because the input might be in cartesian representation)
            pointing_spherical = observatory_pointing.spherical
            pointing_lon = pointing_spherical.lon
            pointing_lat = pointing_spherical.lat

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
            pbar.update(1)
        pbar.close()
        return proximity_flags


def main():
    """Demonstrate usage of the EphemerisFlagger class with example data.

    Creates an EphemerisFlagger instance, defines example celestial bodies and
    observatory parameters, and generates proximity flags for the bodies.
    """
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


def main2():
    import h5py
    from pixell import enmap, utils
    import matplotlib.pyplot as plt
    from astropy.coordinates import get_body, solar_system_ephemeris, EarthLocation

    ephemeris_flagger = EphemerisFlagger()

    celestial_bodies = [
        "Jupiter",
    ]
    proximity_thresholds = [Angle(0.5, unit=u.degree)] * len(celestial_bodies)

    observatory_time = Time(58755.028530304335, format="mjd", scale="tdb")
    observatory_pointing = SkyCoord(
        ra=[256.82296228613416, 0, 90],
        dec=[-22.609898786867273, 76, 43],
        unit=u.deg,
        frame="icrs",
    )

    loc = EarthLocation.of_site("ovro")
    with solar_system_ephemeris.set("builtin"):
        earth = get_body("earth", observatory_time, loc)
    observatory_position = earth.transform_to("heliocentricmeanecliptic")

    flags = ephemeris_flagger.generate_tod_flag(
        celestial_body_names=celestial_bodies,
        observatory_pointing=observatory_pointing,
        observatory_position=observatory_position,
        observatory_time=observatory_time,
        proximity_threshold=proximity_thresholds,
    )

    flags = np.array(flags)
    print("Proximity Flags:", flags)


if __name__ == "__main__":
    main()
    main2()
