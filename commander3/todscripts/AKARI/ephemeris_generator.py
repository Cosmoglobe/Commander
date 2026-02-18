from dataclasses import dataclass
import os
from typing import Optional
import warnings
import argparse

import tqdm
import h5py
import numpy as np
import astropy.units as u
from astropy.coordinates import (
    HeliocentricMeanEcliptic,
    get_body,
    get_body_barycentric,
    solar_system_ephemeris,
    SkyCoord,
)
from astroquery.jplhorizons import Horizons
from astropy.time import Time

warnings.filterwarnings("ignore")


@dataclass
class CelestialBody:
    """Data class for celestial body data."""

    name: str
    time: np.ndarray
    position: np.ndarray


class EphemerisGenerator:
    def __init__(
        self,
        body_names: list[str],
        start_time: Optional[u.quantity.Quantity] = Time(
            "1980-01-01 00:00:00.000", format="iso", scale="utc"
        ).mjd
        * u.day,
        end_time: Optional[u.quantity.Quantity] = Time(
            "2049-12-31 23:59:59.999", format="iso", scale="utc"
        ).mjd
        * u.day,
        time_step: Optional[u.quantity.Quantity] = 3600 * u.s,
    ) -> None:
        """Initialize ephemeris generation parameters and time grid.

        Parameters
        ----------
        body_names : list of str
            List of solar system body names understood by Astropy.
        start_mjd : astropy.units.Quantity, optional
            Start time in MJD with units (default 1980-01-01 UTC).
        end_mjd : astropy.units.Quantity, optional
            End time in MJD with units (default 2049-12-31 UTC).
        time_step : astropy.units.Quantity, optional
            Time step as an Astropy quantity (default 3600 s).
        """
        self.body_names = body_names
        self.start_time = start_time
        self.end_time = end_time
        self.time_step = time_step
        self.time_arr = (
            Time(self.start_time, format="mjd", scale="tai")
            + np.arange(
                0,
                (self.end_time - self.start_time).to(u.min).value,
                int(np.round(self.time_step.to(u.min).value)),
            )
            * u.min
        )

        self.jpl_horizons_query_limit = 90024
        self.number_of_query_batches = int(
            np.ceil(self.time_arr.size / self.jpl_horizons_query_limit)
        )
        self.batch_size = int(
            np.ceil(self.time_arr.size / self.number_of_query_batches)
        )

        print(
            f"\nTotal time steps: {self.time_arr.size}, which will require {self.number_of_query_batches} batch(es) of queries to JPL Horizons."
        )

    def define_bodies(self, time_arrs: Optional[list[Time]] = None) -> None:
        """Create internal `CelestialBody` entries with associated time arrays.

        Parameters
        ----------
        time_arrs : list of astropy.time.Time, optional
            Optional list of time arrays to assign per body. If None, uses the
            generator-wide `self.time_arr` for all bodies.
        """

        self.bodies = []

        print("Defining celestial bodies and associated time arrays...")
        if time_arrs is not None:
            for body, time_arr in zip(self.body_names, time_arrs):
                self.bodies.append(
                    CelestialBody(name=body, time=time_arr, position=None)
                )
        else:
            for body in self.body_names:
                self.bodies.append(
                    CelestialBody(name=body, time=self.time_arr, position=None)
                )

    def _define_major_and_minor_bodies(self):
        """Build lookup tables for major and minor body identifiers.

        Returns
        -------
        major_bodies : dict
            Mapping of major body names (casefolded) to their JPL IDs.
        minor_bodies : dict
            Mapping of minor body names (casefolded) to their identifiers as
            defined in the asteroid orbital elements file.
        """
        major_bodies = {
            "sun": 10,
            "moon": 301,
            "mercury": 199,
            "venus": 299,
            "earth": 399,
            "mars": 499,
            "jupiter": 599,
            "saturn": 699,
            "uranus": 799,
            "neptune": 899,
            "pluto": 999,
        }

        orbital_element_path = "/mn/stornext/d23/cmbco/globe/common/aux/"
        path2asteroid_def = os.path.join(
            orbital_element_path, "asteroid_orbital_elements.txt"
        )

        minor_bodies = {}
        with open(path2asteroid_def, "r") as f:
            asteroid_lines = f.readlines()
            column_widths = []
            for column in asteroid_lines[1].split():
                column_widths.append(len(column))
            for line in asteroid_lines[2:]:
                previous_width = 0
                columns = []
                for width in column_widths:
                    columns.append(
                        line[previous_width : previous_width + width + 1].strip()
                    )
                    previous_width = previous_width + width + 1
                minor_bodies[columns[1].casefold()] = columns[0]

        return major_bodies, minor_bodies

    def generate_ephemeris(self):
        """Populate body positions in the heliocentric mean ecliptic frame."""
        print("Generating ephemeris for celestial bodies...")

        major_bodies, minor_bodies = self._define_major_and_minor_bodies()

        pbar = tqdm.tqdm(
            total=len(self.bodies),
            desc="Processing",
            colour="blue",
            ncols=80,
            position=0,
            leave=True,
        )
        pbar_batch = tqdm.tqdm(
            total=self.number_of_query_batches,
            desc="Querying JPL Horizons",
            colour="#cd44dd",
            leave=False,
            ncols=50,
            position=1,
        )
        for body in self.bodies:
            pbar.set_description(f"\033[92mProcessing {body.name}\033[00m")

            if body.name.casefold() in major_bodies:
                id = f"{major_bodies[body.name.casefold()]}"
            elif body.name.casefold() in minor_bodies:
                id = f"{body.name}"
            else:
                print(
                    f"\n\033[91mWarning:\033[00m Body name '{body.name}' not found in major or minor body lists. Attempting to query JPL Horizons with the name directly.\n"
                )
                id = f"{body.name}"

            body.position = np.zeros((3, body.time.size))
            body.time_mjd = np.zeros(body.time.size)

            pbar_batch.reset()
            for batch in range(self.number_of_query_batches):
                pbar_batch.set_description(f"\033[93mBatch\033[00m")
                batch_start = batch * self.batch_size
                batch_end = min((batch + 1) * self.batch_size, body.time.size)

                # Query JPL Horizons for the current batch of times for the current body
                obj = Horizons(
                    id=f"{id}",  # Celestial body ID
                    location="@10",  # At center of Sun
                    epochs={
                        "start": f"{body.time[batch_start].iso}",
                        "stop": f"{body.time[batch_end-1].iso}",
                        "step": f"{int(np.round(self.time_step.to(u.min).value))}m",
                    },
                )

                vec = obj.vectors()

                batch_time = body.time[batch_start:batch_end]
                batch_time_jd = Time(vec["datetime_jd"].data, format="jd", scale="tai")

                # Define the input coordinates in the original frame (assumed to be
                # heliocentric ecliptic) and transform to Heliocentric Mean Ecliptic frame
                coordinates = SkyCoord(
                    vec["x"].quantity,
                    vec["y"].quantity,
                    vec["z"].quantity,
                    representation_type="cartesian",
                    frame="heliocentriceclipticiau76",
                    obstime=batch_time,
                ).transform_to(HeliocentricMeanEcliptic)

                # Writing to buffer arrays for the current batch
                body.position[:, batch_start:batch_end] = coordinates.cartesian.xyz.to(
                    u.au
                ).value
                body.time_mjd[batch_start:batch_end] = batch_time_jd.mjd

                pbar_batch.update(1)

            # Overwriting body.time with the MJD times from JPL Horizons
            # to ensure consistency, since the original time grid may not
            # align perfectly with the queried times due to rounding of
            # the time step.
            body.time = Time(body.time_mjd, format="mjd", scale="tai")

            pbar.update(1)
        pbar_batch.close()
        pbar.close()

    def save_ephemeris(self, outpath: str, save_format: str = "hdf5") -> None:
        """Save ephemeris to disk in text or HDF5 format.

        Parameters
        ----------
        outpath : str
            Output directory path.
        save_format : str, optional
            Output format, either "txt" or "hdf5".
        """
        if save_format == "txt":
            print(f"Saving ephemeris to {outpath} in txt format...")
            for body in self.bodies:
                body_file_data = np.vstack((body.time.mjd, body.position)).T
                header = f"{len(body.time)}\nMJD\tX[AU]\tY[AU]\tZ[AU] \t Frame: Heliocentric Mean Ecliptic"
                np.savetxt(
                    os.path.join(outpath, f"{body.name}_ephemeris.txt"),
                    body_file_data,
                    header=header,
                )

        if save_format == "hdf5":
            print(f"Saving ephemeris to {outpath} in HDF5 format...")
            with h5py.File(os.path.join(outpath, "ephemeris.h5"), "w") as f:
                for body in self.bodies:
                    grp = f.create_group(body.name)
                    grp.create_dataset("time_mjd", data=body.time.mjd)
                    grp.create_dataset("position_au", data=body.position)
                    grp["position_au"].attrs["units"] = "au"
                    grp["position_au"].attrs["frame"] = "heliocentric_mean_ecliptic"


def parse_commandline_args():
    """Parse command-line arguments for ephemeris generation."""

    parser = argparse.ArgumentParser(
        description="Generate ephemeris for solar system bodies and save to disk."
    )
    parser.add_argument(
        "--outpath",
        type=str,
        help="Output directory for ephemeris files (required)",
        required=True,
    )
    parser.add_argument(
        "--format",
        type=str,
        choices=["txt", "h5"],
        default="hdf5",
        help="Output format for ephemeris files (default: hdf5)",
    )
    parser.add_argument(
        "--time_step",
        type=float,
        default=3600,
        help="Time step in seconds for ephemeris generation (default: 3600 seconds)",
    )
    parser.add_argument(
        "--start_time",
        type=str,
        default="1980-01-01 00:00:00",
        help="Start time for ephemeris generation in ISO format (default: 1980-01-01 00:00:00)",
    )
    parser.add_argument(
        "--end_time",
        type=str,
        default="2049-12-31 23:59:59",
        help="End time for ephemeris generation in ISO format (default: 2049-12-31 23:59:59)",
    )
    parser.add_argument(
        "--bodies",
        nargs="+",
        default=[
            "sun",
            "moon",
            "mercury",
            "venus",
            "earth",
            "mars",
            "jupiter",
            "saturn",
            "uranus",
            "neptune",
        ],
        help="List of solar system bodies to include (default: sun moon mercury venus earth mars jupiter saturn uranus neptune)",
    )

    return parser.parse_args()


def main() -> None:

    args = parse_commandline_args()

    start_time = Time(args.start_time, format="iso", scale="utc").mjd * u.day
    end_time = Time(args.end_time, format="iso", scale="utc").mjd * u.day

    ephem_gen = EphemerisGenerator(
        body_names=args.bodies,
        time_step=args.time_step * u.s,
        start_time=start_time,
        end_time=end_time,
    )

    ephem_gen.define_bodies()
    ephem_gen.generate_ephemeris()
    ephem_gen.save_ephemeris(outpath=args.outpath, save_format=args.format)


if __name__ == "__main__":
    main()
