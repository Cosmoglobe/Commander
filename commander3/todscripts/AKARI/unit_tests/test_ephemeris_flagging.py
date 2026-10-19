"""Test suite for ephemeris_flagging.py module.

This test suite includes unit tests for the EphemerisFlagger class and its methods.

Usage Examples
--------------

Run all tests:
    pytest test_ephemeris_flagging.py -v

Run specific test class:
    pytest test_ephemeris_flagging.py::TestEphemerisFlaggerInitialization -v

Run specific test:
    pytest test_ephemeris_flagging.py::TestEphemerisFlaggerInitialization::test_valid_directory -v

Run tests with detailed output:
    pytest test_ephemeris_flagging.py -vv -s

Run with coverage report:
    pytest test_ephemeris_flagging.py --cov=ephemeris_flagging --cov-report=html
"""

import sys
import os
import pytest
import tempfile
import numpy as np
from unittest.mock import Mock, patch, MagicMock
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, Angle, BarycentricMeanEcliptic

# Add parent directory to path so we can import ephemeris_flagging
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ephemeris_flagging import EphemerisFlagger


# ============================================================================
# Initialization Tests
# ============================================================================


class TestEphemerisFlaggerInitialization:
    """Test EphemerisFlagger initialization and setup."""

    def test_valid_directory_initialization(self):
        """Test initialization with a valid directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create some dummy ephemeris files
            open(os.path.join(tmpdir, "jupiter_ephemeris.txt"), "w").close()
            open(os.path.join(tmpdir, "saturn_ephemeris.txt"), "w").close()

            flagger = EphemerisFlagger(ephemeris_file_dir=tmpdir)

            assert flagger.ephemeris_file_dir == tmpdir
            assert len(flagger.available_ephemeris_files) == 2

    def test_invalid_directory_raises_error(self):
        """Test that invalid directory path raises ValueError."""
        with pytest.raises(ValueError, match="is not a directory"):
            EphemerisFlagger(ephemeris_file_dir="/nonexistent/path/to/directory")

    def test_empty_directory_initialization(self):
        """Test initialization with an empty directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            flagger = EphemerisFlagger(ephemeris_file_dir=tmpdir)

            assert flagger.ephemeris_file_dir == tmpdir
            assert len(flagger.available_ephemeris_files) == 0

    def test_multiple_ephemeris_files_detected(self):
        """Test that multiple ephemeris files are correctly detected."""
        with tempfile.TemporaryDirectory() as tmpdir:
            bodies = ["Jupiter", "Saturn", "Mars", "Mercury", "Venus"]
            for body in bodies:
                open(os.path.join(tmpdir, f"{body.lower()}_ephemeris.txt"), "w").close()

            flagger = EphemerisFlagger(ephemeris_file_dir=tmpdir)

            assert len(flagger.available_ephemeris_files) == 5


# ============================================================================
# Ephemeris Loading Tests
# ============================================================================


class TestEphemerisLoading:
    """Test ephemeris file loading and parsing."""

    def create_sample_ephemeris_file(self, filepath, num_points=10):
        """Create a sample ephemeris file for testing.

        Args:
            filepath: Path where to create the file.
            num_points: Number of data points to create.
        """
        mjd_values = np.linspace(51544.0, 51553.0, num_points)  # Jan 1 to Jan 10, 2000
        x_values = np.linspace(0.5, 1.5, num_points)  # AU
        y_values = np.linspace(-0.2, 0.2, num_points)  # AU
        z_values = np.linspace(-0.1, 0.1, num_points)  # AU

        data = np.column_stack((mjd_values, x_values, y_values, z_values))

        # Write with 2-line header (as expected by the loader)
        with open(filepath, "w") as f:
            f.write("Header line 1\n")
            f.write("Header line 2\n")
            np.savetxt(f, data)

    def test_load_ephemeris_valid_file(self):
        """Test loading a valid ephemeris file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ephemeris_file = os.path.join(tmpdir, "jupiter_ephemeris.txt")
            self.create_sample_ephemeris_file(ephemeris_file)

            flagger = EphemerisFlagger(ephemeris_file_dir=tmpdir)
            mjd, pos = flagger.load_ephemeris(ephemeris_file)

            # Check returned types
            assert isinstance(mjd, Time)
            assert isinstance(pos, SkyCoord)

            # Check data integrity
            assert len(mjd) == 10
            assert len(pos) == 10

            # Check coordinate frame
            assert pos.frame.name == "heliocentricmeanecliptic"

    def test_load_ephemeris_time_format(self):
        """Test that ephemeris time is in correct format (MJD, TDB scale)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ephemeris_file = os.path.join(tmpdir, "mars_ephemeris.txt")
            self.create_sample_ephemeris_file(ephemeris_file, num_points=5)

            flagger = EphemerisFlagger(ephemeris_file_dir=tmpdir)
            mjd, _ = flagger.load_ephemeris(ephemeris_file)

            assert mjd.scale == "tdb"
            # Check approximate values
            assert np.allclose(mjd.mjd[0], 51544.0, atol=0.01)
            assert np.allclose(mjd.mjd[-1], 51553.0, atol=0.01)

    def test_load_ephemeris_position_units(self):
        """Test that ephemeris positions are in correct units (AU)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ephemeris_file = os.path.join(tmpdir, "venus_ephemeris.txt")
            self.create_sample_ephemeris_file(ephemeris_file, num_points=3)

            flagger = EphemerisFlagger(ephemeris_file_dir=tmpdir)
            _, pos = flagger.load_ephemeris(ephemeris_file)

            # Cartesian representation should be in AU
            x_au = pos.cartesian.x[0].to(u.au).value
            assert 0.4 < x_au < 0.6  # Should be ~0.5 AU


# ============================================================================
# Interpolation Tests
# ============================================================================


class TestEphemerisInterpolation:
    """Test ephemeris position interpolation."""

    def create_flagger_with_sample_data(self, tmpdir, num_points=5):
        """Create a flagger and return it with sample ephemeris data."""
        ephemeris_file = os.path.join(tmpdir, "jupiter_ephemeris.txt")

        # Create predictable data for interpolation
        mjd_values = np.array([51544.0, 51545.0, 51546.0, 51547.0, 51548.0])
        x_values = np.array([1.0, 1.5, 2.0, 2.5, 3.0])  # Linear increase
        y_values = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
        z_values = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

        data = np.column_stack((mjd_values, x_values, y_values, z_values))

        with open(ephemeris_file, "w") as f:
            f.write("Header line 1\n")
            f.write("Header line 2\n")
            np.savetxt(f, data)

        flagger = EphemerisFlagger(ephemeris_file_dir=tmpdir)
        return flagger, ephemeris_file

    def test_interpolation_at_exact_time_point(self):
        """Test interpolation at times that match ephemeris data points."""
        with tempfile.TemporaryDirectory() as tmpdir:
            flagger, ephemeris_file = self.create_flagger_with_sample_data(tmpdir)
            mjd, pos = flagger.load_ephemeris(ephemeris_file)

            # Interpolate at an exact time point
            target_time = Time(51545.0, format="mjd", scale="tdb")
            interpolated = flagger.linearly_interpolate_ephemeris(mjd, pos, target_time)

            assert isinstance(interpolated, SkyCoord)
            # Should match the exact value
            assert np.allclose(
                interpolated.cartesian.x.to(u.km).value, 1.5 * u.au.to(u.km)
            )

    def test_interpolation_between_time_points(self):
        """Test interpolation at times between ephemeris data points."""
        with tempfile.TemporaryDirectory() as tmpdir:
            flagger, ephemeris_file = self.create_flagger_with_sample_data(tmpdir)
            mjd, pos = flagger.load_ephemeris(ephemeris_file)

            # Interpolate at midpoint between two data points
            target_time = Time(51544.5, format="mjd", scale="tdb")
            interpolated = flagger.linearly_interpolate_ephemeris(mjd, pos, target_time)

            # Should be approximately halfway between 1.0 and 1.5
            x_value = interpolated.cartesian.x.to(u.km).value
            expected = 1.25 * u.au.to(u.km)
            assert np.allclose(x_value, expected, rtol=0.01)

    def test_interpolation_out_of_bounds_raises_error(self):
        """Test that interpolation outside data bounds raises ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            flagger, ephemeris_file = self.create_flagger_with_sample_data(tmpdir)
            mjd, pos = flagger.load_ephemeris(ephemeris_file)

            # Try to interpolate before the data range
            target_time = Time(51543.0, format="mjd", scale="tdb")
            with pytest.raises(ValueError, match="out of bounds"):
                flagger.linearly_interpolate_ephemeris(mjd, pos, target_time)

            # Try to interpolate after the data range
            target_time = Time(51549.0, format="mjd", scale="tdb")
            with pytest.raises(ValueError, match="out of bounds"):
                flagger.linearly_interpolate_ephemeris(mjd, pos, target_time)

    def test_interpolation_preserves_frame(self):
        """Test that interpolation preserves the coordinate frame."""
        with tempfile.TemporaryDirectory() as tmpdir:
            flagger, ephemeris_file = self.create_flagger_with_sample_data(tmpdir)
            mjd, pos = flagger.load_ephemeris(ephemeris_file)

            target_time = Time(51545.5, format="mjd", scale="tdb")
            interpolated = flagger.linearly_interpolate_ephemeris(mjd, pos, target_time)

            assert interpolated.frame.name == pos.frame.name
            assert interpolated.frame.name == "heliocentricmeanecliptic"


# ============================================================================
# Flag Generation Tests
# ============================================================================


class TestFlagGeneration:
    """Test proximity flag generation."""

    def create_test_setup(self, tmpdir):
        """Create test files and objects for flag generation tests."""
        # Create dummy ephemeris file for Jupiter
        ephemeris_file = os.path.join(tmpdir, "jupiter_ephemeris.txt")

        # Simple constant position data
        mjd_values = np.linspace(51544.0, 51545.0, 5)
        x_values = np.ones(5) * 1.0  # Constant 1 AU
        y_values = np.zeros(5)
        z_values = np.zeros(5)

        data = np.column_stack((mjd_values, x_values, y_values, z_values))

        with open(ephemeris_file, "w") as f:
            f.write("Header line 1\n")
            f.write("Header line 2\n")
            np.savetxt(f, data)

        flagger = EphemerisFlagger(ephemeris_file_dir=tmpdir)
        return flagger

    def test_generate_tod_flag_basic(self):
        """Test basic flag generation with simple inputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            flagger = self.create_test_setup(tmpdir)

            # Simple test case - use spherical coordinates (ra, dec) for pointing
            observatory_pointing = SkyCoord(
                ra=0 * u.deg,
                dec=0 * u.deg,
                frame="icrs",
            )
            observatory_position = SkyCoord(
                x=0 * u.au,
                y=0 * u.au,
                z=0 * u.au,
                representation_type="cartesian",
                frame="heliocentricmeanecliptic",
            )
            observatory_time = Time(51544.5, format="mjd", scale="tdb")
            proximity_threshold = [Angle(90, unit=u.degree)]

            flags = flagger.generate_tod_flag(
                celestial_body_names=["jupiter"],
                observatory_pointing=observatory_pointing,
                observatory_position=observatory_position,
                observatory_time=observatory_time,
                proximity_threshold=proximity_threshold,
            )

            assert isinstance(flags, list)
            assert len(flags) == 1
            assert isinstance(flags[0], (bool, np.bool_))

    def test_generate_tod_flag_invalid_pointing_type(self):
        """Test that invalid pointing object type raises ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            flagger = self.create_test_setup(tmpdir)

            with pytest.raises(ValueError, match="Observatory pointing"):
                flagger.generate_tod_flag(
                    celestial_body_names=["jupiter"],
                    observatory_pointing="invalid",  # Invalid type
                    observatory_position=SkyCoord(
                        x=0 * u.au,
                        y=0 * u.au,
                        z=0 * u.au,
                        representation_type="cartesian",
                        frame="heliocentricmeanecliptic",
                    ),
                    observatory_time=Time(51544.5, format="mjd", scale="tdb"),
                    proximity_threshold=[Angle(90, unit=u.degree)],
                )

    def test_generate_tod_flag_invalid_position_type(self):
        """Test that invalid position object type raises ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            flagger = self.create_test_setup(tmpdir)

            with pytest.raises(ValueError, match="Observatory position"):
                flagger.generate_tod_flag(
                    celestial_body_names=["jupiter"],
                    observatory_pointing=SkyCoord(
                        ra=0 * u.deg,
                        dec=0 * u.deg,
                        frame="icrs",
                    ),
                    observatory_position="invalid",  # Invalid type
                    observatory_time=Time(51544.5, format="mjd", scale="tdb"),
                    proximity_threshold=[Angle(90, unit=u.degree)],
                )

    def test_generate_tod_flag_invalid_time_type(self):
        """Test that invalid time object type raises ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            flagger = self.create_test_setup(tmpdir)

            with pytest.raises(ValueError, match="Observatory time"):
                flagger.generate_tod_flag(
                    celestial_body_names=["jupiter"],
                    observatory_pointing=SkyCoord(
                        ra=0 * u.deg,
                        dec=0 * u.deg,
                        frame="icrs",
                    ),
                    observatory_position=SkyCoord(
                        x=0 * u.au,
                        y=0 * u.au,
                        z=0 * u.au,
                        representation_type="cartesian",
                        frame="heliocentricmeanecliptic",
                    ),
                    observatory_time="2024-01-01",  # Invalid type
                    proximity_threshold=[Angle(90, unit=u.degree)],
                )

    def test_generate_tod_flag_missing_ephemeris_file(self):
        """Test that missing ephemeris file raises ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            flagger = EphemerisFlagger(ephemeris_file_dir=tmpdir)

            with pytest.raises(ValueError, match="Ephemeris file"):
                flagger.generate_tod_flag(
                    celestial_body_names=["nonexistent_body"],
                    observatory_pointing=SkyCoord(
                        ra=0 * u.deg,
                        dec=0 * u.deg,
                        frame="icrs",
                    ),
                    observatory_position=SkyCoord(
                        x=0 * u.au,
                        y=0 * u.au,
                        z=0 * u.au,
                        representation_type="cartesian",
                        frame="heliocentricmeanecliptic",
                    ),
                    observatory_time=Time(51544.5, format="mjd", scale="tdb"),
                    proximity_threshold=[Angle(90, unit=u.degree)],
                )

    def test_generate_tod_flag_case_insensitive_body_names(self):
        """Test that body names are case-insensitive."""
        with tempfile.TemporaryDirectory() as tmpdir:
            flagger = self.create_test_setup(tmpdir)

            # Should work with uppercase
            flags = flagger.generate_tod_flag(
                celestial_body_names=["JUPITER"],
                observatory_pointing=SkyCoord(
                    ra=0 * u.deg,
                    dec=0 * u.deg,
                    frame="icrs",
                ),
                observatory_position=SkyCoord(
                    x=0 * u.au,
                    y=0 * u.au,
                    z=0 * u.au,
                    representation_type="cartesian",
                    frame="heliocentricmeanecliptic",
                ),
                observatory_time=Time(51544.5, format="mjd", scale="tdb"),
                proximity_threshold=[Angle(90, unit=u.degree)],
            )

            assert len(flags) == 1


# ============================================================================
# Integration Tests
# ============================================================================


class TestIntegration:
    """Integration tests combining multiple components."""

    def test_full_workflow(self):
        """Test complete workflow from initialization to flag generation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create ephemeris files
            for body_name in ["jupiter", "saturn"]:
                ephemeris_file = os.path.join(tmpdir, f"{body_name}_ephemeris.txt")

                mjd_values = np.linspace(51544.0, 51545.0, 10)
                x_values = np.linspace(1.0, 1.5, 10)
                y_values = np.linspace(-0.2, 0.2, 10)
                z_values = np.linspace(-0.1, 0.1, 10)

                data = np.column_stack((mjd_values, x_values, y_values, z_values))

                with open(ephemeris_file, "w") as f:
                    f.write("Header line 1\n")
                    f.write("Header line 2\n")
                    np.savetxt(f, data)

            # Initialize flagger
            flagger = EphemerisFlagger(ephemeris_file_dir=tmpdir)

            # Load data
            jupiter_time, jupiter_pos = flagger.load_ephemeris(
                os.path.join(tmpdir, "jupiter_ephemeris.txt")
            )
            saturn_time, saturn_pos = flagger.load_ephemeris(
                os.path.join(tmpdir, "saturn_ephemeris.txt")
            )

            # Interpolate
            target_time = Time(51544.5, format="mjd", scale="tdb")
            jupiter_interp = flagger.linearly_interpolate_ephemeris(
                jupiter_time, jupiter_pos, target_time
            )
            saturn_interp = flagger.linearly_interpolate_ephemeris(
                saturn_time, saturn_pos, target_time
            )

            assert isinstance(jupiter_interp, SkyCoord)
            assert isinstance(saturn_interp, SkyCoord)

            # Generate flags
            flags = flagger.generate_tod_flag(
                celestial_body_names=["jupiter", "saturn"],
                observatory_pointing=SkyCoord(
                    ra=0 * u.deg,
                    dec=0 * u.deg,
                    frame="icrs",
                ),
                observatory_position=SkyCoord(
                    x=0 * u.au,
                    y=0 * u.au,
                    z=0 * u.au,
                    representation_type="cartesian",
                    frame="heliocentricmeanecliptic",
                ),
                observatory_time=target_time,
                proximity_threshold=[
                    Angle(90, unit=u.degree),
                    Angle(90, unit=u.degree),
                ],
            )

            assert len(flags) == 2

    def test_coordinate_frame_invariance(self):
        """Testing that known pointing and time yield correct flags
        and that input coordinate system do not change the flags.
        """
        import h5py
        from pixell import enmap, utils
        import matplotlib.pyplot as plt
        from astropy.coordinates import get_body, solar_system_ephemeris, EarthLocation

        ephemeris_flagger = EphemerisFlagger()

        celestial_bodies = [
            "Jupiter",
        ]
        proximity_thresholds = [Angle(0.5, unit=u.degree)] * len(celestial_bodies)

        # Using a known time and pointing where Jupiter should be within the proximity threshold
        # as well as two pointings that should not
        observatory_time = Time(58755.028530304335, format="mjd", scale="tdb")
        observatory_pointing = SkyCoord(
            ra=[256.82296228613416, 0, 90],
            dec=[-22.609898786867273, 76, 43],
            unit=u.deg,
            frame="icrs",
        )

        # Get the observatory (COMAP telescope) position in heliocentric mean ecliptic frame
        loc = EarthLocation.of_site("ovro")
        with solar_system_ephemeris.set("builtin"):
            earth = get_body("earth", observatory_time, loc)
        observatory_position = earth.transform_to("heliocentricmeanecliptic")

        # Generate test flags, only the first should be True
        flags = ephemeris_flagger.generate_tod_flag(
            celestial_body_names=celestial_bodies,
            observatory_pointing=observatory_pointing,
            observatory_position=observatory_position,
            observatory_time=observatory_time,
            proximity_threshold=proximity_thresholds,
        )

        # Assert that the first pointing is flagged as True (Jupiter within threshold) and the others are False
        flags = np.array(flags)
        assert flags[
            0, 0
        ], "Jupiter should be within the proximity threshold at this time and pointing."
        assert not np.any(
            flags[0, 1:]
        ), "Only Jupiter should be within the proximity threshold at this time and pointing."

        # Check that different input frames yeild the same flags (since the method should handle transformations internally)
        frames = [
            "galactic",
            "fk4",
            "fk5",
            "barycentricmeanecliptic",
        ]
        reference_flags = flags.copy()  # Store the original flags for comparison
        for frame in frames:
            transformed_pointing = observatory_pointing.transform_to(frame)
            flags = ephemeris_flagger.generate_tod_flag(
                celestial_body_names=celestial_bodies,
                observatory_pointing=transformed_pointing,
                observatory_position=observatory_position,
                observatory_time=observatory_time,
                proximity_threshold=proximity_thresholds,
            )
            assert np.array_equal(
                flags, reference_flags
            ), f"Flags should be the same regardless of the frame used for the observatory pointing. Discrepancy found in {frame} frame."
