"""Test suite for ephemeris_generator.py module.

This test suite includes both fast mocked tests and slow integration tests that
perform actual JPL Horizons queries.

Usage Examples
--------------

Run all tests (fast mocked tests only):
    pytest test_ephemeris_generator.py -v

Run all tests including slow integration tests:
    pytest test_ephemeris_generator.py -v -m slow

Run only real JPL Horizons query tests:
    pytest test_ephemeris_generator.py -m real_query -v

Run specific test class:
    pytest test_ephemeris_generator.py::TestArgumentParsing -v

Run specific test:
    pytest test_ephemeris_generator.py::TestIntegration::test_sun_position_at_origin -v

Run tests and show detailed output:
    pytest test_ephemeris_generator.py -vv -s

Skip slow tests explicitly:
    pytest test_ephemeris_generator.py -v -m "not slow"

Run with coverage report:
    pytest test_ephemeris_generator.py --cov=ephemeris_generator --cov-report=html

Available markers:
    - slow: Tests that take longer to run (integration tests)
    - real_query: Tests that make actual JPL Horizons API calls
"""

import sys
import os
import pytest
import argparse
import tempfile
import numpy as np
import h5py
from unittest.mock import Mock, patch, MagicMock
import astropy.units as u
from astropy.time import Time

# Add parent directory to path so we can import ephemeris_generator
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ephemeris_generator import (
    EphemerisGenerator,
    CelestialBody,
    parse_commandline_args,
    main,
)


# ============================================================================
# Argument Parsing Tests
# ============================================================================


class TestArgumentParsing:
    """Test command-line argument parsing."""

    def test_mutually_exclusive_start_time_args(self):
        """Test that ISO and MJD start time formats cannot be used simultaneously."""
        with pytest.raises(SystemExit):
            with patch(
                "sys.argv",
                [
                    "ephemeris_generator.py",
                    "--outpath",
                    "/tmp/test",
                    "--start_time",
                    "2000-01-01 00:00:00",
                    "--start_mjd",
                    "51544.0",
                ],
            ):
                parse_commandline_args()

    def test_mutually_exclusive_end_time_args(self):
        """Test that ISO and MJD end time formats cannot be used simultaneously."""
        with pytest.raises(SystemExit):
            with patch(
                "sys.argv",
                [
                    "ephemeris_generator.py",
                    "--outpath",
                    "/tmp/test",
                    "--end_time",
                    "2000-01-01 00:00:00",
                    "--end_mjd",
                    "51544.0",
                ],
            ):
                parse_commandline_args()

    def test_required_outpath_argument(self):
        """Test that outpath argument is required."""
        with pytest.raises(SystemExit):
            with patch("sys.argv", ["ephemeris_generator.py"]):
                parse_commandline_args()

    def test_valid_argument_parsing(self):
        """Test valid argument combinations are accepted."""
        with patch(
            "sys.argv",
            [
                "ephemeris_generator.py",
                "--outpath",
                "/tmp/test",
                "--start_time",
                "2000-01-01 00:00:00",
                "--end_time",
                "2000-01-02 00:00:00",
                "--bodies",
                "sun",
                "earth",
            ],
        ):
            args = parse_commandline_args()
            assert args.outpath == "/tmp/test"
            assert args.start_time == "2000-01-01 00:00:00"
            assert args.end_time == "2000-01-02 00:00:00"
            assert args.bodies == ["sun", "earth"]

    def test_mjd_time_arguments(self):
        """Test MJD time arguments are parsed correctly."""
        with patch(
            "sys.argv",
            [
                "ephemeris_generator.py",
                "--outpath",
                "/tmp/test",
                "--start_mjd",
                "51544.0",
                "--end_mjd",
                "51545.0",
            ],
        ):
            args = parse_commandline_args()
            assert args.start_mjd == 51544.0
            assert args.end_mjd == 51545.0


# ============================================================================
# Time Handling Tests
# ============================================================================


class TestTimeHandling:
    """Test time array generation and batch calculations."""

    def test_time_array_generation_hourly(self):
        """Test time array generation with 1-hour steps."""
        start = Time("2000-01-01 00:00:00", format="iso", scale="tdb").mjd * u.day
        end = Time("2000-01-02 00:00:00", format="iso", scale="tdb").mjd * u.day
        time_step = 3600 * u.s

        gen = EphemerisGenerator(
            body_names=["sun"], start_time=start, end_time=end, time_step=time_step
        )

        # Should have 24 time steps (0-23 hours inclusive = 25 points, but depends on implementation)
        # Let's check the array exists and has reasonable length
        assert len(gen.time_arr) > 0
        assert len(gen.time_arr) <= 25  # At most 25 for 24 hours + endpoint

    def test_time_array_generation_daily(self):
        """Test time array generation with 1-day steps."""
        start = Time("2000-01-01 00:00:00", format="iso", scale="tdb").mjd * u.day
        end = Time("2000-01-10 00:00:00", format="iso", scale="tdb").mjd * u.day
        time_step = 86400 * u.s

        gen = EphemerisGenerator(
            body_names=["sun"], start_time=start, end_time=end, time_step=time_step
        )

        assert len(gen.time_arr) >= 9
        assert len(gen.time_arr) <= 11

    def test_batch_calculation_small_array(self):
        """Test batch calculation for arrays smaller than JPL limit."""
        start = Time("2000-01-01 00:00:00", format="iso", scale="tdb").mjd * u.day
        end = Time("2000-01-02 00:00:00", format="iso", scale="tdb").mjd * u.day
        time_step = 3600 * u.s

        gen = EphemerisGenerator(
            body_names=["sun"], start_time=start, end_time=end, time_step=time_step
        )

        # Should only need 1 batch for small time range
        assert gen.number_of_query_batches == 1
        assert gen.batch_size >= len(gen.time_arr)

    def test_batch_calculation_large_array(self):
        """Test batch calculation for arrays requiring multiple batches."""
        # Create generator with more than 90024 time steps
        start = Time("2000-01-01 00:00:00", format="iso", scale="tdb").mjd * u.day
        end = Time("2050-01-01 00:00:00", format="iso", scale="tdb").mjd * u.day
        time_step = 3600 * u.s

        gen = EphemerisGenerator(
            body_names=["sun"], start_time=start, end_time=end, time_step=time_step
        )

        # Should require multiple batches
        assert gen.number_of_query_batches > 1
        # Verify batches cover all time steps
        total_coverage = gen.batch_size * gen.number_of_query_batches
        assert total_coverage >= len(gen.time_arr)


# ============================================================================
# Body Definition Tests
# ============================================================================


class TestBodyDefinition:
    """Test celestial body recognition and definition."""

    def test_major_body_identification(self):
        """Test major solar system bodies are correctly identified."""
        gen = EphemerisGenerator(
            body_names=["sun", "earth", "mars", "jupiter"],
            start_time=51544.0 * u.day,
            end_time=51545.0 * u.day,
        )

        major_bodies, _ = gen._define_major_and_minor_bodies()

        assert "sun" in major_bodies
        assert major_bodies["sun"] == 10
        assert "earth" in major_bodies
        assert major_bodies["earth"] == 399
        assert "mars" in major_bodies
        assert major_bodies["mars"] == 499
        assert "jupiter" in major_bodies
        assert major_bodies["jupiter"] == 599

    def test_casefold_body_matching(self):
        """Test that body name matching is case-insensitive."""
        gen = EphemerisGenerator(
            body_names=["SUN", "Earth", "MARS"],
            start_time=51544.0 * u.day,
            end_time=51545.0 * u.day,
        )

        major_bodies, _ = gen._define_major_and_minor_bodies()

        # All should match despite different cases
        assert "sun" in major_bodies
        assert "earth" in major_bodies
        assert "mars" in major_bodies

    def test_define_bodies_creates_list(self):
        """Test define_bodies creates body list."""
        gen = EphemerisGenerator(
            body_names=["sun", "earth"],
            start_time=51544.0 * u.day,
            end_time=51545.0 * u.day,
        )
        gen.define_bodies()

        assert len(gen.bodies) == 2
        assert all(isinstance(body, CelestialBody) for body in gen.bodies)
        assert gen.bodies[0].name == "sun"
        assert gen.bodies[1].name == "earth"

    def test_custom_time_arrays(self):
        """Test bodies can be defined with custom time arrays."""
        custom_times = [
            Time("2000-01-01 00:00:00", format="iso", scale="tdb")
            + np.arange(0, 10) * u.hour,
            Time("2000-01-01 00:00:00", format="iso", scale="tdb")
            + np.arange(0, 5) * u.hour,
        ]

        gen = EphemerisGenerator(
            body_names=["sun", "earth"],
            start_time=51544.0 * u.day,
            end_time=51545.0 * u.day,
        )
        gen.define_bodies(time_arrs=custom_times)

        assert len(gen.bodies[0].time) == 10
        assert len(gen.bodies[1].time) == 5


# ============================================================================
# CelestialBody Dataclass Tests
# ============================================================================


class TestCelestialBodyDataclass:
    """Test CelestialBody dataclass."""

    def test_celestial_body_creation(self):
        """Test CelestialBody can be created with correct attributes."""
        time_arr = (
            Time("2000-01-01 00:00:00", format="iso", scale="tdb")
            + np.arange(0, 10) * u.hour
        )
        position = np.random.rand(3, 10)

        body = CelestialBody(name="earth", time=time_arr, position=position)

        assert body.name == "earth"
        assert len(body.time) == 10
        assert body.position.shape == (3, 10)

    def test_celestial_body_none_position(self):
        """Test CelestialBody can be created with None position."""
        time_arr = (
            Time("2000-01-01 00:00:00", format="iso", scale="tdb")
            + np.arange(0, 10) * u.hour
        )

        body = CelestialBody(name="mars", time=time_arr, position=None)

        assert body.name == "mars"
        assert body.position is None


# ============================================================================
# File I/O Tests
# ============================================================================


class TestFileIO:
    """Test ephemeris file output."""

    def test_hdf5_output_structure(self):
        """Test HDF5 file has correct structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            gen = EphemerisGenerator(
                body_names=["sun"],
                start_time=51544.0 * u.day,
                end_time=51544.1 * u.day,
                time_step=3600 * u.s,
            )
            gen.define_bodies()

            # Create mock position data
            gen.bodies[0].position = np.random.rand(3, 3)
            gen.bodies[0].time_mjd = np.array([51544.0, 51544.04, 51544.08])

            gen.save_ephemeris(outpath=tmpdir, save_format="hdf5")

            # Verify file exists and has correct structure
            hdf5_path = os.path.join(tmpdir, "ephemeris.h5")
            assert os.path.exists(hdf5_path)

            with h5py.File(hdf5_path, "r") as f:
                assert "sun" in f
                assert "time_mjd" in f["sun"]
                assert "position_au" in f["sun"]
                assert f["sun"]["position_au"].attrs["units"] == "au"
                assert (
                    f["sun"]["position_au"].attrs["frame"]
                    == "heliocentric_mean_ecliptic"
                )

    def test_text_output_format(self):
        """Test text file format matches specification."""
        with tempfile.TemporaryDirectory() as tmpdir:
            gen = EphemerisGenerator(
                body_names=["earth"],
                start_time=51544.0 * u.day,
                end_time=51544.1 * u.day,
                time_step=3600 * u.s,
            )
            gen.define_bodies()

            # Create mock position data
            gen.bodies[0].position = np.random.rand(3, 3)
            gen.bodies[0].time = Time([51544.0, 51544.04, 51544.08], format="mjd")

            gen.save_ephemeris(outpath=tmpdir, save_format="txt")

            # Verify file exists
            txt_path = os.path.join(tmpdir, "earth_ephemeris.txt")
            assert os.path.exists(txt_path)

            # Read and verify format
            with open(txt_path, "r") as f:
                lines = f.readlines()
                # Should have header lines and data lines
                assert len(lines) > 3
                # Check header contains expected info
                assert "MJD" in lines[1]
                assert "AU" in lines[1]

    def test_multiple_bodies_hdf5(self):
        """Test HDF5 output with multiple bodies."""
        with tempfile.TemporaryDirectory() as tmpdir:
            gen = EphemerisGenerator(
                body_names=["sun", "earth", "mars"],
                start_time=51544.0 * u.day,
                end_time=51544.1 * u.day,
                time_step=3600 * u.s,
            )
            gen.define_bodies()

            # Create mock data for all bodies
            for body in gen.bodies:
                body.position = np.random.rand(3, 3)
                body.time_mjd = np.array([51544.0, 51544.04, 51544.08])

            gen.save_ephemeris(outpath=tmpdir, save_format="hdf5")

            hdf5_path = os.path.join(tmpdir, "ephemeris.h5")
            with h5py.File(hdf5_path, "r") as f:
                assert "sun" in f
                assert "earth" in f
                assert "mars" in f

    @pytest.mark.slow
    @pytest.mark.real_query
    def test_coordinates_have_physical_meaning(self):
        """Test that real ephemeris coordinates have physically meaningful distances.

        This test queries real JPL Horizons data for inner planets and validates
        that the heliocentric distances match their known orbital parameters.
        """
        gen = EphemerisGenerator(
            body_names=["sun", "earth", "mars"],
            start_time=Time("2000-01-01 00:00:00", format="iso", scale="tdb").mjd
            * u.day,
            end_time=Time("2000-06-30 00:00:00", format="iso", scale="tdb").mjd * u.day,
            time_step=2592000 * u.s,  # ~30 days
        )
        gen.define_bodies()
        gen.generate_ephemeris()

        # Verify Sun is at origin
        sun_pos = gen.bodies[0].position
        sun_distances = np.sqrt(np.sum(sun_pos**2, axis=0))
        assert np.all(
            sun_distances < 1e-6
        ), "Sun should be at origin in heliocentric coordinates"

        # Verify Earth is at ~1 AU
        earth_pos = gen.bodies[1].position
        earth_distances = np.sqrt(np.sum(earth_pos**2, axis=0))
        assert np.all(
            (earth_distances > 0.98) & (earth_distances < 1.02)
        ), "Earth distances should be approximately 1 AU"

        # Verify Mars is at ~1.5 AU
        mars_pos = gen.bodies[2].position
        mars_distances = np.sqrt(np.sum(mars_pos**2, axis=0))
        assert np.all(
            (mars_distances > 1.38) & (mars_distances < 1.70)
        ), "Mars distances should be approximately 1.5 AU"


# ============================================================================
# Edge Case Tests
# ============================================================================


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_single_time_step(self):
        """Test with only one time point."""
        gen = EphemerisGenerator(
            body_names=["sun"],
            start_time=51544.0 * u.day,
            end_time=51544.1 * u.day,  # 2.4 hours
            time_step=3600 * u.s,  # 1 hour
        )

        # Should have at least 1 time step (for the start time)
        assert len(gen.time_arr) >= 1
        # Should have at most 3 time steps (0h, 1h, 2h)
        assert len(gen.time_arr) <= 3

    def test_very_large_time_step(self):
        """Test with time step larger than time range."""
        gen = EphemerisGenerator(
            body_names=["sun"],
            start_time=51544.0 * u.day,
            end_time=51545.0 * u.day,
            time_step=172800 * u.s,  # 2 days
        )

        # Should have minimal time steps
        assert len(gen.time_arr) >= 1

    def test_empty_body_list(self):
        """Test behavior with no bodies specified."""
        gen = EphemerisGenerator(
            body_names=[],
            start_time=51544.0 * u.day,
            end_time=51545.0 * u.day,
        )
        gen.define_bodies()

        assert len(gen.bodies) == 0


# ============================================================================
# Integration Tests
# ============================================================================


class TestIntegration:
    """Integration tests for complete workflows."""

    @pytest.mark.slow
    @patch("ephemeris_generator.Horizons")
    def test_full_workflow_mocked(self, mock_horizons):
        """Test complete workflow with mocked JPL Horizons."""
        # Mock JPL Horizons response
        mock_vec = MagicMock()
        mock_vec.__getitem__.side_effect = lambda key: {
            "x": MagicMock(quantity=np.array([1.0, 1.1, 1.2]) * u.au),
            "y": MagicMock(quantity=np.array([0.0, 0.1, 0.2]) * u.au),
            "z": MagicMock(quantity=np.array([0.0, 0.0, 0.1]) * u.au),
            "datetime_jd": MagicMock(
                data=np.array([2451544.5, 2451544.54, 2451544.58])
            ),
        }[key]

        mock_horizons_instance = MagicMock()
        mock_horizons_instance.vectors.return_value = mock_vec
        mock_horizons.return_value = mock_horizons_instance

        with tempfile.TemporaryDirectory() as tmpdir:
            gen = EphemerisGenerator(
                body_names=["sun"],
                start_time=51544.0 * u.day,
                end_time=51544.1 * u.day,
                time_step=3600 * u.s,
            )
            gen.define_bodies()
            gen.generate_ephemeris()
            gen.save_ephemeris(outpath=tmpdir, save_format="hdf5")

            # Verify output exists
            assert os.path.exists(os.path.join(tmpdir, "ephemeris.h5"))

    @pytest.mark.slow
    @pytest.mark.real_query
    def test_sun_position_at_origin(self):
        """Test that Sun's position is at (0, 0, 0) in heliocentric mean ecliptic coordinates.

        This test performs an actual JPL Horizons query (not mocked).
        Mark with @pytest.mark.real_query to run only when needed.
        """
        # Query a short time range to keep test fast
        gen = EphemerisGenerator(
            body_names=["sun"],
            start_time=Time("2000-01-01 00:00:00", format="iso", scale="tdb").mjd
            * u.day,
            end_time=Time("2000-01-01 12:00:00", format="iso", scale="tdb").mjd * u.day,
            time_step=3600 * u.s,  # 1 hour
        )
        gen.define_bodies()
        gen.generate_ephemeris()

        # Sun should be at origin in heliocentric coordinates
        sun_position = gen.bodies[0].position

        # Check all positions are at or very close to (0, 0, 0)
        # Allow small tolerance for numerical precision
        tolerance = 1e-10  # AU
        assert np.all(
            np.abs(sun_position[0, :]) < tolerance
        ), "Sun X coordinate not at origin"
        assert np.all(
            np.abs(sun_position[1, :]) < tolerance
        ), "Sun Y coordinate not at origin"
        assert np.all(
            np.abs(sun_position[2, :]) < tolerance
        ), "Sun Z coordinate not at origin"

        # Check the maximum distance from origin
        distances = np.sqrt(np.sum(sun_position**2, axis=0))
        max_distance = np.max(distances)
        assert (
            max_distance < tolerance
        ), f"Sun maximum distance from origin {max_distance} AU exceeds tolerance"

    @pytest.mark.slow
    @pytest.mark.real_query
    def test_sun_earth_distance_approximately_1au(self):
        """Test that distance between Sun and Earth is approximately 1 AU.

        This test performs actual JPL Horizons queries (not mocked).
        Mark with @pytest.mark.real_query to run only when needed.
        """
        # Query a time range covering different parts of Earth's orbit
        gen = EphemerisGenerator(
            body_names=["sun", "earth"],
            start_time=Time("2000-01-01 00:00:00", format="iso", scale="tdb").mjd
            * u.day,
            end_time=Time("2000-07-01 00:00:00", format="iso", scale="tdb").mjd * u.day,
            time_step=86400
            * 30
            * u.s,  # ~30 days to sample different orbital positions
        )
        gen.define_bodies()
        gen.generate_ephemeris()

        # Get positions
        sun_position = gen.bodies[0].position
        earth_position = gen.bodies[1].position

        # Calculate distances between Sun and Earth at all time points
        distances = np.sqrt(np.sum((earth_position - sun_position) ** 2, axis=0))

        # Earth's orbit is elliptical with perihelion ~0.983 AU and aphelion ~1.017 AU
        # Mean distance should be close to 1 AU
        mean_distance = np.mean(distances)
        min_distance = np.min(distances)
        max_distance = np.max(distances)

        # Check mean distance is close to 1 AU (within 5%)
        assert (
            0.95 < mean_distance < 1.05
        ), f"Mean Sun-Earth distance {mean_distance} AU not close to 1 AU"

        # Check min/max are within expected range for Earth's elliptical orbit
        assert (
            0.95 < min_distance < 1.0
        ), f"Minimum Sun-Earth distance {min_distance} AU outside expected range"
        assert (
            1.0 < max_distance < 1.05
        ), f"Maximum Sun-Earth distance {max_distance} AU outside expected range"

        # Verify all distances are physically reasonable (between 0.9 and 1.1 AU)
        assert np.all(distances > 0.9), "Some distances unreasonably small"
        assert np.all(distances < 1.1), "Some distances unreasonably large"

    @pytest.mark.slow
    @pytest.mark.real_query
    def test_queried_time_matches_requested_time(self):
        """Test that JPL Horizons returns times matching the requested time grid.

        This test verifies that the time array generated by the EphemerisGenerator
        matches the actual times returned by JPL Horizons queries. This is important
        because time step rounding or conversion between different time formats
        could introduce discrepancies.

        This test performs actual JPL Horizons queries (not mocked).
        """
        # Use a short time range with specific time step
        gen = EphemerisGenerator(
            body_names=["earth"],
            start_time=Time("2000-01-01 00:00:00", format="iso", scale="tdb").mjd
            * u.day,
            end_time=Time("2000-01-02 00:00:00", format="iso", scale="tdb").mjd * u.day,
            time_step=3600 * u.s,  # 1 hour
        )
        gen.define_bodies()

        # Store the original requested time array
        requested_time_mjd = gen.bodies[0].time.mjd

        # Query JPL Horizons
        gen.generate_ephemeris()

        # Get the actual time array returned from JPL Horizons
        actual_time_mjd = gen.bodies[0].time_mjd

        # Both arrays should have the same length
        assert len(requested_time_mjd) == len(
            actual_time_mjd
        ), f"Time array length mismatch: requested {len(requested_time_mjd)}, got {len(actual_time_mjd)}"

        # Calculate time differences in seconds
        time_diff_seconds = (
            actual_time_mjd - requested_time_mjd
        ) * 86400  # Convert days to seconds

        # Maximum difference should be very small (tolerance: 1 second)
        # JPL Horizons may have slight rounding in time representation
        max_diff = np.max(np.abs(time_diff_seconds))
        assert (
            max_diff < 1.0
        ), f"Maximum time difference {max_diff:.3f} seconds exceeds 1 second tolerance"

        # Mean difference should be even smaller
        mean_diff = np.mean(np.abs(time_diff_seconds))
        assert (
            mean_diff < 0.5
        ), f"Mean time difference {mean_diff:.3f} seconds exceeds 0.5 second tolerance"

        # Verify time spacing is consistent (should match requested time step)
        requested_spacing = np.diff(requested_time_mjd)
        actual_spacing = np.diff(actual_time_mjd)

        # Check that spacing is consistent within tolerance
        spacing_diff = (
            np.abs(actual_spacing - requested_spacing) * 86400
        )  # Convert to seconds
        max_spacing_diff = np.max(spacing_diff)
        assert (
            max_spacing_diff < 1.0
        ), f"Time step consistency error: max difference {max_spacing_diff:.3f} seconds"

    @pytest.mark.slow
    @pytest.mark.real_query
    def test_heliocentric_coordinates_physical_validity(self):
        """Test that heliocentric coordinates have physically reasonable distances.

        This test validates that:
        - Inner planets (Mercury, Venus, Earth, Mars) are at ~0.4-1.5 AU
        - Outer planets (Jupiter, Saturn, Uranus, Neptune) are at their known distances
        - All bodies maintain approximately their known orbital semi-major axes

        This test performs actual JPL Horizons queries (not mocked).
        """
        # Expected approximate semi-major axes in AU
        expected_distances = {
            "mercury": (0.35, 0.45),  # perihelion ~0.307, aphelion ~0.467
            "venus": (0.70, 0.73),  # perihelion ~0.718, aphelion ~0.728
            "earth": (0.98, 1.02),  # perihelion ~0.983, aphelion ~1.017
            "mars": (1.38, 1.67),  # perihelion ~1.381, aphelion ~1.666
            "jupiter": (4.95, 5.45),  # perihelion ~4.950, aphelion ~5.458
            "saturn": (9.00, 10.05),  # perihelion ~9.020, aphelion ~10.055
        }

        # Query a time range covering different parts of orbits
        bodies_to_test = list(expected_distances.keys())
        gen = EphemerisGenerator(
            body_names=bodies_to_test,
            start_time=Time("2000-01-01 00:00:00", format="iso", scale="tdb").mjd
            * u.day,
            end_time=Time("2000-12-31 00:00:00", format="iso", scale="tdb").mjd * u.day,
            time_step=2592000 * u.s,  # ~30 days
        )
        gen.define_bodies()
        gen.generate_ephemeris()

        # Validate each body's distance from the Sun
        for i, body_name in enumerate(bodies_to_test):
            position = gen.bodies[i].position

            # Calculate distances from origin (Sun in heliocentric coordinates)
            distances = np.sqrt(np.sum(position**2, axis=0))

            min_dist = np.min(distances)
            max_dist = np.max(distances)
            mean_dist = np.mean(distances)

            expected_min, expected_max = expected_distances[body_name]

            # Check that measured distances fall within expected ranges
            assert min_dist >= expected_min * 0.95, (
                f"{body_name.capitalize()}: minimum distance {min_dist:.3f} AU "
                f"below expected range [{expected_min:.3f}, {expected_max:.3f}] AU"
            )
            assert max_dist <= expected_max * 1.05, (
                f"{body_name.capitalize()}: maximum distance {max_dist:.3f} AU "
                f"above expected range [{expected_min:.3f}, {expected_max:.3f}] AU"
            )

            # Verify mean distance is reasonable
            assert expected_min * 0.95 <= mean_dist <= expected_max * 1.05, (
                f"{body_name.capitalize()}: mean distance {mean_dist:.3f} AU "
                f"outside expected range [{expected_min:.3f}, {expected_max:.3f}] AU"
            )

    @pytest.mark.slow
    @pytest.mark.real_query
    def test_earth_approximately_1au_from_sun(self):
        """Test that Earth's distance from Sun is approximately 1 AU throughout its orbit.

        This test specifically validates that the Earth's orbital parameters are
        correctly represented in heliocentric mean ecliptic coordinates by checking
        that its distance from the Sun (at the origin) is consistently close to 1 AU.

        This test performs actual JPL Horizons queries (not mocked).
        """
        # Query Earth over a full year to sample entire orbit
        gen = EphemerisGenerator(
            body_names=["sun", "earth"],
            start_time=Time("2000-01-01 00:00:00", format="iso", scale="tdb").mjd
            * u.day,
            end_time=Time("2000-12-31 23:59:59", format="iso", scale="tdb").mjd * u.day,
            time_step=86400 * u.s,  # 1 day
        )
        gen.define_bodies()
        gen.generate_ephemeris()

        # Get Sun at origin and Earth position
        sun_pos = gen.bodies[0].position
        earth_pos = gen.bodies[1].position

        # Calculate heliocentric distance of Earth
        # In heliocentric coordinates, Sun should be at origin anyway
        earth_distances = np.sqrt(np.sum(earth_pos**2, axis=0))

        # Statistical validation
        mean_distance = np.mean(earth_distances)
        min_distance = np.min(earth_distances)
        max_distance = np.max(earth_distances)
        std_distance = np.std(earth_distances)

        # Earth's orbit: perihelion ~0.9833 AU, aphelion ~1.0167 AU
        # Mean should be very close to 1 AU
        assert (
            0.98 < mean_distance < 1.02
        ), f"Earth mean distance {mean_distance:.6f} AU not close to 1 AU"

        # Check perihelion and aphelion are within expected ranges
        assert (
            0.98 < min_distance < 0.99
        ), f"Earth perihelion {min_distance:.6f} AU outside expected range [0.98, 0.99]"
        assert (
            1.01 < max_distance < 1.02
        ), f"Earth aphelion {max_distance:.6f} AU outside expected range [1.01, 1.02]"

        # Standard deviation should be small (orbit is nearly circular)
        # Earth's orbital variation is (aphelion - perihelion) / 2 ≈ 0.0167 AU
        # so std should be roughly 0.0167 / sqrt(12) ≈ 0.0048 AU for uniform distribution
        assert (
            std_distance < 0.015
        ), f"Earth orbital variation {std_distance:.6f} AU is too large"

        # Verify Sun position is at origin (heliocentric)
        sun_distances = np.sqrt(np.sum(sun_pos**2, axis=0))
        max_sun_distance = np.max(sun_distances)
        assert (
            max_sun_distance < 1e-6
        ), f"Sun position {max_sun_distance} AU not at origin"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
