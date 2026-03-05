import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from astropy.coordinates import SkyCoord
import astropy.units as u
############################################################
# IRAS Focal Plane Geometry (from ircc_const.c, yloc/zloc)
# yloc and zloc are in ARCMINUTES, indexed by IRAS detector number (1-based)
yloc = [
    28.00, 28.01, 28.01, 23.97, 23.97, 23.98, 23.98,
    19.66, 19.67, 19.67, 19.67, 17.14, 17.15, 17.15, 17.15,
    14.03, 14.03, 14.03, 12.26, 12.26, 12.27, 12.27,
    9.47, 9.48, 9.49, 9.49, 7.71, 7.71, 7.72, 7.72,
    4.54, 4.53, 4.53, 4.53, 2.02, 2.01, 2.01, 2.01,
    -1.17, -1.17, -1.17, -1.17, -2.93, -2.93, -2.93, -2.93,
    -5.68, -5.68, -5.67, -5.67, -7.44, -7.44, -7.43, -7.43,
    -11.33, -11.33, -11.33, -11.32, -15.36, -15.36, -15.36, -15.36,
    99999, 99999, 99999, 99999, 99999, 99999, 99999, 99999,
    -26.35, -26.40, -26.40, -26.50, -26.45
]
zloc = [
    8.71, 0.04, -8.62, 12.86, 4.37, -4.29, -12.77,
    9.80, 1.14, -7.53, -14.46, 13.49, 5.47, -3.20, -11.86,
    8.71, 0.04, -8.62, 12.96, 4.37, -4.29, -12.88,
    9.81, 1.14, -7.52, -14.50, 13.55, 5.47, -3.19, -11.86,
    14.55, 7.61, -1.06, -9.73, 11.94, 3.27, -5.40, -13.41,
    14.05, 6.55, -2.12, -10.78, 10.88, 2.22, -6.45, -13.95,
    14.64, 7.65, -1.02, -9.68, 11.98, 3.32, -5.35, -13.41,
    13.95, 6.55, -2.12, -10.79, 10.88, 2.21, -6.46, -13.85,
    99999, 99999, 99999, 99999, 99999, 99999, 99999, 99999,
    8.15, 3.00, -1.85, 7.35, -0.15
]

# Conversion factor: arcmin to radians
ARC_MIN_TO_RAD = np.pi / (180.0 * 60.0)

# Coordinate transformation functions using astropy
def equatorial_to_galactic(ra_deg, dec_deg):
    """
    Convert J2000 equatorial coordinates (RA, Dec) to Galactic coordinates (l, b).
    
    Args:
        ra_deg: Right Ascension in degrees (float or array)
        dec_deg: Declination in degrees (float or array)
        
    Returns:
        (l, b): Galactic longitude and latitude in degrees
    """
    coord = SkyCoord(ra=u.Quantity(ra_deg, u.deg), dec=u.Quantity(dec_deg, u.deg), frame='icrs')
    gal = coord.transform_to('galactic')
    return gal.l.deg, gal.b.deg


def ecliptic_to_equatorial(lon_ecl, lat_ecl):
    """
    Convert ecliptic coordinates to J2000 equatorial coordinates.
    
    Args:
        lon_ecl: Ecliptic longitude in degrees (float or array)
        lat_ecl: Ecliptic latitude in degrees (float or array)
        
    Returns:
        (ra, dec): Right Ascension and Declination in degrees
    """
    coord = SkyCoord(lon=u.Quantity(lon_ecl, u.deg), lat=u.Quantity(lat_ecl, u.deg), frame='barycentricmeanecliptic')
    icrs = coord.transform_to('icrs')
    return icrs.ra.deg, icrs.dec.deg


def sun_ecliptic_longitude(satcal):
    """
    Calculate the sun's ecliptic longitude for a given SATCAL (Julian seconds since 1983-01-01 00:00).
    
    The sun's ecliptic longitude is approximately its mean longitude in the ecliptic frame,
    which increases by ~1 degree per day as Earth orbits the sun.
    
    Args:
        satcal: Julian seconds since 1983-01-01 (IRAS convention)
        
    Returns:
        Sun's ecliptic longitude in degrees (0-360)
    """
    # SATCAL is seconds since 1983-01-01 00:00 UT
    # Convert to days
    days_since_1983 = satcal / 86400.0
    
    # Mean solar longitude calculation
    # On 1983-01-01, the sun was at approximately 280.46 degrees ecliptic longitude
    # The sun moves approximately 0.98565 degrees per day (360 deg / 365.25 days)
    sun_lon_1983_jan_01 = 280.46  # degrees
    sun_motion_rate = 0.98565  # degrees per day
    
    # Calculate sun's ecliptic longitude
    sun_lon = sun_lon_1983_jan_01 + sun_motion_rate * days_since_1983
    sun_lon = np.mod(sun_lon, 360.0)
    
    return sun_lon


def ecliptic_to_galactic(lon_ecl, lat_ecl):
    """
    Convert ecliptic coordinates to Galactic coordinates.
    
    Args:
        lon_ecl: Ecliptic longitude in degrees (float or array)
        lat_ecl: Ecliptic latitude in degrees (float or array)
        
    Returns:
        (l, b): Galactic longitude and latitude in degrees
    """
    coord = SkyCoord(lon=u.Quantity(lon_ecl, u.deg), lat=u.Quantity(lat_ecl, u.deg), frame='barycentricmeanecliptic')
    gal = coord.transform_to('galactic')
    return gal.l.deg, gal.b.deg


def iras_detector_offset(detnum):
    """
    Given IRAS detector number (1-based), return (dy, dz) in radians.
    If detnum is invalid, returns (None, None).
    Note: yloc/zloc are in ARCMINUTES.
    """
    if detnum is None or detnum < 1 or detnum > len(yloc):
        return (None, None)
    dy = yloc[detnum-1]
    dz = zloc[detnum-1]
    if abs(dy) > 1e4 or abs(dz) > 1e4:
        return (None, None)
    # Convert from ARCMINUTES to radians
    return (dy * ARC_MIN_TO_RAD, dz * ARC_MIN_TO_RAD)

def iras_focalplane_to_sky(lon0, lat0, twist_deg, dy_rad, dz_rad):
    """
    Given boresight lon0, lat0 (deg), twist (deg), and detector offsets dy, dz (radians),
    compute the detector sky position (lon, lat) in degrees.
    This mimics the logic in det_xt() from irds_rd_detpos.c.
    """
    # Convert to radians
    lon0 = np.deg2rad(lon0)
    lat0 = np.deg2rad(lat0)
    twist = np.deg2rad(twist_deg)
    # 1. Compute r (projected radius)
    r = 1.0 / np.sqrt(1.0 + dy_rad**2 + dz_rad**2)
    # 2. Detector offset in tangent plane (u, v)
    u = dy_rad * np.cos(twist) + dz_rad * np.sin(twist)
    v = -dy_rad * np.sin(twist) + dz_rad * np.cos(twist)
    # 3. Projected xyz in telescope (plate) system
    x = r
    y = u * r
    z = v * r
    # 4. Rotate to sky (sun-referenced) system
    cosl = np.cos(lon0)
    sinl = np.sin(lon0)
    cosb = np.cos(lat0)
    sinb = np.sin(lat0)
    xx = x * cosl * cosb - y * sinl - z * cosl * sinb
    yy = x * sinl * cosb + y * cosl - z * sinl * sinb
    zz = x * sinb + z * cosb
    # 5. Convert back to lon/lat
    lat = np.arcsin(zz)
    lon = np.arctan2(yy, xx)
    # Return in degrees
    return (np.rad2deg(lon), np.rad2deg(lat))

# This function attempts to mimic the logic in lrstimes.c for sop files
# It assumes the sop file has a header with keys like NAXIS1, SATCAL, INSTRUME, CRVAL1, CDELT1, CRPIX1
# and a data section with timestream values

def parse_sop_header(f):
    header = {}
    header_bytes = b""
    # Read header blocks (assume 2880 bytes per block, like FITS)
    while True:
        block = f.read(2880)
        if not block:
            break
        header_bytes += block
        if b"END" in block:
            break
    header_text = header_bytes.decode("ascii", errors="ignore")
    for i in range(0, len(header_text), 80):
        card = header_text[i:i+80]
        key = card[:8].strip()
        if key == "END":
            break
        if "=" in card:
            value = card[10:30].strip()
            header[key] = value
    return header

def generate_table():
    """
    Generate the IRAS difference table using the logarithmic compression scheme.
    
    The compression encodes differences as:
        0bs_cccc_xxx
    where s is the sign bit, cccc is the shift code (14 patterns), 
    and xxx is the significand (3 bits).
    
    Returns:
        np.ndarray: Sorted array of positive difference values.
    """
    rows = [
        # Each row: bit pattern string with xxx placeholders
        # From smallest magnitude band to largest
        "1xxx000000000000",      # shift_code 0
        "01xxx00000000000",      # shift_code 1
        "001xxx0000000000",      # shift_code 2
        "0001xxx000000000",      # shift_code 3
        "00001xxx00000000",      # shift_code 4
        "000001xxx0000000",      # shift_code 5
        "0000001xxx000000",      # shift_code 6
        "00000001xxx00000",      # shift_code 7
        "000000001xxx0000",      # shift_code 8
        "0000000001xxx000",      # shift_code 9
        "000000000011xxx0",      # shift_code 10
        "000000000010xxx0",      # shift_code 11
        "0000000000011xxx",      # shift_code 12
        "0000000000010xxx",      # shift_code 13
        "0000000000001xxx",      # shift_code 14
        "0000000000000xxx",      # shift_code 15
    ]

    values = []

    # Iterate from smallest magnitude band to largest
    for pattern in reversed(rows):
        for xxx in range(8):
            bits = pattern.replace("xxx", format(xxx, "03b"))
            value = int(bits, 2)
            values.append(value)

    # Remove zero if present (will be added separately)
    values = [v for v in values if v != 0]

    return sorted(values)


def build_diftab():
    """
    Build the full 256-entry IRAS difference table.
    
    The table structure is:
        D[0] = 0         (blank)
        D[1:128] = negative differences (sorted)
        D[128] = 0       (zero)
        D[129:256] = positive differences (sorted)
    
    Returns:
        np.ndarray: 256-entry difference table.
    """
    table = np.array(generate_table())
    D = np.zeros(256, dtype=np.int32)
    D[0] = 0
    D[128] = 0
    D[1:128] = -table[::-1]
    D[129:] = table
    return D


def extract_timestream(filepath, use_decompression=True):
    with open(filepath, "rb") as f:
        header = parse_sop_header(f)
        naxis = int(header.get("NAXIS", 1))
        naxis1 = int(header.get("NAXIS1", 0))
        naxis2 = int(header.get("NAXIS2", 1))
        satcal = int(header.get("SATCAL", 0))
        crval1 = float(header.get("CRVAL1", 0.0))
        cdelt1 = float(header.get("CDELT1", 1.0))
        crpix1 = float(header.get("CRPIX1", 1.0))
        # Find the end of the FITS header (look for 'END' and pad to 2880 bytes)
        f.seek(0)
        header_bytes = b""
        while True:
            block = f.read(2880)
            if not block:
                break
            header_bytes += block
            if b"END" in block:
                break
        # Find where 'END' occurs
        end_idx = header_bytes.find(b'END')
        if end_idx == -1:
            raise RuntimeError("FITS header END not found")
        # Data starts at next 2880-byte boundary after 'END'
        header_end = ((end_idx + 80 - 1) // 2880 + 1) * 2880
        f.seek(header_end)
        # Read all data as flat array
        total_samples = naxis1 * naxis2
        if use_decompression:
            # Read as uint8 compressed data
            raw_data = np.fromfile(f, dtype='<u1', count=total_samples)
        else:
            # Read as uint16 raw data
            raw_data = np.fromfile(f, dtype='<u2', count=total_samples)
    if naxis2 > 1:
        # Try detector-major first
        if raw_data.size == total_samples:
            data = np.empty((naxis1, naxis2), dtype=raw_data.dtype)
            for det in range(naxis2):
                start = det * naxis1
                end = start + naxis1
                data[:, det] = raw_data[start:end]
        else:
            # Try alternate shapes
            if raw_data.size % naxis2 == 0:
                alt_naxis1 = raw_data.size // naxis2
                data = raw_data[:alt_naxis1 * naxis2].reshape((alt_naxis1, naxis2))
            elif raw_data.size % naxis1 == 0:
                alt_naxis2 = raw_data.size // naxis1
                data = raw_data[:naxis1 * alt_naxis2].reshape((naxis1, alt_naxis2))
            else:
                data = raw_data
    else:
        data = raw_data
    
    # If using decompression, apply the difference table
    if use_decompression:
        diftab = build_diftab()
        if data.ndim == 2:
            # Decompress each detector separately
            decompressed = np.zeros_like(data, dtype=np.int32)
            for det in range(data.shape[1]):
                init = header.get(f"INIT00{det+1:02}", 0)
                decompressed[0,det] = init
                for i in range(1, data.shape[0]):
                    diff_idx = data[i, det]
                    decompressed[i, det] = decompressed[i-1, det] + diftab[diff_idx]
            data = decompressed
        else:
            # Decompress 1D data
            decompressed = np.zeros_like(data, dtype=np.int32)
            for i in range(1, data.shape[0]):
                diff_idx = data[i]
                decompressed[i] = decompressed[i-1] + diftab[diff_idx]
            data = decompressed
    # Adjust xdat to match actual data length (IRPL_EXTEND-style robustness)
    actual_naxis1 = data.shape[0] if data.ndim > 1 else data.size
    if actual_naxis1 < naxis1:
        # Optionally, warn or handle truncated data here
        pass
    xdat = np.arange(actual_naxis1)
    xdat = (xdat + 1.0 - crpix1) * cdelt1
    xdat += satcal
    xdat += crval1
    # Print physical detector labels for this file (must be before return)
    instrume = header.get("INSTRUME", "").strip().upper()
    iras_maps = {
        # 100 micron
        '100B': list(range(1, 8)),         # 1-7
        '100A': list(range(55, 63)),       # 55-62
        # 60 micron
        '60B': list(range(8, 16)),         # 8-15
        '60A': list(range(31, 39)),        # 31-38
        # 25 micron
        '25B': list(range(16, 23)),        # 16-22
        '25A': list(range(39, 47)),        # 39-46
        # 12 micron
        '12B': list(range(23, 31)),        # 23-30
        '12A': list(range(47, 55)),        # 47-54
    }
    band = None
    mode = None
    if "B4" in instrume or "100" in instrume:
        band = "100"
    elif "B3" in instrume or "65" in instrume:
        band = "65"
    elif "B2" in instrume or "25" in instrume:
        band = "25"
    elif "B1" in instrume or "12" in instrume:
        band = "12"
    if "SURVEY" in instrume:
        mode = "A"
    elif "AO" in instrume or "SPLINE" in instrume:
        mode = "B"
    map_key = None
    if band and mode:
        map_key = band + mode
    elif band:
        map_key = band + "A"
    if map_key in iras_maps:
        det_labels = iras_maps[map_key]
    return xdat, data  # data is 1D or 2D depending on file

# --- New function to extract a single detector's timestream ---
def extract_single_detector_timestream(filepath, detnum=1):
    """
    Extract the timestream for a specific physical detector number (1-based).
    Returns (xdat, ydat) for the selected detector, or (xdat, None) if not found.
    """
    xdat, data = extract_timestream(filepath)
    with open(filepath, "rb") as f:
        header = parse_sop_header(f)
    naxis2 = int(header.get("NAXIS2", 1))
    instrume = header.get("INSTRUME", "").strip().upper()
    # Mapping logic (same as in main)
    iras_maps = {
        # 100 micron
        '100B': list(range(1, 8)),         # 1-7
        '100A': list(range(55, 63)),       # 55-62
        # 60 micron
        '60B': list(range(8, 16)),         # 8-15
        '60A': list(range(31, 39)),        # 31-38
        # 25 micron
        '25B': list(range(16, 23)),        # 16-22
        '25A': list(range(39, 47)),        # 39-46
        # 12 micron
        '12B': list(range(23, 31)),        # 23-30
        '12A': list(range(47, 55)),        # 47-54
    }
    band = None
    mode = None
    if "B4" in instrume or "100" in instrume:
        band = "100"
    elif "B3" in instrume or "65" in instrume:
        band = "65"
    elif "B2" in instrume or "25" in instrume:
        band = "25"
    elif "B1" in instrume or "12" in instrume:
        band = "12"
    if "SURVEY" in instrume:
        mode = "A"
    elif "AO" in instrume or "SPLINE" in instrume:
        mode = "B"
    map_key = None
    if band and mode:
        map_key = band + mode
    elif band:
        map_key = band + "A"
    # Find the detector index for the requested detnum
    if map_key in iras_maps:
        try:
            det_idx = iras_maps[map_key].index(detnum)
        except ValueError:
            return xdat, None
        # If data is 2D, select the column for this detector
        if naxis2 > 1 and data.ndim == 2:
            ydat = data[:, det_idx]
        elif naxis2 == 1 and data.ndim == 1:
            ydat = data
        else:
            ydat = None
        return xdat, ydat
    else:
        return xdat, None

if __name__ == "__main__":
    # Example usage
    from glob import glob
    sopfiles = glob('/home/dwatts/IRAS/*/*/survey.b1/*/*')
    sopfiles.sort()
    
    # STRIDE control for file sampling:
    # Each SOP file covers ~80° (half-scan), so you need multiple files for full coverage
    # - stride=1: Use ALL files (best coverage, slow)
    # - stride=5: Use every 5th file (good coverage, moderate speed)
    # - stride=10: Use every 10th file (reasonable coverage, faster)
    # - stride=50: Use every 50th file (sparse, very fast, poor coverage) ← original
    # 
    # To get full sky coverage: Use stride=1 or stride=5
    STRIDE = 5  # Change this to 1 for best coverage, or 50 for fast testing
    sopfiles = sopfiles[::STRIDE]
    print(f"Using {len(sopfiles)} files (stride={STRIDE})")
    
    plt.figure(figsize=(10,6))
    NSIDE = 128
    hits_map = np.zeros(12*NSIDE**2)
    bmap = np.zeros(12*NSIDE**2)
    from tqdm import tqdm
    for sop_file in tqdm(sopfiles):
        # Print all header keys and values
        with open(sop_file, "rb") as f:
            header = parse_sop_header(f)
        #print("Header keys and values:")
        #for k, v in header.items():
        #    print(f"{k}: {v}")

        # Compute derived pointing
        psi = float(header.get("PSI", 0))
        psirate = float(header.get("PSIRATE", 0))
        theta = float(header.get("THETA", 0))
        satcal = int(header.get("SATCAL", 0))
        sop = int(header.get("SOP", 0))
        att = int(header.get("ATT", 0))
        obs = int(header.get("OBS", 0))
        crval1 = float(header.get("CRVAL1", 0))
        crval2 = float(header.get("CRVAL2", 0))
        cdelt1 = float(header.get("CDELT1", 1))
        cdelt2 = float(header.get("CDELT2", 1))
        crpix1 = float(header.get("CRPIX1", 1))
        crpix2 = float(header.get("CRPIX2", 1))
        #print(f"\nDerived pointing parameters:")
        #print(f"PSI: {psi} deg")
        #print(f"PSIRATE: {psirate} deg/sec")
        #print(f"THETA: {theta} deg")
        #print(f"SATCAL: {satcal}")
        #print(f"SOP: {sop}")
        #print(f"ATT: {att}")
        #print(f"OBS: {obs}")
        #print(f"CRVAL1: {crval1}")
        #print(f"CRVAL2: {crval2}")
        #print(f"CDELT1: {cdelt1}")
        #print(f"CDELT2: {cdelt2}")
        #print(f"CRPIX1: {crpix1}")
        #print(f"CRPIX2: {crpix2}")

        xdat, ydat = extract_timestream(sop_file)


        # IRAS detector mappings from legacy code
        iras_maps = {
            # The following mappings are from the C array 'detnrs' in:
            # /home/dwatts/IRAS/diskrog10androg11/rog10/gipsy/sub/ircc_const.c lines 80-100
            # 100B (B4, survey): detnrs[4][0..15]
            '100B': [1, 2, 3, 4, 5, 6, 7, 55, 56, 57, 58, 59, 60, 61, 62, -1],
            # 100A (B4, AO/SPLINE): detnrs[3][0..15]
            '100A': [8, 9, 10, 11, 12, 13, 14, 15, 31, 32, 33, 34, 35, 36, 37, 38],
            # 65A (B3, survey): detnrs[2][0..15]
            '65A': [16, 17, 18, 19, 20, 21, 22, 39, 40, 41, 42, 43, 44, 45, 46, -1],
            # 65B (B3, AO/SPLINE): detnrs[3][0..15] (same as 100A)
            '65B': [8, 9, 10, 11, 12, 13, 14, 15, 31, 32, 33, 34, 35, 36, 37, 38],
            # 25A (B2, survey): detnrs[1][0..15]
            '25A': [23, 24, 25, 26, 27, 28, 29, 30, 47, 48, 49, 50, 51, 52, 53, 54],
            # 25B (B2, AO/SPLINE): detnrs[2][0..15] (first 15 only)
            '25B': [16, 17, 18, 19, 20, 21, 22, 39, 40, 41, 42, 43, 44, 45, 46],
            # 12A (B1, survey): detnrs[4][0..15] (same as 100B)
            '12A': [1, 2, 3, 4, 5, 6, 7, 55, 56, 57, 58, 59, 60, 61, 62, -1],
            # 12B (B1, AO/SPLINE): detnrs[1][0..15] (same as 25A)
            '12B': [23, 24, 25, 26, 27, 28, 29, 30, 47, 48, 49, 50, 51, 52, 53, 54],
        }

        naxis2 = int(header.get("NAXIS2", 1))
        detector_indices = np.arange(naxis2)
        instrume = header.get("INSTRUME", "").strip().upper()
        band = None
        mode = None
        # Try to determine band and mode from INSTRUME
        if "B4" in instrume or "100" in instrume:
            band = "100"
        elif "B3" in instrume or "65" in instrume:
            band = "65"
        elif "B2" in instrume or "25" in instrume:
            band = "25"
        elif "B1" in instrume or "12" in instrume:
            band = "12"
        # Mode: survey vs AO/SPLINE
        if "SURVEY" in instrume:
            mode = "A"
        elif "AO" in instrume or "SPLINE" in instrume:
            mode = "B"

        # Compose mapping key
        map_key = None
        if band and mode:
            map_key = band + mode
        elif band:
            # fallback: if only band is known, prefer survey
            map_key = band + "A"



        #print("\nPlotting all valid physical detectors (1-based, using yloc/zloc):")
        sample_indices = np.arange(len(xdat))
        sky_lons = []
        sky_lats = []
        labels = []
        # Precompute arrays for all samples
        # For meridional pole-to-pole scans:
        # - Longitude is CONSTANT = THETA (boresight longitude direction)
        # - Latitude varies with time = PSI + PSIRATE × t (boresight sweeps poles)
        # CDELT1 = seconds per sample
        # PSIRATE = rate of latitude change per second
        # 
        # Boresight position:
        # lon = THETA (constant)
        # lat(t) = PSI + PSIRATE × t_seconds (varies from pole to pole)
        # 
        # Time in seconds from the reference time SATCAL:
        time_seconds = (sample_indices - crpix1) * cdelt1
        
        # Boresight ecliptic longitude and co-latitude (sun-referenced frame)
        lon0 = theta  # ecliptic LONGITUDE (constant)
        lat0_arr = psi + psirate * time_seconds  # ecliptic co-latitude (varies pole-to-pole)
        
        # Detector offset rotation:
        # For IRAS, the focal plane coordinate system may have its own twist.
        # Based on legacy code, detector offsets (dy, dz) are rotated by twist angle.
        # However, the raw FITS header doesn't specify a time-varying twist angle.
        # For now, assume zero twist (focal plane aligned with scan direction).
        twist_arr = np.zeros_like(time_seconds)  # No twist variation mentioned in header
        
        lon0_rad = np.deg2rad(lon0)  # scalar
        lat0_rad = np.deg2rad(lat0_arr)  # array
        twist_rad = np.deg2rad(twist_arr)
        cosl = np.cos(lon0_rad)  # scalar
        sinl = np.sin(lon0_rad)  # scalar
        cosb_arr = np.cos(lat0_rad)  # array
        sinb_arr = np.sin(lat0_rad)  # array
        cos_twist = np.cos(twist_rad)
        sin_twist = np.sin(twist_rad)
        for detnum in range(1, len(yloc)+1):
            dy = yloc[detnum-1]
            dz = zloc[detnum-1]
            if abs(dy) > 1e4 or abs(dz) > 1e4:
                continue  # skip invalid detectors
            # Convert detector offsets from ARCMINUTES to radians
            dy_rad = dy * ARC_MIN_TO_RAD
            dz_rad = dz * ARC_MIN_TO_RAD
            # Vectorized calculation for all samples
            r = 1.0 / np.sqrt(1.0 + dy_rad**2 + dz_rad**2)
            # Rotate detector offsets by twist angle (constant in this case)
            det_u = dy_rad * cos_twist + dz_rad * sin_twist
            det_v = -dy_rad * sin_twist + dz_rad * cos_twist
            x = np.ones_like(lat0_arr) * r
            y = det_u * r
            z = det_v * r
            # Projected xyz for all samples (lat0_arr varies, lon0 is constant)
            # Using transformation from legacy irds_rd_detpos.c det_xt():
            # x' = x * cosl * cosb - y * sinl - z * cosl * sinb
            # y' = x * sinl * cosb + y * cosl - z * sinl * sinb
            # z' = x * sinb + z * cosb
            # With lon0 scalar and lat0_arr array:
            xx = x * cosl * cosb_arr - y * sinl - z * cosl * sinb_arr
            yy = x * sinl * cosb_arr + y * cosl - z * sinl * sinb_arr
            zz = x * sinb_arr + z * cosb_arr
            # Clamp zz to [-1, 1] to avoid numerical precision issues with arcsin
            zz = np.clip(zz, -1, 1)
            lat_arr = np.arcsin(zz)
            lon_arr = np.arctan2(yy, xx)
            sky_lons.append(np.rad2deg(lon_arr))
            sky_lats.append(np.rad2deg(lat_arr))
            labels.append(f"Det {detnum}")

        #print("uint16: min=", np.min(ydat), "max=", np.max(ydat))
        #plt.subplot(2,1,1)
        #for i in range(naxis2):
        #    plt.plot(xdat, ydat[:,i], '.', label='Signal', ms=1, color=plt.cm.viridis(i/naxis2))
        #plt.xlabel('Time (s)')
        #plt.ylabel('Signal')
        #plt.title('Timestream from SOP file')
        #plt.legend()
        plt.subplot(2,1,2)
        # Plot only the detectors corresponding to the current band/mode
        # Use iras_maps to get the relevant physical detector numbers
        relevant_detnums = []
        if map_key in iras_maps:
            relevant_detnums = [d for d in iras_maps[map_key] if isinstance(d, int) and d > 0]
        # Find indices in the sky_lons/skylats/labels arrays that match relevant_detnums
        indices_to_plot = []
        for detnum in relevant_detnums:
            # labels are 'Det {detnum}'
            try:
                idx = labels.index(f"Det {detnum}")
                indices_to_plot.append(idx)
            except ValueError:
                continue
        for i in indices_to_plot:
            # Sky coordinates from focal plane geometry are in SUN-CENTERED ECLIPTIC
            # (as per legacy irds_rd_detpos.c det_xt() transformation)
            lon_sr = np.array(sky_lons[i], dtype=np.float64)  # sun-referenced longitude
            lat_sr = np.array(sky_lats[i], dtype=np.float64)  # sun-referenced latitude
            
            lat_sr = np.clip(lat_sr, -90, 90)  # Clamp to [-90, 90]
            lon_sr = np.mod(lon_sr + 180, 360) - 180  # Wrap to [-180, 180]
            
            # Convert from sun-referenced ecliptic to ecliptic coordinates
            # by adding the sun's ecliptic longitude
            sun_lon = sun_ecliptic_longitude(satcal)
            lon_ecl = lon_sr + sun_lon
            lon_ecl = np.mod(lon_ecl + 180, 360) - 180
            lat_ecl = lat_sr
            
            # Convert from ecliptic to Galactic coordinates
            lons_gal, lats_gal = ecliptic_to_galactic(lon_ecl, lat_ecl)
            lats_gal = np.clip(lats_gal, -90, 90)
            lons_gal = np.mod(lons_gal + 180, 360) - 180
            
            # Plot the Galactic coordinates
            plt.plot(lons_gal, lats_gal, 'k.', label=labels[i], ms=1)
            
            if i < 2:
                # Remove NaN entries before passing to HEALPix
                valid_mask = ~(np.isnan(lats_gal) | np.isnan(lons_gal))
                lats_gal = lats_gal[valid_mask].astype(np.float64)
                lons_gal = lons_gal[valid_mask].astype(np.float64)
                ydat_filtered = ydat[valid_mask, i]
                
                #print(f"Debug - After filtering: {np.sum(valid_mask)} valid points")
                
                # hp.ang2pix with lonlat=True expects (longitude, latitude) in degrees
                pix = hp.ang2pix(NSIDE, lons_gal, lats_gal, lonlat=True)
                for j in range(len(pix)):
                    hits_map[pix[j]] += 1
                    bmap[pix[j]] += ydat_filtered[j]


        plt.xlabel('Galactic Longitude (deg)')
        plt.ylabel('Galactic Latitude (deg)')
        plt.title(f'Galactic Coordinates for Detectors in {map_key}')
        #plt.legend(fontsize='small', ncol=2)
        plt.tight_layout()
    
    # Save sky coordinates to file for analysis
    #print(f"\nSaving {len(sky_lons)} detectors to timestream_coords.npz...")
    # Convert lists to arrays with consistent shapes by padding
    max_samples = max(len(arr) for arr in sky_lons) if sky_lons else 0
    lon_arr = np.full((len(sky_lons), max_samples), np.nan)
    lat_arr = np.full((len(sky_lons), max_samples), np.nan)
    for i, (lons, lats) in enumerate(zip(sky_lons, sky_lats)):
        lon_arr[i, :len(lons)] = lons
        lat_arr[i, :len(lats)] = lats
    #np.savez('/home/dwatts/IRAS/python/timestream_coords.npz', 
    #         lon_sr=lon_arr, lat_sr=lat_arr, labels=labels)
    #print(f"Saved: timestream_coords.npz with shape {lon_arr.shape}")
    
    hp.mollview(hits_map)
    hp.mollview(bmap/hits_map)
    plt.savefig('/home/dwatts/IRAS/python/healpix_maps.png')
    print("Saved: healpix_maps.png")
    plt.show()
