import argparse
import csv
import functools
import pickle
from pathlib import Path

import fitsio
import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

DEG2RAD = np.pi / 180.0

# Obliquity of the ecliptic at MJD_J2000 = 51544.5 (jdcent = 0)
# From coordinate.cc eps(): epsilon = 23.439292 deg  (IAU 1976 constant term)
ECL_OBLIQUITY_RAD = 23.439292 * DEG2RAD


def wrap_360(deg):
    return np.mod(deg, 360.0)


def wrap_pm180(deg):
    return (deg + 180.0) % 360.0 - 180.0


def normalize_quat(quat_xyzw):
    qnorm = np.linalg.norm(quat_xyzw, axis=1)
    valid = qnorm > 0.0
    if not np.any(valid):
        raise RuntimeError("All attitude quaternions have zero norm.")
    q = quat_xyzw[valid] / qnorm[valid, None]
    return q, valid


def rotate_vectors_by_quat(v, q_xyzw):
    """Rotate vector(s) by quaternion(s), both using xyzw convention.

    v: (N,3)
    q_xyzw: (N,4)
    """
    qv = q_xyzw[:, :3]
    qw = q_xyzw[:, 3]
    t = 2.0 * np.cross(qv, v)
    return v + qw[:, None] * t + np.cross(qv, t)


def rotate_unit_z_by_quat(q_xyzw):
    """Rotate boresight unit vector [0,0,1] by quaternion(s), xyzw convention."""
    x = q_xyzw[:, 0]
    y = q_xyzw[:, 1]
    z = q_xyzw[:, 2]
    w = q_xyzw[:, 3]

    out = np.empty((len(q_xyzw), 3), dtype=np.float64)
    out[:, 0] = 2.0 * (x * z + y * w)
    out[:, 1] = 2.0 * (y * z - x * w)
    out[:, 2] = 1.0 - 2.0 * (x * x + y * y)
    return out


def quat_to_roll_euler(quat_xyzw):
    """Roll extraction matching Quaternion::GetEuler in C++.

    m23 = 2*q1*q2 + 2*q0*q3
    m33 = 2*q3*q3 + 2*q2*q2 - 1
    roll = atan2(m23, m33)
    where [q0,q1,q2,q3] = [x,y,z,w]
    """
    q0 = quat_xyzw[:, 0]
    q1 = quat_xyzw[:, 1]
    q2 = quat_xyzw[:, 2]
    q3 = quat_xyzw[:, 3]

    m23 = 2.0 * q1 * q2 + 2.0 * q0 * q3
    m33 = 2.0 * q3 * q3 + 2.0 * q2 * q2 - 1.0
    return np.arctan2(m23, m33)


def quat_multiply_xyzw(q_left, q_right):
    """Quaternion product q_left * q_right in xyzw convention.

    q_left: (N, 4)
    q_right: (N, 4) or (4,)
    returns: (N, 4)
    """
    qx = q_left[:, 0]
    qy = q_left[:, 1]
    qz = q_left[:, 2]
    qw = q_left[:, 3]

    if q_right.ndim == 1:
        rx, ry, rz, rw = q_right
    else:
        rx = q_right[:, 0]
        ry = q_right[:, 1]
        rz = q_right[:, 2]
        rw = q_right[:, 3]

    ox = qw * rx + qx * rw + qy * rz - qz * ry
    oy = qw * ry - qx * rz + qy * rw + qz * rx
    oz = qw * rz + qx * ry - qy * rx + qz * rw
    ow = qw * rw - qx * rx - qy * ry - qz * rz
    return np.stack((ox, oy, oz, ow), axis=1)


def quat_to_matrix_xyzw(quat_xyzw):
    x = quat_xyzw[:, 0]
    y = quat_xyzw[:, 1]
    z = quat_xyzw[:, 2]
    w = quat_xyzw[:, 3]

    r = np.empty((len(quat_xyzw), 3, 3), dtype=np.float64)
    r[:, 0, 0] = 1.0 - 2.0 * (y * y + z * z)
    r[:, 0, 1] = 2.0 * (x * y - z * w)
    r[:, 0, 2] = 2.0 * (x * z + y * w)
    r[:, 1, 0] = 2.0 * (x * y + z * w)
    r[:, 1, 1] = 1.0 - 2.0 * (x * x + z * z)
    r[:, 1, 2] = 2.0 * (y * z - x * w)
    r[:, 2, 0] = 2.0 * (x * z - y * w)
    r[:, 2, 1] = 2.0 * (y * z + x * w)
    r[:, 2, 2] = 1.0 - 2.0 * (x * x + y * y)
    return r


def quat_to_roll_like_pr_rad(quat_xyzw, ra_rad, dec_rad):
    """PR-compatible roll from quaternion and boresight (radians input)."""
    x = quat_xyzw[:, 0]
    y = quat_xyzw[:, 1]
    z = quat_xyzw[:, 2]
    w = quat_xyzw[:, 3]

    # Body X-axis in sky frame (first column of rotation matrix).
    bx0 = 1.0 - 2.0 * (y * y + z * z)
    bx1 = 2.0 * (x * y + z * w)
    bx2 = 2.0 * (x * z - y * w)

    sinra = np.sin(ra_rad)
    cosra = np.cos(ra_rad)
    sindec = np.sin(dec_rad)
    cosdec = np.cos(dec_rad)

    # Dot products with local East and North basis at boresight.
    dot_east = bx0 * (-sinra) + bx1 * cosra
    dot_north = bx0 * (-sindec * cosra) + bx1 * (-sindec * sinra) + bx2 * cosdec

    return np.degrees(np.arctan2(dot_east, dot_north))


def quat_to_roll_like_pr(quat_xyzw, ra_deg, dec_deg):
    """PR-compatible roll as sky position angle of body X-axis around boresight."""
    return quat_to_roll_like_pr_rad(quat_xyzw, np.radians(ra_deg), np.radians(dec_deg))


def equatorial_to_ecliptic_lon_lat(ra_deg, dec_deg):
    """Ecliptic longitude and latitude [deg] from J2000 equatorial coordinates [deg].

    Matches Coordinate::Equatorial2Ecliptic(MJD_J2000, ...) in coordinate.cc.
    At MJD_J2000 jdcent = 0, so eps() returns exactly 23.439292 deg.
    
    Returns (lon, lat) both in [deg], lon wrapped to [0, 360).
    """
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    
    # Ecliptic latitude
    sin_beta = np.sin(dec) * np.cos(ECL_OBLIQUITY_RAD) - np.cos(dec) * np.sin(ECL_OBLIQUITY_RAD) * np.sin(ra)
    lat = np.degrees(np.arcsin(np.clip(sin_beta, -1.0, 1.0)))
    
    # Ecliptic longitude
    cos_eps = np.cos(ECL_OBLIQUITY_RAD)
    sin_eps = np.sin(ECL_OBLIQUITY_RAD)
    lon_rad = np.arctan2(np.sin(ra) * cos_eps + np.tan(dec) * sin_eps, np.cos(ra))
    lon = np.degrees(lon_rad)
    lon = wrap_360(lon)
    
    return lon, lat


def equatorial_to_ecliptic_lat(ra_deg, dec_deg):
    """Ecliptic latitude [deg] only (for backwards compatibility in aberration calculation)."""
    _, lat = equatorial_to_ecliptic_lon_lat(ra_deg, dec_deg)
    return lat


def xy_to_rtheta_deg(x_arcmin, y_arcmin):
    """Replicate Tel2SS::XYtoRT (returns r [arcmin], theta [deg])."""
    r_rad = np.arccos(np.cos(x_arcmin * DEG2RAD / 60.0) * np.cos(y_arcmin * DEG2RAD / 60.0))

    with np.errstate(divide="ignore", invalid="ignore"):
        sintheta = np.sin(y_arcmin * DEG2RAD / 60.0) / np.sin(r_rad)
    sintheta = np.clip(sintheta, -1.0, 1.0)
    theta_deg = np.degrees(np.arcsin(sintheta))

    r_arcmin = np.degrees(r_rad) * 60.0
    theta_deg = np.where(x_arcmin < 0.0, -180.0 - theta_deg, theta_deg)

    return r_arcmin, theta_deg


def quat_from_xy_arcmin(x_arcmin, y_arcmin):
    """Replicate Tel2SS::QuatXY for rotation=0, output xyzw."""
    r_arcmin, theta_deg = xy_to_rtheta_deg(x_arcmin, y_arcmin)

    r = r_arcmin * DEG2RAD / 60.0
    th = theta_deg * DEG2RAD

    qx = np.sin(-th) * np.sin(r / 2.0)
    qy = np.cos(-th) * np.sin(r / 2.0)
    qz = np.zeros_like(qx)
    qw = np.cos(r / 2.0)
    return np.stack((qx, qy, qz, qw), axis=1)


@functools.lru_cache(maxsize=16)
def _load_pixel_positions_cached(detector, pixelpos_dir_resolved):
    pixelpos_dir = Path(pixelpos_dir_resolved)
    if detector == "SW":
        fn = pixelpos_dir / "ofsts_001_002.txt"
        n_det = 100
    else:
        fn = pixelpos_dir / "ofstl_001_003.txt"
        n_det = 75

    arr = np.loadtxt(fn, comments="#", dtype=np.float64)
    if arr.shape[0] < n_det:
        raise RuntimeError(f"Expected at least {n_det} rows in {fn}, found {arr.shape[0]}")

    arr = arr[:n_det]

    # Match fis_pixpos.cc: first column -> y, second -> x, both negated.
    y = -arr[:, 0]
    x = -arr[:, 1]
    # Keep cached arrays immutable by convention.
    x.setflags(write=False)
    y.setflags(write=False)
    return x, y


def load_pixel_positions(detector, pixelpos_dir):
    """Load detector pixel offsets, using a process-local cache."""
    resolved = str(Path(pixelpos_dir).expanduser().resolve())
    x, y = _load_pixel_positions_cached(detector, resolved)
    return x, y


def infer_detector_type_from_filename(filename):
    """Infer detector type from AKARI file naming convention.

    Returns "SW"/"LW" when possible, otherwise None.
    """
    name = Path(filename).name.upper()
    if "_SW_" in name:
        return "SW"
    if "_LW_" in name:
        return "LW"
    return None


def apply_longcorrect_offsets(x, y, detector, packetid0):
    """Apply LongCorrect x/y offsets (no aberration) for all channels."""
    x = x.copy()
    y = y.copy()

    if detector == "SW":
        wide_end = 60
        in_scan_offset = 130.515 / 60.0
        cross_scan_offset = 7.603 / 60.0
        in_scan_cds_offset = 20.0 / 60.0
    else:
        wide_end = 45
        in_scan_offset = 186.338 / 60.0
        in_scan_offset += -43.0 / 60.0
        cross_scan_offset = 20.920 / 60.0
        in_scan_cds_offset = 10.0 / 60.0

    in_scan_offset_band = np.array([-0.009 * 60.0, -0.003 * 60.0, 0.0115 * 60.0, 0.022 * 60.0])
    in_scan_cds_offset_band = np.array([0.0185 * 60.0, 0.011 * 60.0, 0.002 * 60.0, 0.0 * 60.0])

    # enum Band {N60, WideS, WideL, N160}
    ch = np.arange(len(x))
    if detector == "SW":
        # SW: wide channels first -> WideS, then N60
        band_idx = np.where(ch < wide_end, 1, 0)
    else:
        # LW: wide channels first -> WideL, then N160
        band_idx = np.where(ch < wide_end, 2, 3)

    x += in_scan_offset
    x += in_scan_offset_band[band_idx]
    y += cross_scan_offset

    cds_active = packetid0 == 0x41
    if cds_active:
        x += in_scan_cds_offset
        x += in_scan_cds_offset_band[band_idx]
        # cross_scan_cds_offset is 0.0 in LongCorrect

    return x, y, cds_active


def get_attitude_dependent_alignment_offsets_arcsec(detector, attitude_source):
    """Return attitude-dependent alignment offsets [crossW, inW, crossN, inN] in arcsec.

    Values follow idl/ASTRO-F/reduction/position/pro/pos_bs2detpix.pro
    for default (non-FTS) operation.
    """
    src = attitude_source.upper()
    if src == "PR":
        if detector == "SW":
            return np.array([-9.5, -121.0, -8.0, -126.0], dtype=np.float64)
        return np.array([-19.2, -163.0, -19.2, -163.0], dtype=np.float64)
    # Treat GADS (and any non-PR source) as the IDL AOCU/GADS branch.
    if detector == "SW":
        return np.array([22.8, -75.6, 22.8, -75.6], dtype=np.float64)
    return np.array([12.5, -112.5, 12.5, -112.5], dtype=np.float64)


def apply_attitude_dependent_alignment_offsets(x, y, detector, attitude_source):
    """Apply attitude-dependent alignment offsets to detector positions.

    x, y are arcmin detector offsets as loaded from CCF proxy tables.
    In pos_bs2detpix.pro these map to dx (in-scan) and dy (cross-scan) and are
    offset separately for wide and narrow detector groups.
    """
    dx = x.copy()
    dy = y.copy()

    off = get_attitude_dependent_alignment_offsets_arcsec(detector, attitude_source)
    if detector == "SW":
        iw = np.arange(60)
        inr = np.arange(60, 100)
    else:
        iw = np.arange(45)
        inr = np.arange(45, 75)

    # IDL: dx -= in_scan/60, dy -= cross_scan/60 for wide and narrow separately.
    dx[iw] -= off[1] / 60.0
    dy[iw] -= off[0] / 60.0
    dx[inr] -= off[3] / 60.0
    dy[inr] -= off[2] / 60.0

    return dx, dy


def compute_seasonal_aberration_arcsec(aftime, boresight_ra_deg, boresight_dec_deg):
    """Compute the seasonal/AFTIME aberration model, returning (abx, aby) in arcsec.

    Follows the equations around aberration_x / aberration_y in:
    idl/ASTRO-F/reduction/position/pro/pos_bs2detpix.pro
    """
    elon, elat = equatorial_to_ecliptic_lon_lat(boresight_ra_deg, boresight_dec_deg)

    yphase = (aftime - 6912000.0) / 24.0 / 3600.0 / 365.25
    yphase = (yphase - np.floor(yphase)) * 360.0

    lmdp = yphase - 90.0
    lmdp = np.where(lmdp < 0.0, lmdp + 360.0, lmdp)

    dlmd = elon - lmdp
    dlmd = np.where(dlmd >= 180.0, dlmd - 360.0, dlmd)
    dlmd = np.where(dlmd <= -180.0, dlmd + 360.0, dlmd)

    sdir = np.ones_like(elon, dtype=np.float64)
    sdir[np.abs(dlmd) > 90.0] = -1.0

    yphase_rad = np.radians(yphase)
    oangle = (23.4 + 8.24 * np.sin(yphase_rad)) * DEG2RAD
    sc_y = np.sin(oangle) * np.cos(np.radians(elat))

    aberration_x = 5.0 * np.sqrt(np.clip(1.0 - sc_y**2, 0.0, None)) + 20.5 * np.sin(np.radians(elat)) + 2.0 * sdir
    aberration_y = -5.0 * sc_y * sdir + 2.0 * sdir
    return aberration_x, aberration_y


def infer_detector_type(data):
    flux = data["FIS_OBS"]["FLUX"].read()
    n_det = flux.shape[1]
    if n_det == 100:
        return "SW"
    if n_det == 75:
        return "LW"
    raise RuntimeError(f"Unable to infer detector type from FLUX shape {flux.shape}")


def boresight_solution(quat):
    vrot = rotate_unit_z_by_quat(quat)
    al = np.arctan2(vrot[:, 1], vrot[:, 0])
    dl = np.arcsin(np.clip(vrot[:, 2], -1.0, 1.0))
    return wrap_360(np.degrees(al)), np.degrees(dl)


def _focal_plane_solution_no_aberration_qfov(quat, q_fov, compute_roll=True):
    """Internal no-aberration kernel with precomputed detector quaternions."""
    nsamp = quat.shape[0]
    ndet = q_fov.shape[0]

    out_ra = np.empty((nsamp, ndet), dtype=np.float64)
    out_dec = np.empty((nsamp, ndet), dtype=np.float64)
    out_roll = np.empty((nsamp, ndet), dtype=np.float64) if compute_roll else None

    # Rotate unit boresight by each fixed detector quaternion once.
    v_det = rotate_unit_z_by_quat(q_fov)  # (ndet, 3)

    # Apply time-dependent attitude rotation to all detectors at once.
    r = quat_to_matrix_xyzw(quat)  # (nsamp, 3, 3)
    vrot = np.einsum("sij,dj->sdi", r, v_det, optimize=True)  # (nsamp, ndet, 3)

    al = np.arctan2(vrot[:, :, 1], vrot[:, :, 0])
    dl = np.arcsin(np.clip(vrot[:, :, 2], -1.0, 1.0))
    out_ra[:, :] = wrap_360(np.degrees(al))
    out_dec[:, :] = np.degrees(dl)

    if not compute_roll:
        return out_ra, out_dec, None

    for j in range(ndet):
        # Compose: q_full = q_ads * q_fov  (attitude first, then detector offset).
        # This is the rotation mapping body frame -> sky frame: correct combined roll.
        qj = q_fov[j]
        qx = quat[:, 3] * qj[0] + quat[:, 0] * qj[3] + quat[:, 1] * qj[2] - quat[:, 2] * qj[1]
        qy = quat[:, 3] * qj[1] - quat[:, 0] * qj[2] + quat[:, 1] * qj[3] + quat[:, 2] * qj[0]
        qz = quat[:, 3] * qj[2] + quat[:, 0] * qj[1] - quat[:, 1] * qj[0] + quat[:, 2] * qj[3]
        qw = quat[:, 3] * qj[3] - quat[:, 0] * qj[0] - quat[:, 1] * qj[1] - quat[:, 2] * qj[2]
        q_full = np.stack((qx, qy, qz, qw), axis=1)
        out_roll[:, j] = quat_to_roll_like_pr(q_full, out_ra[:, j], out_dec[:, j])

    return out_ra, out_dec, out_roll


def parse_detector_spec(spec, n_det):
    """Parse detector selection like '1,5,10-12' into 0-based unique indices."""
    parts = [p.strip() for p in spec.split(",") if p.strip()]
    out = []

    for part in parts:
        if "-" in part:
            a_str, b_str = part.split("-", 1)
            a = int(a_str)
            b = int(b_str)
            if a <= b:
                rng = range(a, b + 1)
            else:
                rng = range(a, b - 1, -1)
            for det in rng:
                if 1 <= det <= n_det:
                    out.append(det - 1)
        else:
            det = int(part)
            if 1 <= det <= n_det:
                out.append(det - 1)

    # Keep order, remove duplicates.
    seen = set()
    uniq = []
    for idx in out:
        if idx not in seen:
            seen.add(idx)
            uniq.append(idx)
    return uniq


def _focal_plane_solution_with_aberration(quat, x_arcmin, y_arcmin, compute_roll=True):
    """Internal focal-plane kernel matching LongCorrect's two-pass aberration correction.

    Replicates the LongCorrect.cc pointing loop:
      1. Preliminary pointing with QfovDet = QuatXY(x, y) (no aberration)
      2. Convert preliminary ecliptic latitude beta_ecl
      3. aberration = (5.1 + 20.5 * sin(beta_ecl)) / 60  [arcmin]
      4. Final pointing with q_ab = QuatXY(x - aberration, y)
      5. Roll from the composed quaternion at the final pointing

    x_arcmin, y_arcmin: (ndet,) arrays with all LongCorrect offsets applied
    Returns ra, dec, roll each with shape (nsamp, ndet) [degrees]
    """
    nsamp = quat.shape[0]
    ndet = len(x_arcmin)
    # Precompute attitude rotation matrices once; reused for all detectors.
    r_att = quat_to_matrix_xyzw(quat)  # (nsamp, 3, 3)

    out_ra = np.empty((nsamp, ndet), dtype=np.float64)
    out_dec = np.empty((nsamp, ndet), dtype=np.float64)
    out_roll = np.empty((nsamp, ndet), dtype=np.float64) if compute_roll else None

    for j in range(ndet):
        # ---- Step 1: preliminary no-aberration pointing ----
        q_fov_j = quat_from_xy_arcmin(
            np.array([x_arcmin[j]]), np.array([y_arcmin[j]])
        )  # (1, 4)
        v_det = rotate_unit_z_by_quat(q_fov_j)[0]  # (3,)
        v_pre = np.einsum("sij,j->si", r_att, v_det, optimize=True)  # (nsamp, 3)
        al_pre = np.arctan2(v_pre[:, 1], v_pre[:, 0])
        dl_pre = np.arcsin(np.clip(v_pre[:, 2], -1.0, 1.0))

        # ---- Step 2: ecliptic latitude and per-sample aberration ----
        ecl_beta = equatorial_to_ecliptic_lat(
            np.degrees(al_pre), np.degrees(dl_pre)
        )  # (nsamp,)
        aberration = (5.1 + 20.5 * np.sin(np.radians(ecl_beta))) / 60.0  # arcmin

        # ---- Step 3: aberration-corrected per-sample detector quaternion ----
        q_ab_j = quat_from_xy_arcmin(
            x_arcmin[j] - aberration, np.full(nsamp, y_arcmin[j])
        )  # (nsamp, 4)

        # ---- Step 4: final pointing ----
        v_step1 = rotate_unit_z_by_quat(q_ab_j)   # (nsamp, 3): apply q_ab to boresight
        v_final = np.einsum("sij,sj->si", r_att, v_step1, optimize=True)  # (nsamp, 3)
        al = np.arctan2(v_final[:, 1], v_final[:, 0])
        dl = np.arcsin(np.clip(v_final[:, 2], -1.0, 1.0))
        out_ra[:, j] = wrap_360(np.degrees(al))
        out_dec[:, j] = np.degrees(dl)

        if compute_roll:
            # ---- Step 5: roll from composed quaternion q_ads * q_ab ----
            # Correct body-to-sky rotation is q_ads * q_ab (not q_ab * q_ads).
            # q_ab_j is (nsamp, 4) here (per-sample aberration-adjusted detector quaternion).
            q_full = quat_multiply_xyzw(quat, q_ab_j)
            out_roll[:, j] = quat_to_roll_like_pr_rad(q_full, al, dl)

    return out_ra, out_dec, out_roll


def normalize_offset_style(offset_style):
    """Normalize public offset-style labels to internal canonical values."""
    if offset_style is None:
        return None
    style = offset_style.lower()
    if style in {"none", "raw"}:
        return None
    if style in {"long", "longcorrect", "static"}:
        return "long"
    if style in {"attitude", "seasonal-attitude", "idl"}:
        return "attitude"
    raise ValueError(f"offset_style must be one of 'long', 'attitude', or None; got {offset_style!r}")


def normalize_aberration_style(aberration_style):
    """Normalize public aberration-style labels to internal canonical values."""
    if aberration_style is None:
        return "long"
    style = aberration_style.lower()
    if style in {"long", "longcorrect", "ecliptic"}:
        return "long"
    if style in {"seasonal", "aftime", "idl"}:
        return "seasonal"
    raise ValueError(
        f"aberration_style must be one of 'long' or 'seasonal'; got {aberration_style!r}"
    )


def focal_plane_solution(
    quat,
    x_arcmin,
    y_arcmin,
    aberration=True,
    compute_roll=True,
    offset_style=None,
    detector=None,
    packetid0=None,
    attitude_source=None,
    aberration_style=None,
    aftime=None,
):
    """Unified focal-plane pointing solution.

    Both the C++ (LongCorrect) and IDL (pos_bs2detpix) pipelines reduce to
    this function with different flag combinations.  The results are
    numerically identical when the same flags are used — all differences
    between the two historical pipelines come entirely from the choice of
    offset formula and/or aberration model.

    Pipeline recipes
    ----------------
    **C++ / LongCorrect pipeline** (attitude-agnostic offsets, ecliptic-lat aberration)::

        focal_plane_solution(
            quat, x_arcmin, y_arcmin,
            offset_style="long",  detector=detector, packetid0=packetid0,
            aberration_style="long",
        )

    **Attitude-dependent / seasonal pipeline** (attitude-dependent offsets, seasonal AFTIME aberration)::

        focal_plane_solution(
            quat, x_arcmin, y_arcmin,
            offset_style="attitude",  detector=detector, attitude_source=attitude_source,
            aberration_style="seasonal", aftime=aftime,
        )

    Args
    ----
    quat : (nsamp, 4) attitude quaternion array in xyzw convention.
    x_arcmin, y_arcmin : (ndet,) detector positions in arcmin.
        When offset_style is None the positions must already have LongCorrect
        offsets applied by the caller (legacy behaviour).  When offset_style
        is set the positions must be raw CCF base positions (before any
        alignment offset) and offsets are applied internally.
    aberration : bool
        If True, apply an aberration correction (style selected by
        aberration_style).  If False, no aberration is applied regardless.
    compute_roll : bool
        If False, skip roll computation for speed.
    offset_style : None | "long" | "attitude"
        None   – positions used as-is (caller pre-applied Long offsets).
        "long" – apply LongCorrect attitude-agnostic offsets; requires
                 `detector` and `packetid0`.
        "attitude" – apply pos_bs2detpix attitude-dependent offsets; requires
                     `detector` and `attitude_source`.
    detector : "SW" | "LW"
        Required when offset_style is not None.
    packetid0 : int
        First PACKETID value.  Required when offset_style="long".
    attitude_source : "PR" | "GADS"
        Required when offset_style="attitude".
    aberration_style : None | "long" | "seasonal"
        None or "long" – LongCorrect two-pass ecliptic-latitude model (default).
        "seasonal"     – seasonal/AFTIME model (both x and y components);
                         requires `aftime`.
    aftime : (nsamp,) array
        AKARI mission time in seconds.  Required when aberration_style="seasonal".
    """
    offset_style = normalize_offset_style(offset_style)
    aberration_style = normalize_aberration_style(aberration_style)

    if offset_style is not None:
        if offset_style == "long":
            if detector is None or packetid0 is None:
                raise ValueError("offset_style='long' requires detector and packetid0")
            x_arcmin, y_arcmin, _ = apply_longcorrect_offsets(x_arcmin, y_arcmin, detector, packetid0)
        elif offset_style == "attitude":
            if detector is None or attitude_source is None:
                raise ValueError("offset_style='attitude' requires detector and attitude_source")
            x_arcmin, y_arcmin = apply_attitude_dependent_alignment_offsets(x_arcmin, y_arcmin, detector, attitude_source)

    use_seasonal_aberration = aberration_style == "seasonal"

    if not aberration:
        q_fov = quat_from_xy_arcmin(x_arcmin, y_arcmin)
        return _focal_plane_solution_no_aberration_qfov(quat, q_fov, compute_roll=compute_roll)

    if use_seasonal_aberration:
        if aftime is None:
            raise ValueError("aberration_style='seasonal' requires aftime")
        # Delegate to the IDL aberration path using already-offset positions
        # (pass offset_style=None so no further offsets are applied).
        ra_bore, dec_bore = boresight_solution(quat)
        aberr_x_arcsec, aberr_y_arcsec = compute_seasonal_aberration_arcsec(aftime, ra_bore, dec_bore)
        nsamp = quat.shape[0]
        ndet = len(x_arcmin)
        # Precompute attitude rotation matrices once; reused for all detectors.
        r_att = quat_to_matrix_xyzw(quat)  # (nsamp, 3, 3)
        out_ra = np.empty((nsamp, ndet), dtype=np.float64)
        out_dec = np.empty((nsamp, ndet), dtype=np.float64)
        out_roll = np.empty((nsamp, ndet), dtype=np.float64) if compute_roll else None
        for j in range(ndet):
            x_samp = x_arcmin[j] - aberr_x_arcsec / 60.0
            y_samp = y_arcmin[j] - aberr_y_arcsec / 60.0
            q_det = quat_from_xy_arcmin(x_samp, y_samp)
            v_step1 = rotate_unit_z_by_quat(q_det)
            v_final = np.einsum("sij,sj->si", r_att, v_step1, optimize=True)
            al = np.arctan2(v_final[:, 1], v_final[:, 0])
            dl = np.arcsin(np.clip(v_final[:, 2], -1.0, 1.0))
            out_ra[:, j] = wrap_360(np.degrees(al))
            out_dec[:, j] = np.degrees(dl)
            if compute_roll:
                q_full = quat_multiply_xyzw(quat, q_det)
                out_roll[:, j] = quat_to_roll_like_pr_rad(q_full, al, dl)
        return out_ra, out_dec, out_roll

    # Default: LongCorrect two-pass aberration.
    return _focal_plane_solution_with_aberration(quat, x_arcmin, y_arcmin, compute_roll=compute_roll)


def run_sanity_checks(quat, ra_fp, dec_fp, roll_fp, plot_samples=5000):
    """Three orientation sanity checks for per-pixel roll.

    1. Identity test: zero-offset detector roll must equal boresight roll to float precision.
    2. Roll spread across the focal plane at each sample must be << 1 deg.
    3. Roll vs scan-trajectory PA: body X-axis is cross-scan, so roll ≈ scan_PA ± 90 deg.
    """
    nsamp, ndet = roll_fp.shape

    # ------------------------------------------------------------------ #
    # Boresight reference
    # ------------------------------------------------------------------ #
    ra_bore, dec_bore = boresight_solution(quat)
    roll_bore = quat_to_roll_like_pr(quat, ra_bore, dec_bore)

    # ------------------------------------------------------------------ #
    # Check 1: identity quaternion → roll == boresight roll
    # ------------------------------------------------------------------ #
    q_identity = np.array([[0.0, 0.0, 0.0, 1.0]], dtype=np.float64)
    _, _, roll_id = _focal_plane_solution_no_aberration_qfov(quat, q_identity)
    max_id_diff = np.max(np.abs(wrap_pm180(roll_id[:, 0] - roll_bore)))

    print("\n=== Sanity check 1: identity (zero-offset) detector roll == boresight roll ===")
    print(f"  max |roll_ident - roll_bore| = {max_id_diff:.2e} deg  [should be ~0]")

    # ------------------------------------------------------------------ #
    # Check 2: roll spread across focal plane each sample
    # ------------------------------------------------------------------ #
    roll_spread = roll_fp.max(axis=1) - roll_fp.min(axis=1)

    print(f"\n=== Sanity check 2: roll spread across {ndet} detectors per sample ===")
    print(f"  mean spread = {roll_spread.mean():.4f} deg")
    print(f"  max  spread = {roll_spread.max():.4f} deg")
    print(f"  [expected << 1 deg for the ~10-arcmin FIS focal plane]")

    # ------------------------------------------------------------------ #
    # Check 3: roll vs scan-trajectory position angle
    # ------------------------------------------------------------------ #
    ra_r = np.radians(ra_bore)
    dec_r = np.radians(dec_bore)
    dra = np.diff(ra_r)
    dra = (dra + np.pi) % (2.0 * np.pi) - np.pi      # wrap to [-pi, pi]
    ddec = np.diff(dec_r)

    scan_pa = np.degrees(np.arctan2(np.cos(dec_r[:-1]) * dra, ddec))  # east of north
    roll_at_scan = roll_bore[:-1]

    diff_p90 = wrap_pm180(roll_at_scan - scan_pa - 90.0)
    diff_m90 = wrap_pm180(roll_at_scan - scan_pa + 90.0)
    residual = np.where(np.abs(diff_p90) < np.abs(diff_m90), diff_p90, diff_m90)
    frac_p90 = np.mean(np.abs(diff_p90) < np.abs(diff_m90))

    # Also check whether body X is in-scan rather than cross-scan.
    diff_inscan = wrap_pm180(roll_at_scan - scan_pa)

    print("\n=== Sanity check 3: roll vs scan-trajectory PA ===")
    print(f"  mean |roll - scan_PA|        = {np.abs(diff_inscan).mean():.4f} ± {np.abs(diff_inscan).std():.4f} deg  [in-scan body X]")
    print(f"  mean |roll - (scan_PA ± 90)| = {np.abs(residual).mean():.4f} ± {np.abs(residual).std():.4f} deg  [cross-scan body X]")
    print(f"  fraction of samples closer to scan_PA+90 (vs scan_PA-90): {frac_p90:.2f}")
    if np.abs(diff_inscan).mean() < np.abs(residual).mean():
        print("  Conclusion: body X is IN-SCAN  → roll = scan PA  (beam major axis PA = roll)")
    else:
        print("  Conclusion: body X is CROSS-SCAN → roll = scan PA ± 90  (beam major axis PA = roll ± 90)")

    # ------------------------------------------------------------------ #
    # Plot: roll vs scan PA
    # ------------------------------------------------------------------ #
    nplot = min(plot_samples, len(scan_pa))
    _, axes = plt.subplots(3, 1, sharex=True, figsize=(10, 8))

    axes[0].plot(roll_bore[:nplot], label="Boresight roll", alpha=0.9)
    axes[0].plot(scan_pa[:nplot] + 90.0, label="Scan PA + 90°", ls="--", alpha=0.8)
    axes[0].plot(scan_pa[:nplot], label="Scan PA", ls=":", alpha=0.6)
    axes[0].set_ylabel("Angle [deg]")
    axes[0].legend(loc="best")

    axes[1].plot(residual[:nplot], label="roll − (scan_PA ± 90°)")
    axes[1].set_ylabel("Residual [deg]")
    axes[1].legend(loc="best")

    axes[2].plot(roll_spread[:nplot], label="Roll spread across focal plane")
    axes[2].set_ylabel("Spread [deg]")
    axes[2].set_xlabel("Sample")
    axes[2].legend(loc="best")

    plt.suptitle("Sanity checks: roll vs scan PA and focal-plane spread")
    plt.tight_layout()
    plt.show()


def load_pickle_reference(fits_filename, detector_type, fits_dir_path):
    """Load reference lon/lat pickle files from master's thesis work.
    
    Expects input FITS filename like 'FIS_SW_20061231020000_gb.fits.gz' and 
    looks for pickles: FIS_{SW/LW}_20061231020000_gb_{lon/lat}_gads.pkl
    
    Returns (lon, lat) as (nsamp, ndet) numpy arrays.
    """
    import pickle
    
    fits_dir = Path(fits_dir_path).expanduser()
    
    # Extract the date/time and mode part from the FITS filename
    # e.g., 'FIS_SW_20061231020000_gb.fits.gz' -> '20061231020000_gb'
    fits_name = Path(fits_filename).name
    # Remove 'FIS_SW_' or 'FIS_LW_' prefix and '.fits.gz' suffix
    core = fits_name.replace('FIS_SW_', '').replace('FIS_LW_', '').replace('.fits.gz', '')
    
    lon_file = fits_dir / f"FIS_{detector_type}_{core}_lon_gads.pkl"
    lat_file = fits_dir / f"FIS_{detector_type}_{core}_lat_gads.pkl"
    
    if not lon_file.exists() or not lat_file.exists():
        print(f"WARNING: Pickle reference files not found:")
        print(f"  {lon_file}")
        print(f"  {lat_file}")
        return None, None
    
    with open(lon_file, 'rb') as f:
        lon_list = pickle.load(f)  # list[detector][sample]
    with open(lat_file, 'rb') as f:
        lat_list = pickle.load(f)
    
    # Convert list[det][samp] to numpy (nsamp, ndet)
    ndet = len(lon_list)
    nsamp = len(lon_list[0])
    lon = np.zeros((nsamp, ndet), dtype=np.float64)
    lat = np.zeros((nsamp, ndet), dtype=np.float64)
    
    for d in range(ndet):
        lon[:, d] = np.array(lon_list[d], dtype=np.float64)
        lat[:, d] = np.array(lat_list[d], dtype=np.float64)
    
    return lon, lat


def compare_to_reference(ra_fp, dec_fp, lon_ref, lat_ref, det_indices=None, plot_samples=5000):
    """Compare focal-plane pointing to reference pickle data.
    
    Reference pickle files actually contain equatorial (RA, Dec) coordinates,
    despite being labeled with 'lon'/'lat' naming (legacy artifact).
    """
    # Direct equatorial coordinate comparison
    lon_py = ra_fp  # Reference "lon" is actually RA
    lat_py = dec_fp  # Reference "lat" is actually Dec
    
    # Compute differences (wrap RA to [-180, 180] for proper circular mean)
    dlon = wrap_pm180(lon_py - lon_ref)
    dlat = lat_py - lat_ref
    
    print(f"\n=== Comparison to reference pickle ===")
    print(f"Lon difference [deg]:")
    print(f"  mean: {np.nanmean(dlon):.6e}, std: {np.nanstd(dlon):.6e}")
    print(f"  min:  {np.nanmin(dlon):.6e}, max: {np.nanmax(dlon):.6e}")
    print(f"Lat difference [deg]:")
    print(f"  mean: {np.nanmean(dlat):.6e}, std: {np.nanstd(dlat):.6e}")
    print(f"  min:  {np.nanmin(dlat):.6e}, max: {np.nanmax(dlat):.6e}")
    
    # Per-detector stats
    if det_indices is None:
        det_indices = [0]
    
    for idx in det_indices:
        dlon_det = dlon[:, idx]
        dlat_det = dlat[:, idx]
        print(f"Detector {idx+1}:")
        print(f"  dLon: mean={np.nanmean(dlon_det):.6e}±{np.nanstd(dlon_det):.6e} deg")
        print(f"  dLat: mean={np.nanmean(dlat_det):.6e}±{np.nanstd(dlat_det):.6e} deg")
    
    # Plot comparison for first few detectors
    if det_indices:
        nplot = min(plot_samples, dlon.shape[0])
        _, axes = plt.subplots(2, 1, sharex=True, figsize=(10, 6))
        
        for idx in det_indices:
            axes[0].plot(wrap_pm180(dlon[:nplot, idx]), label=f"Det {idx+1}", alpha=0.8)
            axes[1].plot(dlat[:nplot, idx], label=f"Det {idx+1}", alpha=0.8)
        
        axes[0].set_ylabel("dRA [deg]")
        axes[1].set_ylabel("dDec [deg]")
        axes[1].set_xlabel("Sample")
        axes[0].legend(loc="best")
        axes[1].legend(loc="best")
        
        plt.suptitle("Comparison to reference: Equatorial coordinate differences")
        plt.tight_layout()
        plt.show()


def compute_diff_stats(dra_deg, ddec_deg):
    """Return scalar summary statistics for RA/Dec differences."""
    return {
        "dra_mean": float(np.nanmean(dra_deg)),
        "dra_std": float(np.nanstd(dra_deg)),
        "dra_p95_abs": float(np.nanpercentile(np.abs(dra_deg), 95.0)),
        "dra_max_abs": float(np.nanmax(np.abs(dra_deg))),
        "ddec_mean": float(np.nanmean(ddec_deg)),
        "ddec_std": float(np.nanstd(ddec_deg)),
        "ddec_p95_abs": float(np.nanpercentile(np.abs(ddec_deg), 95.0)),
        "ddec_max_abs": float(np.nanmax(np.abs(ddec_deg))),
    }


def run_pointing_model_for_config(
    data,
    det_type,
    pixelpos_dir,
    attitude,
    offset_style,
    aberration_style,
    use_aberration=True,
):
    """Run the unified focal-plane model for one configuration."""
    quat_raw = data[attitude]["QUATERNION"].read().astype(np.float64)
    quat, valid = normalize_quat(quat_raw)

    x0, y0 = load_pixel_positions(det_type, pixelpos_dir)

    packetid = data["FIS_OBS"]["PACKETID"].read()
    packetid0 = int(packetid[0])
    aftime = data["FIS_OBS"]["AFTIME"].read().astype(np.float64)[valid]

    ra, dec, roll = focal_plane_solution(
        quat,
        x0,
        y0,
        aberration=use_aberration,
        compute_roll=False,
        offset_style=offset_style,
        detector=det_type,
        packetid0=packetid0,
        attitude_source=attitude,
        aberration_style=aberration_style,
        aftime=aftime,
    )

    return ra, dec, roll, valid


def subset_model_to_mask(model_arr, valid_mask, target_mask):
    """Subset model array (already compressed by valid_mask) to target original-mask rows."""
    orig_to_comp = -np.ones(len(valid_mask), dtype=np.int64)
    orig_to_comp[np.where(valid_mask)[0]] = np.arange(np.sum(valid_mask), dtype=np.int64)
    comp_idx = orig_to_comp[np.where(target_mask)[0]]
    return model_arr[comp_idx]


def run_comparison_suite(data, filename, det_type, pixelpos_dir, suite_csv=None):
    """Run a 6-configuration pointing comparison suite against reference pickle data.

    Reference pickles are treated as equatorial (RA/Dec) despite lon/lat naming.
    """
    lon_ref, lat_ref = load_pickle_reference(filename, det_type, Path(filename).parent)
    if lon_ref is None or lat_ref is None:
        raise RuntimeError("Cannot run comparison suite without reference pickle files.")

    available_ext = {h.get_extname() for h in data}
    configs = [
        {"name": "PR_long_offset_long_ab", "attitude": "PR", "offset_style": "long", "aberration_style": "long", "use_aberration": True},
        {"name": "PR_long_offset_noab", "attitude": "PR", "offset_style": "long", "aberration_style": "long", "use_aberration": False},
        {"name": "PR_raw_noab", "attitude": "PR", "offset_style": None, "aberration_style": "long", "use_aberration": False},
        {"name": "GADS_long_offset_long_ab", "attitude": "GADS", "offset_style": "long", "aberration_style": "long", "use_aberration": True},
        {"name": "GADS_long_offset_noab", "attitude": "GADS", "offset_style": "long", "aberration_style": "long", "use_aberration": False},
        {"name": "GADS_raw_noab", "attitude": "GADS", "offset_style": None, "aberration_style": "long", "use_aberration": False},
        {"name": "PR_attitude_offset_seasonal_ab", "attitude": "PR", "offset_style": "attitude", "aberration_style": "seasonal", "use_aberration": True},
        {"name": "GADS_attitude_offset_seasonal_ab", "attitude": "GADS", "offset_style": "attitude", "aberration_style": "seasonal", "use_aberration": True},
    ]

    qpr_raw = data["PR"]["QUATERNION"].read().astype(np.float64)
    pr_valid = np.linalg.norm(qpr_raw, axis=1) > 0.0

    rows = []
    det_rows = []
    print("\n=== Automated Comparison Suite (Reference as RA/Dec) ===")
    for cfg in configs:
        if cfg["attitude"] not in available_ext:
            print(f"SKIP {cfg['name']}: extension {cfg['attitude']} not present")
            continue

        ra_fp, dec_fp, _roll_fp, valid = run_pointing_model_for_config(
            data,
            det_type,
            pixelpos_dir,
            attitude=cfg["attitude"],
            offset_style=cfg["offset_style"],
            aberration_style=cfg["aberration_style"],
            use_aberration=cfg["use_aberration"],
        )

        if cfg["attitude"] == "PR":
            compare_mask = valid
        else:
            # Reference pickle alignment is PR-based; compare GADS on common timeline.
            compare_mask = valid & pr_valid

        ref_ra = lon_ref[compare_mask, :]
        ref_dec = lat_ref[compare_mask, :]
        ra_cmp = subset_model_to_mask(ra_fp, valid, compare_mask)
        dec_cmp = subset_model_to_mask(dec_fp, valid, compare_mask)

        ndet = min(ra_cmp.shape[1], ref_ra.shape[1])

        dra = wrap_pm180(ra_cmp[:, :ndet] - ref_ra[:, :ndet])
        ddec = dec_cmp[:, :ndet] - ref_dec[:, :ndet]
        s = compute_diff_stats(dra, ddec)

        d1_ra = wrap_pm180(ra_cmp[:, 0] - ref_ra[:, 0])
        d1_dec = dec_cmp[:, 0] - ref_dec[:, 0]
        s1 = compute_diff_stats(d1_ra, d1_dec)

        row = {
            "config": cfg["name"],
            "attitude": cfg["attitude"],
            "offset_style": cfg["offset_style"] or "none",
            "aberration_style": cfg["aberration_style"],
            "use_aberration": int(cfg["use_aberration"]),
            "nsamp_valid": int(ra_fp.shape[0]),
            "nsamp_compared": int(np.sum(compare_mask)),
            "ndet_compared": int(ndet),
            **s,
            "det1_dra_mean": s1["dra_mean"],
            "det1_dra_std": s1["dra_std"],
            "det1_ddec_mean": s1["ddec_mean"],
            "det1_ddec_std": s1["ddec_std"],
        }
        rows.append(row)

        print(
            f"{cfg['name']:20s}  nsamp={row['nsamp_compared']:6d}  "
            f"dRA(mean/std)={row['dra_mean']:.3e}/{row['dra_std']:.3e} deg  "
            f"dDec(mean/std)={row['ddec_mean']:.3e}/{row['ddec_std']:.3e} deg"
        )

        # Per-detector ranking metrics
        dra_mean_det = np.nanmean(dra, axis=0)
        dra_std_det = np.nanstd(dra, axis=0)
        ddec_mean_det = np.nanmean(ddec, axis=0)
        ddec_std_det = np.nanstd(ddec, axis=0)
        # Combined RMS metric for ranking detector agreement quality.
        det_rms = np.sqrt(dra_std_det**2 + ddec_std_det**2)

        for j in range(ndet):
            det_rows.append(
                {
                    "config": cfg["name"],
                    "detector": int(j + 1),
                    "dra_mean": float(dra_mean_det[j]),
                    "dra_std": float(dra_std_det[j]),
                    "ddec_mean": float(ddec_mean_det[j]),
                    "ddec_std": float(ddec_std_det[j]),
                    "combined_rms": float(det_rms[j]),
                }
            )

        nshow = min(5, ndet)
        order = np.argsort(det_rms)
        best_idx = order[:nshow]
        worst_idx = order[-nshow:][::-1]

        print(f"  Best {nshow} detectors by combined RMS [deg]:")
        for j in best_idx:
            print(
                f"    det {j+1:3d}: rms={det_rms[j]:.3e}, "
                f"dRA std={dra_std_det[j]:.3e}, dDec std={ddec_std_det[j]:.3e}"
            )
        print(f"  Worst {nshow} detectors by combined RMS [deg]:")
        for j in worst_idx:
            print(
                f"    det {j+1:3d}: rms={det_rms[j]:.3e}, "
                f"dRA std={dra_std_det[j]:.3e}, dDec std={ddec_std_det[j]:.3e}"
            )

    if suite_csv and rows:
        out_csv = Path(suite_csv).expanduser()
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        with open(out_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nSaved suite summary CSV: {out_csv}")

        if det_rows:
            det_csv = out_csv.with_name(out_csv.stem + "_per_detector" + out_csv.suffix)
            with open(det_csv, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(det_rows[0].keys()))
                writer.writeheader()
                writer.writerows(det_rows)
            print(f"Saved detector ranking CSV: {det_csv}")

    return rows


def describe_recipe(attitude, offset_style, aberration_style, use_aberration):
    """Return a short user-facing label for the selected pointing recipe."""
    offset_label = offset_style or "none"
    aberration_label = "none" if not use_aberration else aberration_style

    if offset_label == "long" and aberration_label == "long":
        return "LongCorrect-equivalent recipe"
    if attitude == "PR" and offset_label == "attitude" and aberration_label == "seasonal":
        return "IDL/pos_bs2detpix-equivalent recipe"
    if attitude == "GADS" and offset_label == "attitude" and aberration_label == "seasonal":
        return "GADS + attitude-offset + seasonal-aberration recipe"
    return f"Custom recipe ({attitude}, offset={offset_label}, aberration={aberration_label})"


def maybe_print_pickle_reference_warning(args):
    """Warn when selected flags are not expected to reproduce the pickle reference."""
    if not args.compare_pickle:
        return

    selected_offset = args.offset_style
    selected_aberration = "none" if args.mode != "focalplane" else args.aberration_style

    if args.attitude != "PR" or selected_offset != "attitude" or selected_aberration != "seasonal":
        print("WARNING: The downloaded pickle reference is PR-aligned and is reproduced by")
        print("         --attitude PR --offset-style attitude --aberration-style seasonal")
        if args.attitude == "GADS":
            print("         The current GADS selection is not the recipe used to generate that reference.")


def main():
    parser = argparse.ArgumentParser(
        description="Get AKARI pointing from quaternion (boresight or focal-plane offsets only).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Recipe equivalents:\n"
            "  LongCorrect: --mode focalplane --offset-style long --aberration-style long\n"
            "  IDL/pos_bs2detpix: --attitude PR --mode focalplane --offset-style attitude --aberration-style seasonal\n"
            "\n"
            "Pickle reference note:\n"
            "  The downloaded lon/lat pickle files match the PR + attitude + seasonal recipe,\n"
            "  even though the filenames end in _gads.pkl.\n"
        ),
    )
    parser.add_argument(
        "fits_file",
        nargs="?",
        default="~/FISimage/FIS_SW_20061231020000_gb.fits.gz",
        help="Input GB FITS file",
    )
    parser.add_argument(
        "--attitude",
        choices=["PR", "GADS"],
        default="PR",
        help="Quaternion source extension",
    )
    parser.add_argument(
        "--mode",
        choices=["boresight", "focalplane-noab", "focalplane"],
        default="boresight",
        help=(
            "boresight: quaternion boresight only; "
            "focalplane-noab: detector offsets, no aberration; "
            "focalplane: focal-plane solution using the selected offset/aberration models"
        ),
    )
    parser.add_argument(
        "--offset-style",
        choices=["long", "attitude", "none"],
        default="long",
        help="Detector alignment model: long=LongCorrect, attitude=attitude-dependent pos_bs2detpix, none=raw detector positions",
    )
    parser.add_argument(
        "--aberration-style",
        choices=["long", "seasonal"],
        default="long",
        help="Aberration model: long=two-pass ecliptic-latitude, seasonal=seasonal/AFTIME",
    )
    parser.add_argument(
        "--detector-type",
        choices=["SW", "LW"],
        default=None,
        help="Force detector type for focal-plane mode (default: inferred from FLUX shape)",
    )
    parser.add_argument(
        "--pixelpos-dir",
        default="~/FISimage/PixelPos",
        help="Directory containing ofsts_001_002.txt/ofstl_001_003.txt",
    )
    parser.add_argument(
        "--plot-samples",
        type=int,
        default=5000,
        help="Number of samples to plot",
    )
    parser.add_argument(
        "--plot-det",
        type=int,
        default=1,
        help="1-based detector index to plot in focal-plane mode (legacy single-detector option)",
    )
    parser.add_argument(
        "--plot-dets",
        default=None,
        help="Comma/range list of 1-based detectors to plot, e.g. '1,10,20-25'",
    )
    parser.add_argument(
        "--sanity",
        action="store_true",
        default=False,
        help="Run roll-angle sanity checks (identity test, focal-plane spread, scan-PA comparison)",
    )
    parser.add_argument(
        "--compare-pickle",
        action="store_true",
        default=False,
        help=(
            "Compare focal-plane results to reference pickle files. "
            "The downloaded reference is reproduced by --attitude PR --offset-style attitude "
            "--aberration-style seasonal."
        ),
    )
    parser.add_argument(
        "--comparison-suite",
        action="store_true",
        default=False,
        help="Run 6-configuration automated comparison suite against reference pickle files",
    )
    parser.add_argument(
        "--suite-csv",
        default=None,
        help="Optional CSV path for --comparison-suite summary output",
    )
    args = parser.parse_args()

    filename = str(Path(args.fits_file).expanduser())
    data = fitsio.FITS(filename)

    quat_raw = data[args.attitude]["QUATERNION"].read().astype(np.float64)
    quat, valid = normalize_quat(quat_raw)

    print(f"Attitude source: {args.attitude}")
    print(f"Samples (total): {len(quat_raw)}")
    print(f"Samples (valid quaternion norm>0): {len(quat)}")
    maybe_print_pickle_reference_warning(args)

    if args.comparison_suite:
        det_type_suite = args.detector_type or infer_detector_type(data)
        rows = run_comparison_suite(
            data,
            filename,
            det_type_suite,
            args.pixelpos_dir,
            suite_csv=args.suite_csv,
        )
        if not rows:
            raise RuntimeError("Comparison suite completed with no valid configurations.")
        return

    if args.mode == "boresight":
        ra_deg, dec_deg = boresight_solution(quat)
        roll_cpp_deg = np.degrees(quat_to_roll_euler(quat))
        roll_prlike_deg = quat_to_roll_like_pr(quat, ra_deg, dec_deg)

        print(f"RA range [deg]:   {ra_deg.min():.6f} .. {ra_deg.max():.6f}")
        print(f"Dec range [deg]:  {dec_deg.min():.6f} .. {dec_deg.max():.6f}")
        print(f"Roll(PR-like) range [deg]: {roll_prlike_deg.min():.6f} .. {roll_prlike_deg.max():.6f}")
        print(f"Roll(C++ Euler) range [deg]: {roll_cpp_deg.min():.6f} .. {roll_cpp_deg.max():.6f}")

        has_pr = "PR" in [h.get_extname() for h in data]
        if has_pr:
            pr_ra = data["PR"]["RA"].read().astype(np.float64)[valid]
            pr_dec = data["PR"]["DEC"].read().astype(np.float64)[valid]
            pr_roll = data["PR"]["ROLL"].read().astype(np.float64)[valid]

            dra = wrap_pm180(ra_deg - wrap_360(pr_ra))
            ddec = dec_deg - pr_dec
            droll_pr = wrap_pm180(roll_prlike_deg - pr_roll)
            droll_cpp = wrap_pm180(roll_cpp_deg - pr_roll)

            print("\nComparison against PR columns:")
            print(f"RA   diff [deg] mean/std: {np.nanmean(dra):.6e} / {np.nanstd(dra):.6e}")
            print(f"Dec  diff [deg] mean/std: {np.nanmean(ddec):.6e} / {np.nanstd(ddec):.6e}")
            print(
                f"Roll(PR-like) diff [deg] mean/std: "
                f"{np.nanmean(droll_pr):.6e} / {np.nanstd(droll_pr):.6e}"
            )
            print(
                f"Roll(C++ Euler) diff [deg] mean/std: "
                f"{np.nanmean(droll_cpp):.6e} / {np.nanstd(droll_cpp):.6e}"
            )

            nplot = min(args.plot_samples, len(pr_ra))
            _, axes = plt.subplots(3, 1, sharex=True, figsize=(10, 8))

            axes[0].plot(pr_ra[:nplot], label="PR RA", alpha=0.8)
            axes[0].plot(ra_deg[:nplot], label=f"{args.attitude} quat to RA", alpha=0.8)
            axes[0].set_ylabel("RA [deg]")
            axes[0].legend(loc="best")

            axes[1].plot(pr_dec[:nplot], label="PR Dec", alpha=0.8)
            axes[1].plot(dec_deg[:nplot], label=f"{args.attitude} quat to Dec", alpha=0.8)
            axes[1].set_ylabel("Dec [deg]")
            axes[1].legend(loc="best")

            axes[2].plot(pr_roll[:nplot], label="PR Roll", alpha=0.8)
            axes[2].plot(
                roll_prlike_deg[:nplot], label=f"{args.attitude} quat to Roll (PR-like)", alpha=0.8
            )
            axes[2].set_ylabel("Roll [deg]")
            axes[2].set_xlabel("Sample")
            axes[2].legend(loc="best")

            plt.tight_layout()
            plt.show()

    else:
        det_type = args.detector_type or infer_detector_type_from_filename(filename) or infer_detector_type(data)
        x0, y0 = load_pixel_positions(det_type, args.pixelpos_dir)
        packetid = data["FIS_OBS"]["PACKETID"].read()
        packetid0 = int(packetid[0])
        aftime = data["FIS_OBS"]["AFTIME"].read().astype(np.float64)[valid]
        offset_style = None if args.offset_style == "none" else args.offset_style
        use_ab = args.mode == "focalplane"
        ra_fp, dec_fp, roll_fp = focal_plane_solution(
            quat,
            x0,
            y0,
            aberration=use_ab,
            compute_roll=True,
            offset_style=offset_style,
            detector=det_type,
            packetid0=packetid0,
            attitude_source=args.attitude,
            aberration_style=args.aberration_style,
            aftime=aftime,
        )

        cds_active = bool(offset_style == "long" and packetid0 == 0x41)

        ndet = ra_fp.shape[1]
        if args.plot_dets:
            det_indices = parse_detector_spec(args.plot_dets, ndet)
            if not det_indices:
                raise RuntimeError("--plot-dets did not contain any valid detector IDs.")
        else:
            det_indices = [max(1, min(args.plot_det, ndet)) - 1]
        det_idx = det_indices[0]

        print(f"\nMode: {args.mode}")
        print(f"Recipe: {describe_recipe(args.attitude, args.offset_style, args.aberration_style, use_ab)}")
        print(f"Detector type: {det_type}")
        print(f"Offset model: {args.offset_style}")
        print(f"Aberration model: {'none' if not use_ab else args.aberration_style}")
        print(f"Channels: {ndet}")
        print(f"CDS offsets active: {cds_active} (packetid[0]={packetid0:#x})")
        print(f"RA(full) range [deg]:  {ra_fp.min():.6f} .. {ra_fp.max():.6f}")
        print(f"Dec(full) range [deg]: {dec_fp.min():.6f} .. {dec_fp.max():.6f}")
        print(f"Roll(full) range [deg]: {roll_fp.min():.6f} .. {roll_fp.max():.6f}")
        print(
            f"Detector {det_idx+1} RA/Dec/Roll range [deg]: "
            f"{ra_fp[:,det_idx].min():.6f}..{ra_fp[:,det_idx].max():.6f} / "
            f"{dec_fp[:,det_idx].min():.6f}..{dec_fp[:,det_idx].max():.6f} / "
            f"{roll_fp[:,det_idx].min():.6f}..{roll_fp[:,det_idx].max():.6f}"
        )
        print("Plotted detectors (1-based): " + ", ".join(str(i + 1) for i in det_indices))

        # Compare one detector track to boresight for sanity.
        ra_bore, dec_bore = boresight_solution(quat)
        roll_bore = quat_to_roll_like_pr(quat, ra_bore, dec_bore)
        dra_det = wrap_pm180(ra_fp[:, det_idx] - ra_bore)
        ddec_det = dec_fp[:, det_idx] - dec_bore
        droll_det = wrap_pm180(roll_fp[:, det_idx] - roll_bore)
        print(
            f"Detector {det_idx+1} - boresight diff [deg] mean/std: "
            f"dRA={np.nanmean(dra_det):.6e}/{np.nanstd(dra_det):.6e}, "
            f"dDec={np.nanmean(ddec_det):.6e}/{np.nanstd(ddec_det):.6e}, "
            f"dRoll={np.nanmean(droll_det):.6e}/{np.nanstd(droll_det):.6e}"
        )

        nplot = min(args.plot_samples, len(ra_bore))
        if args.sanity:
            run_sanity_checks(quat, ra_fp, dec_fp, roll_fp, plot_samples=nplot)

        _, axes = plt.subplots(3, 1, sharex=True, figsize=(10, 8))

        axes[0].plot(ra_bore[:nplot], label="Boresight RA", alpha=0.8)
        for idx in det_indices:
            axes[0].plot(ra_fp[:nplot, idx], label=f"Det {idx+1} RA", alpha=0.8)
        axes[0].set_ylabel("RA [deg]")
        axes[0].legend(loc="best")

        axes[1].plot(dec_bore[:nplot], label="Boresight Dec", alpha=0.8)
        for idx in det_indices:
            axes[1].plot(dec_fp[:nplot, idx], label=f"Det {idx+1} Dec", alpha=0.8)
        axes[1].set_ylabel("Dec [deg]")
        axes[1].legend(loc="best")

        axes[2].plot(roll_bore[:nplot], label="Boresight Roll", alpha=0.8)
        for idx in det_indices:
            axes[2].plot(roll_fp[:nplot, idx], label=f"Det {idx+1} Roll", alpha=0.8)
        axes[2].set_ylabel("Roll [deg]")
        axes[2].set_xlabel("Sample")
        axes[2].legend(loc="best")

        plt.tight_layout()
        plt.show()

        # Compare to reference pickle if requested
        if args.compare_pickle:
            lon_ref, lat_ref = load_pickle_reference(filename, det_type, Path(filename).parent)
            if lon_ref is not None and lat_ref is not None:
                # Subset reference data to match valid quaternion samples
                lon_ref_valid = lon_ref[valid, :]
                lat_ref_valid = lat_ref[valid, :]
                compare_to_reference(ra_fp, dec_fp, lon_ref_valid, lat_ref_valid, det_indices=det_indices, plot_samples=nplot)


# ---------------------------------------------------------------------------
# Embedded unit tests
# Run with: pytest -q get_pointing_from_quaternion.py
# ---------------------------------------------------------------------------


def _test_identity_quat(nsamp):
    q = np.zeros((nsamp, 4), dtype=np.float64)
    q[:, 3] = 1.0
    return q


def test_style_normalization_aliases():
    assert normalize_offset_style(None) is None
    assert normalize_offset_style("none") is None
    assert normalize_offset_style("raw") is None
    assert normalize_offset_style("longcorrect") == "long"
    assert normalize_offset_style("idl") == "attitude"
    assert normalize_offset_style("seasonal-attitude") == "attitude"

    assert normalize_aberration_style(None) == "long"
    assert normalize_aberration_style("longcorrect") == "long"
    assert normalize_aberration_style("ecliptic") == "long"
    assert normalize_aberration_style("idl") == "seasonal"
    assert normalize_aberration_style("aftime") == "seasonal"


def test_style_normalization_rejects_invalid_values():
    import pytest

    with pytest.raises(ValueError):
        normalize_offset_style("bad-style")
    with pytest.raises(ValueError):
        normalize_aberration_style("bad-style")


def test_required_arguments_for_styles():
    import pytest

    quat = _test_identity_quat(4)
    x = np.array([0.0], dtype=np.float64)
    y = np.array([0.0], dtype=np.float64)

    with pytest.raises(ValueError, match="packetid0"):
        focal_plane_solution(
            quat,
            x,
            y,
            offset_style="long",
            detector="SW",
            compute_roll=False,
        )

    with pytest.raises(ValueError, match="attitude_source"):
        focal_plane_solution(
            quat,
            x,
            y,
            offset_style="attitude",
            detector="SW",
            compute_roll=False,
        )

    with pytest.raises(ValueError, match="aftime"):
        focal_plane_solution(
            quat,
            x,
            y,
            aberration=True,
            aberration_style="seasonal",
            compute_roll=False,
        )


def test_boresight_identity_quaternion():
    quat = _test_identity_quat(3)
    ra, dec = boresight_solution(quat)

    assert np.allclose(ra, 0.0, atol=1e-12)
    assert np.allclose(dec, 90.0, atol=1e-12)


def test_noab_matches_internal_kernel():
    quat = _test_identity_quat(5)
    x = np.array([0.5, 2.5], dtype=np.float64)
    y = np.array([-0.25, -1.0], dtype=np.float64)

    ra_a, dec_a, _ = focal_plane_solution(
        quat,
        x,
        y,
        aberration=False,
        compute_roll=False,
    )

    q_fov = quat_from_xy_arcmin(x, y)
    ra_b, dec_b, _ = _focal_plane_solution_no_aberration_qfov(
        quat,
        q_fov,
        compute_roll=False,
    )

    assert np.allclose(ra_a, ra_b, atol=1e-12, equal_nan=True)
    assert np.allclose(dec_a, dec_b, atol=1e-12, equal_nan=True)


def test_seasonal_aberration_shapes_and_finiteness():
    n = 8
    aftime = np.linspace(6.0e6, 7.0e6, n)
    ra = np.linspace(0.0, 359.0, n)
    dec = np.linspace(-45.0, 45.0, n)

    abx, aby = compute_seasonal_aberration_arcsec(aftime, ra, dec)

    assert abx.shape == (n,)
    assert aby.shape == (n,)
    assert np.all(np.isfinite(abx))
    assert np.all(np.isfinite(aby))


def test_matched_recipe_equivalence_zero_residual_synthetic():
    """Matched recipes should be numerically identical on synthetic inputs."""
    nsamp = 7
    ndet = 75  # LW channel count; required by offset helper indexing.

    # Synthetic attitude timeline with small smooth motion.
    quat = quat_from_xy_arcmin(
        np.linspace(-4.0, 4.0, nsamp),
        np.linspace(1.5, -1.5, nsamp),
    )
    quat /= np.linalg.norm(quat, axis=1, keepdims=True)

    # Synthetic detector geometry and time.
    x = np.linspace(-3.0, 3.0, ndet, dtype=np.float64)
    y = np.linspace(-1.5, 1.5, ndet, dtype=np.float64)
    aftime = np.linspace(6.1e6, 6.2e6, nsamp, dtype=np.float64)

    # Long recipe: canonical labels vs aliases.
    ra_l1, dec_l1, roll_l1 = focal_plane_solution(
        quat,
        x,
        y,
        aberration=True,
        compute_roll=True,
        offset_style="long",
        detector="LW",
        packetid0=0x42,
        attitude_source="PR",
        aberration_style="long",
        aftime=aftime,
    )
    ra_l2, dec_l2, roll_l2 = focal_plane_solution(
        quat,
        x,
        y,
        aberration=True,
        compute_roll=True,
        offset_style="longcorrect",
        detector="LW",
        packetid0=0x42,
        attitude_source="PR",
        aberration_style="ecliptic",
        aftime=aftime,
    )

    # Seasonal recipe: canonical labels vs historical aliases.
    ra_s1, dec_s1, roll_s1 = focal_plane_solution(
        quat,
        x,
        y,
        aberration=True,
        compute_roll=True,
        offset_style="attitude",
        detector="LW",
        packetid0=0x42,
        attitude_source="PR",
        aberration_style="seasonal",
        aftime=aftime,
    )
    ra_s2, dec_s2, roll_s2 = focal_plane_solution(
        quat,
        x,
        y,
        aberration=True,
        compute_roll=True,
        offset_style="idl",
        detector="LW",
        packetid0=0x42,
        attitude_source="PR",
        aberration_style="idl",
        aftime=aftime,
    )

    assert np.allclose(ra_l1, ra_l2, atol=1e-12, equal_nan=True)
    assert np.allclose(dec_l1, dec_l2, atol=1e-12, equal_nan=True)
    assert np.allclose(roll_l1, roll_l2, atol=1e-12, equal_nan=True)

    assert np.allclose(ra_s1, ra_s2, atol=1e-12, equal_nan=True)
    assert np.allclose(dec_s1, dec_s2, atol=1e-12, equal_nan=True)
    assert np.allclose(roll_s1, roll_s2, atol=1e-12, equal_nan=True)


if __name__ == "__main__":
    main()
