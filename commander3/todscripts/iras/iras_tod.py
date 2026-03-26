#!/usr/bin/env python3
"""
iras_tod.py  —  IRAS time-ordered data (TOD) access layer.

Provides a Python interface to the raw IRAS archive, suitable both for
full-survey map-making and for feeding calibration-ready TOD to external
pipelines such as Commander3.

Quick start
-----------
    from iras_tod import load_cal_tables, iter_obs

    cal = load_cal_tables()          # load all calibration files once

    for obs in iter_obs(band=4, nside=512, cal=cal):
        # obs.dn          -- raw decompressed ADU per detector, dict[int] -> (NOBS,)
        # obs.pix         -- HEALPix RING pixel per sample,   dict[int] -> (NOBS,)
        # obs.utcs        -- UTC-1981 seconds per sample,     shape (NOBS,)
        # obs.utmbb       -- seconds since last bias boost,   shape (NOBS,)
        # obs.apl/bpl/tpl -- SRHF responsivity arrays,        shape (62,)
        # obs.flash_trim  -- B2 ICF edge mask,                bool (NOBS,)
        # obs.corrupt     -- A2 corrupt-times mask,           bool (NOBS,)
        process(obs)

Commander3 / Fortran interface
------------------------------
Each IrasObs contains all inputs that comm_tod_iras_mod.f90 needs to
reproduce the Python DN → V → A → W/m² → MJy sr⁻¹ calibration chain:

  dn            raw decompressed ADU (re-calibrate in Fortran)
  utcs          UTC-1981 seconds per sample (timing-dependent steps)
  utmbb         seconds since last ICF bias boost
  pix           HEALPix RING pixel per sample (-1 = no valid pointing)
  apl, bpl, tpl SRHF flash-model responsivity arrays, shape (62,)
  prhf_deltut   PRHF time tags, shape (M,), or None if unavailable
  prhf_dtr      PRHF responsivity corrections, shape (M, 31), or None
  psi0_deg      PSI angle at satcal[0] (degrees)
  psirate       deg / SATCAL-tick (scan orientation)
  theta         ecliptic co-latitude of boresight (degrees)
  lambda_obs_rad   solar ecliptic longitude at satcal[0] (radians)
  cdelt1        SATCAL ticks per sample
  scan_speed    normalised speed (1.0 = survey; < 1 = AO slow scan)
  flash_trim    B2 scan-edge ICF zone mask (bool)
  corrupt       A2 corrupt-times window mask (bool)

Detector metadata accessible via BAND_CFG[band]:
  yloc, zloc    focal-plane offsets in arcmin (per detector)
  ircc_times    readout-time offsets in SATCAL ticks
  rate          sub-samples per SATCAL tick
  ftail, decay  tail-flagging parameters (irim_flags)
  dead          frozenset of dead detector numbers
"""

import os
import sys
import glob
import math
import re
import warnings
from dataclasses import dataclass
from typing import Iterator, Optional, List

import numpy as np
import healpy as hp
import fitsio
from astropy.io import fits as afits
from astropy.coordinates import SkyCoord, BarycentricMeanEcliptic
from astropy.time import Time
import astropy.units as u

# ── ipaccal (calibration chain) ───────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ipaccal as ip

# ── Default file paths ────────────────────────────────────────────────────────
ROOT     = os.path.dirname('/mn/stornext/d23/cmbco/globe/orig/iras/kester_rawdb/')
IPAC_DIR = os.path.join(ROOT, 'diskrog10androg11/rog11/IPAC')

DEFAULT_PATHS = dict(
    ctype2        = os.path.join(IPAC_DIR, 'tables/ctype2'),
    vfet          = os.path.join(IPAC_DIR, 'tables/vfetsg'),
    bbtimes       = os.path.join(IPAC_DIR, 'tables/bbtimes'),
    srhf_dir      = os.path.join(IPAC_DIR, 'S'),
    prhf_dir      = os.path.join(IPAC_DIR, 'P'),
    utcs_refs     = os.path.join(IPAC_DIR, 'tables/utcs_satc.refs'),
    corrupt_times = os.path.join(ROOT, 'diskrog10androg11/rog11/outage/corrupt_times'),
    bphf_glob     = os.path.join(ROOT,
                        'diskrog10androg11/rog11/bphf/*/sop.{seg}_/[0-9]*'),
)

# ── Physical / instrument constants ──────────────────────────────────────────
M2R              = np.pi / (180.0 * 60.0)   # arcmin → radians
_SURVEY_SPEED_PT = -0.06416667              # nominal survey PSIRATE, deg/tick

# B2: FLASHTIMES scan-edge trim in SATCAL ticks.
FLASH_PRE  = 32   # ticks to discard from the start of each survey observation
FLASH_POST = 20   # ticks to discard from the end   of each survey observation

# D2: planet guard radii (degrees).
_UTC1981_UNIX = Time('1981-01-01T00:00:00', scale='utc').unix
_PLANET_GUARDS = [
    ('jupiter', 60 / 60.0),   # bright + stray-light halo
    ('saturn',  60 / 60.0),   # rings extend > 30'; stray light
    ('mars',    20 / 60.0),   # can saturate Band IV
    ('uranus',  15 / 60.0),   # bright at 60 µm
    ('neptune', 10 / 60.0),   # moderately bright
    ('moon',   600 / 60.0),   # extreme stray light
]

# ── Per-band focal-plane constants ────────────────────────────────────────────
# Sources: ircc_const.c, Faint Source Survey Expl. Supplement Appendix II.2.
# All position values in arcmin.  Detector numbering matches detnrs[band].
# ftail/decay: tail-flagging parameters from irscan.src irim_flags().
# dead: detector numbers with alive[det]==0 in ircc_const.c.
BAND_CFG = {
    1: dict(
        name='Band I (12 \u00b5m)', survey_tag='b1', rate=16,
        edelay=0.0369, wavelength_um=12, ftail=0.04, decay=0.4,
        dets  =[23,  24,  25,  26,  27,  28,  29,  30,
                47,  48,  49,  50,  51,  52,  53,  54],
        yloc  =[ 9.47,  9.48,  9.49,  9.49,  7.71,  7.71,  7.72,  7.72,
                -5.68, -5.68, -5.67, -5.67, -7.44, -7.44, -7.43, -7.43],
        zloc  =[ 9.81,  1.14, -7.52,-14.50, 13.55,  5.47, -3.19,-11.86,
                14.64,  7.65, -1.02, -9.68, 11.98,  3.32, -5.35,-13.41],
        times =[0.01709,0.01904,0.02100,0.02295,0.02490,0.02686,0.02881,0.03076,
                0.01611,0.01807,0.02002,0.02197,0.02393,0.02588,0.02783,0.02979],
        dead  =frozenset(),
        row1  =frozenset([23,24,25,26,27,28,29,30]),
        row2  =frozenset([47,48,49,50,51,52,53,54]),
    ),
    2: dict(
        name='Band II (25 \u00b5m)', survey_tag='b2', rate=16,
        edelay=0.0333, wavelength_um=25, ftail=0.04, decay=0.5,
        dets  =[16,  17,  18,  19,  20,  21,  22,
                39,  40,  41,  42,  43,  44,  45,  46],
        yloc  =[14.03, 14.03, 14.03, 12.26, 12.26, 12.27, 12.27,
                -1.17, -1.17, -1.17, -1.17, -2.93, -2.93, -2.93, -2.93],
        zloc  =[ 8.71,  0.04, -8.62, 12.96,  4.37, -4.29,-12.88,
                14.05,  6.55, -2.12,-10.78, 10.88,  2.22, -6.45,-13.95],
        times =[0.00244,0.00439,0.00635,0.00830,0.01025,0.01221,0.01416,
                0.00146,0.00342,0.00537,0.00732,0.00928,0.01123,0.01318,0.01514],
        dead  =frozenset([17, 20]),
        row1  =frozenset([16,17,18,19,20,21,22]),
        row2  =frozenset([39,40,41,42,43,44,45,46]),
    ),
    3: dict(
        name='Band III (60 \u00b5m)', survey_tag='b3', rate=8,
        edelay=0.0804, wavelength_um=60, ftail=0.0, decay=0.0,
        dets  =[ 8,  9, 10, 11, 12, 13, 14, 15,
                31, 32, 33, 34, 35, 36, 37, 38],
        yloc  =[19.66, 19.67, 19.67, 19.67, 17.14, 17.15, 17.15, 17.15,
                 4.54,  4.53,  4.53,  4.53,  2.02,  2.01,  2.01,  2.01],
        zloc  =[ 9.80,  1.14, -7.53,-14.46, 13.49,  5.47, -3.20,-11.86,
                14.55,  7.61, -1.06, -9.73, 11.94,  3.27, -5.40,-13.41],
        times =[0.03271,0.03467,0.03662,0.03857,0.09521,0.09717,0.09912,0.10107,
                0.03174,0.03369,0.03564,0.03760,0.09424,0.09619,0.09814,0.10010],
        dead  =frozenset([36]),
        row1  =frozenset([ 8, 9,10,11,12,13,14,15]),
        row2  =frozenset([31,32,33,34,35,36,37,38]),
    ),
    4: dict(
        name='Band IV (100 \u00b5m)', survey_tag='b4', rate=4,
        edelay=0.1478, wavelength_um=100, ftail=0.0, decay=0.0,
        dets  =[ 1,  2,  3,  4,  5,  6,  7,
                55, 56, 57, 58, 59, 60, 61, 62],
        yloc  =[28.00, 28.01, 28.01, 23.97, 23.97, 23.98, 23.98,
               -11.33,-11.33,-11.33,-11.32,-15.36,-15.36,-15.36,-15.36],
        zloc  =[ 8.71,  0.04, -8.62, 12.86,  4.37, -4.29,-12.77,
                13.95,  6.55, -2.12,-10.79, 10.88,  2.21, -6.46,-13.85],
        times =[0.04150,0.10205,0.10400,0.16455,0.16650,0.22705,0.22900,
                0.04053,0.04248,0.10303,0.10498,0.16553,0.16748,0.22803,0.22998],
        dead  =frozenset(),
        row1  =frozenset([ 1, 2, 3, 4, 5, 6, 7]),
        row2  =frozenset([55,56,57,58,59,60,61,62]),
    ),
}

# ── B1950 ecliptic → galactic rotation matrix (computed once at import) ───────
_b1950  = Time(1950.0, format='byear')
_basis  = SkyCoord(lon=[0., 90.,  0.]*u.deg,
                   lat=[0.,  0., 90.]*u.deg,
                   frame=BarycentricMeanEcliptic(equinox=_b1950))
R_ECL2GAL = np.column_stack([b.galactic.cartesian.xyz.value for b in _basis])

# ── Delta-decompression lookup table ─────────────────────────────────────────
def _build_D():
    rows = [
        "1xxx000000000000", "01xxx00000000000", "001xxx0000000000",
        "0001xxx000000000", "00001xxx00000000", "000001xxx0000000",
        "0000001xxx000000", "00000001xxx00000", "000000001xxx0000",
        "0000000001xxx000", "000000000011xxx0", "000000000010xxx0",
        "0000000000011xxx", "0000000000010xxx", "0000000000001xxx",
        "0000000000000xxx",
    ]
    values = []
    for pattern in reversed(rows):
        for xxx in range(8):
            b = pattern.replace("xxx", format(xxx, "03b"))
            values.append(int(b, 2))
    values = sorted(v for v in values if v != 0)
    table  = np.array(values)
    D      = np.zeros(256, dtype=np.int32)
    D[1:128] = -table[::-1]
    D[129:]  =  table
    return D

_D_TABLE = _build_D()


# ── Calibration tables container ──────────────────────────────────────────────
@dataclass
class CalTables:
    """Pre-loaded calibration tables, shared across all TOD observations.

    Obtain via :func:`load_cal_tables`; do not construct directly.

    Attributes
    ----------
    a2dc          ADC step size V/DN·gain (scalar, from ctype2).
    gains         Detector gain matrix, shape (62, 3), from ctype2.
    offsets       Detector offset matrix, shape (62, 8), from ctype2.
    utc_tbl       UTC epoch timestamps for vfet baseline, shape (M,).
    vfets_tbl     Vfet baseline table, shape (M, 62).
    bbt           ICF bias-boost event table, shape (N, 2).
    utcs_refs     {sop_int: (sat_off, utc_ref, rate_corr)} SATCAL→UTC.
    corrupt_table {(sop, obs, band): [(satcal_ref, pre, post, frozenset(dets))]}.
    all_srhfs     {sop_int: {obs_int: (apl[62], bpl[62], tpl[62])}}.
    all_prhfs     {sop_int: {obs_int: (deltut[N], dtr[N,31])}}.
    """
    a2dc:          float
    gains:         np.ndarray
    offsets:       np.ndarray
    utc_tbl:       np.ndarray
    vfets_tbl:     np.ndarray
    bbt:           np.ndarray
    utcs_refs:     dict
    corrupt_table: dict
    all_srhfs:     dict
    all_prhfs:     dict


# ── TOD observation container ────────────────────────────────────────────────
@dataclass(eq=False)
class IrasObs:
    """Complete time-ordered data for one IRAS SOP+OBS (all active detectors).

    A SOP+OBS is the longest continuous stretch of data between sky repointings.
    Raw files sharing the same (SOP, OBS) header pair are concatenated and
    sorted by SATCAL time to form a single, contiguous observation object.

    Shared time fields apply uniformly to all detectors; per-detector fields
    are stored in dicts keyed by 1-based detector number.  Only active
    (non-dead) detectors with at least one valid pointing sample are included
    in the per-detector dicts.

    Quality masks follow the convention ``True`` = bad (exclude from
    calibration or map accumulation).
    Invalid pointing is indicated by ``pix[det][i] == -1``.

    Optional true-pointing fields ``glon`` / ``glat`` may be populated when
    ``iter_obs(..., with_angles=True)`` is used. These store per-sample
    galactic longitude and latitude in degrees and use NaN for unmatched
    samples.

    Commander3 / Fortran interface
    --------------------------------
    Everything needed to calibrate DN → MJy sr⁻¹ is present:
      dn[det]         raw ADU;  NaN = hardware blank or ADC saturation
      utcs            UTC-1981 seconds per sample
      utmbb           seconds since last ICF bias boost
      pix[det]        HEALPix RING pixel per sample (-1 = unmatched)
      apl, bpl, tpl   SRHF flash-model responsivity arrays, shape (62,)
      prhf_deltut/dtr PRHF photon-history arrays, or None
      scan_speed, psirate, theta, psi0_deg, lambda_obs_rad — scan orientation
    """
    # -- identification --
    sop:  int   # IRAS Survey Operations Plan number
    obs:  int   # Observation number within the SOP
    band: int   # Photometric band (1=12µm, 2=25µm, 3=60µm, 4=100µm)

    # -- shared time axis (concatenated & sorted across all raw files) --
    satcal:     np.ndarray   # float64 (NOBS,) SATCAL clock ticks
    utcs:       np.ndarray   # float64 (NOBS,) UTC-1981 seconds
    utmbb:      np.ndarray   # float64 (NOBS,) seconds since last ICF bias boost
    flash_trim: np.ndarray   # bool    (NOBS,) B2 ICF scan-edge mask

    # -- per-detector arrays (keys = 1-based detector numbers) --
    dn:      dict   # det -> float64 (NOBS,): raw ADU; NaN = blanked/saturated
    pix:     dict   # det -> int32   (NOBS,): HEALPix RING pixel; -1 = unmatched
    corrupt: dict   # det -> bool    (NOBS,): A2 corrupt-times mask

    # -- scan orientation (from first raw file; constant within one sopobs) --
    scan_speed:     float   # normalised (1.0 = survey; < 1 = AO slow scan)
    psirate:        float   # deg / SATCAL-tick
    theta:          float   # ecliptic co-latitude of boresight (degrees)
    psi0_deg:       float   # PSI angle at first sample (degrees)
    lambda_obs_rad: float   # solar ecliptic longitude at first sample (radians)
    cdelt1:         float   # SATCAL ticks per sample

    # -- SRHF flash-calibration (shared across all detectors) --
    apl: np.ndarray   # Flash amplitudes,      shape (62,)
    bpl: np.ndarray   # Flash decay constants, shape (62,)
    tpl: np.ndarray   # Flash timestamps,      shape (62,)

    # -- PRHF photon-response history (None if unavailable) --
    prhf_deltut: Optional[np.ndarray]   # shape (M,)
    prhf_dtr:    Optional[np.ndarray]   # shape (M, 31)

    # -- optional provenance / true-pointing metadata --
    source_kind: str = 'unknown'          # survey, ao, spline, unknown, mixed
    glon: Optional[dict] = None           # det -> float64 (NOBS,) gal lon [deg]
    glat: Optional[dict] = None           # det -> float64 (NOBS,) gal lat [deg]


# ── Table-loading helpers ─────────────────────────────────────────────────────

def load_utcs_refs(path):
    """
    Parse utcs_satc.refs → dict: sop_int -> (sat_off, utc_ref, rate_corr).

    File format (space-separated): SOP  satcal_offset  utc_ref  scale_times_1e6
    UTC = utc_ref + (SATCAL - satcal_offset) * scale_times_1e6 / 1e6
    """
    refs = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    sop = int(parts[0])
                    refs[sop] = (float(parts[1]),
                                 float(parts[2]),
                                 float(parts[3]) / 1e6)
                except ValueError:
                    pass
    return refs


def satcal_to_utc(satcal, sop, refs):
    """Convert a SATCAL tick value to UTC-1981 seconds."""
    sat_off, utc_ref, rate_corr = refs[sop]
    return utc_ref + (satcal - sat_off) * rate_corr


def load_corrupt_times(path):
    """
    Parse rog11/outage/corrupt_times → lookup dict.

    Format (space-separated): SOP  ATT  BAND  SATCAL  PRE  POST  det1 [det2 ...]
    The corrupt window covers SATCAL ticks in [SATCAL+PRE, SATCAL+POST] for
    the listed detectors in that (SOP, ATT, BAND) observation.

    Returns
    -------
    dict : {(sop, obs, band): [(satcal_ref, pre, post, frozenset(dets)), ...]}
    Empty dict if the file is missing.
    """
    table = {}
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 7:
                    continue
                sop        = int(parts[0])
                att        = int(parts[1])
                band       = int(parts[2])
                satcal_ref = int(parts[3])
                pre        = int(parts[4])
                post       = int(parts[5])
                dets       = frozenset(int(x) for x in parts[6:])
                table.setdefault((sop, att, band), []).append(
                    (satcal_ref, pre, post, dets))
    except FileNotFoundError:
        pass
    return table


def load_all_srhfs(srhf_dir):
    """
    Pre-load all SRHF files → {sop_int: {obs_int: (apl[62], bpl[62], tpl[62])}}.
    Inaccessible files are silently skipped.
    """
    all_srhfs = {}
    for fpath in sorted(glob.glob(os.path.join(srhf_dir, 'IR.SRHF.D0*'))):
        m = re.search(r'D0(\d{3})(\d{3})$', os.path.basename(fpath))
        if not m:
            continue
        sop, obs = int(m.group(1)), int(m.group(2))
        try:
            _, _, apl, bpl, tpl = ip.read_srhf(srhf_dir, sop, obs)
            all_srhfs.setdefault(sop, {})[obs] = (apl, bpl, tpl)
        except Exception:
            pass
    return all_srhfs


def _closest_srhf(srhf_by_obs, obs_int):
    """Return (apl, bpl, tpl) from the closest-OBS SRHF entry."""
    if not srhf_by_obs:
        return None
    best = min(srhf_by_obs, key=lambda o: abs(o - obs_int))
    return srhf_by_obs[best]


def _lookup_srhf(all_srhfs, sop_int, obs_int):
    """
    Look up SRHF for (sop_int, obs_int).

    Primary: exact SOP, closest OBS.  Falls back to nearest-SOP if the exact
    SOP has no SRHF files (SRHF changes slowly over the mission).

    Returns (entry, fallback_sop) where fallback_sop is None when own SOP matched.
    """
    entry = _closest_srhf(all_srhfs.get(sop_int, {}), obs_int)
    if entry is not None:
        return entry, None
    if not all_srhfs:
        return None, None
    nearest_sop = min(all_srhfs, key=lambda s: abs(s - sop_int))
    return _closest_srhf(all_srhfs[nearest_sop], obs_int), nearest_sop


def load_all_prhfs(prhf_dir):
    """
    Pre-load all PRHF files → {sop_int: {obs_int: (deltut[N], dtr[N,31])}}.
    Inaccessible or empty files are silently skipped.
    """
    all_prhfs = {}
    for fpath in sorted(glob.glob(os.path.join(prhf_dir, 'IR.PRHF.D0*'))):
        m = re.search(r'D0(\d{3})(\d{3})$', os.path.basename(fpath))
        if not m:
            continue
        sop, obs = int(m.group(1)), int(m.group(2))
        try:
            deltut, dtr = ip.read_prhf(prhf_dir, sop, obs)
            if len(deltut) > 0:
                all_prhfs.setdefault(sop, {})[obs] = (deltut, dtr)
        except Exception:
            pass
    return all_prhfs


def _lookup_prhf(all_prhfs, sop_int, obs_int):
    """
    Look up PRHF for (sop_int, obs_int), closest-OBS within same SOP only.

    No cross-SOP fallback: PRHF records the actual photon history of that
    specific scan, so a different SOP's data would be physically incorrect.
    Returns (None, None) when unavailable (correction is simply skipped).
    """
    sop_data = all_prhfs.get(sop_int)
    if not sop_data:
        return None, None
    best = min(sop_data, key=lambda o: abs(o - obs_int))
    return sop_data[best]


def utmbb_at_utc(utcs_start, bbt):
    """
    Return seconds elapsed since the last ICF bias-boost event before utcs_start.
    Returns 1e9 if no prior event is found (detector treated as fully relaxed).

    The bbtimes file stores (utbb1, utbb2) = (start, end) of each bias-boost
    window.  The Fortran cigbbi.shl stores bbt = utbb2 (the END) and returns
    lastbb = utbb2, so the Fortran computes utmbb = snip_utc - utbb2.
    We must use column 1 (the end time) to match.
    """
    bbt_ends = bbt[:, 1]   # utbb2 = end of bias-boost window (Fortran: lastbb)
    idx = np.searchsorted(bbt_ends, utcs_start, side='right') - 1
    if idx < 0:
        return 1.0e9
    return float(utcs_start - bbt[idx, 1])


def load_cal_tables(ipac_dir=None, root=None, ctype2=None, vfet=None,
                    bbtimes=None, srhf_dir=None, prhf_dir=None,
                    utcs_refs=None, corrupt_times=None):
    """
    Load all calibration files and return a :class:`CalTables` instance.

    All arguments are optional; unspecified paths fall back to the module-level
    defaults in ``DEFAULT_PATHS``.  Use the keyword arguments to override
    individual paths, e.g. for testing or alternative data sets.

    Parameters
    ----------
    ipac_dir      : str, optional.  Override the entire IPAC directory.
    root          : str, optional.  Override archive root (affects corrupt_times).
    ctype2 … corrupt_times : str, optional.  Per-file path overrides.
    """
    p = dict(DEFAULT_PATHS)
    if ipac_dir is not None:
        p['ctype2']        = os.path.join(ipac_dir, 'tables/ctype2')
        p['vfet']          = os.path.join(ipac_dir, 'tables/vfetsg')
        p['bbtimes']       = os.path.join(ipac_dir, 'tables/bbtimes')
        p['srhf_dir']      = os.path.join(ipac_dir, 'S')
        p['prhf_dir']      = os.path.join(ipac_dir, 'P')
        p['utcs_refs']     = os.path.join(ipac_dir, 'tables/utcs_satc.refs')
    if root is not None:
        p['corrupt_times'] = os.path.join(
            root, 'diskrog10androg11/rog11/outage/corrupt_times')
    for k, v in [('ctype2', ctype2), ('vfet', vfet), ('bbtimes', bbtimes),
                 ('srhf_dir', srhf_dir), ('prhf_dir', prhf_dir),
                 ('utcs_refs', utcs_refs), ('corrupt_times', corrupt_times)]:
        if v is not None:
            p[k] = v

    a2dc, gains, offsets = ip.load_ctype2(p['ctype2'])
    utc_tbl, vfets_tbl   = ip.load_vfet(p['vfet'])
    bbt                   = ip.load_bbtimes(p['bbtimes'])
    _utcs_refs            = load_utcs_refs(p['utcs_refs'])
    _corrupt              = load_corrupt_times(p['corrupt_times'])
    _all_srhfs            = load_all_srhfs(p['srhf_dir'])
    _all_prhfs            = load_all_prhfs(p['prhf_dir'])

    return CalTables(
        a2dc=a2dc, gains=gains, offsets=offsets,
        utc_tbl=utc_tbl, vfets_tbl=vfets_tbl, bbt=bbt,
        utcs_refs=_utcs_refs, corrupt_table=_corrupt,
        all_srhfs=_all_srhfs, all_prhfs=_all_prhfs,
    )


# ── Decompression and pointing ────────────────────────────────────────────────

def decompress_det(data_int8, det_row, h):
    """
    Delta-decode one detector row with periodic INIT absolute resets.

    Parameters
    ----------
    data_int8 : int8 ndarray, shape (N_detectors, NOBS)
    det_row   : row index (0-based) for this detector in data_int8
    h         : FITS header dict (provides INITSTEP and INIT## keys)

    Returns
    -------
    signal : float64 ndarray, shape (NOBS,).
        NaN at A1 (telemetry outage) and A4 (ADC saturation) positions.
    """
    NOBS     = data_int8.shape[1]
    INITSTEP = int(h['INITSTEP'])
    raw_u8   = (data_int8[det_row].astype(np.int16) + 128).astype(np.uint8)
    deltas   = _D_TABLE[raw_u8]
    signal   = np.empty(NOBS, dtype=np.float64)
    step, i0 = 0, 0
    while i0 < NOBS:
        i1       = min(i0 + INITSTEP, NOBS)
        init_key = f'INIT{step:02}{det_row + 1:02}'
        init_val = float(h.get(init_key, 0)) + 32768.0  # header stores DN−32768
        seg      = np.empty(i1 - i0)
        seg[0]   = init_val
        if i1 - i0 > 1:
            seg[1:] = init_val + np.cumsum(deltas[i0 + 1: i1])
        signal[i0:i1] = seg
        i0 = i1; step += 1
    # A1: byte value 0 is the BLANK_C sentinel (telemetry outage).
    signal[raw_u8 == 0] = np.nan
    # A4: bytes 1 and 255 are D_TABLE extremes (ADC saturation / clipped delta).
    signal[(raw_u8 == 1) | (raw_u8 == 255)] = np.nan
    return signal


def pointing_det(fpath, yloc, zloc, ircc_times, rate, nside):
    """
    Compute (satcal_hz, healpix_pix) for all sub-samples in a BPHF ATT file.

    Parameters
    ----------
    fpath      : path to a BPHF FITS attitude file
    yloc, zloc : detector focal-plane offsets in arcmin
    ircc_times : readout time offset in SATCAL ticks (negative for Band IV)
    rate       : sub-samples per tick (e.g. 4 for Band IV survey)
    nside      : HEALPix NSIDE for pixel assignment

    Returns
    -------
    satcal_hz : float64 array, length n*rate — sub-sample SATCAL timestamps
    pix       : int32 array,   length n*rate — RING HEALPix pixels (galactic)
    """
    with afits.open(fpath) as ff:
        h   = ff[0].header
        raw = ff[0].data

    n      = raw.shape[0]
    crval2 = int(round(h['CRVAL2']))

    psis       = np.empty(n); thetas = np.empty(n)
    psis[0]    = float(h['PSI-0']);   thetas[0] = float(h['THETA-0'])
    psis[1:]   = psis[0]   + np.cumsum(raw[1:, 0].astype(np.float64))
    thetas[1:] = thetas[0] + np.cumsum(raw[1:, 1].astype(np.float64))
    twists     = np.radians(raw[:, 2].astype(np.float64))
    sunlon     = float(h['LNGSUN']) + np.arange(n) * float(h['SUNRATE'])

    P  = np.radians(psis); T = np.radians(thetas); L = np.radians(sunlon)
    dy = yloc * M2R; dz = zloc * M2R
    r  = 1.0 / np.sqrt(1 + dy**2 + dz**2)
    tw = twists
    u_ =  dy*np.cos(tw) + dz*np.sin(tw)
    v_ = -dy*np.sin(tw) + dz*np.cos(tw)

    xx = (np.cos(P)*np.sin(T) - u_*np.sin(P) - v_*np.cos(P)*np.cos(T)) * r
    yy = (np.sin(P)*np.sin(T) + u_*np.cos(P) - v_*np.sin(P)*np.cos(T)) * r
    zz = (np.cos(T)                            + v_*np.sin(T)           ) * r

    sL = np.sin(L); cL = np.cos(L)
    ex =  sL*yy + cL*zz
    ey = -cL*yy + sL*zz
    ez =  xx

    k   = np.arange(rate, dtype=np.float64)
    p   = 2.0*ircc_times - 1.0 + 2.0*k/rate
    p2  = p*p
    pf1 = 0.5*(p2-p); pf2 = 1.0-p2; pf3 = 0.5*(p2+p)

    def ext(v):
        ve = np.empty(n + 2)
        ve[1:-1] = v; ve[0] = 2*v[0]-v[1]; ve[-1] = 2*v[-1]-v[-2]
        return ve

    ex_e = ext(ex); ey_e = ext(ey); ez_e = ext(ez)
    idx  = np.arange(n)
    ix = (pf1[np.newaxis,:]*ex_e[idx,np.newaxis] +
          pf2[np.newaxis,:]*ex_e[idx+1,np.newaxis] +
          pf3[np.newaxis,:]*ex_e[idx+2,np.newaxis]).ravel()
    iy = (pf1[np.newaxis,:]*ey_e[idx,np.newaxis] +
          pf2[np.newaxis,:]*ey_e[idx+1,np.newaxis] +
          pf3[np.newaxis,:]*ey_e[idx+2,np.newaxis]).ravel()
    iz = (pf1[np.newaxis,:]*ez_e[idx,np.newaxis] +
          pf2[np.newaxis,:]*ez_e[idx+1,np.newaxis] +
          pf3[np.newaxis,:]*ez_e[idx+2,np.newaxis]).ravel()

    norm = np.sqrt(ix**2 + iy**2 + iz**2)
    norm = np.where(norm > 0, norm, 1.0)
    ix /= norm; iy /= norm; iz /= norm

    gx = R_ECL2GAL[0,0]*ix + R_ECL2GAL[0,1]*iy + R_ECL2GAL[0,2]*iz
    gy = R_ECL2GAL[1,0]*ix + R_ECL2GAL[1,1]*iy + R_ECL2GAL[1,2]*iz
    gz = R_ECL2GAL[2,0]*ix + R_ECL2GAL[2,1]*iy + R_ECL2GAL[2,2]*iz

    theta_hp = np.arccos(np.clip(gz, -1.0, 1.0))
    phi_hp   = np.arctan2(gy, gx) % (2*np.pi)
    pix      = hp.ang2pix(nside, theta_hp, phi_hp, nest=False)

    satcal_hz = (crval2 + np.arange(n)[:,np.newaxis] +
                 ircc_times + k[np.newaxis,:]/rate).ravel().astype(np.float64)
    return satcal_hz, pix.astype(np.int32)


def pointing_det_angles(fpath, yloc, zloc, ircc_times, rate):
    """
    Compute (satcal_hz, glon_deg, glat_deg) for all sub-samples in a BPHF ATT file.

    Identical geometry to pointing_det() but returns the continuous floating-point
    galactic longitude and latitude (degrees) rather than HEALPix pixel indices.

    Returns
    -------
    satcal_hz : float64 array, length n*rate
    glon_deg  : float64 array, length n*rate — galactic longitude [0, 360)
    glat_deg  : float64 array, length n*rate — galactic latitude  (-90, 90]
    """
    with afits.open(fpath) as ff:
        h   = ff[0].header
        raw = ff[0].data

    n      = raw.shape[0]
    crval2 = int(round(h['CRVAL2']))

    psis       = np.empty(n); thetas = np.empty(n)
    psis[0]    = float(h['PSI-0']);   thetas[0] = float(h['THETA-0'])
    psis[1:]   = psis[0]   + np.cumsum(raw[1:, 0].astype(np.float64))
    thetas[1:] = thetas[0] + np.cumsum(raw[1:, 1].astype(np.float64))
    twists     = np.radians(raw[:, 2].astype(np.float64))
    sunlon     = float(h['LNGSUN']) + np.arange(n) * float(h['SUNRATE'])

    P  = np.radians(psis); T = np.radians(thetas); L = np.radians(sunlon)
    dy = yloc * M2R; dz = zloc * M2R
    r  = 1.0 / np.sqrt(1 + dy**2 + dz**2)
    tw = twists
    u_ =  dy*np.cos(tw) + dz*np.sin(tw)
    v_ = -dy*np.sin(tw) + dz*np.cos(tw)

    xx = (np.cos(P)*np.sin(T) - u_*np.sin(P) - v_*np.cos(P)*np.cos(T)) * r
    yy = (np.sin(P)*np.sin(T) + u_*np.cos(P) - v_*np.sin(P)*np.cos(T)) * r
    zz = (np.cos(T)                            + v_*np.sin(T)           ) * r

    sL = np.sin(L); cL = np.cos(L)
    ex =  sL*yy + cL*zz
    ey = -cL*yy + sL*zz
    ez =  xx

    k   = np.arange(rate, dtype=np.float64)
    p   = 2.0*ircc_times - 1.0 + 2.0*k/rate
    p2  = p*p
    pf1 = 0.5*(p2-p); pf2 = 1.0-p2; pf3 = 0.5*(p2+p)

    def ext(v):
        ve = np.empty(n + 2)
        ve[1:-1] = v; ve[0] = 2*v[0]-v[1]; ve[-1] = 2*v[-1]-v[-2]
        return ve

    ex_e = ext(ex); ey_e = ext(ey); ez_e = ext(ez)
    idx  = np.arange(n)
    ix = (pf1[np.newaxis,:]*ex_e[idx,np.newaxis] +
          pf2[np.newaxis,:]*ex_e[idx+1,np.newaxis] +
          pf3[np.newaxis,:]*ex_e[idx+2,np.newaxis]).ravel()
    iy = (pf1[np.newaxis,:]*ey_e[idx,np.newaxis] +
          pf2[np.newaxis,:]*ey_e[idx+1,np.newaxis] +
          pf3[np.newaxis,:]*ey_e[idx+2,np.newaxis]).ravel()
    iz = (pf1[np.newaxis,:]*ez_e[idx,np.newaxis] +
          pf2[np.newaxis,:]*ez_e[idx+1,np.newaxis] +
          pf3[np.newaxis,:]*ez_e[idx+2,np.newaxis]).ravel()

    norm = np.sqrt(ix**2 + iy**2 + iz**2)
    norm = np.where(norm > 0, norm, 1.0)
    ix /= norm; iy /= norm; iz /= norm

    gx = R_ECL2GAL[0,0]*ix + R_ECL2GAL[0,1]*iy + R_ECL2GAL[0,2]*iz
    gy = R_ECL2GAL[1,0]*ix + R_ECL2GAL[1,1]*iy + R_ECL2GAL[1,2]*iz
    gz = R_ECL2GAL[2,0]*ix + R_ECL2GAL[2,1]*iy + R_ECL2GAL[2,2]*iz

    theta_hp = np.arccos(np.clip(gz, -1.0, 1.0))
    phi_hp   = np.arctan2(gy, gx) % (2*np.pi)

    satcal_hz = (crval2 + np.arange(n)[:,np.newaxis] +
                 ircc_times + k[np.newaxis,:]/rate).ravel().astype(np.float64)
    return satcal_hz, np.degrees(phi_hp), 90.0 - np.degrees(theta_hp)


def planet_pixel_set(utc_1981, nside):
    """
    Return a sorted int64 array of HEALPix RING pixels (galactic coordinates)
    within the guard radius of any tracked solar-system body at the given time.

    Parameters
    ----------
    utc_1981 : float — seconds since 1981-01-01 00:00 UTC
    nside    : int   — HEALPix NSIDE

    Returns
    -------
    pix : int64 ndarray — unique pixel indices (may be empty)
    """
    from astropy.coordinates import get_body
    t = Time(_UTC1981_UNIX + float(utc_1981), format='unix', scale='utc')
    pieces = []
    for body, r_deg in _PLANET_GUARDS:
        try:
            pos = get_body(body, t)
            gal = pos.galactic
            vec = hp.ang2vec(np.pi / 2.0 - float(gal.b.rad),
                             float(gal.l.rad))
            pix = hp.query_disc(nside, vec, np.radians(r_deg), inclusive=True)
            pieces.append(pix)
        except Exception:
            pass
    return np.unique(np.concatenate(pieces)) if pieces else np.array([], dtype=np.int64)


# ── Raw-file glob helpers ─────────────────────────────────────────────────────

def _normalize_source_kinds(include_ao, source_kinds):
    if source_kinds is None:
        return ('survey', 'ao', 'spline', 'unknown') if include_ao else ('survey',)

    norm = []
    for kind in source_kinds:
        k = str(kind).strip().lower()
        if k in {'survey', 'ao', 'spline', 'unknown'} and k not in norm:
            norm.append(k)
    return tuple(norm)


def _raw_kind_dirname(kind):
    return 'survey' if kind == 'survey' else kind.upper()


def _classify_raw_path(fpath, band):
    lower = fpath.lower()
    tag = BAND_CFG[band]['survey_tag']
    if f'/survey.{tag}/' in lower:
        return 'survey'
    if f'/ao.{tag}/' in lower:
        return 'ao'
    if f'/spline.{tag}/' in lower:
        return 'spline'
    if f'/unknown.{tag}/' in lower:
        return 'unknown'
    return 'unknown'


def _raw_globs_for_band(band, root, include_ao, source_kinds=None):
    """Return list of glob pattern strings (with {seg} placeholder) for raw TOD."""
    tag = BAND_CFG[band]['survey_tag']
    kinds = _normalize_source_kinds(include_ao, source_kinds)
    patterns = []

    if band == 4:
        # disk1to15A is a byte-identical duplicate of disk1to15B for B4 (SOPs 3-25).
        # Use only B-side trees so we don't double-count the first half of the sky.
        for kind in kinds:
            dirname = _raw_kind_dirname(kind)
            patterns.extend([
                os.path.join(root, f'disk16to27?/*/{dirname}.{tag}/sop.{{seg}}_/[0-9]*'),
                os.path.join(root, f'disk1to15B/*/{dirname}.{tag}/sop.{{seg}}_/[0-9]*'),
                os.path.join(root, f'diskrog10androg11/*/*/{dirname}.{tag}/sop.{{seg}}_/[0-9]*'),
            ])
    else:
        for kind in kinds:
            dirname = _raw_kind_dirname(kind)
            patterns.append(os.path.join(root, f'disk*/*/{dirname}.{tag}/sop.{{seg}}_/[0-9]*'))

    return patterns


def _raw_files_for_seg(band, seg_str, root, include_ao, source_kinds=None):
    """Expand raw-file globs for one segment string, returning a sorted unique list."""
    import itertools
    globs = _raw_globs_for_band(band, root, include_ao, source_kinds)
    return sorted(set(itertools.chain.from_iterable(
        glob.glob(g.format(seg=seg_str)) for g in globs)))


# ── Calibration table loading ────────────────────────────────────────────────

def load_cal_tables(ipac_dir=None, root=None, ctype2=None, vfet=None,
                    bbtimes=None, srhf_dir=None, prhf_dir=None,
                    utcs_refs=None, corrupt_times=None):
    """
    Load all calibration files and return a :class:`CalTables` instance.

    All arguments are optional; unspecified paths fall back to the module-level
    defaults in ``DEFAULT_PATHS``.  Use the keyword arguments to override
    individual paths, e.g. for testing or alternative data sets.

    Parameters
    ----------
    ipac_dir      : str, optional.  Override the entire IPAC directory.
    root          : str, optional.  Override archive root (affects corrupt_times).
    ctype2 … corrupt_times : str, optional.  Per-file path overrides.
    """
    p = dict(DEFAULT_PATHS)
    if ipac_dir is not None:
        p['ctype2']        = os.path.join(ipac_dir, 'tables/ctype2')
        p['vfet']          = os.path.join(ipac_dir, 'tables/vfetsg')
        p['bbtimes']       = os.path.join(ipac_dir, 'tables/bbtimes')
        p['srhf_dir']      = os.path.join(ipac_dir, 'S')
        p['prhf_dir']      = os.path.join(ipac_dir, 'P')
        p['utcs_refs']     = os.path.join(ipac_dir, 'tables/utcs_satc.refs')
    if root is not None:
        p['corrupt_times'] = os.path.join(
            root, 'diskrog10androg11/rog11/outage/corrupt_times')
    for k, v in [('ctype2', ctype2), ('vfet', vfet), ('bbtimes', bbtimes),
                 ('srhf_dir', srhf_dir), ('prhf_dir', prhf_dir),
                 ('utcs_refs', utcs_refs), ('corrupt_times', corrupt_times)]:
        if v is not None:
            p[k] = v

    a2dc, gains, offsets = ip.load_ctype2(p['ctype2'])
    utc_tbl, vfets_tbl   = ip.load_vfet(p['vfet'])
    bbt                   = ip.load_bbtimes(p['bbtimes'])
    _utcs_refs            = load_utcs_refs(p['utcs_refs'])
    _corrupt              = load_corrupt_times(p['corrupt_times'])
    _all_srhfs            = load_all_srhfs(p['srhf_dir'])
    _all_prhfs            = load_all_prhfs(p['prhf_dir'])

    return CalTables(
        a2dc=a2dc, gains=gains, offsets=offsets,
        utc_tbl=utc_tbl, vfets_tbl=vfets_tbl, bbt=bbt,
        utcs_refs=_utcs_refs, corrupt_table=_corrupt,
        all_srhfs=_all_srhfs, all_prhfs=_all_prhfs,
    )


# ── Main TOD iterator ─────────────────────────────────────────────────────────

def iter_obs(band=4, *, nside=512, cal, root=None, bphf_glob=None,
             include_ao=False, source_kinds=None, detectors=None, with_angles=False,
             max_seg=None, segments=None,
             flash_pre=FLASH_PRE, flash_post=FLASH_POST):
    """
    Iterate over complete SOP+OBS observations for one IRAS photometric band.

    Yields one :class:`IrasObs` per unique (SOP, OBS) pair found within each
    disk segment.  Raw files sharing the same (SOP, OBS) are concatenated and
    sorted by SATCAL time.  If a sopobs is split across disk-segment boundaries
    (rare), its two halves appear as separate IrasObs objects.

    BPHF pointing is computed once per disk segment and covers all sopobs
    within that segment.

    Observations are silently skipped when:
      - the SOP is absent from ``cal.utcs_refs`` (no UTC reference)
      - no SRHF data are available anywhere (all_srhfs is empty)
      - no detector has any valid pointing match to the BPHF table

    Parameters
    ----------
    band       : int (1-4) — photometric band.
    nside      : int — HEALPix NSIDE for pixel assignment (default 512).
    cal        : CalTables — pre-loaded calibration tables.
    root       : str — archive root directory (default: module ROOT).
    bphf_glob  : str — glob pattern for BPHF files with {seg} placeholder.
    include_ao : bool — include AO/SPLINE/UNKNOWN observation directories.
                 Default False keeps survey-only input selection.
    source_kinds : iterable[str] or None — restrict raw-file classes to any of
                 {'survey', 'ao', 'spline', 'unknown'}. Default follows
                 include_ao.
    detectors : iterable[int] or None — restrict processing to this subset of
                 detector numbers. Invalid or dead detectors are ignored.
    with_angles : bool — also populate per-sample galactic longitude / latitude
                 arrays in obs.glon / obs.glat.
    max_seg    : str or None — zero-padded segment ceiling, e.g. '03'.
    segments   : iterable of str or None — restrict to these segment strings.
    flash_pre  : float — SATCAL ticks to trim at each raw-file start (B2).
    flash_post : float — SATCAL ticks to trim at each raw-file end   (B2).

    Yields
    ------
    IrasObs
    """
    if root is None:
        root = ROOT
    if bphf_glob is None:
        bphf_glob = DEFAULT_PATHS['bphf_glob']

    bcfg      = BAND_CFG[band]
    rate      = bcfg['rate']
    edelay    = bcfg['edelay']
    all_dets  = bcfg['dets']
    dead      = bcfg['dead']
    live_dets = [d for d in all_dets if d not in dead]
    if detectors is not None:
        det_filter = set(detectors)
        live_dets = [d for d in live_dets if d in det_filter]
    if not live_dets:
        return
    det_to_row = {d: i for i, d in enumerate(all_dets)}
    yloc = dict(zip(all_dets, bcfg['yloc']))
    zloc = dict(zip(all_dets, bcfg['zloc']))
    ircc = {d: t - edelay for d, t in zip(all_dets, bcfg['times'])}
    max_gap   = 1.0 / rate + 0.01   # pointing match tolerance in SATCAL ticks

    # --- Discover available path segments ---
    bphf_segs = sorted(set(
        m.group(1)
        for f in glob.glob(bphf_glob.format(seg='*'))
        for m in [re.search(r'sop\.(\d+)_', f)] if m))
    raw_segs = sorted(set(
        m.group(1)
        for g in _raw_globs_for_band(band, root, include_ao, source_kinds)
        for f in glob.glob(g.format(seg='*'))
        for m in [re.search(r'sop\.(\d+)_', f)] if m))

    all_segs = sorted(set(bphf_segs) & set(raw_segs))
    if max_seg is not None:
        all_segs = [s for s in all_segs if s <= max_seg]
    if segments is not None:
        seg_set  = set(segments)
        all_segs = [s for s in all_segs if s in seg_set]

    for seg_str in all_segs:
        bphf_files = sorted(glob.glob(bphf_glob.format(seg=seg_str)))
        raw_files  = _raw_files_for_seg(band, seg_str, root, include_ao, source_kinds)
        if not bphf_files or not raw_files:
            continue

        # -- Per-segment BPHF sky context (for zodycal / orientation interpolation) --
        _bphf_sky = []   # [(crval2, psi0_deg, lngsun_deg, sunrate_deg/tick)]
        for fp in bphf_files:
            try:
                with afits.open(fp) as ff:
                    hh = ff[0].header
                _bphf_sky.append((int(round(float(hh['CRVAL2']))),
                                  float(hh['PSI-0']),
                                  float(hh['LNGSUN']),
                                  float(hh['SUNRATE'])))
            except Exception:
                pass
        _bphf_sky.sort(key=lambda x: x[0])
        _bphf_crval2s = np.array([d[0] for d in _bphf_sky], dtype=np.float64)

        # -- Per-segment pointing table (one entry per active detector) --
        # det -> dict(sc, pix[, glon, glat])
        pointing_dets = {}
        for det in live_dets:
            sc_list, pix_list = [], []
            glon_list, glat_list = [], []
            for fp in bphf_files:
                try:
                    if with_angles:
                        sc, glon, glat = pointing_det_angles(
                            fp, yloc[det], zloc[det], ircc[det], rate)
                        px = hp.ang2pix(
                            nside,
                            np.radians(90.0 - glat),
                            np.radians(glon),
                            nest=False,
                        ).astype(np.int32)
                        glon_list.append(glon)
                        glat_list.append(glat)
                    else:
                        sc, px = pointing_det(fp, yloc[det], zloc[det],
                                              ircc[det], rate, nside)
                    sc_list.append(sc); pix_list.append(px)
                except Exception:
                    pass
            if sc_list:
                bsc = np.concatenate(sc_list)
                bpx = np.concatenate(pix_list)
                order = np.argsort(bsc)
                entry = dict(sc=bsc[order], pix=bpx[order])
                if with_angles:
                    entry['glon'] = np.concatenate(glon_list)[order]
                    entry['glat'] = np.concatenate(glat_list)[order]
                pointing_dets[det] = entry

        if not pointing_dets:
            continue

        # -- First pass: group raw files by (SOP, OBS) --
        # Peek at headers only; actual data reading happens per-observation.
        obs_files = {}   # (sop_int, obs_int) -> [(satcal0, fpath)] sorted later
        for fpath in raw_files:
            try:
                with fitsio.FITS(fpath) as ff:
                    h = ff[0].read_header()
                key = (int(h['SOP']), int(h['OBS']))
                obs_files.setdefault(key, []).append(
                    (float(h['SATCAL']), fpath))
            except Exception:
                pass

        # -- Second pass: process each (SOP, OBS) group ---
        for (sop_int, obs_int), tagged in sorted(obs_files.items()):
            tagged.sort(key=lambda x: x[0])   # sort by SATCAL
            fpaths = [fp for _, fp in tagged]

            if sop_int not in cal.utcs_refs:
                continue

            srhf, _ = _lookup_srhf(cal.all_srhfs, sop_int, obs_int)
            if srhf is None:
                continue
            apl, bpl, tpl = srhf

            prhf_deltut, prhf_dtr = _lookup_prhf(cal.all_prhfs, sop_int, obs_int)
            ct_windows = cal.corrupt_table.get((sop_int, obs_int, band), [])

            # Parts to concatenate; populated per raw file
            sc_parts    = []
            utcs_parts  = []
            utmbb_parts = []
            ft_parts    = []
            dn_parts    = {det: [] for det in live_dets}
            pix_parts   = {det: [] for det in live_dets}
            if with_angles:
                glon_parts = {det: [] for det in live_dets}
                glat_parts = {det: [] for det in live_dets}
            cor_parts   = {det: [] for det in live_dets}
            valid_dets  = set()   # dets with ≥1 matched sample
            first_meta  = None
            obs_kinds   = set()

            for fpath in fpaths:
                try:
                    raw_kind = _classify_raw_path(fpath, band)
                    obs_kinds.add(raw_kind)
                    with fitsio.FITS(fpath) as ff:
                        raw_h = ff[0].read_header()
                        data  = ff[0].read()   # (N_det_rows, NOBS) int8

                    satcal0  = float(raw_h['SATCAL'])
                    cdelt1   = float(raw_h['CDELT1'])
                    NOBS     = data.shape[1]

                    _psirate = float(raw_h.get('PSIRATE', _SURVEY_SPEED_PT))
                    _theta   = float(raw_h.get('THETA',   90.0))
                    _cos_lat = math.cos(math.radians(90.0 - _theta))
                    scan_speed = max(abs(_psirate * _cos_lat / _SURVEY_SPEED_PT), 0.05)

                    psi0_deg = 0.0; lambda_rad = 0.0
                    if _bphf_sky:
                        bi = int(np.searchsorted(_bphf_crval2s, satcal0))
                        bi = max(0, min(bi, len(_bphf_sky) - 1))
                        if (bi > 0 and
                                abs(_bphf_crval2s[bi - 1] - satcal0) <
                                abs(_bphf_crval2s[bi]     - satcal0)):
                            bi -= 1
                        bc2, bp0, bln, bsr = _bphf_sky[bi]
                        dt         = satcal0 - float(bc2)
                        psi0_deg   = bp0 + _psirate * dt
                        lambda_rad = np.radians(bln + bsr * dt)

                    if first_meta is None:
                        first_meta = dict(
                            scan_speed=scan_speed, psirate=_psirate,
                            theta=_theta, psi0_deg=psi0_deg,
                            lambda_obs_rad=lambda_rad, cdelt1=cdelt1)

                    tod_sc      = satcal0 + np.arange(NOBS, dtype=np.float64) * cdelt1
                    utcs_start  = satcal_to_utc(satcal0, sop_int, cal.utcs_refs)
                    utcs_inc    = cdelt1 * cal.utcs_refs[sop_int][2]
                    utcs_arr    = utcs_start + np.arange(NOBS, dtype=np.float64) * utcs_inc
                    utmbb_start = utmbb_at_utc(utcs_start, cal.bbt)
                    utmbb_arr   = utmbb_start + np.arange(NOBS, dtype=np.float64) * utcs_inc

                    ft = np.zeros(NOBS, dtype=bool)
                    # FLASHTIMES trimming applies to survey-mode scan edges.
                    # Pointed AO/SPLINE/UNKNOWN observations are untrimmed.
                    if raw_kind == 'survey':
                        n_pre  = int(round(flash_pre  / cdelt1))
                        n_post = int(round(flash_post / cdelt1))
                        if 0 < n_pre  < NOBS: ft[:n_pre]   = True
                        if 0 < n_post < NOBS: ft[-n_post:] = True

                    sc_parts.append(tod_sc)
                    utcs_parts.append(utcs_arr)
                    utmbb_parts.append(utmbb_arr)
                    ft_parts.append(ft)

                    for det in live_dets:
                        if det not in pointing_dets:
                            pix_parts[det].append(
                                np.full(NOBS, -1, dtype=np.int32))
                            if with_angles:
                                glon_parts[det].append(np.full(NOBS, np.nan))
                                glat_parts[det].append(np.full(NOBS, np.nan))
                            dn_parts[det].append(
                                np.full(NOBS, np.nan))
                            cor_parts[det].append(
                                np.zeros(NOBS, dtype=bool))
                            continue

                        pdet = pointing_dets[det]
                        bsc = pdet['sc']
                        bpx = pdet['pix']
                        idx    = np.searchsorted(bsc, tod_sc)
                        idx    = np.clip(idx, 0, len(bsc) - 1)
                        idx_m1 = np.maximum(idx - 1, 0)
                        closer = (np.abs(bsc[idx_m1] - tod_sc) <
                                  np.abs(bsc[idx]    - tod_sc))
                        idx    = np.where(closer, idx_m1, idx)
                        valid  = np.abs(bsc[idx] - tod_sc) <= max_gap

                        dn  = decompress_det(data, det_to_row[det], raw_h)
                        px  = np.where(valid, bpx[idx],
                                       np.int32(-1)).astype(np.int32)
                        if with_angles:
                            glon = np.where(valid, pdet['glon'][idx], np.nan)
                            glat = np.where(valid, pdet['glat'][idx], np.nan)
                        corrupt = np.zeros(NOBS, dtype=bool)
                        for sat_ref, pre, post, dets in ct_windows:
                            if det in dets:
                                corrupt |= ((tod_sc >= sat_ref + pre) &
                                            (tod_sc <= sat_ref + post))

                        dn_parts[det].append(dn)
                        pix_parts[det].append(px)
                        if with_angles:
                            glon_parts[det].append(glon)
                            glat_parts[det].append(glat)
                        cor_parts[det].append(corrupt)
                        if valid.any():
                            valid_dets.add(det)

                except Exception:
                    pass

            if not sc_parts or not valid_dets:
                continue

            if first_meta is None:
                first_meta = dict(
                    scan_speed=1.0, psirate=_SURVEY_SPEED_PT,
                    theta=90.0, psi0_deg=0.0,
                    lambda_obs_rad=0.0, cdelt1=1.0 / rate)

            if not obs_kinds:
                source_kind = 'unknown'
            elif len(obs_kinds) == 1:
                source_kind = next(iter(obs_kinds))
            else:
                source_kind = 'mixed'

            yield IrasObs(
                sop=sop_int, obs=obs_int, band=band,
                satcal=np.concatenate(sc_parts),
                utcs=np.concatenate(utcs_parts),
                utmbb=np.concatenate(utmbb_parts),
                flash_trim=np.concatenate(ft_parts),
                dn={d: np.concatenate(dn_parts[d])  for d in valid_dets},
                pix={d: np.concatenate(pix_parts[d]) for d in valid_dets},
                corrupt={d: np.concatenate(cor_parts[d]) for d in valid_dets},
                apl=apl, bpl=bpl, tpl=tpl,
                prhf_deltut=prhf_deltut, prhf_dtr=prhf_dtr,
                    source_kind=source_kind,
                    glon=({d: np.concatenate(glon_parts[d]) for d in valid_dets}
                        if with_angles else None),
                    glat=({d: np.concatenate(glat_parts[d]) for d in valid_dets}
                        if with_angles else None),
                **first_meta,
            )
