#!/usr/bin/env python3
"""
fullsurvey_b4_alldet_cal_map.py

Full-survey calibrated HEALPix map combining all 15 Band IV (100 µm) detectors.

Strategy: path-segment-by-segment (groups of raw/BPHF files sharing the same
sop.XX_ directory name, which is a disk-segment index, not the IRAS SOP number).

For each segment:
  1. Compute pointing (HEALPix pixel, SATCAL timestamp) for all 15 band-IV
     detectors using the BPHF attitude files.
  2. Load each raw Band IV TOD file → extract actual SOP/OBS from its FITS
     header → decompress all 15 detector rows.
  3. Convert SATCAL → UTC via the utcs_satc.refs table; derive utmbb (time
     since last bias boost) from the bbtimes table.
  4. Calibrate each detector's timeline: DN → MJy/sr (full IPACCAL pipeline).
  5. Nearest-neighbour match calibrated samples to HEALPix pixels via SATCAL
     timestamp; accumulate sum_map / hits_map (arithmetic mean = equal weights).

Output
------
  python/fullsurvey_b4_alldet_cal_map_nside512.fits  — calibrated MJy/sr map
  python/fullsurvey_b4_alldet_cal_hits_nside512.fits — hits map
  python/fullsurvey_b4_alldet_cal_map_nside512.png   — Mollweide plot
"""

import os
import sys
import glob
import re
import warnings
import multiprocessing
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import fitsio
from scipy.ndimage import convolve1d
from astropy.io import fits as afits
from astropy.coordinates import SkyCoord, BarycentricMeanEcliptic
from astropy.time import Time
import astropy.units as u
from tqdm import tqdm

# ── Path to ipaccal ──────────────────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ipaccal as ip

# ── Configuration ─────────────────────────────────────────────────────────────
NSIDE = 512        # HEALPix resolution (~6.9 arcmin pixels)
M2R   = np.pi / (180.0 * 60.0)   # arcmin → radians

ROOT     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IPAC_DIR = os.path.join(ROOT, 'diskrog10androg11/rog11/IPAC')
CTYPE2   = os.path.join(IPAC_DIR, 'tables/ctype2')
VFET     = os.path.join(IPAC_DIR, 'tables/vfetsg')
BBT_FILE = os.path.join(IPAC_DIR, 'tables/bbtimes')
SRHF_DIR = os.path.join(IPAC_DIR, 'S')
UTCS_REF      = os.path.join(IPAC_DIR, 'tables/utcs_satc.refs')
CORRUPT_TIMES = os.path.join(ROOT, 'diskrog10androg11/rog11/outage/corrupt_times')

BPHF_GLOB = os.path.join(ROOT, 'diskrog10androg11/rog11/bphf/*/sop.{seg}_/[0-9]*')
OUT_DIR   = os.path.dirname(os.path.abspath(__file__))

# Optional segment cap for test runs.  Set MAX_SEG='03' (zero-padded string)
# to process only sop.00_ … sop.03_.  None = no limit (full survey).
MAX_SEG   = os.environ.get('MAX_SEG', None)   # e.g. MAX_SEG=03 python3 ...

# Parallelism: number of worker processes.  NPROC=1 → serial (default when
# NPROC not set).  NPROC=0 → use all available CPUs.
NPROC     = int(os.environ.get('NPROC', '1'))

# B2: FLASHTIMES scan-edge trim (SATCAL ticks; same defaults as PLATE v3.2).
# ICF flashes are scheduled 32 ticks after scan start and 20 ticks before end.
FLASH_PRE  = 32   # ticks to discard from the start of every survey snip
FLASH_POST = 20   # ticks to discard from the end   of every survey snip

# C1: irim_flags glitch/tail filter constants (irscan.src, irim_flags())
# fiGli = {-3, 6, -3}  MG=3 (3-tap Mexican-hat, 2nd derivative)
# fiPnt = {-1,-1,-1, 2,2,2, -1,-1,-1}  MP=9 (point-source filter)
# SUMF  = 12 = sum(|fiGli|) = sum(|fiPnt|);  thres = 5 for all bands
_FGLI        = np.array([-3.0,  6.0, -3.0])
_FPNT        = np.array([-1.0, -1.0, -1.0,  2.0, 2.0, 2.0,  -1.0, -1.0, -1.0])
_IRIM_SUMF   = 12.0
_IRIM_THRES  =  5.0

# D2: planet / bright moving-object guard radii (arcmin → degrees).
# Guard radii follow the IRAS Explanatory Supplement recommendations.
_UTC1981_UNIX = Time('1981-01-01T00:00:00', scale='utc').unix
_PLANET_GUARDS = [
    ('jupiter', 60 / 60.0),   # degrees: bright + stray-light halo
    ('saturn',  60 / 60.0),   # rings extend >30'; stray light
    ('mars',    20 / 60.0),   # can saturate Band IV
    ('uranus',  15 / 60.0),   # bright at 60 µm
    ('neptune', 10 / 60.0),   # moderately bright
    ('moon',   600 / 60.0),   # extreme stray light
]

def planet_pixel_set(utc_1981, nside):
    """
    Return a (sorted) int64 numpy array of HEALPix RING pixels in galactic
    coordinates that fall within the guard radius of any tracked solar-system
    body at the time given by utc_1981 (seconds since 1981-01-01 00:00 UTC).
    Planet positions are geocentric (parallax negligible for guard radii used).
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


# ── Per-band focal-plane constants (ircc_const.c + Faint Source Survey supp.) ─
# All position values in arcmin, from Appendix II.2 of the FSS Expl. Supplement.
# Detector number ordering matches detnrs[band] in ircc_const.c.
# ftail / decay: tail-flagging parameters from irscan.src irim_flags().
# dead: detector numbers with alive[det]=0 in ircc_const.c (skip in loop).
_BAND_CFG = {
    1: dict(
        name='Band I (12 \u00b5m)', survey_tag='b1', rate=16, edelay=0.0369, wavelength_um=12,
        ftail=0.04, decay=0.4,
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
        name='Band II (25 \u00b5m)', survey_tag='b2', rate=16, edelay=0.0333, wavelength_um=25,
        ftail=0.04, decay=0.5,
        dets  =[16,  17,  18,  19,  20,  21,  22,
                39,  40,  41,  42,  43,  44,  45,  46],
        yloc  =[14.03, 14.03, 14.03, 12.26, 12.26, 12.27, 12.27,
                -1.17, -1.17, -1.17, -1.17, -2.93, -2.93, -2.93, -2.93],
        zloc  =[ 8.71,  0.04, -8.62, 12.96,  4.37, -4.29,-12.88,
                14.05,  6.55, -2.12,-10.78, 10.88,  2.22, -6.45,-13.95],
        times =[0.00244,0.00439,0.00635,0.00830,0.01025,0.01221,0.01416,
                0.00146,0.00342,0.00537,0.00732,0.00928,0.01123,0.01318,0.01514],
        dead  =frozenset([17, 20]),          # alive[17]=0, alive[20]=0
        row1  =frozenset([16,17,18,19,20,21,22]),
        row2  =frozenset([39,40,41,42,43,44,45,46]),
    ),
    3: dict(
        name='Band III (60 \u00b5m)', survey_tag='b3', rate=8, edelay=0.0804, wavelength_um=60,
        ftail=0.0, decay=0.0,
        dets  =[ 8,  9, 10, 11, 12, 13, 14, 15,
                31, 32, 33, 34, 35, 36, 37, 38],
        yloc  =[19.66, 19.67, 19.67, 19.67, 17.14, 17.15, 17.15, 17.15,
                 4.54,  4.53,  4.53,  4.53,  2.02,  2.01,  2.01,  2.01],
        zloc  =[ 9.80,  1.14, -7.53,-14.46, 13.49,  5.47, -3.20,-11.86,
                14.55,  7.61, -1.06, -9.73, 11.94,  3.27, -5.40,-13.41],
        times =[0.03271,0.03467,0.03662,0.03857,0.09521,0.09717,0.09912,0.10107,
                0.03174,0.03369,0.03564,0.03760,0.09424,0.09619,0.09814,0.10010],
        dead  =frozenset([36]),              # alive[36]=0
        row1  =frozenset([ 8, 9,10,11,12,13,14,15]),
        row2  =frozenset([31,32,33,34,35,36,37,38]),
    ),
    4: dict(
        name='Band IV (100 \u00b5m)', survey_tag='b4', rate=4, edelay=0.1478, wavelength_um=100,
        ftail=0.0, decay=0.0,
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

# Select band via BAND env var (default 4).
BAND  = int(os.environ.get('BAND', '4'))

# Diagnostic mode: MASK=0 disables all QC flags (calibrated but unmasked).
# Useful for diagnosing survey coverage gaps.
MASK  = int(os.environ.get('MASK', '1'))
_bcfg = _BAND_CFG[BAND]
RATE  = _bcfg['rate']           # samples per SATCAL tick

# Detector numbers for this band, excluding known dead detectors.
# _DET_TO_ROW maps det number → row index in the raw FITS (NDET, NOBS) array,
# using the full unfiltered detnrs ordering (needed even for dead dets).
_ALL_DETS  = _bcfg['dets']
_DEAD_DETS = _bcfg['dead']
BAND_DETS  = [d for d in _ALL_DETS if d not in _DEAD_DETS]
_DET_TO_ROW = {d: i for i, d in enumerate(_ALL_DETS)}

BAND_ROW1 = _bcfg['row1'] - _DEAD_DETS
BAND_ROW2 = _bcfg['row2'] - _DEAD_DETS

_ircc_list = [t - _bcfg['edelay'] for t in _bcfg['times']]
BAND_YLOC  = dict(zip(_ALL_DETS, _bcfg['yloc']))
BAND_ZLOC  = dict(zip(_ALL_DETS, _bcfg['zloc']))
BAND_IRCC  = dict(zip(_ALL_DETS, _ircc_list))

RAW_GLOB = os.path.join(ROOT,
    f'disk*/*/survey.{_bcfg["survey_tag"]}/sop.{{seg}}_/[0-9]*')

def _band_row(det):
    """Return focal-plane row index (1 or 2) for a detector number."""
    return 1 if det in BAND_ROW1 else 2

# ── B1950 ecliptic → galactic rotation matrix ─────────────────────────────────
_b1950 = Time(1950.0, format='byear')
_basis  = SkyCoord(lon=[0., 90.,  0.]*u.deg,
                   lat=[0.,  0., 90.]*u.deg,
                   frame=BarycentricMeanEcliptic(equinox=_b1950))
R_ECL2GAL = np.column_stack([b.galactic.cartesian.xyz.value for b in _basis])

# ── Delta-decompression table ─────────────────────────────────────────────────
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

D_TABLE = _build_D()


def decompress_det(data_int8, det_row, h):
    """Delta-decode one detector row with periodic INIT absolute resets."""
    NOBS     = data_int8.shape[1]
    INITSTEP = int(h['INITSTEP'])
    raw_u8   = (data_int8[det_row].astype(np.int16) + 128).astype(np.uint8)
    deltas   = D_TABLE[raw_u8]
    signal   = np.empty(NOBS, dtype=np.float64)
    step, i0 = 0, 0
    while i0 < NOBS:
        i1       = min(i0 + INITSTEP, NOBS)
        init_key = f'INIT{step:02}{det_row + 1:02}'
        init_val = float(h.get(init_key, 0)) + 32768.0  # header stores DN-32768
        seg      = np.empty(i1 - i0)
        seg[0]   = init_val
        if i1 - i0 > 1:
            seg[1:] = init_val + np.cumsum(deltas[i0 + 1: i1])
        signal[i0:i1] = seg
        i0 = i1; step += 1
    # A1: mark telemetry outage positions (BLANK_C = 0 compressed byte) as NaN.
    # D_TABLE[128] = 0 for a genuine zero-delta sample, so byte 0 is exclusively
    # the BLANK_C sentinel.  NaN propagates through calibrate_tod and is caught
    # by the ~np.isnan(mjy_sr) guard before map accumulation.
    signal[raw_u8 == 0] = np.nan
    # A4: bytes 1 and 255 are the D_TABLE extremes (delta = -61440 and +61440 DN
    # respectively, = 94% of the full 16-bit range).  A clipped delta means the
    # true DN step exceeded the encodable range (ADC saturation).
    signal[(raw_u8 == 1) | (raw_u8 == 255)] = np.nan
    return signal


def irim_glitch(mjy_sr, ftail=0.0, decay=0.0):
    """
    Detect glitch / spike and (for bands 1 & 2) tail samples following the
    GIPSY irim_flags algorithm (irscan.src, irim_flags()).

    Applies the 3-tap Mexican-hat filter (fiGli = {-3,6,-3}) and the 9-tap
    point-source filter (fiPnt) with circular boundary conditions.

    A sample is assigned SPIKEFLAG when its glitch-filter response exceeds
    5 sigma AND exceeds the concurrent point-source filter response (so
    genuine point sources are not misclassified as glitches).

    For bands 1 & 2 (ftail > 0), TAILFLAG is propagated forward after each
    SOURCEFLAG event using the vectorised state-machine equivalence:
        pstail_k = max_{j≤k}(bps_j + psdec·j) − psdec·k
    where bps_j = min(ftail/sn · ps_j, 80) if source[j] else −∞.
    This is exact modulo floating-point for standard survey speed (=1).

    Parameters
    ----------
    mjy_sr : ndarray, shape (N,)
    ftail  : float  — irscan.src ftail[band] (0 for bands 3 & 4)
    decay  : float  — irscan.src decay[band] (0 for bands 3 & 4)

    Returns
    -------
    spike  : bool ndarray  — SPIKEFLAG samples (excluded from co-add).
    source : bool ndarray  — SOURCEFLAG samples (used by C2; not excluded
                             from co-add unless reclassified by C2).
    tail   : bool ndarray  — TAILFLAG samples (excluded from co-add).
                             Always all-False for bands 3 & 4.
    """
    N = len(mjy_sr)
    spike  = np.zeros(N, dtype=bool)
    source = np.zeros(N, dtype=bool)
    tail   = np.zeros(N, dtype=bool)

    finite = np.isfinite(mjy_sr)
    if not finite.any():
        return spike, source, tail

    # Legacy step: replace outage/blank samples with the snip mean before
    # filtering ("set the outages temporarily to the average values").
    data = mjy_sr.copy()
    data[~finite] = data[finite].mean()

    # Noise estimate: replaces scan->sample.noise from IRDS.
    # Use MAD / 0.6745 so that the glitches we are trying to find don't
    # inflate the threshold that is used to detect them.
    med = np.median(data[finite])
    noise = np.median(np.abs(data[finite] - med)) / 0.6745
    if noise <= 0.0:
        return spike, source, tail

    thold = _IRIM_THRES * noise * _IRIM_SUMF

    gl = convolve1d(data, _FGLI, mode='wrap')
    ps = convolve1d(data, _FPNT, mode='wrap')

    src    = (ps > thold) & (ps > gl)
    glitch = (~src) & (gl > thold)

    # Never newly flag positions that are already NaN.
    glitch[~finite] = False
    src[~finite]    = False

    spike[:]  = glitch
    source[:] = src

    # D1 tail flagging (bands 1 & 2 only; ftail=0 for bands 3 & 4).
    if ftail > 0.0 and src.any():
        MAXPS = 80.0
        sn    = noise * _IRIM_SUMF
        psdec = decay   # speed = 1 for standard survey scan
        btail = ftail / sn
        # Vectorised state-machine (see docstring): bps_j stored at source
        # positions; -inf elsewhere so they don't influence the max-accumulate.
        bps = np.where(src, np.minimum(btail * ps, MAXPS), -np.inf)
        launch  = bps + psdec * np.arange(N, dtype=np.float64)
        pstail  = np.maximum.accumulate(launch) - psdec * np.arange(N, dtype=np.float64)
        tail[:] = (pstail > 0.0) & ~src & finite

    return spike, source, tail


# ── Pointing computation ───────────────────────────────────────────────────────
def pointing_det(fpath, yloc, zloc, ircc_times, rate):
    """
    Compute (satcal_hz, healpix_pix) for all sub-samples in a BPHF ATT file.

    Parameters
    ----------
    yloc, zloc    : detector focal-plane offsets in arcmin
    ircc_times    : readout time offset in SATCAL ticks (negative for B4)
    rate          : sub-samples per tick (4 for Band IV survey)

    Returns
    -------
    satcal_hz : float64 array, length n*rate
    pix       : int32 array,   length n*rate  (RING-ordered HEALPix pixel)
    """
    with afits.open(fpath) as ff:
        h   = ff[0].header
        raw = ff[0].data

    n      = raw.shape[0]
    crval2 = int(round(h['CRVAL2']))

    psis   = np.empty(n); thetas = np.empty(n)
    psis[0]   = float(h['PSI-0']);   thetas[0] = float(h['THETA-0'])
    psis[1:]   = psis[0]   + np.cumsum(raw[1:, 0].astype(np.float64))
    thetas[1:] = thetas[0] + np.cumsum(raw[1:, 1].astype(np.float64))
    twists = np.radians(raw[:, 2].astype(np.float64))
    sunlon = float(h['LNGSUN']) + np.arange(n) * float(h['SUNRATE'])

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
    pix      = hp.ang2pix(NSIDE, theta_hp, phi_hp, nest=False)

    satcal_hz = (crval2 + np.arange(n)[:,np.newaxis] +
                 ircc_times + k[np.newaxis,:]).ravel().astype(np.float64)
    return satcal_hz, pix.astype(np.int32)


# ── SATCAL → UTC conversion ────────────────────────────────────────────────────
def load_utcs_refs(path):
    """
    Parse utcs_satc.refs.  Format (space-separated):
        SOP  satcal_offset  utc_ref  scale_times_1e6

    UTC = utc_ref + (SATCAL - satcal_offset) * scale_times_1e6 / 1e6

    Returns dict: sop_int -> (sat_off, utc_ref, rate_corr)
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
    sat_off, utc_ref, rate_corr = refs[sop]
    return utc_ref + (satcal - sat_off) * rate_corr


def load_corrupt_times(path):
    """
    Parse rog11/outage/corrupt_times into a lookup dict.

    Format (space-separated, one entry per line):
        SOP  ATT  BAND  SATCAL  PRE  POST  det1 [det2 ...]

    The corrupt window is SATCAL ticks in [SATCAL+PRE, SATCAL+POST] (inclusive)
    for the listed detectors in that (SOP, ATT, BAND) observation.

    Returns
    -------
    dict : {(sop, att, band): [(satcal_ref, pre, post, frozenset(dets)), ...]}
        Empty dict if the file is missing or unreadable.
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
                sop  = int(parts[0])
                att  = int(parts[1])
                band = int(parts[2])
                satcal_ref = int(parts[3])
                pre  = int(parts[4])
                post = int(parts[5])
                dets = frozenset(int(x) for x in parts[6:])
                table.setdefault((sop, att, band), []).append(
                    (satcal_ref, pre, post, dets))
    except FileNotFoundError:
        pass
    return table


# ── SRHF loading per SOP ───────────────────────────────────────────────────────
def load_srhf_for_sop(srhf_dir, sop_int):
    """
    Load all SRHF files for a given SOP.
    Returns dict: obs_int -> (apl[62], bpl[62], tpl[62])
    """
    srhf_by_obs = {}
    pattern = os.path.join(srhf_dir, f'IR.SRHF.D0{sop_int:03d}*')
    for fpath in glob.glob(pattern):
        m = re.search(r'D0\d{3}(\d{3})$', os.path.basename(fpath))
        if not m:
            continue
        obs = int(m.group(1))
        try:
            _, _, apl, bpl, tpl = ip.read_srhf(srhf_dir, sop_int, obs)
            srhf_by_obs[obs] = (apl, bpl, tpl)
        except Exception:
            pass
    return srhf_by_obs


def closest_srhf(srhf_by_obs, obs_int):
    """Return (apl, bpl, tpl) from the SRHF entry whose OBS is closest."""
    if not srhf_by_obs:
        return None
    best = min(srhf_by_obs, key=lambda o: abs(o - obs_int))
    return srhf_by_obs[best]


def load_all_srhfs(srhf_dir):
    """
    Pre-load every SRHF file found on disk into a nested dict.
    Returns: {sop_int: {obs_int: (apl[62], bpl[62], tpl[62])}}
    Entries for inaccessible files are silently skipped.
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


def _lookup_srhf(all_srhfs, sop_int, obs_int):
    """
    Look up SRHF calibration data for a given SOP/OBS.

    Primary path: exact SOP, closest OBS.  If the SOP has no SRHF files,
    fall back to the nearest available SOP by number.  The SRHF flash
    amplitudes change slowly over the mission, so a neighbour-SOP estimate
    is a good approximation.

    Returns
    -------
    (entry, fallback_sop)
      entry        : (apl, bpl, tpl) or None if all_srhfs is empty
      fallback_sop : the SOP actually used, or None when own SOP matched
    """
    entry = closest_srhf(all_srhfs.get(sop_int, {}), obs_int)
    if entry is not None:
        return entry, None
    if not all_srhfs:
        return None, None
    nearest_sop = min(all_srhfs, key=lambda s: abs(s - sop_int))
    return closest_srhf(all_srhfs[nearest_sop], obs_int), nearest_sop


# ── utmbb: time since last bias boost ──────────────────────────────────────────
def utmbb_at_utc(utcs_start, bbt):
    """
    Return seconds elapsed since the last bias-boost event before utcs_start.
    Returns 1e9 if no prior event is found (detector treated as fully relaxed).
    """
    bbt_starts = bbt[:, 0]
    idx = np.searchsorted(bbt_starts, utcs_start, side='right') - 1
    if idx < 0:
        return 1.0e9
    return float(utcs_start - bbt[idx, 0])


# ── Inline calibration ─────────────────────────────────────────────────────────
def calibrate_tod(dn, det, utcs, utmbb,
                  a2dc, gains, offsets, utc_tbl, vfets_tbl,
                  apl, bpl, tpl, bbt):
    """
    Calibrate one detector's decompressed timeline to MJy/sr.

    Parameters
    ----------
    dn           : float64 array (NOBS,)  decompressed ADU
    det          : int,  1-based detector number
    utcs         : float64 array (NOBS,)  UTC seconds for each sample
    utmbb        : float64 array (NOBS,)  seconds since last bias boost
    a2dc, gains, offsets : from load_ctype2()
    utc_tbl, vfets_tbl   : from load_vfet()
    apl, bpl, tpl        : SRHF flash-amplitude arrays (62 elements)
    bbt                  : from load_bbtimes()

    Returns
    -------
    mjy_sr : float64 array (NOBS,).  ICF-flagged samples are set to NaN.
    """
    flags  = ip.icf_flag(utcs, bbt)
    volts, _ = ip.dn_to_volts(dn, det, a2dc, gains, offsets, 'standard')
    volts  = ip.subtract_baseline(volts, det, utcs, utc_tbl, vfets_tbl, utmbb)
    amps   = ip.volts_to_amps(volts)
    wm2    = ip.amps_to_wm2(amps, det, utmbb, apl, bpl, tpl)
    mjy_sr = ip.wm2_to_mjy_sr(wm2, det)
    mjy_sr[flags] = np.nan
    return mjy_sr


# ── Process one path segment ───────────────────────────────────────────────────
def process_seg(seg_str, sum_map, hits_map,
                a2dc, gains, offsets, utc_tbl, vfets_tbl, bbt, utcs_refs,
                all_srhfs, corrupt_table=None):
    """
    Process all BPHF + raw Band IV files in one sop.{seg}_ directory pair.
    Calibrated MJy/sr samples are accumulated into sum_map / hits_map in-place.

    Returns (n_matched, stats_dict).
    """
    _zero = dict(n_raw=0, n_exc=0, n_no_utc=0, n_no_srhf=0, n_srhf_fallback=0)
    bphf_files = sorted(glob.glob(BPHF_GLOB.format(seg=seg_str)))
    raw_files  = sorted(glob.glob(RAW_GLOB.format(seg=seg_str)))
    if not bphf_files or not raw_files:
        return 0, dict(_zero, n_raw=len(raw_files))

    # ── Build per-detector pointing table (sorted by SATCAL) ────────────────
    pointing_dets = {}   # det -> (bsc_sorted, bpx_sorted)
    for det in BAND_DETS:
        sc_list, pix_list = [], []
        for fpath in bphf_files:
            try:
                sc, pix = pointing_det(fpath, BAND_YLOC[det], BAND_ZLOC[det],
                                       BAND_IRCC[det], RATE)
                sc_list.append(sc); pix_list.append(pix)
            except Exception:
                pass
        if not sc_list:
            continue
        bsc = np.concatenate(sc_list)
        bpx = np.concatenate(pix_list)
        order = np.argsort(bsc)
        pointing_dets[det] = (bsc[order], bpx[order])

    if not pointing_dets:
        return 0, dict(_zero, n_raw=len(raw_files))

    # ── Process raw TOD files ─────────────────────────────────────────────────
    n_matched       = 0
    n_exc           = 0   # raw files dropped due to unhandled exceptions
    n_no_utc        = 0   # raw files dropped: SOP absent from utcs_refs
    n_no_srhf       = 0   # raw files dropped: no SRHF anywhere (all_srhfs empty)
    n_srhf_fallback = 0   # raw files calibrated using nearest-SOP SRHF
    bad_files       = []  # files that produced non-finite calibrated values
    max_gap         = 1.0 / RATE + 0.01   # matching tolerance in SATCAL ticks

    for fpath in raw_files:
        try:
            with fitsio.FITS(fpath) as ff:
                raw_h = ff[0].read_header()
                data  = ff[0].read()   # (15, NOBS) int8

            sop_int = int(raw_h['SOP'])
            obs_int = int(raw_h['OBS'])
            # A2: fetch corrupt-times windows for this SOP/OBS/band snip.
            ct_windows = corrupt_table.get((sop_int, obs_int, BAND), []) if corrupt_table else []
            satcal0 = float(raw_h['SATCAL'])
            cdelt1  = float(raw_h['CDELT1'])   # 0.25 ticks/sample for B4
            NOBS    = data.shape[1]

            if sop_int not in utcs_refs:
                n_no_utc += 1
                continue

            # SATCAL axis for pointing match
            tod_sc = satcal0 + np.arange(NOBS, dtype=np.float64) * cdelt1

            # B2: FLASHTIMES scan-edge trim.
            # FLASH_PRE / FLASH_POST are in SATCAL ticks; cdelt1 = 0.25 ticks/sample.
            n_pre  = int(round(FLASH_PRE  / cdelt1))   # = 128 samples
            n_post = int(round(FLASH_POST / cdelt1))   # =  80 samples
            flash_trim = np.zeros(NOBS, dtype=bool)
            if n_pre  > 0 and n_pre  < NOBS: flash_trim[:n_pre]   = True
            if n_post > 0 and n_post < NOBS: flash_trim[-n_post:] = True

            _, utc_ref, rate_corr = utcs_refs[sop_int]
            utcs_start = satcal_to_utc(satcal0, sop_int, utcs_refs)
            utcs_inc   = cdelt1 * rate_corr
            utcs       = utcs_start + np.arange(NOBS, dtype=np.float64) * utcs_inc

            # D2: planet guard — pixel set for the midpoint time of this file
            planet_pix = planet_pixel_set(utcs[NOBS // 2], NSIDE)

            # Time since last bias boost
            utmbb_start = utmbb_at_utc(utcs_start, bbt)
            utmbb       = utmbb_start + np.arange(NOBS, dtype=np.float64) * utcs_inc

            # SRHF flash-calibration parameters (with cross-SOP fallback)
            srhf, fallback_sop = _lookup_srhf(all_srhfs, sop_int, obs_int)
            if srhf is None:
                n_no_srhf += 1
                continue
            if fallback_sop is not None:
                n_srhf_fallback += 1
            apl, bpl, tpl = srhf

            # ── Pass 1: calibrate all detectors + compute C1 spike/source/tail flags
            det_cache = {}   # det -> dict(mjy_sr, idx, valid, corrupt, spike, source, tail)
            for det in BAND_DETS:
                if det not in pointing_dets:
                    continue
                bsc, bpx = pointing_dets[det]

                # Nearest-neighbour match in SATCAL time
                idx    = np.searchsorted(bsc, tod_sc)
                row_idx = _DET_TO_ROW[det]
                idx    = np.clip(idx, 0, len(bsc) - 1)
                idx_m1 = np.maximum(idx - 1, 0)
                closer = np.abs(bsc[idx_m1] - tod_sc) < np.abs(bsc[idx] - tod_sc)
                idx    = np.where(closer, idx_m1, idx)
                valid  = np.abs(bsc[idx] - tod_sc) <= max_gap
                if not valid.any():
                    continue

                dn = decompress_det(data, row_idx, raw_h)

                # A2: corrupt-times rejection mask for this detector.
                corrupt = np.zeros(NOBS, dtype=bool)
                for sat_ref, pre, post, dets in ct_windows:
                    if det in dets:
                        corrupt |= (tod_sc >= sat_ref + pre) & (tod_sc <= sat_ref + post)

                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    mjy_sr = calibrate_tod(
                        dn, det, utcs, utmbb,
                        a2dc, gains, offsets, utc_tbl, vfets_tbl,
                        apl, bpl, tpl, bbt)

                # C1: glitch / spike / tail detection (irim_flags).
                spike, source, tail = irim_glitch(mjy_sr,
                                                  ftail=_bcfg['ftail'],
                                                  decay=_bcfg['decay'])

                det_cache[det] = dict(mjy_sr=mjy_sr, idx=idx, valid=valid,
                                      corrupt=corrupt, spike=spike,
                                      source=source, tail=tail)

            # ── C2: large-CR persistence — multi-detector coherence check ─────
            # A large CR causes slow detector hysteresis; the resulting broad
            # profile is too wide for the 3-tap glitch filter and gets
            # misclassified as SOURCEFLAG.  In the legacy pipeline FLAGS=p,t
            # means SOURCEFLAG survives the co-add.  We correct this by checking
            # spatial coherence across same-row detectors: a genuine bright source
            # illuminates adjacent in-scan detectors sequentially; a CR hits
            # exactly one detector.  Lone-detector SOURCEFLAG events are
            # reclassified as spike (excluded from the co-add).
            src_pix = {}   # det -> set of HEALPix pixels carrying SOURCEFLAG
            for det, d in det_cache.items():
                usable = (d['source'] & d['valid']
                          & ~np.isnan(d['mjy_sr']) & ~d['corrupt'] & ~flash_trim)
                if usable.any():
                    bpx_d = pointing_dets[det][1]
                    src_pix[det] = set(bpx_d[d['idx'][usable]].tolist())
                else:
                    src_pix[det] = set()

            for det, d in det_cache.items():
                if not src_pix[det]:
                    continue
                row = _band_row(det)
                # Build union of source pixels + 1-ring HEALPix neighbours from
                # all OTHER detectors in the same focal-plane row.
                other_pix = set()
                for det2, pix2 in src_pix.items():
                    if det2 == det or _band_row(det2) != row or not pix2:
                        continue
                    p2_arr = np.fromiter(pix2, dtype=np.int64, count=len(pix2))
                    nb = hp.get_all_neighbours(NSIDE, p2_arr)   # (8, N), RING
                    nb_flat = nb.ravel()
                    other_pix.update(p2_arr.tolist())
                    other_pix.update(nb_flat[nb_flat >= 0].tolist())
                if not other_pix:
                    # No other same-row detector has SOURCEFLAG this snip;
                    # cannot discriminate CR from real source — leave unchanged.
                    continue
                usable = (d['source'] & d['valid']
                          & ~np.isnan(d['mjy_sr']) & ~d['corrupt'] & ~flash_trim)
                bpx_d      = pointing_dets[det][1]
                sample_pix = bpx_d[d['idx']]
                other_arr  = np.fromiter(other_pix, dtype=np.int64, count=len(other_pix))
                lone = usable & ~np.isin(sample_pix, other_arr)
                d['spike'] |= lone   # reclassify large-CR SOURCEFLAG as spike

            # ── Pass 2: accumulate ─────────────────────────────────────────────
            for det, d in det_cache.items():
                bpx_d = pointing_dets[det][1]
                planet = (np.isin(bpx_d[d['idx']], planet_pix)
                          if planet_pix.size else np.zeros(NOBS, dtype=bool))
                good  = (d['valid'] & np.isfinite(d['mjy_sr'])
                         & ~d['corrupt'] & ~flash_trim & ~d['spike'] & ~d['tail']
                         & ~planet)
                if not MASK:
                    # Diagnostic: skip all QC flags, keep only pointing validity
                    # and physical blanks (NaN/Inf).  Produces a calibrated but
                    # unmasked map useful for diagnosing coverage gaps.
                    good = d['valid'] & np.isfinite(d['mjy_sr'])
                pix_v = bpx_d[d['idx'][good]]
                val_v = d['mjy_sr'][good]
                # Paranoia check: the isfinite mask above should make this
                # impossible, but catch it here so we know exactly which file
                # produced non-finite values.
                if not np.all(np.isfinite(val_v)):
                    bad_idx  = np.where(~np.isfinite(val_v))[0]
                    bad_vals = val_v[bad_idx]
                    bad_files.append(
                        f'{fpath}  det={det}  '
                        f'{len(bad_idx)} non-finite sample(s): '
                        f'{bad_vals[:5].tolist()}{" ..." if len(bad_idx) > 5 else ""}')
                    # Strip them so they never enter sum_map.
                    keep  = np.isfinite(val_v)
                    pix_v = pix_v[keep]
                    val_v = val_v[keep]
                np.add.at(sum_map,  pix_v, val_v)
                np.add.at(hits_map, pix_v, 1)
                n_matched += int(len(pix_v))

        except Exception as exc:
            n_exc += 1
            import traceback
            print(f'  [WARN] seg={seg_str} file={os.path.basename(fpath)}: '
                  f'{type(exc).__name__}: {exc}', file=sys.stderr)
            if os.environ.get('TRACEBACK'):
                traceback.print_exc()

    stats = dict(
        n_raw=len(raw_files),
        n_exc=n_exc,
        n_no_utc=n_no_utc,
        n_no_srhf=n_no_srhf,
        n_srhf_fallback=n_srhf_fallback,
        bad_files=bad_files,
    )
    return n_matched, stats


# ── Multiprocessing worker ─────────────────────────────────────────────────────
def _seg_worker(args):
    """
    Top-level worker for multiprocessing.Pool.  Creates local accumulation
    arrays, delegates to process_seg, and returns the partial results.
    Defined at module level so it can be pickled by the worker processes.

    args = (seg_str, npix, a2dc, gains, offsets,
            utc_tbl, vfets_tbl, bbt, utcs_refs, all_srhfs, corrupt_table)
    """
    (seg_str, npix, a2dc, gains, offsets,
     utc_tbl, vfets_tbl, bbt, utcs_refs, all_srhfs, corrupt_table) = args
    local_sum  = np.zeros(npix, dtype=np.float64)
    local_hits = np.zeros(npix, dtype=np.int32)
    n, stats = process_seg(seg_str, local_sum, local_hits,
                           a2dc, gains, offsets,
                           utc_tbl, vfets_tbl, bbt, utcs_refs,
                           all_srhfs, corrupt_table)
    return n, stats, local_sum, local_hits


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    print(f'IRAS {_bcfg["name"]} — all {len(BAND_DETS)} detectors — calibrated HEALPix map')
    print(f'NSIDE={NSIDE}  ({hp.nside2resol(NSIDE, arcmin=True):.1f} arcmin pixels)')
    print()

    # Preload calibration tables (done once for the entire run)
    print('Loading calibration tables...')
    a2dc, gains, offsets = ip.load_ctype2(CTYPE2)
    utc_tbl, vfets_tbl   = ip.load_vfet(VFET)
    bbt                   = ip.load_bbtimes(BBT_FILE)
    utcs_refs             = load_utcs_refs(UTCS_REF)
    corrupt_table         = load_corrupt_times(CORRUPT_TIMES)
    print('Loading SRHF flash-calibration files (all SOPs)...')
    all_srhfs             = load_all_srhfs(SRHF_DIR)
    n_srhf_sops = len(all_srhfs)
    n_srhf_obs  = sum(len(v) for v in all_srhfs.values())
    print(f'  ctype2 : a2dc = {a2dc:.4e} V/DN·gain')
    print(f'  vfet   : {len(utc_tbl)} epochs,  UTC {utc_tbl[0]:.0f} – {utc_tbl[-1]:.0f}')
    print(f'  bbtimes: {len(bbt)} events')
    print(f'  UTCS   : {len(utcs_refs)} SOP entries')
    n_ct = sum(len(v) for v in corrupt_table.values())
    print(f'  corrupt: {n_ct} windows across {len(corrupt_table)} (SOP,ATT,band) triplets')
    print(f'  SRHF   : {n_srhf_obs} observations across {n_srhf_sops} SOPs')
    print()

    # Find all path segments with both BPHF and raw data for this band
    bphf_segs = sorted(set(
        m.group(1)
        for f in glob.glob(BPHF_GLOB.format(seg='*'))
        for m in [re.search(r'sop\.(\d+)_', f)] if m))
    raw_segs  = sorted(set(
        m.group(1)
        for f in glob.glob(RAW_GLOB.format(seg='*'))
        for m in [re.search(r'sop\.(\d+)_', f)] if m))
    segs = sorted(set(bphf_segs) & set(raw_segs))
    if MAX_SEG is not None:
        segs = [s for s in segs if s <= MAX_SEG]
    print(f'Path segments (BPHF + {_bcfg["survey_tag"].upper()} raw both present): {len(segs)}')
    print(f'Detectors: {BAND_DETS}')
    print()

    npix     = hp.nside2npix(NSIDE)
    sum_map  = np.zeros(npix, dtype=np.float64)
    hits_map = np.zeros(npix, dtype=np.int32)
    total    = 0

    nproc = multiprocessing.cpu_count() if NPROC == 0 else NPROC
    print(f'Processing {len(segs)} segments with {nproc} worker(s)...')

    worker_args = [
        (seg, npix, a2dc, gains, offsets,
         utc_tbl, vfets_tbl, bbt, utcs_refs, all_srhfs, corrupt_table)
        for seg in segs
    ]

    tot_raw = tot_exc = tot_no_utc = tot_no_srhf = tot_srhf_fallback = 0
    all_bad_files = []

    if nproc == 1:
        # Serial path — keeps tqdm progress bar and avoids fork overhead.
        for args in tqdm(worker_args, desc='segments'):
            n, st, s, h = _seg_worker(args)
            total              += n
            sum_map            += s
            hits_map           += h
            tot_raw            += st['n_raw']
            tot_exc            += st['n_exc']
            tot_no_utc         += st['n_no_utc']
            tot_no_srhf        += st['n_no_srhf']
            tot_srhf_fallback  += st['n_srhf_fallback']
            all_bad_files      += st['bad_files']
    else:
        # Parallel path — each worker processes one segment independently;
        # partial maps are summed in the main process.
        with multiprocessing.Pool(processes=nproc) as pool:
            for n, st, s, h in tqdm(
                    pool.imap_unordered(_seg_worker, worker_args),
                    total=len(worker_args), desc='segments'):
                total              += n
                sum_map            += s
                hits_map           += h
                tot_raw            += st['n_raw']
                tot_exc            += st['n_exc']
                tot_no_utc         += st['n_no_utc']
                tot_no_srhf        += st['n_no_srhf']
                tot_srhf_fallback  += st['n_srhf_fallback']
                all_bad_files      += st['bad_files']

    tot_skipped = tot_exc + tot_no_utc + tot_no_srhf
    tot_ok      = tot_raw - tot_skipped
    print()
    print(f'Raw file summary:')
    print(f'  Total found                   : {tot_raw:,}')
    print(f'  Processed OK                  : {tot_ok:,}')
    if tot_srhf_fallback:
        print(f'    of which, fallback SRHF     : {tot_srhf_fallback:,}  (nearest-SOP flash data)')
    print(f'  Skipped (no UTC ref)          : {tot_no_utc:,}')
    print(f'  Skipped (no SRHF anywhere)    : {tot_no_srhf:,}')
    print(f'  Dropped (exceptions)          : {tot_exc:,}')
    if tot_skipped > 0:
        print(f'  *** {tot_skipped:,} raw files contributed no data ***')
        if tot_exc > 0:
            print(f'      Re-run with TRACEBACK=1 to see full stack traces.')

    hit_pix = np.count_nonzero(hits_map)
    print()
    print(f'Total matched+valid samples : {total:,}')
    print(f'Pixels hit                  : {hit_pix:,}  '
          f'({hit_pix / npix * 100:.1f}% of sky at nside={NSIDE})')
    if hit_pix > 0:
        print(f'Max hits per pixel          : {hits_map.max()}')
        print(f'Median hits (non-zero pix)  : {np.median(hits_map[hits_map > 0]):.0f}')

    if all_bad_files:
        print()
        print(f'ERROR: {len(all_bad_files)} file/detector combination(s) produced '
              f'non-finite calibrated values (stripped before accumulation):')
        for entry in all_bad_files:
            print(f'  {entry}')
        raise RuntimeError(
            f'{len(all_bad_files)} non-finite calibration result(s) — '
            f'see list above.  Fix the calibration or add those files to '
            f'the corrupt_times table.')

    # Build signal map
    signal_map = np.full(npix, hp.UNSEEN, dtype=np.float64)
    mask = hits_map > 0
    signal_map[mask] = sum_map[mask] / hits_map[mask]
    # Guard: the isfinite accumulation check above should make this impossible;
    # if it fires there is a bug in the accumulation logic.
    bad = mask & ~np.isfinite(signal_map)
    if bad.any():
        raise RuntimeError(
            f'{bad.sum()} pixel(s) have non-finite signal after co-add despite '
            f'per-file isfinite guard.  This is a bug — please report it with '
            f'the pixel indices: {np.where(bad)[0][:20].tolist()}')

    # Save FITS
    btag     = _bcfg['survey_tag']
    mask_sfx = '' if MASK else '_nomask'
    out_sig  = os.path.join(OUT_DIR, f'fullsurvey_{btag}_alldet_cal_map_nside{NSIDE}{mask_sfx}.fits')
    out_hits = os.path.join(OUT_DIR, f'fullsurvey_{btag}_alldet_cal_hits_nside{NSIDE}{mask_sfx}.fits')
    hp.write_map(out_sig, signal_map, nest=False, dtype=np.float64, coord='G',
                 extra_header=[
                     ('NSIDE',    NSIDE,              'HEALPix nside'),
                     ('ORDERING', 'RING',              'HEALPix ordering'),
                     ('COORDSYS', 'G',                 'G=galactic'),
                     ('BUNIT',    'MJy/sr',            'IRAS IPACCAL calibrated flux'),
                     ('BANDNAME', btag.upper(),        f'IRAS Band {BAND} ({_bcfg["wavelength_um"]}um)'),
                     ('WEIGHT',   'EQUAL',             'Arithmetic mean; w_i = 1'),
                 ], overwrite=True)
    print(f'\nSaved: {out_sig}')

    hp.write_map(out_hits, hits_map.astype(float), nest=False, coord='G',
                 overwrite=True)
    print(f'Saved: {out_hits}')

    # Plot
    if hit_pix > 100:
        vals   = signal_map[mask]
        finite = vals[np.isfinite(vals)]
        vmin, vmax = np.percentile(finite, 1), np.percentile(finite, 99)
        fig, axes = plt.subplots(1, 2, figsize=(16, 5))
        fig.suptitle(
            f'IRAS {_bcfg["name"]} -- all {len(BAND_DETS)} detectors -- calibrated map (nside={NSIDE})',
            fontsize=13)
        plt.sca(axes[0])
        hp.mollview(signal_map, title='Signal (MJy/sr)', unit='MJy/sr',
                    min=vmin, max=vmax, coord='G', cmap='inferno', hold=True)
        hp.graticule(dpar=30, dmer=60, alpha=0.4)
        plt.sca(axes[1])
        hp.mollview(hits_map.astype(float), title='Hits per pixel',
                    coord='G', cmap='Blues', hold=True,
                    min=0, max=np.percentile(hits_map[mask], 99))
        hp.graticule(dpar=30, dmer=60, alpha=0.4)
        plt.tight_layout()
        out_png = os.path.join(OUT_DIR, f'fullsurvey_{btag}_alldet_cal_map_nside{NSIDE}{mask_sfx}.png')
        plt.savefig(out_png, dpi=150)
        print(f'Saved: {out_png}')
        plt.show()


if __name__ == '__main__':
    main()
