#!/usr/bin/env python3
"""
Full-survey calibrated HEALPix map for any IRAS photometric band.

Imports IRAS TOD access infrastructure from iras_tod (see that module for
data-model details).  This script focuses on the map-making steps:
  1. Calibrate DN → MJy sr⁻¹ (IPACCAL pipeline via ipaccal.py).
  2. Detect and flag glitches / spikes / tails (C1 irim_flags).
  3. Multi-detector coherence check for large cosmic-ray events (C2).
  4. Accumulate valid samples into sum_map / hits_map (equal-weight mean).

Band, masking, and AO-inclusion are controlled by environment variables.

Output
------
  python/IRAS_{tag}_cal_map_nside{N}[_nomask][_noao].fits
  python/IRAS_{tag}_cal_hits_nside{N}[_nomask][_noao].fits
  python/IRAS_{tag}_cal_map_nside{N}[_nomask][_noao].png

Environment variables
---------------------
  BAND=4        Photometric band (1-4, default 4).
  NSIDE=512     HEALPix resolution (default 512).
  MASK=1        Enable all QC flags; MASK=0 = diagnostic unmasked map.
    INCLUDE_AO=0  Include AO/SPLINE/UNKNOWN observation directories.
  NPROC=1       Worker processes; 0 = all CPUs.
  MAX_SEG=None  Process only segments ≤ MAX_SEG (e.g. MAX_SEG=03).
    DESTRIPE=0         E1 per-snip polynomial destripe (1 = enabled).
                                         Fits and subtracts a robust polynomial baseline from each
                                         detector-observation after calibration and glitch detection.
                                         Removes inter-detector DC offsets and slow thermal drifts
                                         that survive the IPACCAL vfet / GXN correction.  Required
                                         for useful AO / SPLINE maps; optional for survey (minor
                                         improvement to faint extended emission).
    DESTRIPE_DEGREE=1  Polynomial degree for DESTRIPE (0 = offset, 1 = linear,
                                         2 = quadratic).  Linear is appropriate for most uses.
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
from scipy.ndimage import convolve1d
from tqdm import tqdm

# ── iras_tod provides all TOD access infrastructure ───────────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import iras_tod as tod
import ipaccal as ip

# ── Configuration ─────────────────────────────────────────────────────────────
NSIDE      = int(os.environ.get('NSIDE',      '512'))
BAND       = int(os.environ.get('BAND',       '4'))
MASK       = int(os.environ.get('MASK',       '1'))
INCLUDE_AO = int(os.environ.get('INCLUDE_AO', '0'))
MAX_SEG    = os.environ.get('MAX_SEG', None)
NPROC      = int(os.environ.get('NPROC', '1'))
DESTRIPE        = int(os.environ.get('DESTRIPE',         '0'))
DESTRIPE_DEGREE = int(os.environ.get('DESTRIPE_DEGREE',  '1'))
LEVEL1_PARITY_MASK = int(os.environ.get('LEVEL1_PARITY_MASK', '0'))
OUTAGE_ALLDET      = int(os.environ.get('OUTAGE_ALLDET', '0'))

# Band-specific derived constants (use BAND_CFG from iras_tod)
_bcfg      = tod.BAND_CFG[BAND]
_ALL_DETS  = _bcfg['dets']
_DEAD_DETS = _bcfg['dead']
BAND_DETS  = [d for d in _ALL_DETS if d not in _DEAD_DETS]
BAND_ROW1  = _bcfg['row1'] - _DEAD_DETS
BAND_ROW2  = _bcfg['row2'] - _DEAD_DETS
_ircc_list = [t - _bcfg['edelay'] for t in _bcfg['times']]
BAND_YLOC  = dict(zip(_ALL_DETS, _bcfg['yloc']))
BAND_ZLOC  = dict(zip(_ALL_DETS, _bcfg['zloc']))

OUT_DIR    = os.path.dirname(os.path.abspath(__file__))
BPHF_GLOB  = tod.DEFAULT_PATHS['bphf_glob']

# ── Zodycal sky map (loaded once at startup) ──────────────────────────────────
# Only bands 3 & 4 have sky maps (file000001.mt / file000002.mt); bands 1 & 2
# receive no zodycal gain correction (irzc_gain is not called for them).
_ZODYCAL_MAPS_DIR = os.path.join(tod.ROOT,
    'diskrog10androg11/rog11/zodycal/MAPS')
if BAND >= 3:
    try:
        _ZODYCAL_SKY_MAP = ip.load_zodycal_map(_ZODYCAL_MAPS_DIR, band=BAND)
    except Exception as _ze:
        import warnings as _wmod
        _wmod.warn(f'zodycal sky map not loaded ({_ze}); gain correction disabled')
        _ZODYCAL_SKY_MAP = None
else:
    _ZODYCAL_SKY_MAP = None   # no zodycal gain correction for B1/B2

# C1: irim_flags glitch/tail filter constants (irscan.src irim_flags())
_FGLI       = np.array([-3.0,  6.0, -3.0])
_FPNT       = np.array([-1.0, -1.0, -1.0,  2.0, 2.0, 2.0,  -1.0, -1.0, -1.0])
_IRIM_SUMF  = 12.0
_IRIM_THRES =  5.0


def _band_row(det):
    """Return focal-plane row (1 or 2) for a detector number."""
    return 1 if det in BAND_ROW1 else 2


def irim_glitch(mjy_sr, ftail=0.0, decay=0.0, speed=1.0):
    """
    Detect glitch / spike and (for bands 1 & 2) tail samples following the
    GIPSY irim_flags algorithm (irscan.src, irim_flags()).

    Applies the 3-tap Mexican-hat filter (fiGli = {-3,6,-3}) and the 9-tap
    point-source filter (fiPnt) with circular boundary conditions.

    A sample is assigned SPIKEFLAG when its glitch-filter response exceeds
    5 sigma AND exceeds the concurrent point-source filter response (so
    genuine point sources are not misclassified as glitches).

    For slow scans (speed < 1), the SOURCEFLAG exclusion zone is extended
    symmetrically using the exponential decay tau = 0.1^(1/slow), matching
    the irim_flags() SLOW= logic.  At survey speed (speed=1), tau=0.1 and
    the flanks are negligible (only 0–1 extra samples for typical sources).

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
    speed  : float  — normalised scan speed (1.0 for survey, <1 for AO slow scans)

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

    # Legacy step: replace outage/blank samples with the observation mean
    # before filtering ("set the outages temporarily to the average values").
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

    # Extend SOURCEFLAG exponentially on both sides of each detected source
    # (mirrors irim_flags C code: "set flag on exponential sides of a ps").
    # slow = 1/speed; tau = 0.1^(1/max(slow,1)).
    # For survey speed: tau=0.1 → flanks are at most 1 sample wide.
    # For 2× slow (AO): tau=0.316 → flanks extend ~2-3 samples per side.
    _slow = 1.0 / max(speed, 1e-3)
    _tau  = 0.1 ** (1.0 / max(_slow, 1.0))
    if _tau > 0.1 and src.any():
        # Only run the loop for slow scans (tau>0.1), since at survey speed
        # the flanks are always negligibly narrow.
        for k in np.where(src)[0]:
            ps_cur = ps[k] * _tau
            j = int(k); i = int(k)
            while ps_cur > thold:
                if j > 0:
                    j -= 1
                    source[j] = True
                if i < N - 1:
                    i += 1
                    source[i] = True
                ps_cur *= _tau

    # D1 tail flagging (bands 1 & 2 only; ftail=0 for bands 3 & 4).
    if ftail > 0.0 and src.any():
        MAXPS = 80.0
        sn    = noise * _IRIM_SUMF
        # irim_flags formula: psdec = decay * max(v, 1).
        # For survey / slow AO scans (v ≤ 1): psdec = decay (unchanged).
        # For fast scans (v > 1): psdec increases so the tail decays faster per sample.
        psdec = decay * max(speed, 1.0)
        btail = ftail / sn
        # Vectorised state-machine (see docstring): bps_j stored at source
        # positions; -inf elsewhere so they don't influence the max-accumulate.
        bps = np.where(src, np.minimum(btail * ps, MAXPS), -np.inf)
        launch  = bps + psdec * np.arange(N, dtype=np.float64)
        pstail  = np.maximum.accumulate(launch) - psdec * np.arange(N, dtype=np.float64)
        tail[:] = (pstail > 0.0) & ~src & finite

    return spike, source, tail



# ── Inline calibration ────────────────────────────────────────────────────────
def calibrate_tod(det, obs, cal):
    """
    Calibrate one detector within an IrasObs to MJy sr⁻¹.

    Parameters
    ----------
    det  : int      — 1-based IRAS detector number (must be a key in obs.dn).
    obs  : IrasObs  — complete time-ordered data for one (SOP, OBS).
    cal  : CalTables — shared calibration tables.

    Returns
    -------
    mjy_sr : float64 array (NOBS,).
        ICF-flagged samples are set to NaN. If LEVEL1_PARITY_MASK=1, also
        applies PLATE-like pre-map masking:
          - flash_trim (survey scan-edge FLASHTIMES)
          - outage/corrupt windows (per-detector by default, or all-detector
            if OUTAGE_ALLDET=1 to emulate legacy IRPL_UNTAB behavior).
    """
    flags    = ip.icf_flag(obs.utcs, cal.bbt)

    if LEVEL1_PARITY_MASK:
        flags = flags | obs.flash_trim
        if OUTAGE_ALLDET:
            all_corrupt = np.zeros_like(flags)
            for dmask in obs.corrupt.values():
                all_corrupt |= dmask
            flags = flags | all_corrupt
        else:
            flags = flags | obs.corrupt.get(det, np.zeros_like(flags))

    volts, _ = ip.dn_to_volts(obs.dn[det], det, cal.a2dc, cal.gains,
                               cal.offsets, 'standard')
    volts    = ip.subtract_baseline(volts, det, obs.utcs,
                                    cal.utc_tbl, cal.vfets_tbl, obs.utmbb)
    amps     = ip.volts_to_amps(volts)
    wm2      = ip.amps_to_wm2(amps, det, obs.utmbb,
                               obs.apl, obs.bpl, obs.tpl,
                               prhf_deltut=obs.prhf_deltut,
                               prhf_dtr=obs.prhf_dtr)
    mjy_sr   = ip.wm2_to_mjy_sr(wm2, det)
    mjy_sr[flags] = np.nan
    return mjy_sr


# ── Process one path segment ───────────────────────────────────────────────────
def process_seg(seg_str, sum_map, hits_map, cal):
    """
    Calibrate and accumulate all TOD in one sop.{seg}_ directory pair.

    Uses iter_obs() to stream one IrasObs per (SOP, OBS) pair.  Each obs
    already has all its raw files concatenated across time, so no further
    file-level buffering is needed.  Applies C1 glitch detection, C2
    coherence check, and accumulates valid samples into sum_map / hits_map.

    Returns (n_matched, stats_dict).
    """
    n_matched  = 0
    n_obs      = 0
    bad_files  = []

    for obs in tod.iter_obs(BAND, nside=NSIDE, cal=cal,
                             root=tod.ROOT, bphf_glob=BPHF_GLOB,
                             include_ao=bool(INCLUDE_AO),
                             segments=[seg_str]):
        n_obs += 1
        obs_label = f'SOP{obs.sop:04d}/OBS{obs.obs:03d}'
        try:
            NOBS       = len(obs.satcal)
            planet_pix = tod.planet_pixel_set(obs.utcs[NOBS // 2], NSIDE)

            # ── Pass 1: calibrate + ZODYCAL correction + C1 glitch detect ────
            det_cache = {}   # det -> dict(mjy_sr, spike, source, tail)
            for det in obs.dn:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    mjy_sr = calibrate_tod(det, obs, cal)

                # ZODYCAL position-dependent gain correction (IRZC_GAIN).
                if _ZODYCAL_SKY_MAP is not None and ip._BAND[det - 1] >= 3:
                    try:
                        _gain_norm = ip.irzc_gain(
                            det,
                            obs.psi0_deg, obs.psirate, obs.theta,
                            obs.lambda_obs_rad,
                            obs.cdelt1, NOBS,
                            _ZODYCAL_SKY_MAP,
                            BAND_YLOC[det], BAND_ZLOC[det])
                        with np.errstate(invalid='ignore', divide='ignore'):
                            _safe = np.where(_gain_norm > 0.0, _gain_norm, 1.0)
                            mjy_sr = mjy_sr / _safe
                    except Exception:
                        pass

                spike, source, tail = irim_glitch(
                    mjy_sr,
                    ftail=_bcfg['ftail'], decay=_bcfg['decay'],
                    speed=obs.scan_speed)

                det_cache[det] = dict(mjy_sr=mjy_sr,
                                      spike=spike, source=source, tail=tail)

            # ── C2: large-CR persistence — multi-detector coherence check ─────
            # A large CR causes slow detector hysteresis; the resulting broad
            # profile is too wide for the 3-tap glitch filter and gets
            # misclassified as SOURCEFLAG.  Lone-detector SOURCEFLAG events
            # are reclassified as spike (excluded from the co-add).
            src_pix = {}   # det -> set of HEALPix pixels carrying SOURCEFLAG
            for det, d in det_cache.items():
                usable = (d['source'] & (obs.pix[det] >= 0)
                          & ~np.isnan(d['mjy_sr'])
                          & ~obs.corrupt[det] & ~obs.flash_trim)
                if usable.any():
                    src_pix[det] = set(obs.pix[det][usable].tolist())
                else:
                    src_pix[det] = set()

            for det, d in det_cache.items():
                if not src_pix[det]:
                    continue
                row       = _band_row(det)
                other_pix = set()
                for det2, pix2 in src_pix.items():
                    if det2 == det or _band_row(det2) != row or not pix2:
                        continue
                    p2_arr = np.fromiter(pix2, dtype=np.int64, count=len(pix2))
                    nb      = hp.get_all_neighbours(NSIDE, p2_arr)
                    nb_flat = nb.ravel()
                    other_pix.update(p2_arr.tolist())
                    other_pix.update(nb_flat[nb_flat >= 0].tolist())
                if not other_pix:
                    continue
                usable    = (d['source'] & (obs.pix[det] >= 0)
                             & ~np.isnan(d['mjy_sr'])
                             & ~obs.corrupt[det] & ~obs.flash_trim)
                other_arr = np.fromiter(other_pix, dtype=np.int64, count=len(other_pix))
                lone = usable & ~np.isin(obs.pix[det], other_arr)
                d['spike'] |= lone

            # ── Pass 2: accumulate ─────────────────────────────────────────────
            for det, d in det_cache.items():
                pix    = obs.pix[det]
                valid  = pix >= 0
                # E1 — per-snip polynomial destripe (optional).
                # Fit a robust polynomial to "clean" background samples and
                # subtract it to remove thermal drifts and inter-detector DC
                # offsets that survive the IPACCAL vfet / GXN correction.
                # Sources and flagged samples are excluded from fitting but
                # still receive the baseline subtraction.
                if DESTRIPE:
                    _bad = (obs.corrupt[det] | obs.flash_trim
                        | d['spike'] | d['tail'] | d['source']
                        | ~np.isfinite(d['mjy_sr']))
                    d['mjy_sr'], _ = ip.poly_destripe(
                    d['mjy_sr'], _bad, degree=DESTRIPE_DEGREE)

                planet = (np.isin(pix, planet_pix)
                          if planet_pix.size else np.zeros(NOBS, dtype=bool))
                good   = (valid & np.isfinite(d['mjy_sr'])
                          & ~obs.corrupt[det] & ~obs.flash_trim
                          & ~d['spike'] & ~d['tail']
                          & ~planet)
                if not MASK:
                    good = valid & np.isfinite(d['mjy_sr'])

                pix_v = pix[good]
                val_v = d['mjy_sr'][good]

                if not np.all(np.isfinite(val_v)):
                    bad_idx  = np.where(~np.isfinite(val_v))[0]
                    bad_vals = val_v[bad_idx]
                    bad_files.append(
                        f'{obs_label}  det={det}  '
                        f'{len(bad_idx)} non-finite sample(s): '
                        f'{bad_vals[:5].tolist()}{" ..." if len(bad_idx) > 5 else ""}')
                    keep  = np.isfinite(val_v)
                    pix_v = pix_v[keep]
                    val_v = val_v[keep]

                np.add.at(sum_map,  pix_v, val_v)
                np.add.at(hits_map, pix_v, 1)
                n_matched += int(len(pix_v))

        except Exception as exc:
            import traceback
            print(f'  [WARN] seg={seg_str} {obs_label}: '
                  f'{type(exc).__name__}: {exc}', file=sys.stderr)
            if os.environ.get('TRACEBACK'):
                traceback.print_exc()

    stats = dict(n_obs=n_obs, bad_files=bad_files)
    return n_matched, stats


# ── Multiprocessing worker ─────────────────────────────────────────────────────
def _seg_worker(args):
    """
    Top-level worker for multiprocessing.Pool.  Defined at module level so
    it can be pickled by worker processes.

    args = (seg_str, npix, cal)
    """
    seg_str, npix, cal = args
    local_sum  = np.zeros(npix, dtype=np.float64)
    local_hits = np.zeros(npix, dtype=np.int32)
    n, stats   = process_seg(seg_str, local_sum, local_hits, cal)
    return n, stats, local_sum, local_hits


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    print(f'IRAS {_bcfg["name"]} — all {len(BAND_DETS)} detectors — calibrated HEALPix map')
    print(f'NSIDE={NSIDE}  ({hp.nside2resol(NSIDE, arcmin=True):.1f} arcmin pixels)')
    print()

    # Preload calibration tables once for the entire run
    print('Loading calibration tables...')
    cal = tod.load_cal_tables()
    n_srhf_sops = len(cal.all_srhfs)
    n_srhf_obs  = sum(len(v) for v in cal.all_srhfs.values())
    n_prhf_sops = len(cal.all_prhfs)
    n_prhf_obs  = sum(len(v) for v in cal.all_prhfs.values())
    print(f'  ctype2 : a2dc = {cal.a2dc:.4e} V/DN·gain')
    print(f'  vfet   : {len(cal.utc_tbl)} epochs,  UTC {cal.utc_tbl[0]:.0f} – {cal.utc_tbl[-1]:.0f}')
    print(f'  bbtimes: {len(cal.bbt)} events')
    print(f'  UTCS   : {len(cal.utcs_refs)} SOP entries')
    n_ct = sum(len(v) for v in cal.corrupt_table.values())
    print(f'  corrupt: {n_ct} windows across {len(cal.corrupt_table)} (SOP,OBS,band) triplets')
    print(f'  SRHF   : {n_srhf_obs} observations across {n_srhf_sops} SOPs')
    print(f'  PRHF   : {n_prhf_obs} observations across {n_prhf_sops} SOPs')
    print()

    # Find all path segments with BPHF + raw data for this band
    bphf_segs = sorted(set(
        m.group(1)
        for f in glob.glob(BPHF_GLOB.format(seg='*'))
        for m in [re.search(r'sop\.(\d+)_', f)] if m))
    raw_segs = sorted(set(
        m.group(1)
        for g in tod._raw_globs_for_band(BAND, tod.ROOT, bool(INCLUDE_AO))
        for f in glob.glob(g.format(seg='*'))
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

    worker_args = [(seg, npix, cal) for seg in segs]

    tot_obs   = 0
    all_bad   = []

    if nproc == 1:
        for args in tqdm(worker_args, desc='segments'):
            n, st, s, h = _seg_worker(args)
            total     += n
            sum_map   += s
            hits_map  += h
            tot_obs   += st['n_obs']
            all_bad   += st['bad_files']
    else:
        with multiprocessing.Pool(processes=nproc) as pool:
            for n, st, s, h in tqdm(
                    pool.imap_unordered(_seg_worker, worker_args),
                    total=len(worker_args), desc='segments'):
                total     += n
                sum_map   += s
                hits_map  += h
                tot_obs   += st['n_obs']
                all_bad   += st['bad_files']

    print()
    print(f'Observation summary:')
    print(f'  SOP+OBS processed (≥1 valid detector) : {tot_obs:,}')
    print(f'  Total matched+valid samples            : {total:,}')

    hit_pix = np.count_nonzero(hits_map)
    print()
    print(f'Pixels hit : {hit_pix:,}  ({hit_pix / npix * 100:.1f}% of sky at nside={NSIDE})')
    if hit_pix > 0:
        print(f'Max hits per pixel    : {hits_map.max()}')
        print(f'Median hits (non-zero): {np.median(hits_map[hits_map > 0]):.0f}')

    if all_bad:
        print()
        print(f'ERROR: {len(all_bad)} file/detector combination(s) produced '
              f'non-finite calibrated values (stripped before accumulation):')
        for entry in all_bad:
            print(f'  {entry}')
        raise RuntimeError(
            f'{len(all_bad)} non-finite calibration result(s) — see list above.')

    # Build signal map
    signal_map = np.full(npix, hp.UNSEEN, dtype=np.float64)
    mask_arr   = hits_map > 0
    signal_map[mask_arr] = sum_map[mask_arr] / hits_map[mask_arr]
    bad = mask_arr & ~np.isfinite(signal_map)
    if bad.any():
        raise RuntimeError(
            f'{bad.sum()} pixel(s) have non-finite signal after co-add — '
            f'this is a bug: {np.where(bad)[0][:20].tolist()}')

    # Save FITS
    btag     = _bcfg['survey_tag']
    mask_sfx = '' if MASK else '_nomask'
    ao_sfx   = '_withao' if INCLUDE_AO else ''
    out_sig  = os.path.join(OUT_DIR,
            f'IRAS_{btag}_alldet_cal_map_n{NSIDE:04}{mask_sfx}{ao_sfx}.fits')
    out_hits = os.path.join(OUT_DIR,
            f'IRAS_{btag}_alldet_cal_hits_n{NSIDE:04}{mask_sfx}{ao_sfx}.fits')
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
        vals   = signal_map[mask_arr]
        finite = vals[np.isfinite(vals)]
        vmin, vmax = np.percentile(finite, 1), np.percentile(finite, 99)
        fig, axes = plt.subplots(1, 2, figsize=(16, 5))
        fig.suptitle(
            f'IRAS {_bcfg["name"]} — all {len(BAND_DETS)} detectors — '
            f'calibrated map (nside={NSIDE})', fontsize=13)
        plt.sca(axes[0])
        hp.mollview(signal_map, title='Signal (MJy/sr)', unit='MJy/sr',
                    min=vmin, max=vmax, coord='G', cmap='inferno', hold=True)
        hp.graticule(dpar=30, dmer=60, alpha=0.4)
        plt.sca(axes[1])
        hp.mollview(hits_map.astype(float), title='Hits per pixel',
                    coord='G', cmap='Blues', hold=True,
                    min=0, max=np.percentile(hits_map[mask_arr], 99))
        hp.graticule(dpar=30, dmer=60, alpha=0.4)
        plt.tight_layout()
        out_png = os.path.join(OUT_DIR,
                f'IRAS_{btag}_map_n{NSIDE:04}{mask_sfx}{ao_sfx}.png')
        plt.savefig(out_png, dpi=150)
        print(f'Saved: {out_png}')


if __name__ == '__main__':
    main()
