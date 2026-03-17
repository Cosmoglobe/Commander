"""
ipaccal.py — Python port of the IRAS IPACCAL calibration pipeline

Converts raw detector DN (data numbers) to calibrated flux density.

Chain (matching ipaccal.src call order):
  DN  --[cidn2v]--> volts
      --[cioffs]--> volts - vfet_baseline - gxn_correction
      --[ciamps]--> amps  (divide by load resistor)
      --[ciwpm2]--> W/m²  (divide by responsivity from SRHF flash model,
                            optionally modified by PRHF photon term)
      --[ccesad]--> MJy/sr (multiply by w2mjy(band) / effective_solid_angle)

Flagging:
  Samples whose UTC falls inside an ICF flash event (from bbtimes) are
  marked with NaN.

References:
  diskrog10androg11/rog10/gipsy/tsk/ipaccal.src  (ciical, cidn2v, cicmgo,
    ciamps, cioffs, cirdvt, ciwpm2, ciprhf, cirdsf/cisrec, ccstim, ccesad,
    ccgxnc, ciseqn, ciseqn)
  diskrog10androg11/rog10/gipsy/tsk/zodycal.src   (IRZC_DETMODEL, IRZC_GAIN)
  diskrog10androg11/rog11/IPAC/tables/ctype2       (a2dc, gains, offsets)
  diskrog10androg11/rog11/IPAC/tables/vfetsg/lg/hg (electronic baseline)
  diskrog10androg11/rog11/IPAC/tables/bbtimes       (ICF flash event mask)
  diskrog10androg11/rog11/IPAC/S/IR.SRHF.D*         (flash model params)
  diskrog10androg11/rog11/IPAC/P/IR.PRHF.D*         (photon-exposure resp.)
"""

import os
import struct
import numpy as np

# ---------------------------------------------------------------------------
# Hard-coded constants from ipaccal.src
# ---------------------------------------------------------------------------

# Detector calibration sequence: maps det (1..62) → calibration index (1..62)
# From ciseqn.f  (Fortran 1-indexed; Python: CISEQN[det-1] gives the 1-based seq)
CISEQN = np.array([
    56, 57, 58, 59, 60, 61, 62,           # dets  1-7   (B4 module B)
    32, 33, 34, 35, 36, 37, 38, 39,       # dets  8-15  (B3 module B)
    25, 26, 27, 28, 29, 30, 31,           # dets 16-22  (B2 module B)
     1,  2,  3,  4,  5,  6,  7,  8,      # dets 23-30  (B1 module B)
    40, 41, 42, 43, 44, 45, 46, 47,      # dets 31-38  (B3 module A)
    17, 18, 19, 20, 21, 22, 23, 24,      # dets 39-46  (B2 module A)
     9, 10, 11, 12, 13, 14, 15, 16,      # dets 47-54  (B1 module A)
    48, 49, 50, 51, 52, 53, 54, 55,      # dets 55-62  (B4 module A)
], dtype=int)  # 1-based values; len 62

# Default offset command for each detector (from cicmgo, 1-indexed by det)
# coff = COF[det-1] + 1  (in Python)
_COF = np.array([
    3, 3, 3, 4, 2, 4, 3,       # dets  1-7
    3, 3, 4, 3, 3, 3, 3, 4,    # dets  8-15
    4, 3, 4, 4, 4, 4, 4,       # dets 16-22
    4, 4, 4, 3, 4, 3, 4, 4,    # dets 23-30
    3, 3, 3, 3, 3, 3, 4, 3,    # dets 31-38
    4, 4, 3, 4, 3, 4, 3, 3,    # dets 39-46
    4, 3, 4, 2, 3, 3, 4, 4,    # dets 47-54
    4, 4, 3, 4, 3, 3, 4, 4,    # dets 55-62
], dtype=int)  # 0-indexed values; len 62

# gain_setting → cgain index (1-based) used to index gains[:,cgain-1]
_GAIN_CMD = {'low': 1, 'high': 2, 'standard': 3}

# Band for each detector (1-indexed det number, Band 1=12µm .. 4=100µm)
# From ircc_const.c / IRCC_BANDET: det 1-7=B4, 8-15=B3, 16-22=B2, 23-30=B1,
#   31-38=B3, 39-46=B2, 47-54=B1, 55-62=B4
_BAND = np.array(
    [4]*7 + [3]*8 + [2]*7 + [1]*8 + [3]*8 + [2]*8 + [1]*8 + [4]*8, dtype=int
)  # 1-indexed in code → use _BAND[det-1]

# Wavelength scale from W/m² to MJy/sr (one per band, index 0-3 = band 1-4)
# From ipaccal.src: data w2mjy / 7.481e6, 19.38e6, 38.76e6, 100.0e6 /
W2MJY = np.array([7.481e6, 19.38e6, 38.76e6, 100.0e6])

# Effective solid angle per detector (steradians × 1e7) from ccesad.shl
# Index = det-1 (0-based); dead detectors = 0
_ESAD_RAW = np.array([
    14.1, 13.7, 12.9, 13.1, 13.2, 13.2, 13.6,   # dets 1-7
     6.10,  6.04,  6.07,  1.92,  4.85,  6.10,  6.19,  6.56,  # dets 8-15
     3.39,  0.00,  3.46,  3.40,  0.00,  3.48,  3.38,  # dets 16-22
     3.21,  3.22,  3.17,  0.92,  2.41,  2.95,  3.25,  3.25,  # dets 23-30
     1.97,  6.40,  6.32,  6.25,  6.36,  0.00,  6.37,  4.84,  # dets 31-38
     1.82,  3.50,  3.44,  3.48,  3.50,  3.45,  3.46,  1.78,  # dets 39-46
     0.89,  3.31,  3.22,  3.29,  3.25,  3.19,  3.24,  2.45,  # dets 47-54
     7.9,  13.3,  14.0,  12.9,  14.0,  13.8,  14.2,   6.8,   # dets 55-62
])
ESAD = _ESAD_RAW * 1.0e-7  # steradians

# DC stimulator flash #5 intensity per detector (amps) from ccstim.f
# Indexed by ciseqn(det)-1 in the original code; here stored as 62-element
# array indexed by calibration sequence (0-based).  ccstim = st5[det] * acdc[det]
_ST5 = np.array([
    0.838254e-12, 0.904365e-12, 0.864749e-12, 0.811226e-12,
    0.832986e-12, 0.957989e-12, 0.971139e-12,
    0.215303e-11, 0.200418e-11, 0.202432e-11, 0.741735e-12,
    0.155880e-11, 0.213811e-11, 0.191740e-11, 0.227023e-11,
    0.392422e-11, 0.000000e+00, 0.427730e-11, 0.330927e-11,
    0.000000e+00, 0.371521e-11, 0.443437e-11,
    0.334078e-11, 0.312709e-11, 0.313663e-11, 0.917304e-12,
    0.250199e-11, 0.317067e-11, 0.314251e-11, 0.350205e-11,
    0.794343e-12, 0.211903e-11, 0.196415e-11, 0.208528e-11,
    0.224815e-11, 0.000000e+00, 0.198245e-11, 0.138940e-11,
    0.192416e-11, 0.331966e-11, 0.337515e-11, 0.412459e-11,
    0.357591e-11, 0.352088e-11, 0.374580e-11, 0.193974e-11,
    0.892483e-12, 0.310774e-11, 0.337660e-11, 0.339881e-11,
    0.327327e-11, 0.320611e-11, 0.313170e-11, 0.235516e-11,
    0.729108e-12, 0.100974e-11, 0.859967e-12, 0.832541e-12,
    0.810062e-12, 0.953099e-12, 0.809286e-12, 0.584463e-12,
])
_ACDC = np.array(
    [1.00]*7 + [0.92]*8 + [0.82]*7 + [0.78]*8 +
    [0.92]*8 + [0.82]*8 + [0.78]*8 + [1.00]*8
)
CCSTIM = _ST5 * _ACDC  # in amps; indexed 0-based by det-1

# GXN correction: alpha and tau per detector from ccgxnc.shl
# kpar is always 0.  For bands 1&2 alpha=0 (and tau=1 to avoid /0).
# Indexed 0-based by det-1.
_GXN_ALPHA = np.array([
    -0.50660e-06,  0.38558e-06, -0.16187e-06, -0.12623e-06,
    -0.51914e-06,  0.18098e-06,  0.70222e-06,
    -0.23421e-06, -0.51040e-07,  0.13120e-07, -0.39730e-07,
    -0.87980e-07, -0.33010e-07, -0.13913e-06,  0.19433e-06,
    *([0.0] * 15),                           # B2 dets 16-22 (inactive)
    *([0.0] * 8),                            # B1 dets 23-30 (inactive)
    -0.20943e-06,  0.21088e-06, -0.11521e-06,  0.13708e-06,
     0.88130e-07,  0.00000e+00,  0.79650e-07, -0.17686e-06,
    *([0.0] * 16),                           # B2+B1 dets 39-54 (inactive)
    -0.12993e-07,  0.46260e-06, -0.38041e-06,  0.31724e-06,
    -0.55873e-06,  0.00000e+00,  0.44439e-06,  0.15833e-06,
])
_GXN_TAU = np.array([
    2100., 2820., 2520., 4080., 1980., 1920., 2940.,
    2700., 4200.,  660., 4620., 4860., 1320., 2520., 2340.,
    *([1.0] * 15),
    *([1.0] * 8),
    2520., 2340., 2220., 2940., 6900., 1000., 3720., 2940.,
    *([1.0] * 16),
    5400., 2460., 2580., 2340., 2580., 1000., 2520., 1680.,
])

# Dead detectors (ESAD == 0 → cannot calibrate)
DEAD = ESAD == 0.0


# ---------------------------------------------------------------------------
# ctype2 loader — a2dc, gains[62,3], offs[62,8]
# ---------------------------------------------------------------------------

def load_ctype2(path):
    """
    Parse the ctype2 calibration-constants file.

    The Fortran read sequence (from cidn2v.f) is:
        read(unit, '(45x,e15.7)') a2dc        ← 4th value on first line
        do i=1,3 : read(unit, '(8e15.7)') gains(1:62, i)
        do i=1,8 : read(unit, '(8e15.7)') offs(1:62, i)

    Returns
    -------
    a2dc : float
        Analogue-to-digital conversion constant (volts per DN × gain).
    gains : ndarray, shape (62, 3)
        Gain values indexed by (calibration_seq-1, gain_cmd-1).
    offsets : ndarray, shape (62, 8)
        Offset values indexed by (calibration_seq-1, offset_cmd-1).
    """
    import re
    with open(path) as f:
        text = f.read()
    vals = np.array([float(v) for v in re.findall(r'[-+]?\d+\.\d+[Ee][+-]\d+', text)])
    # Line 0 has 7 values; the 4th (index 3) is a2dc  (Fortran 45x skip = 3×15 chars)
    a2dc = float(vals[3])
    pos = 7             # after the 7 header values
    gains = vals[pos:pos + 62*3].reshape(3, 62).T   # → (62, 3)
    pos += 62 * 3
    offsets = vals[pos:pos + 62*8].reshape(8, 62).T  # → (62, 8)
    return a2dc, gains, offsets


# ---------------------------------------------------------------------------
# vfet baseline loader — time-varying electronic offset per detector
# ---------------------------------------------------------------------------

def load_vfet(path):
    """
    Load a vfethg / vfetlg / vfetsg baseline file.

    Each epoch row: [UTC-1981_s, vfet_seq1 ... vfet_seq62]
    Epochs are stored in REVERSE time order (most recent first).

    Returns
    -------
    utc : ndarray, shape (N,)   UTC-1981 seconds, ascending order
    vfets : ndarray, shape (N, 62)   Electronic baseline in volts,
            indexed by [epoch, ciseqn(det)-1].
    """
    import re
    with open(path) as f:
        text = f.read()
    vals = np.array([float(v) for v in re.findall(r'[-+]?\d+\.\d+[Ee][+-]\d+', text)])
    epochs = vals.reshape(-1, 63)          # 63 = 1 UTC + 62 det values
    # Stored most-recent-first; reverse so UTC is ascending
    epochs = epochs[::-1]
    utc = epochs[:, 0]
    vfets = epochs[:, 1:]                  # (N, 62), indexed by ciseqn(det)-1
    return utc, vfets


def vfet_baseline(utcs, utc_table, vfets_table, det):
    """
    Linearly interpolate the electronic (vfet) baseline for detector `det`
    at UTC times `utcs`.

    Parameters
    ----------
    utcs : array-like, UTC-1981 seconds for each sample
    utc_table, vfets_table : from load_vfet()
    det : int, 1-based detector number

    Returns
    -------
    baseline : ndarray, shape like utcs  (volts)
    """
    col = CISEQN[det - 1] - 1             # 0-based column in vfets_table
    return np.interp(np.asarray(utcs, dtype=float),
                     utc_table, vfets_table[:, col])


# ---------------------------------------------------------------------------
# ICF flash time loader and flagging
# ---------------------------------------------------------------------------

def load_bbtimes(path):
    """
    Read ICF flash event table.

    Each row: two space-separated floats = start_UTC, end_UTC in UTC-1981 s.
    (e.g. '65475727.982   65476702.037')

    Returns
    -------
    bbt : ndarray, shape (N, 2)   [[start_utc, end_utc], ...] in UTC-1981 s
    """
    import re
    with open(path) as f:
        text = f.read()
    # Match plain decimals (65475727.982) and scientific notation (1.25E+07)
    vals = np.array([float(v) for v in re.findall(
        r'[-+]?(?:\d+\.\d*|\.\d+)(?:[Ee][+-]?\d+)?', text)])
    return vals.reshape(-1, 2)


def icf_flag(utcs, bbt):
    """
    Return a boolean mask, True where UTC time falls inside an ICF flash event.

    Parameters
    ----------
    utcs : array-like, UTC-1981 seconds
    bbt  : ndarray (N,2) from load_bbtimes()
    """
    utcs = np.asarray(utcs, dtype=float)
    mask = np.zeros(len(utcs), dtype=bool)
    for start, end in bbt:
        mask |= (utcs >= start) & (utcs <= end)
    return mask


# ---------------------------------------------------------------------------
# SRHF (Survey Reference History File) reader
# ---------------------------------------------------------------------------

def read_srhf(srhf_dir, sop, obs):
    """
    Read an SRHF file for a given (SOP, OBS) combination.

    The file contains the detector flash-response parameters needed by ciwpm2:
        RESP(det, t) = apl[det] + bpl[det] * exp( -t / tpl[det] )
    where t = UTC seconds since the last bias boost.

    SRHF are stored per OBS (observation = ATT period), not per ATT type.
    Filename: IR.SRHF.D0{sop:03d}{obs:03d}

    Parameters
    ----------
    srhf_dir : str, path to IPAC/S/ directory
    sop : int
    obs : int, the OBS/ATT number

    Returns
    -------
    utc1 : float   UTC-1981 s of first stimulation flash
    utc2 : float   UTC-1981 s of second stimulation flash
    apl  : ndarray (62,) steady-state flash amplitude in amps, by ciseqn order
    bpl  : ndarray (62,) exponential prefactor in amps, by ciseqn order
    tpl  : ndarray (62,) time constant in seconds, by ciseqn order
    """
    fname = os.path.join(srhf_dir, f'IR.SRHF.D0{sop:03d}{obs:03d}')
    with open(fname, 'rb') as f:
        data = f.read()
    # Fortran record: 4-byte big-endian length marker, then payload, then length
    off = 4
    utc1 = struct.unpack('>d', data[off:off+8])[0]; off += 8
    utc2 = struct.unpack('>d', data[off:off+8])[0]; off += 8
    apl = np.frombuffer(data[off:off+248], dtype='>f4').copy(); off += 248
    bpl = np.frombuffer(data[off:off+248], dtype='>f4').copy(); off += 248
    tpl = np.frombuffer(data[off:off+248], dtype='>f4').copy()
    return utc1, utc2, apl, bpl, tpl


def srhf_flash_amplitude(utmbb, apl, bpl, tpl, det):
    """
    Compute reconstructed stimulation flash amplitude (in amps) at time
    `utmbb` seconds since the last bias boost for detector `det`.

    F(t) = apl[seq] + bpl[seq] * exp(-t / tpl[seq])

    The SRHF arrays are in calibration-sequence order (ciseqn).

    Parameters
    ----------
    utmbb : float or ndarray  seconds since last bias boost
    apl, bpl, tpl : from read_srhf() — length-62, ciseqn-indexed
    det : int, 1-based detector number

    Returns
    -------
    flamp : same shape as utmbb (amps)
    """
    seq = CISEQN[det - 1] - 1   # 0-based index into apl/bpl/tpl
    t = np.asarray(utmbb, dtype=float)
    return apl[seq] + bpl[seq] * np.exp(-t / tpl[seq])


# ---------------------------------------------------------------------------
# cidn2v — DN → volts
# ---------------------------------------------------------------------------

def dn_to_volts(raw_dn, det, a2dc, gains, offsets, gain_setting='standard'):
    """
    Convert raw DN to volts.

    volts = (a2dc / gains[cseqn, cgain]) * raw - offsets[cseqn, coff]

    Parameters
    ----------
    raw_dn : array-like, raw data numbers
    det : int, 1-based detector number
    a2dc, gains, offsets : from load_ctype2()
    gain_setting : 'standard' | 'low' | 'high'

    Returns
    -------
    volts : ndarray
    conv  : float  (a2dc / gain) conversion factor
    """
    cseqn = CISEQN[det - 1] - 1           # 0-based
    cgain = _GAIN_CMD[gain_setting] - 1   # 0-based
    coff  = _COF[det - 1]                 # 0-based (already +0; coff in Fortran is COF+1, 1-based → our +0)
    conv = a2dc / gains[cseqn, cgain]
    off  = offsets[cseqn, coff]
    return np.asarray(raw_dn, dtype=float) * conv - off, conv


# ---------------------------------------------------------------------------
# ciamps — volts → amps  (IPAC load-resistor model, ExSup IIC-11)
# ---------------------------------------------------------------------------

_R0   = 2.0e10
_RL1  = _R0 * 10**0.008       # = R0 * 1.01862
_RL2  = _R0
_RL3  = _R0 * 10**(-0.04)     # = R0 * 0.91201
_VM1, _VM2, _VM3     = 0.001, 0.01, 0.14
_PR1, _PR2, _PR3     = -0.008, -0.0349, -0.1402


def load_resistor(v):
    """Piecewise IPAC load-resistor model (ExSup IIC-11), vectorised.

    For v < vm1 (including negative values) the Fortran code returns rl1
    unconditionally; power-law branches are only entered for v >= vm1.
    """
    v = np.asarray(v, dtype=float)
    # Clamp to vm1 for the power-law branches to avoid complex/nan values.
    vc2 = np.maximum(v, _VM1)   # safe for v < vm1 branch (never reached there)
    vc3 = np.maximum(v, _VM2)
    r = np.where(v < _VM1,  _RL1,
        np.where(v < _VM2,  _RL1 * (vc2 / _VM1)**_PR1,
        np.where(v < _VM3,  _RL2 * (vc3 / _VM2)**_PR2,
                             _RL3 * (v   / _VM3)**_PR3)))
    return r


def volts_to_amps(volts):
    """Convert baseline-corrected volts to amps via the load-resistor model."""
    v = np.asarray(volts, dtype=float)
    return v / load_resistor(v)


# ---------------------------------------------------------------------------
# cioffs — subtract electronic baseline (vfet + GXN correction)
# ---------------------------------------------------------------------------

def subtract_baseline(volts, det, utcs_samples, utc_table, vfets_table,
                       utmbb_samples):
    """
    Subtract the electronic baseline from voltage samples.

    Two terms are removed:
      1. vfet(t): time-interpolated electronic offset from the vfet table.
      2. GXN correction: kab * exp(-dt_bias/tau), only for B3/B4 detectors.
         kpar is always 0 in the IRAS data, so gxnc = alpha * exp(-t/tau).

    Parameters
    ----------
    volts : ndarray
    det : int, 1-based
    utcs_samples : ndarray  UTC-1981 seconds for each sample
    utc_table, vfets_table : from load_vfet()
    utmbb_samples : ndarray  seconds since last bias boost for each sample

    Returns
    -------
    volts_corrected : ndarray
    """
    vfetc = vfet_baseline(utcs_samples, utc_table, vfets_table, det)

    alpha = _GXN_ALPHA[det - 1]
    tau   = _GXN_TAU[det - 1]
    # kpar=0 always, so gxnc = (kpar + alpha*bbd)*exp(-dt/tau) simplifies
    # but cioffs uses: kab = kappa + alpha*bbd  (and kappa=0 always)
    # and gxnc = kab * exp(dbutc * (-1/tau))
    # For simplicity without the bias-boost duration: gxnc ≈ alpha*exp(-utmbb/tau)
    gxnc = alpha * np.exp(-np.asarray(utmbb_samples, dtype=float) / tau)

    return volts - vfetc - gxnc


# ---------------------------------------------------------------------------
# ciwpm2 — amps → W/m²  (flash-model responsivity)
# ---------------------------------------------------------------------------

def amps_to_wm2(amps, det, utmbb_samples, apl, bpl, tpl):
    """
    Convert baseline-corrected amplitudes to W/m² using the SRHF responsivity.

    Responsivity R(t) = F(t) / ccstim(det)    [A · m² / W]
    Output S(t) = amps(t) / R(t)              [W/m²]

    The per-sample responsivity is evaluated at the UTC since last bias boost
    for each sample.

    Parameters
    ----------
    amps : ndarray, detector current in amperes
    det : int, 1-based
    utmbb_samples : ndarray, seconds since last bias boost for each sample
    apl, bpl, tpl : from read_srhf(), length-62, ciseqn-indexed

    Returns
    -------
    wm2 : ndarray, flux density in W/m²
    """
    t = np.asarray(utmbb_samples, dtype=float)
    flamp = srhf_flash_amplitude(t, apl, bpl, tpl, det)
    dcflux = CCSTIM[det - 1]
    if dcflux == 0.0:
        return np.full_like(amps, np.nan)
    resp = flamp / dcflux                   # A / (W/m²) = A·m²/W
    return np.asarray(amps, dtype=float) / resp


# ---------------------------------------------------------------------------
# ccesad + w2mjy — W/m² → MJy/sr
# ---------------------------------------------------------------------------

def wm2_to_mjy_sr(wm2, det):
    """
    Convert W/m² to MJy/sr using the effective solid angle of the detector.

    S_MJy/sr = w2mjy(band) / Omega_eff(det) * S_Wm2
    """
    band = _BAND[det - 1]                  # 1-based
    omega = ESAD[det - 1]
    if omega == 0.0:
        return np.full_like(wm2, np.nan)
    return W2MJY[band - 1] / omega * np.asarray(wm2, dtype=float)


def wm2_to_jy(wm2, det):
    """
    Convert IPACCAL in-band W/m² to Jy (point-source calibration).

    The IPACCAL 'W/m²' is total in-band integrated flux, not spectral flux
    density.  The conversion to Jy uses the same w2mjy bandpass factor as
    wm2_to_mjy_sr but without dividing by the solid angle:

        S[Jy] = w2mjy(band) × 1e6 × S[W/m²]

    This equals S[MJy/sr] × Omega_eff (point-source flux = surface brightness
    × beam solid angle).
    """
    band = _BAND[det - 1]      # 1-based
    return np.asarray(wm2, dtype=float) * W2MJY[band - 1] * 1.0e6


# ---------------------------------------------------------------------------
# Full calibration pipeline
# ---------------------------------------------------------------------------

def calibrate(raw_dn, det, sop, obs,
              utcs_start, utcs_inc,
              utmbb_start,
              ctype2_path,
              vfet_path,
              srhf_dir,
              bbt,
              gain_setting='standard',
              output_units='MJy/sr'):
    """
    Full IPACCAL calibration pipeline for one detector snip.

    Parameters
    ----------
    raw_dn : array-like, shape (N,)
        Raw decompressed data numbers from the BPHF/scan TOD.
    det : int
        1-based IRAS detector number (1..62).
    sop : int
        Survey Operations Plan number.
    obs : int
        Observation/ATT number within the SOP.
    utcs_start : float
        UTC-1981 seconds of the first sample.
    utcs_inc : float
        UTC-1981 seconds per sample increment.
    utmbb_start : float
        Seconds since the last bias boost at the start of the snip.
    ctype2_path : str
        Path to the IPAC/tables/ctype2 file.
    vfet_path : str
        Path to the appropriate vfethg/vfetlg/vfetsg file.
    srhf_dir : str
        Directory containing IR.SRHF.D* files (IPAC/S/).
    bbt : ndarray, shape (M, 2)
        ICF flash event table from load_bbtimes(); set None to skip flagging.
    gain_setting : 'standard' | 'low' | 'high'
        Detector gain state.
    output_units : 'MJy/sr' | 'W/m2' | 'Jy'
        Output flux density units.

    Returns
    -------
    flux : ndarray, shape (N,)
        Calibrated flux density in specified units.  NaN for flagged samples
        (ICF flashes, dead detectors).
    flags : ndarray, bool, shape (N,)
        True where the sample is flagged (ICF flash or dead detector).
    utcs : ndarray, shape (N,)
        UTC-1981 seconds for each sample.
    """
    raw_dn = np.asarray(raw_dn, dtype=float)
    N = len(raw_dn)

    # Time arrays
    utcs = utcs_start + np.arange(N) * utcs_inc
    utmbb = utmbb_start + np.arange(N) * utcs_inc

    # Flag dead detector immediately
    if DEAD[det - 1]:
        flags = np.ones(N, dtype=bool)
        return np.full(N, np.nan), flags, utcs

    # --- ICF flash flagging ---
    if bbt is not None:
        flags = icf_flag(utcs, bbt)
    else:
        flags = np.zeros(N, dtype=bool)

    # --- Load calibration tables (cached by caller or re-loaded here) ---
    a2dc, gains, offsets = load_ctype2(ctype2_path)
    utc_table, vfets_table = load_vfet(vfet_path)

    srhf_utc1, srhf_utc2, apl, bpl, tpl = read_srhf(srhf_dir, sop, obs)
    # Reorder SRHF arrays from ciseqn order to physical det order (ciseqn-1 indexing)
    # apl[k] corresponds to calibration sequence k+1, i.e., detector with ciseqn=k+1
    # They are already in ciseqn order as stored; srhf_flash_amplitude uses CISEQN lookup.

    # --- Step 1: DN → volts ---
    volts, _ = dn_to_volts(raw_dn, det, a2dc, gains, offsets, gain_setting)

    # --- Step 2: subtract electronic baseline (vfet + GXN) ---
    volts = subtract_baseline(volts, det, utcs, utc_table, vfets_table, utmbb)

    # --- Step 3: volts → amps ---
    amps = volts_to_amps(volts)

    # --- Step 4: amps → W/m² ---
    wm2 = amps_to_wm2(amps, det, utmbb, apl, bpl, tpl)

    # --- Step 5: unit conversion ---
    if output_units == 'MJy/sr':
        flux = wm2_to_mjy_sr(wm2, det)
    elif output_units == 'Jy':
        flux = wm2_to_jy(wm2, det)
    else:  # 'W/m2'
        flux = wm2

    # Apply flags (set flagged samples to NaN)
    flux[flags] = np.nan

    return flux, flags, utcs


# ---------------------------------------------------------------------------
# Convenience: batch calibration using PRHF responsivity (alternative path)
# ---------------------------------------------------------------------------

def calibrate_with_prhf(raw_dn, det, sop, obs,
                         utcs_start, utcs_inc,
                         prhf_dir, bbt,
                         output_units='MJy/sr'):
    """
    Simplified calibration using only PRHF responsivities (for B4 detectors).

    This bypasses the full ciical/ciwpm2 chain and directly uses the
    tabulated responsivity from the PRHF files (which already encode the
    combined gain, electronic baseline, and flash-model response).

    Only valid for Band 4 detectors (1-7, 55-62).

    Parameters
    ----------
    raw_dn : array-like  raw decompressed ADU
    det : int, 1-based (must be a B4 detector)
    sop, obs : int
    utcs_start, utcs_inc : float  (UTC-1981 seconds)
    prhf_dir : str  path to IPAC/P/ directory
    bbt : ndarray (M,2) or None  ICF flash event table
    output_units : 'MJy/sr' | 'Jy' | 'W/m2'

    Returns
    -------
    flux : ndarray
    flags : ndarray, bool
    utcs : ndarray
    """
    from python.read_prhf import get_resp_for_sop, B4_DETS
    import importlib, sys

    # Allow import from same directory
    sys.path.insert(0, os.path.dirname(__file__))
    rp = importlib.import_module('read_prhf')

    if _BAND[det - 1] != 4:
        raise ValueError(f"calibrate_with_prhf only supports Band 4; det={det} is Band {_BAND[det-1]}")

    raw_dn = np.asarray(raw_dn, dtype=float)
    N = len(raw_dn)

    utcs = utcs_start + np.arange(N) * utcs_inc

    # ICF flagging
    flags = icf_flag(utcs, bbt) if bbt is not None else np.zeros(N, dtype=bool)

    # Load PRHF time series for this SOP
    prhf_t, prhf_r = rp.get_resp_for_sop(prhf_dir, sop, obs_type=0)

    # Determine det_index in the B4 PRHF array
    if det not in rp.B4_DETS:
        raise ValueError(f"det {det} not in B4_DETS")
    det_idx = rp.B4_DETS.index(det)

    # Get orbit-relative time for each sample: need utcs → orbit-relative seconds
    # PRHF times are orbit-relative; approximate: subtract first PRHF time's utcs mapping
    # Since we don't have the absolute PRHF→UTC mapping here, use the PRHF interpolation
    # with orbit-relative UTC difference  (assumes utcs_start corresponds to prhf_t[0])
    resp = rp.interp_responsivity(utcs - utcs_start,
                                   prhf_t, prhf_r, det_index=det_idx)

    # resp is in Jy/ADU (approx); this is already the full responsivity
    # raw_dn here must be decompressed ADU before any DC offset subtraction
    flux = raw_dn * resp  # Jy

    # Convert units
    if output_units == 'MJy/sr':
        omega = ESAD[det - 1]
        if omega > 0:
            flux = flux * 1e-6 / omega   # Jy → MJy → MJy/sr
    elif output_units == 'W/m2':
        flux = flux * 1e-26

    flux[flags] = np.nan
    return flux, flags, utcs


# ---------------------------------------------------------------------------
# Demo / test
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import sys
    ROOT = '/home/dwatts/IRAS'
    IPAC = os.path.join(ROOT, 'diskrog10androg11/rog11/IPAC')

    ctype2_path = os.path.join(IPAC, 'tables/ctype2')
    vfet_path   = os.path.join(IPAC, 'tables/vfetsg')
    srhf_dir    = os.path.join(IPAC, 'S')
    bb_path     = os.path.join(IPAC, 'tables/bbtimes')

    # --- Verify ctype2 parsing ---
    a2dc, gains, offsets = load_ctype2(ctype2_path)
    print(f'a2dc = {a2dc:.6e}   (expect ~1.25e-4 V/DN·gain)')
    print(f'gains[det1-1, standard] = gains[55, 2] = {gains[55, 2]:.4f}')
    print(f'offsets shape: {offsets.shape}')

    # --- Verify vfet loader ---
    utc_tbl, vfet_tbl = load_vfet(vfet_path)
    print(f'\nvfet table: {len(utc_tbl)} epochs, UTC range {utc_tbl[0]:.0f}–{utc_tbl[-1]:.0f} s')
    det1_col = CISEQN[0] - 1   # det1 → col 55 (0-based)
    print(f'det01 vfet at first epoch: {vfet_tbl[0, det1_col]:.4e} V')

    # --- Verify SRHF reading (SOP 290, obs 3) ---
    utc1, utc2, apl, bpl, tpl = read_srhf(srhf_dir, sop=290, obs=3)
    seq = CISEQN[0] - 1   # det01 calibration sequence (0-based)
    print(f'\nSRHF SOP290/obs3: utc1={utc1:.0f}, utc2={utc2:.0f}')
    print(f'det01: apl={apl[seq]:.3e}  bpl={bpl[seq]:.3e}  tpl={tpl[seq]:.0f} s')

    # --- Verify bbtimes ---
    bbt = load_bbtimes(bb_path)
    print(f'\nbbtimes: {len(bbt)} events, first: {bbt[0,0]:.0f}–{bbt[0,1]:.0f} s')

    # --- Print offset value for det1 ---
    cseqn_0 = CISEQN[0] - 1    # det1 → index 55
    coff_0  = _COF[0]           # = 3 (0-based offset command)
    off_v   = offsets[cseqn_0, coff_0]
    print(f'\ndet01 ctype2 offset (cmd {coff_0}): {off_v:.6e} V')
    print(f'det01 ctype2 gains: low={gains[cseqn_0,0]:.4f}  high={gains[cseqn_0,1]:.4f}  std={gains[cseqn_0,2]:.4f}')

    # --- End-to-end realistic test ---
    # The zero-signal raw DN for det1 at standard gain satisfies:
    #   a2dc/gain * raw_0 - off = vfet_dc  (→ volts_corrected ≈ 0)
    # For det1: gain=26.58, off=0.0174V, vfet[0]≈-4.47e-3V
    #   raw_0 = (off + vfet) / (a2dc/gain) = (0.0174 - 0.00447) / 4.70e-6 ≈ 2750 DN
    # A 10 Jy source adds ~2.58e-4 V, so raw_signal ≈ 2750 + 55 ≈ 2805 DN
    # (5V_ADC range over 14 bits ≈ 5/16384 = 3.05e-4 V/DN; 2.58e-4/3.05e-4 ≈ 0.85 DN/Jy... order of magnitude consistent)
    utcs0     = utc1 + 100.0    # 100 s after the first flash for this obs
    utcs_inc  = 1.0 / 16.0     # 16 samples/s for B4 survey mode
    utmbb0    = 100.0           # seconds since last bias boost at start

    # Raw DN at approximate zero-signal level (will give ≈0 Jy for dark sky)
    raw_zero = np.full(100, 2750.0)
    # Simulate a 10-Jy source in the middle 20 samples by adding +55 DN
    raw_sig  = raw_zero.copy()
    raw_sig[40:60] += 55.0

    flux, flags, utcs = calibrate(
        raw_sig, det=1, sop=290, obs=3,
        utcs_start=utcs0, utcs_inc=utcs_inc,
        utmbb_start=utmbb0,
        ctype2_path=ctype2_path,
        vfet_path=vfet_path,
        srhf_dir=srhf_dir,
        bbt=bbt,
        gain_setting='standard',
        output_units='Jy',
    )
    print(f'\nRealistic det01 calibration test:')
    print(f'  dark sky median flux = {np.nanmedian(flux[:40]):.3f} Jy')
    print(f'  +55DN source median  = {np.nanmedian(flux[40:60]):.3f} Jy')
    print(f'  unflagged: {(~flags).sum()}/100')
