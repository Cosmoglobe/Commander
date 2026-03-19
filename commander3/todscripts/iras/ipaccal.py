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
# Stored in PHYSICAL DETECTOR ORDER (indexed 0-based by det-1, NOT by ciseqn).
# In the original Fortran ccstim.f the lookup was ccstim(ciseqn(det)), but
# the values here were rearranged into physical order so the Python lookup is
# simply CCSTIM[det-1].  Dead detectors (17, 20, 36) are 0.
# ccstim = st5[det] * acdc[det]
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
    # The full Fortran cioffs formula is:
    #   kab   = kappa_d + alpha_d * bbd          (bbd = bias-boost firing duration)
    #   gxnc  = kab * exp(-utmbb / tau_d)
    # kpar (kappa_d) is zero for all IRAS detectors, so the kappa term drops.
    # bbd (bias-boost duration) is NOT loaded here; the field comes from cigbbi
    # in the Fortran code.  Setting bbd=0 gives:
    #   gxnc = alpha_d * exp(-utmbb / tau_d)
    # This is an approximation: when bbd > 0 the B3/B4 correction will be
    # slightly underestimated.  The effect is small (bbd is typically a few
    # seconds; alpha_d is ~1e-6 A/s), but it is an undocumented omission.
    gxnc = alpha * np.exp(-np.asarray(utmbb_samples, dtype=float) / tau)

    return volts - vfetc - gxnc


# ---------------------------------------------------------------------------
# PRHF (Photon Response History File) — ciprec.f / ciprhf.f
# ---------------------------------------------------------------------------

def read_prhf(prhf_dir, sop, obs):
    """
    Read one PRHF binary file (IR.PRHF.D0{sop:03d}{obs:03d}).

    Format: Fortran big-endian unformatted sequential.
    Each record: [deltut(f4), dtr[31](f4), lambda(f4), beta(f4)] = 34 floats.
    Records where lambda < -1000 AND beta < -1000 are EOF sentinels.

    Parameters
    ----------
    prhf_dir : str   path to the IPAC/P/ directory
    sop, obs : int

    Returns
    -------
    deltut : ndarray (N,) float64  time since last bias boost (UTC seconds)
    dtr    : ndarray (N, 31) float32  fractional responsivity depression,
             0-indexed by prhf_idx = CISEQN[det-1] - 32  (≥0 only for B3/B4)
    """
    fname = os.path.join(prhf_dir, f'IR.PRHF.D0{sop:03d}{obs:03d}')
    with open(fname, 'rb') as fh:
        raw = fh.read()

    N_FLOATS = 34          # deltut + 31 dtr + lambda + beta
    PAYLOAD  = N_FLOATS * 4
    STRIDE   = PAYLOAD + 8  # 4-byte FRL before + after

    deltut_list = []
    dtr_list    = []
    offset = 0
    while offset + STRIDE <= len(raw):
        # Fortran record-length markers (should both equal PAYLOAD = 136)
        offset += 4                            # skip leading FRL
        vals = np.frombuffer(raw[offset:offset + PAYLOAD], dtype='>f4')
        offset += PAYLOAD
        offset += 4                            # skip trailing FRL

        lambda_ = float(vals[32])
        beta    = float(vals[33])
        if lambda_ < -1000.0 and beta < -1000.0:
            continue                            # dummy/padding record — skip, keep reading

        deltut_list.append(float(vals[0]))
        dtr_list.append(vals[1:32].astype(np.float32))

    if not deltut_list:
        return np.empty(0, dtype=np.float64), np.empty((0, 31), dtype=np.float32)

    return (np.array(deltut_list, dtype=np.float64),
            np.array(dtr_list,   dtype=np.float32))


# PRHF Python index for each detector:  CISEQN[det-1] - 32  (>=0 for B3/B4)
# Negative values → no PRHF correction for that detector.
_PRHF_IDX = (CISEQN - 32).astype(int)   # shape (62,); index 0-30 for B3/B4


def apply_prhf_correction(flamp, det, utmbb_samples, prhf_deltut, prhf_dtr):
    """
    Apply PRHF photon-response-history correction to flash amplitudes.

    Implements ciprhf.f / ciprec.f logic, vectorised over samples.

    For each sample at time utmbb_t (seconds since last bias boost):
        delt = utmbb_t - 4                              # SLW convention
        find record j: prhf_deltut[j] ∈ (delt-4, delt+4]
                                       = (utmbb_t-8,   utmbb_t]
        flamp[sample] /= 1 - prhf_dtr[j, prhf_idx]

    Only applied when prhf_idx = CISEQN[det-1] - 32 ≥ 0  (Bands 3 & 4).
    If no record falls in the 8-second window the correction is skipped
    (matches ciprec ierr=3 "no correction for this case" path).

    Parameters
    ----------
    flamp         : ndarray (N,) float64  modified in-place
    det           : int  1-based detector number
    utmbb_samples : ndarray (N,) float64  seconds since last bias boost
    prhf_deltut   : ndarray (M,) float64  UTC-since-BB time tags (ascending)
    prhf_dtr      : ndarray (M, 31) float32  responsivity depression fractions
    """
    prhf_idx = int(_PRHF_IDX[det - 1])
    if prhf_idx < 0 or len(prhf_deltut) == 0:
        return                                      # not a B3/B4 detector, or no data

    t = np.asarray(utmbb_samples, dtype=np.float64)

    # Largest j such that prhf_deltut[j] ≤ t  (time-ordered records)
    j = np.searchsorted(prhf_deltut, t, side='right') - 1

    # Valid: record must be within 8 s behind current time
    # (mirrors ciprec's acceptance window: curdut ∈ (delt-4, delt+4] with delt=t-4)
    in_bounds = j >= 0
    in_bounds[in_bounds] &= prhf_deltut[j[in_bounds]] > t[in_bounds] - 8.0
    if not in_bounds.any():
        return

    deltr = prhf_dtr[j[in_bounds], prhf_idx].astype(np.float64)
    # Guard: deltr must be < 1 (fractional depression; typically 0–0.05)
    safe  = deltr < 1.0
    flamp[in_bounds] = np.where(safe,
                                flamp[in_bounds] / (1.0 - deltr),
                                flamp[in_bounds])


# ---------------------------------------------------------------------------
# ciwpm2 — amps → W/m²  (flash-model responsivity)
# ---------------------------------------------------------------------------

def amps_to_wm2(amps, det, utmbb_samples, apl, bpl, tpl,
                prhf_deltut=None, prhf_dtr=None):
    """
    Convert baseline-corrected amplitudes to W/m² using the SRHF responsivity,
    optionally applying the PRHF photon-response-history correction (ciprhf.f).

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
    prhf_deltut : ndarray (M,) or None.  If provided together with prhf_dtr,
        the PRHF photon-response correction is applied to the flash amplitude
        before computing responsivity.  Only affects Band 3 & 4 detectors.
    prhf_dtr : ndarray (M, 31) or None.  PRHF correction fractions.

    Returns
    -------
    wm2 : ndarray, flux density in W/m²
    """
    t = np.asarray(utmbb_samples, dtype=float)
    flamp = srhf_flash_amplitude(t, apl, bpl, tpl, det)

    # Optional PRHF correction (bands 3 & 4 only; no-op for bands 1 & 2)
    if prhf_deltut is not None and prhf_dtr is not None and len(prhf_deltut) > 0:
        apply_prhf_correction(flamp, det, t, prhf_deltut, prhf_dtr)

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
# E1 — Per-snip polynomial destripe  (Kester, IPAC IOM 701-039-89(1))
# ---------------------------------------------------------------------------

def poly_destripe(mjy_sr, bad, degree=1, n_iter=3, sigma=5.0):
    """
    Fit and subtract a robust polynomial baseline from one detector-snip.

    Implements the per-snip polynomial destripe step (E1) described in the
    pipeline notes (§ todq_destripe), corresponding to the GIPSY IMAGE task's
    ``destripe / tune = offset drift`` mode with a Tukey M-estimator.

    The IPAC IMAGE task performed this against map residuals (snip − map).
    This implementation fits directly to the TOD, using provided flags to
    exclude outliers (glitches, ICF events, bright sources, corrupted samples)
    so that the polynomial tracks the thermal / vfet drift baseline rather
    than sky signal.  This is appropriate before map-making; Commander3
    handles the map-residual version.

    Algorithm
    ---------
    1. Build an initial fit mask: all *un*-flagged, finite samples.
    2. For *n_iter* iterations:
       a. Fit a degree-*k* polynomial (``numpy.polyfit``) on the masked samples.
       b. Compute residuals for those samples; estimate noise via MAD / 0.6745
          (robust against the point sources that remain in early iterations).
       c. Reject samples with |residual| > *sigma* × rms from the next iteration.
    3. Take a final ``polyfit`` on the surviving samples; subtract the evaluated
       polynomial from **all** samples (flagged and source samples included —
       they receive the baseline correction but play no part in fitting it).

    Parameters
    ----------
    mjy_sr : ndarray, shape (N,)
        Calibrated signal in MJy sr⁻¹ (the destripe step operates in sky units).
    bad : bool ndarray, shape (N,)
        True = exclude from baseline fitting (but the subtraction is still
        applied).  Typical content: ``corrupt | flash_trim | icf_flag | spike |
        tail | source``.  NaN samples are also excluded regardless of this mask.
    degree : int
        Polynomial degree.  0 = DC offset only, 1 = offset + linear drift
        (the recommended choice for most survey and AO data), 2 = quadratic.
    n_iter : int
        Number of sigma-clip iterations.  3 is sufficient for typical data.
    sigma : float
        Rejection threshold in units of the robust RMS (MAD / 0.6745).

    Returns
    -------
    clean    : ndarray, shape (N,) — baseline-subtracted signal; NaN positions
               on input are NaN on output.
    baseline : ndarray, shape (N,) — the subtracted polynomial (zero if the
               fit could not be performed due to insufficient clean samples).
    """
    N = len(mjy_sr)
    # Normalised time axis on [-0.5, 0.5]; avoids numerical ill-conditioning
    # in polyfit for large N.
    t = np.linspace(-0.5, 0.5, N)
    baseline = np.zeros(N, dtype=np.float64)

    # Initial fit mask: flagged samples and NaN are excluded from fitting.
    use = (~np.asarray(bad, dtype=bool)) & np.isfinite(mjy_sr)
    if use.sum() < degree + 1:
        return mjy_sr.copy(), baseline

    for _ in range(n_iter):
        if use.sum() < degree + 1:
            break
        coeffs = np.polyfit(t[use], mjy_sr[use], degree)
        fit = np.polyval(coeffs, t)
        resid = mjy_sr[use] - fit[use]
        # Robust RMS via MAD so that bright sources still inside the fit
        # window do not inflate the noise estimate and widen the rejection band.
        mad = np.median(np.abs(resid - np.median(resid)))
        rms = mad / 0.6745
        if rms <= 0.0:
            break
        use = use & (np.abs(mjy_sr - fit) <= sigma * rms)

    if use.sum() < degree + 1:
        return mjy_sr.copy(), baseline

    # Final polynomial on surviving clean samples.
    coeffs = np.polyfit(t[use], mjy_sr[use], degree)
    baseline = np.polyval(coeffs, t)
    clean = mjy_sr - baseline
    # Preserve NaN at positions that were NaN on entry.
    clean[~np.isfinite(mjy_sr)] = np.nan
    return clean, baseline


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
    Full IPACCAL calibration pipeline for one detector observation segment.

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
        Seconds since the last bias boost at the start of the observation segment.
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

    .. WARNING::
        This function uses ``python/read_prhf.py``, which incorrectly
        interprets PRHF binary records: it reads floats at positions [17:32]
        as direct Jy/ADU responsivities, but those bytes are actually
        ``dtr[16:31]`` — fractional responsivity depressions (0–0.05,
        dimensionless).  Multiplying raw DN by these values does not produce
        flux in Jy.  Use ``calibrate()`` instead, which applies the correct
        ``ciprhf.f`` correction via ``apply_prhf_correction()``.

        This function is preserved for reference but should NOT be used for
        science.

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
# ZODYCAL position-dependent gain correction
# (IRZC_GAIN / IRZC_DETMODEL / IRZC_ZEM / IRZC_INTSKY from zodycal.src)
# ---------------------------------------------------------------------------

# Per-detector photoconductor model parameters (IRZC_DETMODEL, zodycal.src).
# Three values per detector: (Gmin, t=exp(-Δt/τ), b=increment).
# Indexed 0-based by det-1; ordering matches _BAND (B4×7, B3×8, B2×7, B1×8,
# B3×8, B2×8, B1×8, B4×8).  Bands 1 & 2 use gain=1.0 regardless of these
# stored values (t=b=0 is already the no-effect case).
_IRZC_DETPAR = np.array([
    # dets 1-7  (B4 module B)
    [0.4437, 0.9900, 2.612e-4], [0.4670, 0.9900, 2.537e-4],
    [0.9713, 0.9900, 2.779e-4], [0.4256, 0.9900, 2.443e-4],
    [0.4950, 0.9900, 2.055e-4], [0.6701, 0.9900, 2.593e-4],
    [0.7070, 0.9900, 3.548e-4],
    # dets 8-15 (B3 module B)
    [1.024,  0.8824, 4.553e-4], [0.8665, 0.8484, 2.587e-4],
    [0.8497, 0.9102, 5.356e-4], [0.2923, 0.8785, 1.057e-4],
    [0.7586, 0.9149, 5.371e-4], [0.9753, 0.8749, 5.750e-4],
    [0.9756, 0.9084, 5.819e-4], [1.757,  0.9101, 8.486e-4],
    # dets 16-22 (B2 module B — gain=1.0, t and b unused)
    [3.158, 0.0, 0.0], [1.000, 0.0, 0.0], [3.309, 0.0, 0.0],
    [3.480, 0.0, 0.0], [1.000, 0.0, 0.0], [3.220, 0.0, 0.0],
    [2.827, 0.0, 0.0],
    # dets 23-30 (B1 module B — gain=1.0)
    [1.128, 0.0, 0.0], [1.331, 0.0, 0.0], [1.224,  0.0, 0.0],
    [0.3352,0.0, 0.0], [0.8626,0.0, 0.0], [0.7178, 0.0, 0.0],
    [1.180, 0.0, 0.0], [1.294, 0.0, 0.0],
    # dets 31-38 (B3 module A)
    [0.3372, 0.9011, 1.554e-4], [1.061, 0.9003, 7.980e-4],
    [1.037,  0.9006, 5.888e-4], [1.072, 0.9338, 5.256e-4],
    [1.226,  0.9020, 7.758e-4], [1.000, 0.9000, 0.000   ],
    [0.9921, 0.9060, 7.272e-4], [0.6422,0.9254, 3.214e-4],
    # dets 39-46 (B2 module A — gain=1.0)
    [2.155, 0.0, 0.0], [3.870, 0.0, 0.0], [3.914, 0.0, 0.0],
    [4.050, 0.0, 0.0], [3.915, 0.0, 0.0], [3.686, 0.0, 0.0],
    [4.025, 0.0, 0.0], [2.075, 0.0, 0.0],
    # dets 47-54 (B1 module A — gain=1.0)
    [0.4628, 0.0, 0.0], [1.736, 0.0, 0.0], [1.845, 0.0, 0.0],
    [1.880,  0.0, 0.0], [1.625, 0.0, 0.0], [1.708, 0.0, 0.0],
    [1.769,  0.0, 0.0], [1.306, 0.0, 0.0],
    # dets 55-62 (B4 module A)
    [0.2225, 0.9900, 1.004e-4], [0.4359, 0.9900, 2.634e-4],
    [0.4347, 0.9900, 1.921e-4], [0.4376, 0.9900, 2.426e-4],
    [0.4123, 0.9900, 1.893e-4], [0.3103, 0.9900, 1.357e-4],
    [0.4356, 0.9900, 2.461e-4], [0.2744, 0.9900, 1.572e-4],
], dtype=np.float64)   # shape (62, 3)

# Zodiacal emission model parameters per band (IRZC_ZEM, zodycal.src).
# Columns: power_p, ellip_e, offset_alfa_rad, node_rad, scale_MJy_sr.
# Band index 0-based (band-1).
_IRZC_ZEMPAR = np.array([
    [1.779, 3.073, -0.03112, 1.184, 100.0],   # Band 1
    [1.300, 3.198, -0.02939, 1.296, 100.0],   # Band 2
    [1.046, 3.950, -0.02996, 1.275,  81.81],  # Band 3
    [1.000, 4.500, -0.01380, 0.994,  54.36],  # Band 4
], dtype=np.float64)

# Gauss-Legendre nodes and weights mapped to [0, 1] (for ZEM integration).
_GL_X, _GL_W = np.polynomial.legendre.leggauss(20)
_IRZC_GL_T  = (_GL_X + 1.0) / 2.0          # (20,) on [0, 1]
_IRZC_GL_WT = _GL_W / 2.0                   # (20,) weight scaled for [0, 1]

# IRCO_SUNLONG constants (irco_sunlong.shl)
_IRZC_SUNLONG_SUN0   = np.radians(279.10303475)
_IRZC_SUNLONG_SUN1   = np.radians(  0.98564735)
_IRZC_SUNLONG_EARTH0 = np.radians(356.45501990)
_IRZC_SUNLONG_EARTH1 = np.radians(  0.98560026)
_IRZC_SUNLONG_EQOC0  = np.radians(  1.915476)
_IRZC_SUNLONG_EQOC1  = np.radians(  0.020010)


def irco_sunlong(satcal):
    """
    Approximate ecliptic longitude of the sun (radians) from SATCAL tick.

    Port of irco_sunlong.shl.  Accuracy ~30 arcsec.

    Parameters
    ----------
    satcal : int or float — SATCAL clock ticks (seconds since IRAS epoch)
    """
    # Ephemeris time in days since 0 Jan 1983 0hr ET
    et     = (satcal + 53) * 1.0000526 / 86400.0 + 26.0
    slong  = _IRZC_SUNLONG_SUN0   + _IRZC_SUNLONG_SUN1   * et
    anom   = _IRZC_SUNLONG_EARTH0 + _IRZC_SUNLONG_EARTH1 * et
    eqoc   = (_IRZC_SUNLONG_EQOC0 * np.sin(anom)
            + _IRZC_SUNLONG_EQOC1 * np.sin(2.0 * anom))
    slong  = (slong + eqoc) % (2.0 * np.pi)
    return slong


def load_zodycal_map(maps_dir, band):
    """
    Load the galactic background sky map used by IRZC_INTSKY.

    The MAPS directory contains two FITS files (file000001.mt for Band 3,
    file000002.mt for Band 4).  GIPSY reads exactly MX×MY = 361×181 floats
    from each file (first 181 rows of NAXIS2=363).  The undefined value is
    identified as the value found in the first all-blank row and replaced
    with np.nan.

    Parameters
    ----------
    maps_dir : str — path to the MAPS directory (contains file*.mt)
    band     : int, 1-4 — IRAS band number

    Returns
    -------
    sky_map : float64 ndarray, shape (181, 361), values in MJy/sr (or
              NaN for unobserved pixels)
    """
    import fitsio
    skip = band - 3    # 0 for B3, 1 for B4
    fname = os.path.join(maps_dir, f'file{skip + 1:06d}.mt')
    with fitsio.FITS(fname) as ff:
        data = ff[0].read()[:181, :].astype(np.float64)  # (181, 361)
    # Identify blank value: row 181 is entirely blank in the source file.
    # Use the most common value in the extra rows as the sentinel.
    with fitsio.FITS(fname) as ff:
        extra = ff[0].read()[181, :]
    blank_val = float(extra[0])
    data[data == blank_val] = np.nan
    return data


def irzc_zem(psi_rad, theta_rad, lambda_rad, band):
    """
    Zodiacal emission model (IRZC_ZEM, zodycal.src).

    Integrates the ellipsoidal ZEM model along the line of sight for each
    of nz sample positions.  Integration performed via 20-point Gauss-
    Legendre quadrature, which is accurate to machine precision for the
    smooth 1/r^p integrand.

    Parameters
    ----------
    psi_rad   : (nz,) radians — scan position angle
    theta_rad : float radians — solar aspect angle (constant within one observation)
    lambda_rad: float radians — ecliptic longitude of the sun
    band      : int, 1-4

    Returns
    -------
    zem : (nz,) float64, zodiacal emission in MJy/sr
    """
    nz   = len(psi_rad)
    par  = _IRZC_ZEMPAR[band - 1]      # [p, e, alfa, node_rad, scale]
    p, e, alfa, node_rad, scale = par

    ecsun  = 0.016722       # Earth orbit eccentricity (IRZC_ZEM constant)
    perige = -1.374975      # perihelion angle in radians

    orbit  = 1.0 / (1.0 + ecsun * np.cos(lambda_rad - perige))
    orb2   = orbit * orbit
    e2     = e * e
    emin   = 1.0 - e2
    sunode = lambda_rad - node_rad
    sinsn  = np.sin(sunode)
    cossn  = np.cos(sunode)
    snor   = sinsn * orbit

    sinth  = np.sin(theta_rad)
    costh  = np.cos(theta_rad)
    ct2    = costh * costh
    st2    = sinth * sinth

    aest   = alfa * emin * sinth
    ssct   = sinsn * costh      # sin(sunode)*cos(theta)
    csst   = cossn * sinth      # cos(sunode)*sin(theta)
    orct   = orbit * costh
    aestsn = aest * snor

    sinp   = np.sin(psi_rad)    # (nz,)
    cosp   = np.cos(psi_rad)    # (nz,)
    cp2    = cosp * cosp
    sp2    = 1.0 - cp2
    cprot  = 2.0 * cosp * (sinp * csst - ssct)

    # Quadratic coefficients of the ZEM integrand denominator D(x) = a*x²-b*x+c
    # before the boundary shift from (0,∞) to (1,∞).
    a_orig = ct2 + st2 * (sp2 + e2 * cp2) + aest * cprot  # (nz,)
    b_orig = orct - aestsn * cosp                           # (nz,) — but see note
    c_orig = np.full(nz, orb2)                              # (nz,)

    # Boundary shift: x → y+1 (integral from y=0 to ∞ instead of x=1 to ∞).
    # New coefficients after y = x-1 substitution (see IRZC_MIDINF comments):
    a_sh = a_orig
    b_sh = 2.0 * (a_orig + b_orig)
    c_sh = a_orig + 2.0 * b_orig + c_orig

    # 20-point Gauss-Legendre quadrature of MIDINF integrand over [0, 1]:
    #   f(t) = scale * t² / (a_sh*t² - b_sh*t + c_sh)^p
    # with t = 1/y (so this covers y: 1 → ∞, x: 2 → ∞ + boundary shift).
    t  = _IRZC_GL_T[np.newaxis, :]   # (1, 20)
    wt = _IRZC_GL_WT[np.newaxis, :]  # (1, 20)
    a  = a_sh[:, np.newaxis]          # (nz, 1)
    b  = b_sh[:, np.newaxis]
    c  = c_sh[:, np.newaxis]
    denom     = np.abs(a * t**2 - b * t + c) ** p
    integrand = scale * t**2 / np.where(denom > 0, denom, 1e-30)
    zem = np.sum(integrand * wt, axis=1)   # (nz,)
    return np.maximum(zem, 0.0)


def irzc_intsky(psi_rad, theta_rad, lambda_rad, sky_map):
    """
    Interpolate the sky galactic background map at scan positions (IRZC_INTSKY).

    Converts sun-reference (psi, theta) to ecliptic (lon, lat) using the
    solar longitude rotation (IRCO_TRANSFORM SUNREF→ECLIP), then performs
    bilinear interpolation in the 181×361 sinusoidal-projection sky map.

    Parameters
    ----------
    psi_rad   : (nz,) radians
    theta_rad : float radians (constant within one observation)
    lambda_rad: float radians — solar longitude
    sky_map   : (181, 361) float64 — galactic background map, NaN=undefined

    Returns
    -------
    sky : (nz,) float64, sky brightness in MJy/sr (≥ 0)
    """
    nz   = len(psi_rad)
    # Sun-reference unit vector (identical to pointing_det without det offsets)
    xx   = np.cos(psi_rad) * np.sin(theta_rad)  # (nz,)
    yy   = np.sin(psi_rad) * np.sin(theta_rad)
    zz   = np.cos(theta_rad) * np.ones(nz)       # broadcast to (nz,)

    # Rotate to ecliptic (IRCO_TRANSFORM SUNREF→ECLIP):
    # ex = sin(L)*yy + cos(L)*zz, ey = -cos(L)*yy + sin(L)*zz, ez = xx
    sL = np.sin(lambda_rad); cL = np.cos(lambda_rad)
    ex = sL * yy + cL * zz   # ecliptic x
    ey = -cL * yy + sL * zz  # ecliptic y
    ez = xx                   # ecliptic z = sin(eclat)

    eclon_deg = np.degrees(np.arctan2(ey, ex))   # (-180, 180]
    eclat_deg = np.degrees(np.arcsin(np.clip(ez, -1.0, 1.0)))

    # Map coordinates: MY=181 rows (lat -90 to +90), MX=361 cols (lon -180 to +180)
    # Sinusoidal projection: column compressed by cos(lat) at each row.
    HMX, HMY = 181, 91   # Fortran 1-based centres (Python: HMX=180, HMY=90)
    FF = 180.0            # Fortran normalisation for lonf

    # Row index (0-based): bb = 90 + eclat_deg, ib = floor(bb)
    bb   = 90.0 + eclat_deg   # (nz,)
    ib   = np.clip(np.floor(bb).astype(int), 0, 179)

    # Longitude compression factor at row ib (matches Fortran nint rounding)
    lonf = np.round(np.cos(np.radians(ib - 90)) * FF) / FF   # (nz,)

    # Column index (0-based): ll = 180 - eclon_deg * lonf
    ll   = 180.0 - eclon_deg * lonf   # (nz,)
    il   = np.clip(np.floor(ll).astype(int), 0, 359)

    # Bilinear weights
    frac_lon = ll - np.floor(ll)   # along-longitude weight for il+1 col
    frac_lat = bb - ib             # along-latitude weight for ib+1 row
    il1  = np.minimum(il + 1, 360)
    ib1  = np.minimum(ib + 1, 180)

    def _lookup(row, col):
        v = sky_map[row, col]
        return np.where(np.isfinite(v), v, 0.0)   # treat undefined as 0

    v00 = _lookup(ib,  il )   # (nz,) lower-row, left-col
    v01 = _lookup(ib,  il1)   # lower-row, right-col
    v10 = _lookup(ib1, il )   # upper-row, left-col
    v11 = _lookup(ib1, il1)
    sky = ((1.0 - frac_lat) * ((1.0 - frac_lon) * v00 + frac_lon * v01) +
                  frac_lat  * ((1.0 - frac_lon) * v10 + frac_lon * v11))
    return np.maximum(sky, 0.0)


def irzc_detmodel(stim, det):
    """
    IIR photoconductor gain filter (IRZC_DETMODEL, zodycal.src).

    For bands 3 & 4: models the history-dependent gain of stressed Ge:Ga
    photoconductors.  The IIR recursion is:
        bug  = b * stim[k] + Gmin       (equilibrium gain at current stimulus)
        g[k] = bug + (g[k-1] - bug) * t (exponential decay toward equilibrium)
    For bands 1 & 2: gain = 1.0 (no memory effect).

    Parameters
    ----------
    stim : (nz,) float64 — sky stimulus (ZEM + galactic sky, MJy/sr)
    det  : int, 1-based detector number

    Returns
    -------
    g : (nz,) float64 — detector gain at each model time step
    """
    band = _BAND[det - 1]
    nz   = len(stim)
    if band <= 2:
        return np.ones(nz, dtype=np.float64)

    gmin, t, b = _IRZC_DETPAR[det - 1]
    g_arr = np.empty(nz, dtype=np.float64)
    g = gmin + b * stim[0]
    for k in range(nz):
        bug  = b * stim[k] + gmin
        g    = bug + (g - bug) * t
        g_arr[k] = g
    return g_arr


def irzc_gain(det, psi0_deg, psirate_deg_pt, theta_deg, lambda_rad,
              cdelt, ng, sky_map, yloc_arcmin=0.0, zloc_arcmin=0.0):
    """
    Compute per-sample photoconductor gain correction (IRZC_GAIN, zodycal.src).

    Builds a coarse sky-trajectory grid at TIMESTEP=15-tick intervals, looks
    up ZEM + galactic sky stimulus at each point, runs the IRZC_DETMODEL IIR
    filter, then interpolates to the full sample rate.

    The gain values are normalised to g[0] so the correction is purely
    relative: samples at higher gain state (scanning a bright region) are
    divided by a factor > 1, reducing their apparent brightness.

    Parameters
    ----------
    det             : int, 1-based detector number
    psi0_deg        : float — PSI angle at observation start (degrees)
    psirate_deg_pt  : float — PSI rate in deg/SATCAL-tick
    theta_deg       : float — solar aspect angle at observation start (degrees)
    lambda_rad      : float — solar ecliptic longitude (radians)
    cdelt           : float — SATCAL ticks per sample (e.g. 0.25 for B4 survey)
    ng              : int   — number of samples in the observation
    sky_map         : (181, 361) ndarray or None — galactic background map.
                      If None, only the ZEM component is used.
    yloc_arcmin     : float — detector focal-plane cross-scan offset (arcmin)
    zloc_arcmin     : float — detector focal-plane in-scan offset (arcmin)

    Returns
    -------
    gain_norm : (ng,) float64 — gain relative to observation start (1.0 = unchanged).
                Divide calibrated MJy/sr by this to correct for gain variation.
    """
    TIMESTEP = 15        # model time step in SATCAL ticks (zodycal.src MZ grid)
    M2R      = np.pi / (180.0 * 60.0)    # arcminutes to radians

    # Focal-plane offsets: dth = yloc (cross-scan, theta direction)
    #                      dpsi = zloc (in-scan, psi direction)
    dth_rad  = yloc_arcmin * M2R
    dpsi_rad = zloc_arcmin * M2R
    theta_rad = np.radians(theta_deg)
    psi0_rad  = np.radians(psi0_deg)

    # Detector centre angles (mirrors IRZC_GAIN: IRCC_MASK adjustment)
    th0 = theta_rad - dth_rad
    ps0 = psi0_rad  + dpsi_rad / max(np.sin(th0), 1e-6)

    # Number of model time steps (crval=0: observation starts at its own origin)
    nz     = int(round(cdelt * ng)) // TIMESTEP + 1
    nz     = max(2, min(nz, 200))   # clamp to MZ=200

    # Build position trajectory at TIMESTEP intervals, backing up 7 ticks
    # (matches ps0 = ps0 + 7 * psirate in Fortran)
    psdot  = np.radians(psirate_deg_pt) * TIMESTEP
    ps_arr = ps0 + 7.0 * np.radians(psirate_deg_pt) + np.arange(nz) * psdot
    th_arr = np.full(nz, th0)

    # Sky stimulus at each model time step
    zem  = irzc_zem(ps_arr, th0, lambda_rad, _BAND[det - 1])
    if sky_map is not None:
        sky  = irzc_intsky(ps_arr, th_arr, lambda_rad, sky_map)
    else:
        sky  = np.zeros(nz)
    stim = zem + sky

    # IIR photoconductor gain filter
    g = irzc_detmodel(stim, det)

    # Interpolate gain from nz model points to ng samples
    # start = (crval - 7) / TIMESTEP = -7/15 ≈ -0.467 with crval=0
    start_arr = -(7.0 / TIMESTEP) + np.arange(ng) * (cdelt / TIMESTEP)
    # Fortran 1-indexed: ks = int(start); Python: ks_py = int(start) - 1
    ks_float  = start_arr - 1.0          # fractional index into g[0..]
    ks_py     = np.clip(np.floor(ks_float).astype(int), 0, nz - 2)
    frac      = np.clip(ks_float - np.floor(ks_float), 0.0, 1.0)
    gain_ng   = np.where(start_arr < 1.0,
                         g[0],
                         (1.0 - frac) * g[ks_py] + frac * g[np.minimum(ks_py + 1, nz - 1)])

    # Normalise to the first sample so we only correct the RELATIVE gain change
    # (the absolute calibration level is already handled by the SRHF flash model)
    g0 = gain_ng[0] if gain_ng[0] > 0.0 else 1.0
    return gain_ng / g0


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
