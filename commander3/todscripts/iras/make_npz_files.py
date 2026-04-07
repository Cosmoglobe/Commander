#!/usr/bin/env python3
"""Parallel continuous timestream builder for IRAS detectors (all bands or specific detectors).

Writes one merged file per SOP/OBS/detector:
  <out-dir>/sopXXX/obsYYY/detNN_continuous.npz

For each task (sop, obs, det):
- read all matching plate files from data_sopos_*/sopXXX/obsYYY/plate*/detNN.tbl
- compute per-sample UTC from utcs1/utcs2/npts
- mask bad rows (ra == -999 or flux < -1e10)
- sort by time and collapse near-duplicates within band-dependent dedup threshold
"""

import argparse
import os
import re
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

# Detector groups from ircc_const.c / BAND_CFG in iras_tod.py
# Band numbering here follows IRAS convention (I..IV mapped to 1..4).
BAND1_DETS = [23, 24, 25, 26, 27, 28, 29, 30, 47, 48, 49, 50, 51, 52, 53, 54]
BAND2_DETS = [16, 17, 18, 19, 20, 21, 22, 39, 40, 41, 42, 43, 44, 45, 46]
BAND3_DETS = [8, 9, 10, 11, 12, 13, 14, 15, 31, 32, 33, 34, 35, 36, 37, 38]
BAND4_DETS = [1, 2, 3, 4, 5, 6, 7, 55, 56, 57, 58, 59, 60, 61, 62]

ALL_DETS = sorted(BAND1_DETS + BAND2_DETS + BAND3_DETS + BAND4_DETS)

BAND_DEFINITIONS = {
    1: BAND1_DETS,
    2: BAND2_DETS,
    3: BAND3_DETS,
    4: BAND4_DETS,
}

# Survey sample rates in samples per SATCAL tick.
# Dedup threshold is the nominal inter-sample interval: satcal / rate.
DEDUP_SATCAL_DIVISORS = {
    1: 16,
    2: 16,
    3: 8,
    4: 4,
}


def det_to_band(det: int) -> int | None:
    """Map detector number to band (1-4), or None if unknown."""
    for band, dets in BAND_DEFINITIONS.items():
        if det in dets:
            return band
    return None


def get_active_detectors(args: argparse.Namespace) -> list[int]:
    """Determine which detectors to include based on band selection."""
    if args.all_detectors:
        return ALL_DETS
    elif args.band:
        if args.band not in BAND_DEFINITIONS:
            raise SystemExit(f"Invalid band: {args.band}. Must be 1, 2, 3, or 4.")
        return BAND_DEFINITIONS[args.band]
    else:
        return BAND1_DETS  # Default to band-1


def get_env_var_name(args: argparse.Namespace) -> str:
    """Get the appropriate environment variable name based on selection."""
    if args.all_detectors:
        return "IRAS_ALL_TIMESTREAM_ROOT"
    elif args.band:
        return f"IRAS_BAND{args.band}_TIMESTREAM_ROOT"
    else:
        return "IRAS_BAND1_TIMESTREAM_ROOT"


def get_dedup_threshold(det: int, satcal: float, explicit_threshold: float | None) -> float:
    """
    Get the dedup threshold for a detector.
    
    If explicit_threshold is provided, use it.
    Otherwise, use the band-dependent divisor.
    """
    if explicit_threshold is not None:
        return explicit_threshold
    
    band = det_to_band(det)
    if band is None:
        # Unknown detector; use default divisor
        divisor = DEDUP_SATCAL_DIVISORS.get(1, 16)
    else:
        divisor = DEDUP_SATCAL_DIVISORS.get(band, 16)
    
    return satcal / divisor


def resolve_roots(args: argparse.Namespace) -> tuple[Path, Path]:
    """Resolve path inputs with CLI > environment variable precedence."""
    data_root_str = args.data_root or os.getenv("IRAS_DATA_ROOT")
    env_var_name = get_env_var_name(args)
    out_dir_str = args.out_dir or os.getenv(env_var_name)

    if not data_root_str:
        raise SystemExit("Missing data root. Provide --data-root or set IRAS_DATA_ROOT.")
    if not out_dir_str:
        raise SystemExit(
            f"Missing output root. Provide --out-dir or set {env_var_name}."
        )

    data_root = Path(data_root_str).expanduser().resolve()
    out_dir = Path(out_dir_str).expanduser().resolve()

    if not data_root.exists():
        raise SystemExit(f"Data root does not exist: {data_root}")
    if not data_root.is_dir():
        raise SystemExit(f"Data root is not a directory: {data_root}")

    validate_reorganized_layout(data_root)

    return data_root, out_dir


def validate_reorganized_layout(data_root: Path) -> None:
    """
    Validate that data_root follows the reorganized layout:
      data_root/sopXXX/obsYYY/plateZZZZ/detNN.tbl
    """
    sop_dirs = sorted(p for p in data_root.glob("sop*") if p.is_dir())
    if not sop_dirs:
        plate_dirs = sorted(p for p in data_root.glob("plate*") if p.is_dir())
        if plate_dirs:
            raise SystemExit(
                "Detected old plate-root layout (plate*/sop*/obs*/det*.tbl). "
                "This script only supports reorganized layout (sop*/obs*/plate*/det*.tbl). "
                "Point --data-root to reorg_l1_files.py output."
            )
        raise SystemExit(
            "Could not find any sop* directories under data root. "
            "Expected reorganized layout: sop*/obs*/plate*/det*.tbl"
        )

    found_obs = False
    found_plate = False
    for sop_dir in sop_dirs:
        obs_dirs = [p for p in sop_dir.glob("obs*") if p.is_dir()]
        if obs_dirs:
            found_obs = True
        for obs_dir in obs_dirs:
            if any(p.is_dir() for p in obs_dir.glob("plate*")):
                found_plate = True
                break
        if found_plate:
            break

    if not found_obs:
        raise SystemExit(
            "Found sop* directories, but no obs* directories beneath them. "
            "Expected reorganized layout: sop*/obs*/plate*/det*.tbl"
        )
    if not found_plate:
        raise SystemExit(
            "Found sop*/obs* directories, but no plate* directories beneath obs*. "
            "Expected reorganized layout: sop*/obs*/plate*/det*.tbl"
        )


def parse_det_from_name(filename: str) -> int | None:
    """Extract detector number from detXX.tbl filename."""
    if not filename.startswith("det") or not filename.endswith(".tbl"):
        return None
    try:
        return int(filename[3:-4])
    except ValueError:
        return None


def discover_obs_dirs(data_root: Path, sop: int | None, obs: int | None) -> list[Path]:
    """
    Discover observation directories in reorganized layout:
      data_root/sopXXX/obsYYY/plateZZZZ/detNN.tbl
    """
    if sop is not None and obs is not None:
        target = data_root / f"sop{sop:03d}" / f"obs{obs:03d}"
        return [target] if target.is_dir() else []

    if sop is not None:
        sop_dir = data_root / f"sop{sop:03d}"
        if not sop_dir.is_dir():
            return []
        return sorted(p for p in sop_dir.glob("obs*") if p.is_dir())

    obs_dirs: list[Path] = []
    for sop_dir in sorted(p for p in data_root.glob("sop*") if p.is_dir()):
        obs_dirs.extend(sorted(p for p in sop_dir.glob("obs*") if p.is_dir()))
    return obs_dirs


def parse_header(path: Path) -> dict:
    hdr = {}
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s.startswith("\\") and "=" in s:
                k, v = s[1:].split("=", 1)
                hdr[k.strip().lower()] = v.strip()
            elif s.startswith("|"):
                break
    return hdr


def data_start_line(path: Path) -> int:
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f):
            s = line.strip()
            if not s or s.startswith("\\") or s.startswith("|"):
                continue
            return i
    return 0


def load_file(path: Path, bad_flux_floor: float = -1e10) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    hdr = parse_header(path)
    for k in ("utcs1", "utcs2", "npts"):
        if k not in hdr:
            raise ValueError(f"Missing header key '{k}' in {path}")

    utcs1 = float(hdr["utcs1"])
    utcs2 = float(hdr["utcs2"])
    npts = int(hdr["npts"])

    skip = data_start_line(path)
    data = np.loadtxt(path, skiprows=skip)
    if data.ndim == 1 and data.size > 0:
        data = data[None, :]

    if data.shape[0] == 0:
        return (
            np.empty(0, dtype=np.float64),
            np.empty(0, dtype=np.float64),
            np.empty(0, dtype=np.float64),
            np.empty(0, dtype=np.float64),
        )

    bad = (data[:, 0] == -999) | (data[:, 2] < bad_flux_floor)
    data = data.astype(np.float64, copy=True)
    data[bad, :] = np.nan

    n_use = min(npts, data.shape[0])
    if n_use <= 1:
        t = np.array([utcs1], dtype=np.float64)
    else:
        t = np.linspace(utcs1, utcs2, n_use, dtype=np.float64)

    return t, data[:n_use, 0], data[:n_use, 1], data[:n_use, 2]


def merge_near_duplicates(
    t: np.ndarray,
    ra: np.ndarray,
    dec: np.ndarray,
    flux: np.ndarray,
    eps: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if t.size == 0:
        return t, ra, dec, flux, np.empty(0, dtype=np.int32)

    order = np.argsort(t)
    t = t[order]
    ra = ra[order]
    dec = dec[order]
    flux = flux[order]

    out_t, out_ra, out_dec, out_flux, out_n = [], [], [], [], []

    i = 0
    n = t.size
    while i < n:
        j = i + 1
        t0 = t[i]
        while j < n and (t[j] - t0) < eps:
            j += 1

        out_t.append(np.nanmedian(t[i:j]))
        out_ra.append(np.nanmedian(ra[i:j]))
        out_dec.append(np.nanmedian(dec[i:j]))
        out_flux.append(np.nanmedian(flux[i:j]))
        out_n.append(j - i)
        i = j

    return (
        np.array(out_t, dtype=np.float64),
        np.array(out_ra, dtype=np.float64),
        np.array(out_dec, dtype=np.float64),
        np.array(out_flux, dtype=np.float64),
        np.array(out_n, dtype=np.int32),
    )


def merge_file_group(files: list[str], eps: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    t_all, ra_all, dec_all, flux_all = [], [], [], []
    for fp in files:
        path = Path(fp)
        try:
            t, ra, dec, flux = load_file(path)
        except Exception:
            continue
        if t.size == 0:
            continue
        t_all.append(t)
        ra_all.append(ra)
        dec_all.append(dec)
        flux_all.append(flux)

    if not t_all:
        raise ValueError("No usable samples found")

    t = np.concatenate(t_all)
    ra = np.concatenate(ra_all)
    dec = np.concatenate(dec_all)
    flux = np.concatenate(flux_all)

    keep = ~(np.isnan(ra) & np.isnan(dec) & np.isnan(flux))
    t = t[keep]
    ra = ra[keep]
    dec = dec[keep]
    flux = flux[keep]

    t_m, ra_m, dec_m, flux_m, n_merged = merge_near_duplicates(t, ra, dec, flux, eps)
    return t_m, ra_m, dec_m, flux_m, n_merged, int(t.size)


def run_task(task: tuple[str, str, int, list[str], str, float, bool, float]) -> tuple[str, str, int, int, int, int, str]:
    sop, obs, det, files, out_dir, satcal, overwrite, explicit_eps = task
    eps = get_dedup_threshold(det, satcal, explicit_eps)
    out_path = Path(out_dir) / sop / obs / f"det{det:02d}_continuous.npz"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and not overwrite:
        return sop, obs, det, len(files), -1, -1, str(out_path)

    try:
        t_m, ra_m, dec_m, flux_m, n_merged, n_raw = merge_file_group(files, eps)
    except Exception as exc:
        return sop, obs, det, len(files), -2, -2, f"ERROR:{exc}"

    np.savez_compressed(
        out_path,
        t=t_m,
        ra=ra_m,
        dec=dec_m,
        flux=flux_m,
        n_merged=n_merged,
        dedup_seconds=float(eps),
        sop=sop,
        obs=obs,
        det=int(det),
    )

    return sop, obs, det, len(files), n_raw, int(t_m.size), str(out_path)


def discover_tasks_for_obs(
    obs_dir: Path,
    active_dets: set[int],
    out_dir: str,
    satcal: float,
    explicit_eps: float,
    overwrite: bool,
) -> list[tuple[str, str, int, list[str], str, float, bool, float]]:
    """Discover all detector tasks for one sopXXX/obsYYY directory."""
    parts = obs_dir.parts
    try:
        sop_name = next(p for p in parts if re.fullmatch(r"sop\d{3}", p))
        obs_name = next(p for p in parts if re.fullmatch(r"obs\d{3}", p))
    except StopIteration:
        return []

    groups: dict[int, list[str]] = {}
    plate_dirs = sorted(p for p in obs_dir.glob("plate*") if p.is_dir())
    for plate_dir in plate_dirs:
        for tbl in plate_dir.glob("det*.tbl"):
            det = parse_det_from_name(tbl.name)
            if det is None or det not in active_dets:
                continue
            groups.setdefault(det, []).append(str(tbl))

    tasks: list[tuple[str, str, int, list[str], str, float, bool, float]] = []
    for det in sorted(groups):
        tasks.append((sop_name, obs_name, det, groups[det], out_dir, satcal, overwrite, explicit_eps))
    return tasks


def build_tasks(
    data_root: str,
    dets: list[int],
    sop: int | None,
    obs: int | None,
    out_dir: str,
    satcal: float,
    explicit_eps: float,
    overwrite: bool,
) -> list[tuple[str, str, int, list[str], str, float, bool, float]]:
    active_dets = set(dets)
    obs_dirs = discover_obs_dirs(Path(data_root), sop, obs)
    print(f"Discovered {len(obs_dirs)} sop/obs directories")

    tasks: list[tuple[str, str, int, list[str], str, float, bool, float]] = []
    workers = min(max(len(obs_dirs), 1), os.cpu_count() or 4)
    done = 0
    progress_every = max(100, len(obs_dirs) // 20) if obs_dirs else 100
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {
            ex.submit(discover_tasks_for_obs, obs_dir, active_dets, out_dir, satcal, explicit_eps, overwrite): obs_dir
            for obs_dir in obs_dirs
        }
        for fut in as_completed(futs):
            obs_tasks = fut.result()
            tasks.extend(obs_tasks)
            done += 1
            if done % progress_every == 0 or done == len(obs_dirs):
                print(f"  Discovery progress: {done}/{len(obs_dirs)} obs dirs, {len(tasks)} tasks")

    return tasks


def main() -> None:
    p = argparse.ArgumentParser(description="Parallel multi-band IRAS SOP/OBS continuous timestream builder")
    p.add_argument(
        "--data-root",
        default=None,
        help="Root reorganized data directory containing sop*/obs*/plate*/det*.tbl (or IRAS_DATA_ROOT)",
    )
    p.add_argument(
        "--out-dir",
        default=None,
        help="Output root (or band-dependent IRAS_BAND*_TIMESTREAM_ROOT)",
    )
    p.add_argument("--workers", type=int, default=None, help="Worker processes (default: os.cpu_count())")
    p.add_argument("--satcal", type=float, default=1.0, help="SATCAL tick period in seconds (default: 1.0)")
    p.add_argument(
        "--band",
        type=int,
        choices=[1, 2, 3, 4],
        default=None,
        help="Target IRAS band index (1..4). Defaults to band 1. Ignored if --all-detectors is set.",
    )
    p.add_argument(
        "--all-detectors",
        action="store_true",
        help="Include all 62 detectors (overrides --band selection)",
    )
    p.add_argument(
        "--dedup-seconds",
        type=float,
        default=None,
        help="Explicit dedup threshold in seconds (overrides band-dependent defaults)",
    )
    p.add_argument("--sop", type=int, default=None, help="Optional SOP filter")
    p.add_argument("--obs", type=int, default=None, help="Optional OBS filter")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing output files")
    p.add_argument("--progress-every", type=int, default=200, help="Progress print cadence in completed tasks")
    args = p.parse_args()

    # Determine active detectors
    active_dets = get_active_detectors(args)
    
    # Explicit dedup threshold overrides band-dependent defaults
    explicit_eps = args.dedup_seconds if args.dedup_seconds is not None else None
    
    workers = args.workers if args.workers is not None else (os.cpu_count() or 4)

    data_root, out_root = resolve_roots(args)
    out_root.mkdir(parents=True, exist_ok=True)

    selection_info = f"all {len(active_dets)} detectors" if args.all_detectors else f"band {args.band or 1}"
    print(f"Selection: {selection_info}")
    
    print("Building task list...")
    tasks = build_tasks(
        str(data_root),
        active_dets,
        args.sop,
        args.obs,
        str(out_root),
        args.satcal,
        explicit_eps,
        args.overwrite,
    )
    if not tasks:
        raise SystemExit("No tasks found. Check filters and data-root.")

    print(f"Tasks               : {len(tasks)}")
    print(f"Workers             : {workers}")
    print(f"SATCAL              : {args.satcal:.9f} s")
    if explicit_eps is not None:
        print(f"Dedup threshold [s] : {explicit_eps:.9f} (explicit, fixed for all detectors)")
    else:
        print(f"Dedup threshold [s] : band-dependent (satcal / divisor)")
        for band in [1, 2, 3, 4]:
            divisor = DEDUP_SATCAL_DIVISORS.get(band, 16)
            threshold = args.satcal / divisor
            print(f"  Band {band}: rate={divisor} samples/tick -> satcal / {divisor} = {threshold:.9f}")
    print(f"Data root           : {data_root}")
    print(f"Output root         : {out_root}")

    done = 0
    wrote = 0
    skipped = 0
    errors = 0

    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(run_task, t) for t in tasks]

        for fut in as_completed(futs):
            done += 1
            try:
                sop, obs, det, n_files, n_raw, n_merged, out_path = fut.result()
            except Exception as exc:
                errors += 1
                if done % args.progress_every == 0 or done == len(tasks):
                    print(f"[{done}/{len(tasks)}] wrote={wrote} skipped={skipped} errors={errors} EXCEPTION: {exc}")
                continue

            if n_raw == -2:
                errors += 1
                if done % args.progress_every == 0 or done == len(tasks):
                    print(f"[{done}/{len(tasks)}] wrote={wrote} skipped={skipped} errors={errors} ERROR {sop}/{obs}/det{det:02d}: {out_path}")
                continue
            elif n_raw == -1:
                skipped += 1
            else:
                wrote += 1

            if done % args.progress_every == 0 or done == len(tasks):
                print(
                    f"[{done}/{len(tasks)}] wrote={wrote} skipped={skipped} errors={errors} "
                    f"last={sop}/{obs}/det{det:02d} raw={n_raw} merged={n_merged}"
                )

    print("\nDone.")
    print(f"Wrote files         : {wrote}")
    print(f"Skipped existing    : {skipped}")
    print(f"Errors (skipped)    : {errors}")
    print(f"Output root         : {out_root}")


if __name__ == "__main__":
    main()
