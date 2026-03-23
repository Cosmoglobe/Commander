#!/usr/bin/env python3
"""Parallel continuous timestream builder for IRAS band-1 detectors.

Writes one merged file per SOP/OBS/detector:
  <out-dir>/sopXXX/obsYYY/detNN_continuous.npz

For each task (sop, obs, det):
- read all matching plate files from data/plate*/sopXXX/obsYYY/detNN.tbl
- compute per-sample UTC from utcs1/utcs2/npts
- mask bad rows (ra == -999 or flux < -1e10)
- sort by time and collapse near-duplicates within dedup threshold
"""

import argparse
import glob
import os
import re
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

BAND1_DETS = [23, 24, 25, 26, 27, 28, 29, 30, 47, 48, 49, 50, 51, 52, 53, 54]


def resolve_roots(args: argparse.Namespace) -> tuple[Path, Path]:
    """Resolve path inputs with CLI > environment variable precedence."""
    data_root_str = args.data_root or os.getenv("IRAS_DATA_ROOT")
    out_dir_str = args.out_dir or os.getenv("IRAS_BAND1_TIMESTREAM_ROOT")

    if not data_root_str:
        raise SystemExit("Missing data root. Provide --data-root or set IRAS_DATA_ROOT.")
    if not out_dir_str:
        raise SystemExit(
            "Missing output root. Provide --out-dir or set IRAS_BAND1_TIMESTREAM_ROOT."
        )

    data_root = Path(data_root_str).expanduser().resolve()
    out_dir = Path(out_dir_str).expanduser().resolve()

    if not data_root.exists():
        raise SystemExit(f"Data root does not exist: {data_root}")
    if not data_root.is_dir():
        raise SystemExit(f"Data root is not a directory: {data_root}")

    return data_root, out_dir


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


def run_task(task: tuple[str, str, int, list[str], str, float, bool]) -> tuple[str, str, int, int, int, int, str]:
    sop, obs, det, files, out_dir, eps, overwrite = task
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


def discover_tasks_for_det(
    data_root: str,
    sop_pat: str,
    obs_pat: str,
    det: int,
    out_dir: str,
    eps: float,
    overwrite: bool,
) -> list[tuple[str, str, int, list[str], str, float, bool]]:
    """Discover (sop, obs, det) tasks for one detector."""
    pattern = os.path.join(data_root, f"plate*/{sop_pat}/{obs_pat}/det{det:02d}.tbl")
    det_files = [Path(p) for p in glob.glob(pattern)]
    groups: dict[tuple[str, str], list[str]] = {}

    for p in det_files:
        sop_name = None
        obs_name = None
        for part in p.parts:
            if sop_name is None and re.fullmatch(r"sop\d{3}", part):
                sop_name = part
            if obs_name is None and re.fullmatch(r"obs\d{3}", part):
                obs_name = part
        if sop_name is None or obs_name is None:
            continue
        groups.setdefault((sop_name, obs_name), []).append(str(p))

    tasks = []
    for (sop_name, obs_name), files in groups.items():
        tasks.append((sop_name, obs_name, det, files, out_dir, eps, overwrite))
    return tasks


def build_tasks(data_root: str, dets: list[int], sop: int | None, obs: int | None, out_dir: str, eps: float, overwrite: bool) -> list[tuple[str, str, int, list[str], str, float, bool]]:
    sop_pat = f"sop{sop:03d}" if sop is not None else "sop*"
    obs_pat = f"obs{obs:03d}" if obs is not None else "obs*"

    # Discover each detector's task list concurrently.
    per_det: dict[int, list] = {}
    with ThreadPoolExecutor(max_workers=min(len(dets), os.cpu_count() or 4)) as ex:
        futs = {
            ex.submit(discover_tasks_for_det, data_root, sop_pat, obs_pat, det, out_dir, eps, overwrite): det
            for det in dets
        }
        for fut in as_completed(futs):
            det = futs[fut]
            det_tasks = fut.result()
            per_det[det] = det_tasks
            print(f"Discovered det{det:02d}: {len(det_tasks)} sop/obs tasks")

    # Interleave tasks round-robin across detectors so all detectors
    # make progress simultaneously rather than processing one detector
    # at a time.
    tasks = []
    lists = [per_det[d] for d in dets if d in per_det]
    max_len = max((len(l) for l in lists), default=0)
    for i in range(max_len):
        for l in lists:
            if i < len(l):
                tasks.append(l[i])
    return tasks


def main() -> None:
    p = argparse.ArgumentParser(description="Parallel band-1 SOP/OBS continuous timestream builder")
    p.add_argument("--data-root", default=None, help="Root data directory (or IRAS_DATA_ROOT)")
    p.add_argument(
        "--out-dir",
        default=None,
        help="Output root (or IRAS_BAND1_TIMESTREAM_ROOT)",
    )
    p.add_argument("--workers", type=int, default=None, help="Worker processes (default: os.cpu_count())")
    p.add_argument("--satcal", type=float, default=1.0, help="SATCAL seconds (dedup threshold default is satcal/16)")
    p.add_argument("--dedup-seconds", type=float, default=None, help="Explicit dedup threshold in seconds")
    p.add_argument("--dets", type=int, nargs="+", default=BAND1_DETS, help="Detector IDs (default: band-1 detectors)")
    p.add_argument("--sop", type=int, default=None, help="Optional SOP filter")
    p.add_argument("--obs", type=int, default=None, help="Optional OBS filter")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing output files")
    p.add_argument("--progress-every", type=int, default=200, help="Progress print cadence in completed tasks")
    args = p.parse_args()

    eps = args.dedup_seconds if args.dedup_seconds is not None else (args.satcal / 16.0)
    workers = args.workers if args.workers is not None else (os.cpu_count() or 4)

    data_root, out_root = resolve_roots(args)
    out_root.mkdir(parents=True, exist_ok=True)

    print("Building task list...")
    tasks = build_tasks(str(data_root), args.dets, args.sop, args.obs, str(out_root), eps, args.overwrite)
    if not tasks:
        raise SystemExit("No tasks found. Check filters and data-root.")

    print(f"Tasks               : {len(tasks)}")
    print(f"Workers             : {workers}")
    print(f"Dedup threshold [s] : {eps:.9f}")
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
