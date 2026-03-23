#!/usr/bin/env python3
"""
Build band-1 SOP/OBS parallel dataset by scanning obs directories.
Instead of iterating over millions of files, we iterate over obs directories
and filter band-1 .tbl files within each obs.

Reorganizes /data/plate*/sop*/obs*/det*.tbl into /data_sopos_b1/sop*/obs*/plate*/det*.tbl
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set

# Band-1 detector IDs per IRAS Explanatory Supplement
BAND1_DETS = {23, 24, 25, 26, 27, 28, 29, 30, 47, 48, 49, 50, 51, 52, 53, 54}


def resolve_roots(args: argparse.Namespace) -> tuple[Path, Path]:
    """Resolve path inputs with CLI > environment variable precedence."""
    data_root_str = args.data_root or os.getenv("IRAS_DATA_ROOT")
    out_root_str = args.out_root or os.getenv("IRAS_BAND1_REORG_ROOT")

    if not data_root_str:
        raise SystemExit(
            "Missing data root. Provide --data-root or set IRAS_DATA_ROOT."
        )
    if not out_root_str:
        raise SystemExit(
            "Missing output root. Provide --out-root or set IRAS_BAND1_REORG_ROOT."
        )

    data_root = Path(data_root_str).expanduser().resolve()
    out_root = Path(out_root_str).expanduser().resolve()

    if not data_root.exists():
        raise SystemExit(f"Data root does not exist: {data_root}")
    if not data_root.is_dir():
        raise SystemExit(f"Data root is not a directory: {data_root}")

    return data_root, out_root


def parse_det_from_name(filename: str) -> Optional[int]:
    """Extract detector number from 'detXX.tbl' filename."""
    if not filename.startswith("det") or not filename.endswith(".tbl"):
        return None
    try:
        det_str = filename[3:-4]  # "det<XX>.tbl" -> "<XX>"
        return int(det_str)
    except ValueError:
        return None


def parse_tbl_jband(filepath: Path) -> int:
    """Extract \\jband value from .tbl header, or return -1 if not found/error."""
    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if line.startswith("\\jband"):
                    parts = line.split("=")
                    if len(parts) >= 2:
                        try:
                            return int(parts[1].strip())
                        except ValueError:
                            pass
                    break
    except (OSError, IOError):
        pass
    return -1


def find_obs_dirs(data_root: Path, sop: Optional[str] = None, obs: Optional[str] = None) -> List[Path]:
    """Use find to list obs directories efficiently, optionally filtered by SOP/OBS."""
    try:
        result = subprocess.run(
            ["find", str(data_root), "-maxdepth", "3", "-type", "d",
             "-name", "obs*"],
            capture_output=True,
            text=True,
            timeout=None
        )
        if result.returncode == 0:
            obs_dirs = [
                Path(line.strip()) for line in result.stdout.splitlines()
                if line.strip()
            ]
            
            # Apply SOP/OBS filters if provided
            if sop is not None or obs is not None:
                filtered = []
                for obs_dir in obs_dirs:
                    try:
                        parts = obs_dir.relative_to(data_root).parts
                        if len(parts) >= 3:
                            sop_name, obs_name = parts[1], parts[2]
                            if sop is not None and sop_name != f"sop{sop:03d}":
                                continue
                            if obs is not None and obs_name != f"obs{obs:03d}":
                                continue
                            filtered.append(obs_dir)
                    except (ValueError, IndexError):
                        continue
                return sorted(filtered)
            
            return sorted(obs_dirs)
    except subprocess.TimeoutExpired:
        pass
    return []


def is_band1_file(filepath: Path) -> bool:
    """Check if a .tbl file belongs to band 1."""
    det = parse_det_from_name(filepath.name)
    if det is None:
        return False
    
    # Quick check: detector ID membership
    if det in BAND1_DETS:
        return True
    
    # Fallback: read header for jband metadata
    jband = parse_tbl_jband(filepath)
    return jband == 1


def ensure_link_or_copy(src: Path, dst: Path, mode: str) -> bool:
    """Create symlink, hardlink, or copy. Return True if successful."""
    if dst.exists():
        return True  # Already exists
    
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        if mode == "symlink":
            dst.symlink_to(src)
        elif mode == "hardlink":
            os.link(src, dst)
        elif mode == "copy":
            shutil.copy2(src, dst)
        else:
            return False
        return True
    except (OSError, IOError):
        return False


def process_obs_dir(obs_dir: Path, data_root: Path, out_root: Path, mode: str) -> Tuple[List[Dict], Dict[Tuple[str, str], Dict]]:
    """Process a single obs directory. Returns (manifest_rows, coverage_data)."""
    manifest_rows: List[Dict] = []
    coverage_data: Dict[Tuple[str, str], Dict] = {}
    
    # Parse obs_dir path: plate*/sop*/obs*
    try:
        parts = obs_dir.relative_to(data_root).parts
        if len(parts) < 3:
            return manifest_rows, coverage_data
        plate_name = parts[0]
        sop_name = parts[1]
        obs_name = parts[2]
    except (ValueError, IndexError):
        return manifest_rows, coverage_data

    # List .tbl files in this obs directory
    try:
        tbl_files = sorted(obs_dir.glob("det*.tbl"))
    except (OSError, PermissionError):
        return manifest_rows, coverage_data

    sop_obs_key = (sop_name, obs_name)
    coverage_data[sop_obs_key] = {
        "det_present": set(),
        "det_missing": set(BAND1_DETS)
    }

    for tbl_file in tbl_files:
        if not is_band1_file(tbl_file):
            continue

        det = parse_det_from_name(tbl_file.name)
        if det is not None:
            coverage_data[sop_obs_key]["det_present"].add(det)
            coverage_data[sop_obs_key]["det_missing"].discard(det)

        # Create output path: out_root/sop*/obs*/plate*/det*.tbl
        out_subdir = out_root / sop_name / obs_name / plate_name
        out_file = out_subdir / tbl_file.name

        # Create link/copy
        ensure_link_or_copy(tbl_file, out_file, mode)

        # Add to manifest
        jband = parse_tbl_jband(tbl_file)
        manifest_rows.append({
            "plate": plate_name,
            "sop": sop_name,
            "obs": obs_name,
            "det": det if det is not None else -1,
            "jband": jband,
            "source": str(tbl_file),
            "target": str(out_file)
        })
    
    return manifest_rows, coverage_data


def main() -> None:
    p = argparse.ArgumentParser(
        description="Build band-1 SOP/OBS parallel dataset from obs directories (parallelized)"
    )
    p.add_argument(
        "--data-root",
        default=None,
        help="Root directory containing plate*/sop*/obs*/ (or IRAS_DATA_ROOT)"
    )
    p.add_argument(
        "--out-root",
        default=None,
        help="Output root for reorganized band-1 data (or IRAS_BAND1_REORG_ROOT)"
    )
    p.add_argument(
        "--mode",
        choices=["symlink", "hardlink", "copy"],
        default="symlink",
        help="How to organize output files"
    )
    p.add_argument(
        "--overwrite-index",
        action="store_true",
        help="Regenerate manifest files"
    )
    p.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of worker threads (default: auto-detect from CPU count)"
    )
    p.add_argument(
        "--sop",
        type=int,
        default=None,
        help="Filter to specific SOP number (e.g., 340)"
    )
    p.add_argument(
        "--obs",
        type=int,
        default=None,
        help="Filter to specific OBS number (e.g., 20); requires --sop"
    )
    args = p.parse_args()

    data_root, out_root = resolve_roots(args)
    out_root.mkdir(parents=True, exist_ok=True)

    # Auto-detect worker count
    if args.workers is None:
        args.workers = os.cpu_count() or 4
    print(f"Finding obs directories under {data_root}...", file=sys.stderr)
    obs_dirs = find_obs_dirs(data_root, sop=args.sop, obs=args.obs)
    print(f"Found {len(obs_dirs)} obs directories", file=sys.stderr)
    print(f"Processing with {args.workers} workers...", file=sys.stderr)

    manifest_data: List[Dict] = []
    coverage_data: Dict[Tuple[str, str], Dict] = {}

    # Process obs directories in parallel
    total = len(obs_dirs)
    done = 0
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(process_obs_dir, obs_dir, data_root, out_root, args.mode): obs_dir
            for obs_dir in obs_dirs
        }
        
        for future in as_completed(futures):
            manifest_rows, cov_data = future.result()
            manifest_data.extend(manifest_rows)
            done += 1
            if done % 1000 == 0 or done == total:
                print(f"  [{done}/{total}] {len(manifest_data)} band-1 files so far...", file=sys.stderr, flush=True)
            
            # Merge coverage data
            for sop_obs_key, cov in cov_data.items():
                if sop_obs_key not in coverage_data:
                    coverage_data[sop_obs_key] = {
                        "det_present": set(),
                        "det_missing": set(BAND1_DETS)
                    }
                coverage_data[sop_obs_key]["det_present"].update(cov["det_present"])
                coverage_data[sop_obs_key]["det_missing"].intersection_update(cov["det_missing"])

    print(f"\nProcessed {len(obs_dirs)} obs directories, selected {len(manifest_data)} band-1 files",
          file=sys.stderr)
    print(f"Unique SOP/OBS pairs: {len(coverage_data)}", file=sys.stderr)

    # Write manifests
    manifest_file = out_root / "band1_manifest.csv"
    coverage_file = out_root / "band1_sopobs_coverage.csv"

    if args.overwrite_index or not manifest_file.exists():
        with open(manifest_file, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["plate", "sop", "obs", "det", "jband", "source", "target"]
            )
            writer.writeheader()
            writer.writerows(manifest_data)
        print(f"Wrote {len(manifest_data)} entries to {manifest_file}", file=sys.stderr)

    if args.overwrite_index or not coverage_file.exists():
        coverage_rows: List[Dict] = []
        for (sop, obs), cov in sorted(coverage_data.items()):
            coverage_rows.append({
                "sop": sop,
                "obs": obs,
                "n_det_present": len(cov["det_present"]),
                "det_present": ",".join(str(d) for d in sorted(cov["det_present"])),
                "det_missing": ",".join(str(d) for d in sorted(cov["det_missing"]))
            })

        with open(coverage_file, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["sop", "obs", "n_det_present", "det_present", "det_missing"]
            )
            writer.writeheader()
            writer.writerows(coverage_rows)
        print(f"Wrote {len(coverage_rows)} SOP/OBS pairs to {coverage_file}", file=sys.stderr)

    print(f"\nDone! Output directory: {out_root}", file=sys.stderr)
    print(f"Manifests: {manifest_file.name}, {coverage_file.name}", file=sys.stderr)


if __name__ == "__main__":
    main()
