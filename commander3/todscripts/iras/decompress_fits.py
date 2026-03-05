#!/usr/bin/env python3
"""
Decompress raw IRAS FITS files by applying the IRAS difference table.

This script performs true decompression of IRAS survey data by:
1. Reading differenced data (indices into difference table) from FITS files
2. Building the IRAS logarithmic compression difference table
3. Undifferencing the data using cumulative sum
4. Writing fully decompressed detector data in FITS or NPY format

This replaces the legacy 'blowup' tool and uses modern Python libraries
with multiprocessing support for parallel decompression.

Reuses decompression logic from extract_timestream.py for consistency.

Usage:
    python decompress_fits.py [--sop SOP_NUM] [--output-dir DIR] [--format FORMAT] [--workers N]
    
    --sop SOP_NUM:        Only decompress specific SOP (e.g., 29). Default: all
    --output-dir DIR:     Output directory for decompressed files. Default: ./decompressed/
    --format FORMAT:      Output format: 'fits' (FITS with detector values) or 'npy' (numpy binary)
    --raw-root PATH:      Path to raw data root. Default: /home/dwatts/IRAS/disk*
    --workers N:          Number of parallel workers. Default: CPU count

Examples:
    # Decompress SOP029 using all CPU cores
    python decompress_fits.py --sop 29
    
    # Decompress all SOPs with 8 parallel workers to detector values
    python decompress_fits.py --workers 8
    
    # Use compact NPY format with 4 workers
    python decompress_fits.py --format npy --workers 4 --output-dir ~/detector_data/
"""

import argparse
import numpy as np
from pathlib import Path
from tqdm import tqdm
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import os

# Import the decompression functions from extract_timestream.py
sys.path.insert(0, str(Path(__file__).parent))
from extract_timestream import generate_table, build_diftab, parse_sop_header


def decompress_fits_file(input_path, output_path, output_format='fits'):
    """
    Decompress a single FITS file by undifferencing the data.
    
    Uses the exact same file reading and decompression logic as extract_timestream.py
    to ensure consistency and proper handling of raw binary FITS data.
    
    Args:
        input_path: Path to input FITS file
        output_path: Path to output file
        output_format: 'fits' or 'npy'
        
    Returns:
        (success: bool, shape: tuple, input_path: str, output_path: str)
    """
    try:
        # Use the exact same file reading approach as extract_timestream.py
        with open(input_path, "rb") as f:
            header = parse_sop_header(f)
            naxis1 = int(header.get("NAXIS1", 0))
            naxis2 = int(header.get("NAXIS2", 1))
            
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
            
            # Read all data as flat array (uint8 for compressed data)
            total_samples = naxis1 * naxis2
            raw_data = np.fromfile(f, dtype='<u1', count=total_samples)
        
        # Reshape data (detector-major order, matching extract_timestream.py)
        if naxis2 > 1:
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
        
        # Decompress using extract_timestream logic
        diftab = build_diftab()
        if data.ndim == 2:
            # Decompress each detector separately
            decompressed = np.zeros_like(data, dtype=np.int32)
            for det in range(data.shape[1]):
                init = header.get(f"INIT00{det+1:02}", 0)
                decompressed[0, det] = init
                for i in range(1, data.shape[0]):
                    diff_idx = data[i, det]
                    decompressed[i, det] = decompressed[i-1, det] + diftab[diff_idx]
        else:
            # Decompress 1D data
            decompressed = np.zeros_like(data, dtype=np.int32)
            decompressed[0] = header.get("INIT0001", 0)
            for i in range(1, data.shape[0]):
                diff_idx = data[i]
                decompressed[i] = decompressed[i-1] + diftab[diff_idx]
        
        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Write output
        if output_format == 'fits':
            # Write as FITS with decompressed data, preserving metadata and order from input
            import fitsio
            
            # Read input FITS with fitsio to get properly-typed headers
            with fitsio.FITS(str(input_path)) as fits_in:
                input_header = fits_in[0].read_header()
            
            # Keys to preserve with their order from the original file
            # (skip compression-specific keys: BSCALE, BZERO, BLANK)
            preserve_keys_ordered = [
                # Coordinate system (in original order)
                'CRVAL1', 'CRPIX1', 'CDELT1', 'CTYPE1',
                'CRVAL2', 'CRPIX2', 'CDELT2', 'CTYPE2',
                # Data info
                'BUNIT',
                # Source/Observatory
                'ORIGIN', 'INSTRUME', 'TELESCOP',
                # Attitude/Orientation
                'PSI', 'THETA', 'PSIRATE',
                # Observation parameters
                'SOP', 'SATCAL', 'OBS', 'ATT', 'SCANTYPE',
                # History
                'HISTORY'
            ]
            
            # Build output header in the specified order with proper types
            output_fits_header = fitsio.FITSHDR()
            for key in preserve_keys_ordered:
                if key in input_header:
                    value = input_header[key]
                    output_fits_header[key] = value
            
            # Write FITS file with preserved headers in original order
            # fitsio will automatically handle NAXIS, NAXIS1, NAXIS2 based on array shape
            fitsio.write(str(output_path), decompressed, header=output_fits_header, clobber=True)
            
        elif output_format == 'npy':
            # Write as numpy binary (more compact)
            np.save(output_path, decompressed)
        else:
            raise ValueError(f"Unknown format: {output_format}")
        
        return True, decompressed.shape, str(input_path), str(output_path)
    
    except Exception as e:
        return False, str(e), str(input_path), str(output_path)


def decompress_worker(args):
    """
    Worker function for parallel decompression.
    Takes tuple of (input_path, output_path, output_format).
    """
    input_path, output_path, output_format = args
    return decompress_fits_file(input_path, Path(output_path), output_format)


def decompress_sop(sop_num=None, output_dir=None, output_format='fits', raw_root=None, workers=None):
    """
    Decompress all FITS files for a given SOP (or all SOPs).
    
    Args:
        sop_num: SOP number (e.g., 29) or None for all
        output_dir: Output directory path
        output_format: 'fits' or 'npy'
        raw_root: Root directory for raw data
        workers: Number of parallel workers (default: CPU count)
    """
    if output_dir is None:
        output_dir = Path("./decompressed/")
    else:
        output_dir = Path(output_dir)
    
    if raw_root is None:
        # Find all disk directories
        raw_dirs = []
        for disk_dir in Path("/home/dwatts/IRAS").glob("disk*"):
            if disk_dir.is_dir() and not str(disk_dir).endswith(".tgz"):
                raw_dirs.append(disk_dir)
    else:
        raw_dirs = [Path(raw_root)]
    
    # Find all SOP directories
    sop_pattern = f"sop.{sop_num:02d}_" if sop_num is not None else "sop.*"
    
    fits_files = []
    for raw_dir in raw_dirs:
        for survey_dir in raw_dir.glob("d*/*.b?/"):
            for sop_dir in survey_dir.glob(sop_pattern):
                for fits_file in sop_dir.glob("*"):
                    if fits_file.is_file():
                        fits_files.append(fits_file)
    
    fits_files = sorted(fits_files)
    
    if not fits_files:
        print(f"No FITS files found matching pattern: {sop_pattern}")
        return 0
    
    print(f"Found {len(fits_files)} FITS files to decompress")
    
    # Prepare work items: (input_path, output_path, format)
    work_items = []
    for fits_file in fits_files:
        # Construct output path maintaining directory structure
        rel_path = fits_file.relative_to(fits_file.parent.parent.parent.parent)
        output_path = output_dir / rel_path
        
        if output_format == 'fits':
            output_path = output_path.with_suffix('.fits')
        else:
            output_path = output_path.with_suffix('.npy')
        
        work_items.append((fits_file, output_path, output_format))
    
    # Determine number of workers
    if workers is None:
        workers = os.cpu_count() or 4
    
    print(f"Using {workers} parallel workers")
    
    # Decompress with parallel processing
    successful = 0
    failed = 0
    failed_files = []
    
    start_time = time.time()
    
    with ProcessPoolExecutor(max_workers=workers) as executor:
        # Submit all tasks
        futures = {executor.submit(decompress_worker, item): item for item in work_items}
        
        # Process completed tasks with progress bar
        with tqdm(as_completed(futures), total=len(futures), desc="Decompressing") as pbar:
            for future in pbar:
                try:
                    success, result, input_path, output_path = future.result()
                    
                    if success:
                        successful += 1
                    else:
                        failed += 1
                        failed_files.append((input_path, result))
                        
                except Exception as e:
                    failed += 1
                    item = futures[future]
                    failed_files.append((str(item[0]), str(e)))
    
    elapsed = time.time() - start_time
    
    print(f"\nResults:")
    print(f"  Successful: {successful}")
    print(f"  Failed: {failed}")
    print(f"  Time elapsed: {elapsed:.1f} seconds")
    if successful > 0:
        print(f"  Throughput: {successful/elapsed:.1f} files/second")
    print(f"  Output directory: {output_dir}")
    
    if failed_files:
        print(f"\nFailed files:")
        for inp, err in failed_files[:10]:  # Show first 10 failures
            print(f"  {Path(inp).name}: {err}")
        if len(failed_files) > 10:
            print(f"  ... and {len(failed_files) - 10} more")
    
    return successful


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Decompress IRAS FITS files")
    parser.add_argument("--sop", type=int, default=None, help="SOP number (e.g., 29)")
    parser.add_argument("--output-dir", default="./decompressed/", help="Output directory")
    parser.add_argument("--format", choices=['fits', 'npy'], default='fits', help="Output format")
    parser.add_argument("--raw-root", default=None, help="Raw data root directory")
    parser.add_argument("--workers", type=int, default=None, help="Number of parallel workers (default: CPU count)")
    
    args = parser.parse_args()
    
    count = decompress_sop(
        sop_num=args.sop,
        output_dir=args.output_dir,
        output_format=args.format,
        raw_root=args.raw_root,
        workers=args.workers
    )
    
    sys.exit(0 if count > 0 else 1)
