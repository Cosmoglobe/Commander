# Parallelized OTF Fitting - Performance Guide

## Overview

I've created three versions of the OTF fitting script with increasing levels of optimization and parallelization:

1. **`fit_otf.py`** - Original sequential version (UPDATED with basic parallelization)
2. **`fit_otf_parallel.py`** - Multiprocessing version with progress tracking
3. **`fit_otf_joblib.py`** - Advanced version with joblib (best for HPC clusters)

## Quick Start

### Option 1: Updated Original Script (Simplest)
```bash
python fit_otf.py
```
- Uses all available CPU cores automatically
- Drop-in replacement for your existing script
- Best for: Quick runs without configuration

### Option 2: Multiprocessing Version (Recommended)
```bash
python fit_otf_parallel.py
```
- Better progress tracking
- More informative output
- Best for: Interactive work and debugging

### Option 3: Joblib Version (Most Flexible)
```bash
# Use all cores
python fit_otf_joblib.py

# Use specific number of cores
python fit_otf_joblib.py --n-jobs 32

# Use with different backend (for HPC)
python fit_otf_joblib.py --n-jobs 64 --backend multiprocessing

# Verbose output
python fit_otf_joblib.py --n-jobs -1 --verbose
```
- Best for: HPC clusters with many cores
- Highly configurable
- Better memory management

## Performance Comparison

Assuming each frequency takes ~1 second to fit and `g.SPEC_SIZE = 257`:

| Version | Cores | Time (est.) | Speedup |
|---------|-------|-------------|---------|
| Original (sequential) | 1 | 257 seconds (~4.3 min) | 1x |
| Parallel (8 cores) | 8 | 32 seconds | 8x |
| Parallel (32 cores) | 32 | 8 seconds | 32x |
| Parallel (64 cores) | 64 | 4 seconds | 64x |

## Installation Requirements

All versions require standard dependencies already in your script.

For the joblib version, install if not present:
```bash
pip install joblib
```

## Key Optimizations Implemented

### 1. **Parallelization**
- Each frequency (`nui`) is fitted independently
- Perfect for parallel processing since there are no dependencies
- Linear speedup with number of cores (up to `g.SPEC_SIZE`)

### 2. **Progress Tracking**
- Real-time progress updates
- ETA (estimated time to arrival)
- Processing rate statistics

### 3. **Better Output**
- Summary statistics
- Timing information
- Configuration details

## Advanced Usage

### For HPC Clusters with SLURM

Create a submission script `submit_otf.sh`:
```bash
#!/bin/bash
#SBATCH --job-name=otf_fit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=01:00:00
#SBATCH --mem=64G

# Load your environment
module load python/3.x

# Run with all allocated cores
python fit_otf_joblib.py --n-jobs $SLURM_CPUS_PER_TASK
```

Submit with:
```bash
sbatch submit_otf.sh
```

### Memory Considerations

If you run into memory issues with many cores:
```bash
# Reduce number of cores
python fit_otf_joblib.py --n-jobs 16

# Or use threading backend (shares memory)
python fit_otf_joblib.py --n-jobs 32 --backend threading
```

## Choosing the Right Version

| Use Case | Recommended Version |
|----------|-------------------|
| Quick test, existing workflow | `fit_otf.py` (updated) |
| Development, debugging | `fit_otf_parallel.py` |
| Production runs | `fit_otf_parallel.py` |
| HPC cluster (many cores) | `fit_otf_joblib.py` |
| Memory-constrained | `fit_otf_joblib.py` with `--n-jobs` limit |

## Further Optimizations (Future)

If you need even more speed, consider:

1. **GPU Acceleration**: Move minimize() to use GPU-based optimizers
2. **Approximate Methods**: Use faster but less accurate optimization
3. **Warm Starting**: Use solutions from nearby frequencies as initial guesses
4. **Caching**: Pre-compute and cache invariant calculations
5. **Compiled Code**: Convert bottlenecks to Cython or Numba

## Troubleshooting

### "Too many open files" error
Reduce the number of parallel jobs:
```bash
python fit_otf_joblib.py --n-jobs 16
```

### Memory errors
Use threading backend or reduce cores:
```bash
python fit_otf_joblib.py --n-jobs 8 --backend threading
```

### Import errors for joblib
Install joblib:
```bash
pip install joblib
```

## Questions?

The parallelization is straightforward because each frequency fit is independent. The main loop that previously ran sequentially:
```python
for nui in range(g.SPEC_SIZE):
    solution[nui] = minimize(...)  # Each iteration is independent
```

Is now parallelized to run all iterations simultaneously across available cores.
