# Changelog

All notable changes to the Commander3 project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

#### CMake:
- Added `COMM3_BACKEND` variable to identify what BLAS/LAPACK & FFT
  implementation to use. Possible values are `aocl`, `mkl`, `opensrc` and
  `any`. If "any" was chosen, the choice would be made automatically based 
  on the processor Vendor (Intel, AMD, Unknown).
- Added detection of SSE3, SSE4_1, AVX, AVX2 etc. CPU features 
- Added search for, compilation and linking of AMD BLIS 3.1
- Added search for, compilation and linking of AMD FLAME 3.1
- Added compilation and linking of AMD FFTW 3.1 with both CMake and `configure` 
  (the former code is commented out since it is very slow, but can be used if needed)
- Added search for MKL based on the 3d party `FindMKL.cmake` package (there is also 
  option to use `$MKLROOT` and `MKLConfig.cmake` shipped with Intel OneAPI MKL installation)

### Changed

#### CMake:
- Changed OpenBLAS version to 0.3.20.
- Changed HDF5 version to 1.12.2

#### Documentation
- Reformatted the Commander Installation Chapter into more logical structure:
  Installing Prerequisites => Automatic Installation => Manual Installation =>
  Precompiled Binaries (Docker) => Files Downloader => Running Commander =>
  Extending Commander 
- Reformatted FAQ section
- Reformatted README.md to point to GitHub pages instead of BeyondPlanck
- Reformatted File Formats section

### Fixed

#### CMake
- Fixed bug associated with libaec installation while it tried to create a symlink

### Removed 

#### CMake:
- Removed OpenMP from all LAPACK/BLAS libraries because we are not using those
and because it was causing issues when compiling OpenBLAS on owls with more
than 1 thread.
- Removed compilation and linking of FFTW when MKL is present
