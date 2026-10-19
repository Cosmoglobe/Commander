from __future__ import annotations

import h5py

PATH_TO_HDF5_FILES = "/mn/stornext/d16/cmbco/bp/gustavbe/master/dirbe_hdf5_files/"
PATH_TO_ZODI_SIMULATION_FILES = "/mn/stornext/d16/cmbco/bp/metins/dirbe/data/zodi_sims/"
TEMP_OUTPUT_PATH = "/mn/stornext/d16/cmbco/bp/metins/dirbe/data/filelists"
DETECTORS = range(1, 11)

DEFAULT_NSIDE = 128
DEFAULT_WEIGHT = 1
DEFUALT_SPIN_AXIS = 0


def write_filelists(output_path: str | None = None, nside: int = DEFAULT_NSIDE) -> None:
    """Writes the DIRBE filelists for Commander3."""

    for detector in DETECTORS:
        hdf5_file = PATH_TO_HDF5_FILES + f"Phot{detector:02}_{nside}.hdf5"

        with h5py.File(hdf5_file, "r") as file:
            n_scan_ids = len(file)

        filelist_name = f"filelist_{detector:02}_{nside}.txt"
        if output_path is not None:
            filelist_name = f"{output_path}/{filelist_name}"

        with open(filelist_name, "w") as file:
            file.write(f"{n_scan_ids}\n")
            for scan_id in range(1, n_scan_ids + 1):
                file.write(
                    f'{scan_id:>10}  "{hdf5_file:<10}" {DEFAULT_WEIGHT:>10} {DEFUALT_SPIN_AXIS:^10} {DEFUALT_SPIN_AXIS:^10}\n'
                )


def write_zodi_simulation_filelists(output_path: str | None = None, nside: int = DEFAULT_NSIDE) -> None:
    """Writes the DIRBE filelists for Commander3."""

    for detector in DETECTORS:
        hdf5_file = PATH_TO_ZODI_SIMULATION_FILES + f"zodi{detector:02}.hdf5"

        with h5py.File(hdf5_file, "r") as file:
            n_scan_ids = len(file)

        filelist_name = f"filelist_sims_{detector:02}_{nside}.txt"
        if output_path is not None:
            filelist_name = f"{output_path}/{filelist_name}"

        with open(filelist_name, "w") as file:
            file.write(f"{n_scan_ids}\n")
            for scan_id in range(1, n_scan_ids + 1):
                file.write(
                    f'{scan_id:>10}  "{hdf5_file:<10}" {DEFAULT_WEIGHT:>10} {DEFUALT_SPIN_AXIS:^10} {DEFUALT_SPIN_AXIS:^10}\n'
                )



def main() -> None:
    nside = 128
    write_zodi_simulation_filelists(output_path=TEMP_OUTPUT_PATH, nside=nside)


if __name__ == "__main__":
    main()
