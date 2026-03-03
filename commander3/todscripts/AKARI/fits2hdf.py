import akari_commander_data_adapter
from cosmoglobe.tod_tools.commander_hdf_writer import CommanderHDFWriter

# Main script - run to create the AKARI HDF files

AKARI_FITS_DIR = '/mn/stornext/d23/cmbco/akari/akari_TSD/www.ir.isas.jaxa.jp/~yamamura/DR2_decompressed/'

OUTPATH = '/mn/stornext/d23/cmbco/globe/akari/tod/eirik_newhdf/'
NSIDE = 2048
NUM_PROCESSES = 64


akari_comm_data_adapter = akari_commander_data_adapter.AKARICommanderDataAdapter(AKARI_FITS_DIR, NSIDE)

comm_todwriter = CommanderHDFWriter(akari_comm_data_adapter)
comm_todwriter.write_hdf_files(OUTPATH, overwrite=True, num_processes=NUM_PROCESSES)
