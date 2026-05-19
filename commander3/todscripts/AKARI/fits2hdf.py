import akari_commander_data_adapter
from cosmoglobe.tod_tools.commander_hdf_writer import CommanderHDFWriter

# Main script - run to create the AKARI HDF files

AKARI_FITS_DIR = '/mn/stornext/d23/cmbco/akari/akari_TSD/www.ir.isas.jaxa.jp/~yamamura/DR2_decompressed/'

OUTPATH = '/mn/stornext/d23/cmbco/globe/akari/tod/eirik_newhdf/'
OUTPATH = '/mn/stornext/d23/cmbco/globe/akari/tod/duncan_testhdf/'
NSIDE = 2048
NUM_SEGMENT_PROCESSES = 60
# Do this if you want error messages
#NUM_SEGMENT_PROCESSES = None
# You can run all together, but it's typically faster to start four different
# instances, one for each band
#BANDS = ['065', '090', '140', '160']
BANDS = ['065']

def run_single_band_write(band):
    akari_comm_data_adapter = akari_commander_data_adapter.AKARICommanderDataAdapter(
        AKARI_FITS_DIR, NSIDE, bands=[band],
        extend_reset_flag=EXTEND_RESET_FLAG)
    comm_todwriter = CommanderHDFWriter(akari_comm_data_adapter)
    comm_todwriter.write_hdf_files(OUTPATH, overwrite=True,
                                   num_processes=NUM_SEGMENT_PROCESSES,
                                   bands=[band],
                                   include_solarcentric_pixpointing=True)

def main():
    for band in BANDS:
        run_single_band_write(band)

if __name__ == '__main__':
    main()
