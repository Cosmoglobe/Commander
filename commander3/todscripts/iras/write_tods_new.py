from cosmoglobe.tod_tools import CommanderHDFWriter
from cosmoglobe.tod_tools.iras.comm_iras_data_adapter import CommIRASDataAdapter

print('Initializing IRAS Data adapter')
iras_adapter = CommIRASDataAdapter()
print('Initializing comm_writer')
comm_writer = CommanderHDFWriter(iras_adapter)

print('Writing hdf files')
OUTDIR = '/mn/stornext/d5/data/duncanwa/IRAS/hdf_writer_test'
comm_writer.write_hdf_files(OUTDIR,
                            overwrite=True, num_processes=72)
