from cosmoglobe.tod_tools import CommanderHDFWriter
from cosmoglobe.tod_tools.iras.comm_iras_data_adapter import CommIRASDataAdapter

iras_adapter = CommIRASDataAdapter()
comm_writer = CommanderHDFWriter(iras_adapter)

comm_writer.write_hdf_files('/mn/stornext/u3/eirikgje/data/hdfwriter_test/',
                            overwrite=False, num_processes=128)
