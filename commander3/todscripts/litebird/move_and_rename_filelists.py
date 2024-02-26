import pathlib
import litebird_sim as lbs

current_filelist_dir = '/mn/stornext/u3/eirikgje/data/litebird_tods/'
imo_location = '/mn/stornext/u3/eirikgje/src/litebird_imo/IMO/'
imo_version='v1.3'
instrument = ['LFT', 'MFT', 'HFT']

imo_db_interface = lbs.Imo(flatfile_location=imo_location)


