import pathlib
import litebird_sim as lbs

should_only_move_filelist = True
suffix = "_tod_eirik_cmb_fg_wn_ncorr30_dipol_v2"
curr_dir = "/mn/stornext/d16/cmbco/litebird/TOD_analysis_temp_eirik/TODS_eirik_cmb_fg_wn_ncorr30_dipol_v2/"
new_dir = "/mn/stornext/d16/cmbco/litebird/TOD_analysis/data_LB/"

def move_only_filelist(old_dir, new_dir, freqname, new_suffix=None,
                       filelist_in_freq_subdir=True,
                       move_to_freq_subdir=False):
    fname = pathlib.Path(old_dir)
    if filelist_in_freq_subdir:
        fname = fname / freqname
    fname = fname / f"filelist_{freqname}.txt"
    new_fname = pathlib.Path(new_dir)
    if move_to_freq_subdir:
        new_fname = new_fname / freqname
    new_fname = new_fname / f"filelist_{freqname}"
    if new_suffix is not None:
        new_fname = new_fname.parent / (new_fname.name + new_suffix)
    new_fname = new_fname.parent / (new_fname.name + '.txt')
    fname.rename(pathlib.Path(new_fname))

imo_location = '/mn/stornext/u3/eirikgje/src/litebird_imo/IMO/'
imo_version='v1.3'
instrument = ['LFT', 'MFT', 'HFT']

imo_db_interface = lbs.Imo(flatfile_location=imo_location)
imo_db_datapath = f"/releases/{imo_version}/satellite"


for inst in instrument:
    instrument_info = imo_db_interface.query(
        f"{imo_db_datapath}/{inst}/instrument_info")
    freqs = instrument_info.metadata['channel_names']
    freqs = ['LB_' + '_'.join(freq.split('-')[::-1]) for freq in freqs]
    nsides = []
    for freq in freqs:
        if should_only_move_filelist:
            move_only_filelist(curr_dir, new_dir, freq, new_suffix=suffix,
                               filelist_in_freq_subdir=True,
                               move_to_freq_subdir=False)
