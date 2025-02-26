import pathlib
import litebird_sim as lbs
import healpy
import shutil

#current_filelist_dir = '/mn/stornext/u3/eirikgje/data/litebird_tods/'
imo_location = '/mn/stornext/u3/eirikgje/src/litebird_imo/IMO/'
imo_version='v1.3'
instrument = ['LFT', 'MFT', 'HFT']

imo_db_interface = lbs.Imo(flatfile_location=imo_location)

data_dir = '/mn/stornext/d16/cmbco/litebird/TOD_analysis/data_LB/'
tod_dir = '/mn/stornext/d16/cmbco/litebird/TOD_analysis_temp_eirik/TODS_eirik_cmb_fg_wn_ncorr30_dipol/'
imo_db_datapath = f"/releases/{imo_version}/satellite"
defaults_path = "/mn/stornext/u3/eirikgje/src/Commander/commander3/parameter_files/defaults/bands/LiteBIRD/"

generate_map_rms = False
regenerate_defaults = False
detector_list_rename = False
should_rename_beams = False
should_rename_bp_init_files = False
should_rename_tod_h5_files = False

def generate_map_rms_files(old_freqnames, new_freqnames):
    for freq, chan in zip(freqs, chan_names):
        mapfilename_old = f'tod_{freq}_map_c0001_k000002.fits'
        noisefilename_old = f'tod_{freq}_rms_c0001_k000002.fits'
        mapfilename_new_512 = f'tod_{chan}_map_c0001_k000002_n0512.fits'
        noisefilename_new_512 = f'tod_{chan}_rms_c0001_k000002_n0512.fits'
        mapfilename_new_1024 = f'tod_{chan}_map_c0001_k000002_n1024.fits'
        noisefilename_new_1024 = f'tod_{chan}_rms_c0001_k000002_n1024.fits'
        smoothfilename_512 = f'{chan}_rms_90arc_n0512.fits'
        smoothfilename_1024 = f'{chan}_rms_90arc_n1024.fits'
        mapfilepath_old = pathlib.Path(data_dir + mapfilename_old)
        mapfilepath_new_512 = pathlib.Path(data_dir + mapfilename_new_512)
        mapfilepath_new_1024 = pathlib.Path(data_dir + mapfilename_new_1024)
        noisefilepath_old = pathlib.Path(data_dir + noisefilename_old)
        noisefilepath_new_512 = pathlib.Path(data_dir + noisefilename_new_512)
        noisefilepath_new_1024 = pathlib.Path(data_dir + noisefilename_new_1024)
        smoothfilepath_512 = pathlib.Path(data_dir + smoothfilename_512)
        smoothfilepath_1024 = pathlib.Path(data_dir + smoothfilename_1024)
        if mapfilepath_old.exists():
            shutil.copy(mapfilepath_old, data_dir + mapfilename_new_512)
            shutil.copy(noisefilepath_old, data_dir + noisefilename_new_512)
        if not mapfilepath_new_1024.exists():
            mapfile_512 = healpy.read_map(mapfilepath_new_512, field=None)
            mapfile_1024 = healpy.ud_grade(mapfile_512, 1024)
            healpy.write_map(mapfilepath_new_1024, mapfile_1024, overwrite=True)
            print(f"Converted {chan} map to 1024")
        if not noisefilepath_new_1024.exists():
            noisefile_512 = healpy.read_map(noisefilepath_new_512, field=None)
            noisefile_1024 = healpy.ud_grade(noisefile_512, 1024, power=1)
            healpy.write_map(noisefilepath_new_1024, noisefile_1024, overwrite=True)
            print(f"Converted {chan} rms to 1024")
        if not smoothfilepath_1024.exists() and smoothfilepath_512.exists():
            smoothfile_512 = healpy.read_map(smoothfilepath_512, field=None)
            smoothfile_1024 = healpy.ud_grade(smoothfile_512, 1024, power=1)
            healpy.write_map(smoothfilepath_1024, smoothfile_1024, overwrite=True)
            print(f"Converted {chan} smooth rms to 1024")


def rename_detector_lists(inst, old_freqnames, new_freqnames):
    for freq, chan in zip(old_freqnames, new_freqnames):
        fname_old = pathlib.Path(data_dir + f'detectors_{inst}_{freq}_T+B.txt')
        fname_old_4 = pathlib.Path(data_dir + f'detectors_{inst}_{freq}_T+B_4.txt')
        fname_new = pathlib.Path(data_dir + f'detectors_{chan}_T+B.txt')
        fname_new_4 = pathlib.Path(data_dir + f'detectors_{chan}_T+B_4.txt')
        if not fname_new.exists():
            shutil.copy(fname_old, fname_new)
        if not fname_new_4.exists():
            shutil.copy(fname_old_4, fname_new_4)


def regenerate_defaults_files(inst, old_freqnames, new_freqnames, nsides):
    for freq, chan, nside in zip(freqs, chan_names, nsides):
        old_defaultsfile = pathlib.Path(defaults_path + freq + '_TOD.defaults')
        new_defaultsfile = pathlib.Path(defaults_path + chan + '_TOD.defaults')
        with open(old_defaultsfile, 'r') as old_default, open(new_defaultsfile, 'w') as new_default:
            for line in old_default:
                line = line.replace(freq, chan)
                if not line.startswith('BAND_TOD_RIMO') and not line.startswith('BAND_BANDPASSFILE'):
                    line = line.replace(inst+'_', '')
                if nside == 1024:
                    line = line.replace('0512', '1024')
                    line = line.replace('512', '1024')
                if line.startswith('BAND_MAPFILE'):
                    line = line.split('=')[:-1]
                    line.append(f' tod_{chan}_map_c0001_k000002_n{nside:04d}.fits\n')
                    line = "=".join(line)
                if line.startswith('BAND_NOISEFILE'):
                    line = line.split('=')[:-1]
                    line.append(f' tod_{chan}_rms_c0001_k000002_n{nside:04d}.fits\n')
                    line = "=".join(line)
                if line.startswith('BAND_NOISE_RMS') and "SMOOTH01" in line and not "none" in line:
                    line = line.split('=')[:-1]
                    line.append(f' {chan}_rms_90arc_n{nside:04d}.fits\n')
                    line = "=".join(line)
                if line.startswith("BAND_TOD_FILELIST"):
                    line = line.split('=')[:-1]
                    line.append(f' filelist_{chan}.txt\n')
                    line = "=".join(line)
                if line.startswith("BAND_BEAM_B_L_FILE"):
                    line = line.split('=')[:-1]
                    if chan == 'LB_402_H3':
                        line.append(f' {chan}_IMoV1_beam_ext.fits\n')
                    else:
                        line.append(f' {chan}_IMoV1_beam.fits\n')
                    line = "=".join(line)
                new_default.write(line)


def rename_beams(old_freqs, new_freqs):
    for freq, chan in zip(old_freqs, new_freqs):
        fname_old = pathlib.Path(data_dir + f'LB_IMoV1_beam_{freq}.fits')
        fname_new = pathlib.Path(data_dir + f'{chan}_IMoV1_beam.fits')
        if not fname_new.exists():
            shutil.copy(fname_old, fname_new)


def rename_bp_init_files(inst, old_freqs, new_freqs):
    for freq, chan in zip(old_freqs, new_freqs):
        fname_old = pathlib.Path(data_dir + f'bp_init_{inst}_{freq}.dat')
        fname_new = pathlib.Path(data_dir + f'bp_init_{chan}.dat')
        if not fname_new.exists():
            shutil.copy(fname_old, fname_new)


def rename_tod_h5_files(freqs):
    for chan in freqs:
        print(f"Renaming tod files for {chan}")
        k = 1
        decon_chan = chan.split('_')
        print(decon_chan)
        old_chan = f'{decon_chan[0]}_{decon_chan[1]}-{decon_chan[2]}'
        fname_old = pathlib.Path(tod_dir + f'{chan}/{old_chan}_{k:06d}.h5')
        while fname_old.exists():
            fname_old.rename(tod_dir + f'{chan}/{chan}_{k:06d}.h5')
            k += 1
            fname_old = pathlib.Path(tod_dir + f'{chan}/{old_chan}_{k:06d}.h5')
        print(k)



for inst in instrument:
    instrument_info = imo_db_interface.query(
        f"{imo_db_datapath}/{inst}/instrument_info")
    freqs = instrument_info.metadata['channel_names']
    chan_names = ['LB_' + '_'.join(freq.split('-')[::-1]) for freq in freqs]
    nsides = []
    for freq in freqs:
        channel_info = imo_db_interface.query(f"{imo_db_datapath}/{inst}/{freq}/channel_info")
        if channel_info.metadata['fwhm_arcmin'] > 30:
            nsides.append(512)
        else:
            nsides.append(1024)
    for chan, nside in zip(chan_names, nsides):
        print(f"{chan}, {nside}")
    if generate_map_rms:
        generate_map_rms_files(freqs, chan_names)
    if regenerate_defaults:
        regenerate_defaults_files(inst, freqs, chan_names, nsides)
    if detector_list_rename:
        rename_detector_lists(inst, freqs, chan_names)
    if should_rename_beams:
        rename_beams(freqs, chan_names)
    if should_rename_bp_init_files:
        rename_bp_init_files(inst, freqs, chan_names)
    if should_rename_tod_h5_files:
        rename_tod_h5_files(chan_names)
