import litebird_sim as lbs
import healpy
import numpy as np

instruments = ['LFT', 'MFT', 'HFT']

imo_db_interface = lbs.Imo(flatfile_location='/mn/stornext/u3/eirikgje/src/litebird_imo/IMO/')
imo_version = 'v1.3'
imo_db_datapath = f"/releases/{imo_version}/satellite/"

output_dir = '/mn/stornext/u3/eirikgje/data/litebird_sim_analysis/feb_2024_ragnhild_analysis/data/'
lmax = 6000


for instrument in instruments:
    instrument_info = imo_db_interface.query(
        f'{imo_db_datapath}/{instrument}/instrument_info')

    freqs = instrument_info.metadata['channel_names']

    for freq in freqs:
        output_freqname = 'LB_' + '_'.join(freq.split('-')[::-1])
        channel_info = imo_db_interface.query(
            f'{imo_db_datapath}/{instrument}/{freq}/channel_info')
        metadata = channel_info.metadata
        fwhm = metadata['fwhm_arcmin']
        fwhm_rad = fwhm / 60 / 180 * np.pi
        beam = healpy.gauss_beam(fwhm_rad, lmax=lmax)
        healpy.write_cl(f'{output_dir}beam_{output_freqname}_fwhm{fwhm}arcmin.fits', beam)
        print(f"Beam for {output_freqname} generated")
