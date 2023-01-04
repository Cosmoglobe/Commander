import h5py
import numpy as np
import litebird_sim as lbs
from commander_tools.tod_tools.litebird_imo import LitebirdImo

fnames = {
    'LFT': 'LFT_instrument.h5',
    'MFT': 'MFT_instrument.h5',
    'HFT': 'HFT_instrument.h5'
}

imo_version = 'v1.3'
imo_db_interface = lbs.Imo()
for instrument in ('LFT', 'MFT', 'HFT'):
    imo = LitebirdImo(imo_db_interface, imo_version, instrument) 
    f = h5py.File(fnames[instrument], 'w')
    channels = imo.get_channel_names()
#    center_freqs = []
    for channel in channels:
        center_freq = imo.get_detector_frequency(channel)
        bandwidth = imo.get_detector_bandwidth(channel)
        fwhm = imo.get_detector_fwhm(channel)
#        center_freqs.append(imo.get_detector_frequency(channel))
        grp = f.create_group(channel)
        yval = grp.create_dataset("bandpass", (3,), dtype='f')
        xval = grp.create_dataset("bandpassx", (3,), dtype='f')
        #Tophat
        yval[:] = 1.
        xval[:] = np.array([center_freq - bandwidth/2,
                            center_freq,
                            center_freq + bandwidth / 2])
        for detector in imo.get_channel_dets(channel):
            grp = f.create_group(detector)
            yval = grp.create_dataset("bandpass", (3,), dtype='f')
            xval = grp.create_dataset("bandpassx", (3,), dtype='f')
            #Tophat
            yval[:] = 1.
            xval[:] = np.array([center_freq - bandwidth/2,
                                center_freq,
                                center_freq + bandwidth / 2])
            b = grp.create_group("beam")
            bT = b.create_dataset("T", (1,), dtype="f")
            bT[:] = 0.
            bB = b.create_dataset("B", (1,), dtype="f")
            bB[:] = 0
            bE = b.create_dataset("E", (1,), dtype="f")
            bE[:] = 0
            blmax = grp.create_dataset("beamlmax", (1,), dtype='f')
            blmax[0] = 2400
            bmmax = grp.create_dataset('beammax', (1,), dtype='f')
            bmmax[0] = 100
            el = grp.create_dataset("elip", (1,), dtype='f')
            el[0] = 1.
            fw = grp.create_dataset("fwhm", (1,), dtype='f')
            fw[0] = fwhm
            psi = grp.create_dataset("psi_ell", (1,), dtype='f')
            psi[0] = 1.
            sl = grp.create_group("sl")
            slT = sl.create_dataset("T", (1,), dtype="f")
            slT[:] = 0
            slB = sl.create_dataset("B", (1,), dtype="f")
            slB[:] = 0
            slE = sl.create_dataset("E", (1,), dtype="f")
            slE[:] = 0
            sll= grp.create_dataset("sllmax", (1,), dtype='f')
            sll[0]=512
            slm= grp.create_dataset("slmmax", (1,), dtype='f')
            slm[0]=100
            eff= grp.create_dataset("mbeam_eff", (1,), dtype='f')
            eff[0]=1.

    f.close()
