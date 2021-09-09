from commander_tools.tod_tools.hfi import hfi

import h5py 
from astropy.io import fits
import sqlite3
import os
import numpy as np

planck_dir = '/mn/stornext/d16/cmbco/bp/HFI/hfi_miss03_adc'

flags_file = hfi.init_extra_flags('/mn/stornext/d16/cmbco/bp/HFI/aux/hfi_bad_intervals_15s_elephants.txt')

conn = sqlite3.connect('/mn/stornext/d16/cmbco/bp/HFI/aux/hfi_raw_rings_v3.db')
c = conn.cursor()

for od in range(91, 975):

    print('OD:', od)

    pointFile = os.path.join(planck_dir, str(od).zfill(4), 'pointing-' + str(od).zfill(4) + '.fits')
    pointFits = fits.open(pointFile)
    outFile = h5py.File(os.path.join(planck_dir, str(od).zfill(4), 'extra_flags_' + str(od).zfill(4) + '.h5'), 'w')

    times = pointFits[1].data.field('obt')
    
    for dbentry in c.execute("SELECT * FROM ring_times_hfi WHERE stop_time >= '{0}' AND start_time < '{1}'".format(times[0]/1e9, times[-1]/1e9)):

        ring = dbentry[0]
        start_time = dbentry[2]
        end_time = dbentry[3]

        print(ring)

        startIndex = np.where(times/1e9 > start_time)
        endIndex = np.where(times/1e9 > end_time)

        if len(startIndex[0]) > 0:
            pid_start = startIndex[0][0]
        else:#catch days with no pids
            continue
        if len(endIndex[0]) != 0:
            pid_end = endIndex[0][0]
        else:#catch final pid per od
            pid_end = len(exFile[1].data['obt'])
        if pid_start == pid_end:#catch chunks with no data like od 1007
            continue

        subTimes = times[pid_start:pid_end]
    
        outFile.create_dataset(str(ring).zfill(6)+ '/times', data=subTimes, dtype=np.float64)

        extra_flags = hfi.get_extra_flags(subTimes, 'ALL', flags_file)
        outFile.create_dataset(str(ring).zfill(6) + '/flag_extra', data=extra_flags, dtype=np.uint8)
        flags_353 = hfi.get_extra_flags(subTimes, '353-1', flags_file)

        flags_353 = extra_flags + flags_353
        flags_353[flags_353 == 256] = 128
        outFile.create_dataset(str(ring).zfill(6) + '/flags_extra_353-1', data=flags_353, dtype=np.uint8)

    pointFits.close()
    outFile.close()
