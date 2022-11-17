from commander_tools.tod_tools.litebird import litebird
from commander_tools.tod_tools.litebird_imo import LitebirdImo
from commander_tools.tod_tools import commander_tod as tod
import litebird_sim as lbs
from commander_tools.tod_tools.litebird_native_tod_reader import LitebirdTodReader
import argparse
import multiprocessing as mp
import os
import numpy as np
import h5py

def main():
    parser = argparse.ArgumentParser()
    outpath = '/home/eirik/data/litebird_sims/'
    name = 'litebird'
    version = 'v1.3'
    dicts = None
    overwrite = True

#    pool = mp.Pool(processes=num_procs)
#    manager = mp.Manager()
    manager_dicts = {}

    comm_tod = tod.commander_tod(outpath, name, version, manager_dicts, overwrite)
    make_od(comm_tod, 'L1-060', 1, None, version)

    
def make_od(comm_tod, freq, od, args, imo_version):
    if freq[0] == 'L':
        instrument = 'LFT'
    elif freq[0] == 'M':
        instrument = 'MFT'
    elif freq[0] == 'H':
        instrument = 'HFT'
    else:
        raise ValueError('Frequency format not recognized - must start with L, M or H')
    comm_tod.init_file(int(freq[-3:]), od, mode='w')
    imo_db_interface = lbs.Imo()
    imo = LitebirdImo(imo_db_interface, imo_version, instrument) 
    dets = imo.get_channel_dets(freq)
    fsamps = [imo.get_detector_property(freq, det, 'sampling_rate_hz') for det in dets]
    nside = 512
    for fsamp in fsamps:
        assert(fsamp == fsamps[0])

#    chunk_size = 2 ** 22
    chunk_size = 2 ** 16
    scan_time = chunk_size / fsamps[0] # in seconds. About 2.6 days when the sampling frequency is 19 Hz

    if instrument in ('LFT', 'HFT'):
        polang = [0 + 45 * int(it.split('_')[3][0] == 'U') + 90 * int(it.split('_')[3][1] == 'B') for it in dets]
    else:
        polang = [int(it.split('_')[3][:1]) + 90 * int(it.split('_')[3][2] == 'B') for it in dets]

    prefix ='/common'
    comm_tod.add_field(prefix + '/fsamp', fsamps)
    comm_tod.add_field(prefix + '/nside', [nside])
    comm_tod.add_field(prefix + '/det', np.string_(', '.join(dets)))
    comm_tod.add_field(prefix + '/polang', polang)
    comm_tod.add_field(prefix + '/mbang', [0.0]*len(dets))
    for field in ('tod', 'pointings', 'psi'):
        print(f"Field: {field}")
        with LitebirdTodReader('/home/eirik/data/litebird_sims/',
                               'LFT_L1-060_obs_cmb',
                               start_day=100,
                               end_day=105,
                               field=field) as reader:
            for curr_scan_id, chunk in enumerate(reader.iterate_chunks(chunk_size)):
                print(f"Curr scan id: {curr_scan_id}")
                for i, detector in enumerate(dets):
#                    print(f"Detector: {detector}")
                    prefix = f'{curr_scan_id:06}/{detector}' 
                    if field == 'pointings':
                        comm_tod.add_field(f'{prefix}/theta', chunk[i, :, 0])
                        comm_tod.add_field(f'{prefix}/phi', chunk[i, :, 1])
                    elif field == 'tod':
                        comm_tod.add_field(f'{prefix}/tod', chunk[i, :])
                    elif field == 'psi':
                        comm_tod.add_field(f'{prefix}/psi', chunk[i, :])

if __name__ == '__main__':
    main()
