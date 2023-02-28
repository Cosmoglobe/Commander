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
    outpath = '/mn/stornext/u3/eirikgje/data/litebird_tods/'
#    outpath = '/home/eirik/data/litebird_sims/'
    name = 'litebird'
    version = '1.0'
    imo_version = 'v1.3'
    dicts = {60: {}}
    overwrite = True

#    pool = mp.Pool(processes=num_procs)
#    manager = mp.Manager()
#    manager_dicts = {}

#    comm_tod = tod.commander_tod(outpath, name, version, manager_dicts, overwrite)
    comm_tod = tod.commander_tod(outpath, name, version, dicts, overwrite)
    make_ods(comm_tod, 'L1-060', None, imo_version)


def create_new_tod(comm_tod, od, freq, fsamps, nside, dets, polang):
    comm_tod.init_file(int(freq[-3:]), od, mode='w')
    prefix ='/common'
    comm_tod.add_field(prefix + '/fsamp', fsamps)
    comm_tod.add_field(prefix + '/nside', [nside])
    comm_tod.add_field(prefix + '/det', np.string_(', '.join(dets)))
    comm_tod.add_field(prefix + '/polang', polang)
    comm_tod.add_field(prefix + '/mbang', [0.0]*len(dets))
    return comm_tod

    
def make_ods(comm_tod, freq, args, imo_version,
             simpath='/mn/stornext/d16/cmbco/bp/mathew/litebird/sim0000/LFT_L1-060/tods/'):
#             simpath='/home/eirik/data/litebird_sims/'):
    if freq[0] == 'L':
        instrument = 'LFT'
    elif freq[0] == 'M':
        instrument = 'MFT'
    elif freq[0] == 'H':
        instrument = 'HFT'
    else:
        raise ValueError('Frequency format not recognized - must start with L, M or H')
    imo_db_interface = lbs.Imo()
    imo = LitebirdImo(imo_db_interface, imo_version, instrument) 
    dets = imo.get_channel_dets(freq)
    fsamps = [imo.get_detector_property(freq, det, 'sampling_rate_hz') for det in dets]
    nside = 512
    for fsamp in fsamps:
        assert(fsamp == fsamps[0])

    chunk_size = 2 ** 22
#    chunk_size = 2 ** 16
    scan_time = chunk_size / fsamps[0] # in seconds. About 2.6 days when the sampling frequency is 19 Hz

    if instrument in ('LFT', 'HFT'):
        polang = [0 + 45 * int(it.split('_')[3][0] == 'U') + 90 * int(it.split('_')[3][1] == 'B') for it in dets]
    else:
        polang = [int(it.split('_')[3][:1]) + 90 * int(it.split('_')[3][2] == 'B') for it in dets]

    with LitebirdTodReader(simpath,
                           ['LFT_L1-060_obs_cmb',
                            'LFT_L1-060_obs_1_over_f_noise_realistic',
                            'LFT_L1-060_obs_dipole_total',
                            'LFT_L1-060_obs_fg',
                            'LFT_L1-060_obs_w_noise'],
                           start_day=0,
                           end_day=364, 
                           fields=['tod', 'pointings', 'psi']) as reader:
        for curr_scan_id, chunk in enumerate(reader.iterate_chunks(chunk_size)):
            print(f"Curr scan id: {curr_scan_id}")
            comm_tod = create_new_tod(comm_tod, curr_scan_id+1, freq, fsamps, nside, dets, polang)
            for field_idx, field in enumerate(('tod', 'pointings', 'psi')):
                for i, detector in enumerate(dets):
#                    print(f"Detector: {detector}")
                    prefix = f'{curr_scan_id:06}/{detector}' 
                    if field == 'pointings':
                        comm_tod.add_field(f'{prefix}/theta', chunk[field_idx][i, :, 0])
                        comm_tod.add_field(f'{prefix}/phi', chunk[field_idx][i, :, 1])
                    elif field == 'tod':
                        comm_tod.add_field(f'{prefix}/tod', chunk[field_idx][i, :])
                    elif field == 'psi':
                        comm_tod.add_field(f'{prefix}/psi', chunk[field_idx][i, :])
            comm_tod.finalize_chunk(curr_scan_id)
            comm_tod.finalize_file()
    comm_tod.make_filelists()


if __name__ == '__main__':
    main()
