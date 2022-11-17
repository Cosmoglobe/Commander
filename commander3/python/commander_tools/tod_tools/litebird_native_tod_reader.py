import h5py
import numpy as np

class LitebirdTodReader:
    def __init__(self, directory, file_prefix, start_day=0, end_day=364, field='tod'):
        self.file_prefix = file_prefix
        self.directory = directory
        self.curr_day = start_day
        self.curr_day_tod_idx = 0
        self.end_day = end_day
        self.field = field

    def __enter__(self):
        self.currfile = h5py.File(f'{self.directory}/{self.file_prefix}_day{self.curr_day:04}.hdf5')
        return self

    def __exit__(self, type, value, traceback):
        self.currfile.close()

    def iterate_chunks(self, chunk_size):
        while self.curr_day <= self.end_day:
            remaining_chunk_size = chunk_size
            remaining_curr_day_tod_size = self.currfile[self.field].shape[1] - self.curr_day_tod_idx
            if self.field == 'pointings':
                chunk = np.zeros((self.currfile[self.field].shape[0], chunk_size, 2))
            else:
                chunk = np.zeros((self.currfile[self.field].shape[0], chunk_size))
            curr_chunk_idx = 0
            while remaining_curr_day_tod_size <= remaining_chunk_size:
                chunk[:, curr_chunk_idx:remaining_curr_day_tod_size + curr_chunk_idx] = self.currfile[self.field][:, self.curr_day_tod_idx:]
                curr_chunk_idx += remaining_curr_day_tod_size
                remaining_chunk_size -= remaining_curr_day_tod_size
                self.curr_day += 1
                self.currfile.close()
                if self.curr_day > self.end_day:
                    break
                self.currfile = h5py.File(f'{self.directory}/{self.file_prefix}_day{self.curr_day:04}.hdf5')
                self.curr_day_tod_idx = 0
                remaining_curr_day_tod_size = self.currfile[self.field].shape[1]
            if self.curr_day <= self.end_day:
                # Somewhere in the current day is the end of this chunk
                chunk[:, curr_chunk_idx:] = self.currfile[self.field][:, self.curr_day_tod_idx:self.curr_day_tod_idx+remaining_chunk_size]
                self.curr_day_tod_idx += remaining_chunk_size
            yield chunk
