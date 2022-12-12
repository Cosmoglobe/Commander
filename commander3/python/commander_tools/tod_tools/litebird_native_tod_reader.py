import h5py
import numpy as np

class LitebirdTodReader:
    def __init__(self, directory, file_prefixes, start_day=0, end_day=364, fields=['tod', 'pointings', 'psi']):
        # The first file prefix must be the file that contains psi and pointing info
        self.file_prefixes = file_prefixes
        self.directory = directory
        self.curr_day = start_day
        self.curr_day_tod_idx = 0
        self.end_day = end_day
        self.fields = fields

    def __enter__(self):
        self._open_files(self.curr_day)
        return self

    def __exit__(self, type, value, traceback):
        self._close_files()
#        self.currfile.close()

    def _open_files(self, curr_day):
        currfiles = []
        for prefix in self.file_prefixes:
            currfiles.append(h5py.File(f'{self.directory}/{prefix}_day{self.curr_day:04}.hdf5'))
        self.currfiles = currfiles

    def _close_files(self):
        for currfile in self.currfiles:
            currfile.close()

    def iterate_chunks(self, chunk_size):
        while self.curr_day <= self.end_day:
            remaining_chunk_size = chunk_size
            remaining_curr_day_tod_size = self.currfiles[0][self.fields[0]].shape[1] - self.curr_day_tod_idx
            chunks = []
            tod_idx = self.fields.index('tod')

            for field in self.fields:
                if field == 'pointings':
                    chunks.append(np.zeros((self.currfiles[0][field].shape[0], chunk_size, 2)))
                else:
                    chunks.append(np.zeros((self.currfiles[0][field].shape[0], chunk_size)))
            curr_chunk_idx = 0
            while remaining_curr_day_tod_size <= remaining_chunk_size:
                coadded_tod = None
                for i, currfile in enumerate(self.currfiles):
                    if i == 0:
                        for field_idx, field in enumerate(self.fields):
                            if field == 'tod': continue
                            chunks[field_idx][:, curr_chunk_idx:remaining_curr_day_tod_size + curr_chunk_idx] = currfile[field][:, self.curr_day_tod_idx:]
                        coadded_tod = currfile['tod'][:, self.curr_day_tod_idx:]
                    else:
                        coadded_tod += currfile['tod'][:, self.curr_day_tod_idx:]
                chunks[tod_idx][:, curr_chunk_idx:remaining_curr_day_tod_size + curr_chunk_idx] = coadded_tod
                curr_chunk_idx += remaining_curr_day_tod_size
                remaining_chunk_size -= remaining_curr_day_tod_size
                self.curr_day += 1
                self._close_files()
#                self.currfile.close()
                if self.curr_day > self.end_day:
                    break
                self._open_files(self.curr_day)
#                self.currfile = h5py.File(f'{self.directory}/{self.file_prefix}_day{self.curr_day:04}.hdf5')
                self.curr_day_tod_idx = 0
                remaining_curr_day_tod_size = self.currfiles[0][self.fields[0]].shape[1]
            if self.curr_day <= self.end_day:
                # Somewhere in the current day is the end of this chunk
                for field_idx, field in enumerate(self.fields):
                    chunks[field_idx][:, curr_chunk_idx:] = self.currfiles[0][field][:, self.curr_day_tod_idx:self.curr_day_tod_idx+remaining_chunk_size]
                self.curr_day_tod_idx += remaining_chunk_size
            yield chunks
