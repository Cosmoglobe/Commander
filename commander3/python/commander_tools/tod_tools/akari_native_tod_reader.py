import glob
import fitsio
import numpy as np
import healpy as hp
import pickle
import random
import math

"""
A class for reading AKARI TODs within a given range.

The current format is a bit awkward, as it requires the user to know which global TOD index range
they want. It also currently assumes that the user knows the internal workings of the FITS files. It
is possible to rework this, but for now I don't know whether that is worthwhile or not.

It assumes that the user has access to the original .fits.gz (or just .fits) AKARI TODs.

In other words, what this class does is to allow the user to specify an output format, and then
deliver that format for the range the user gives, automatically iterating over the .fits files
involved.
"""

class AKARITODReader:
    """
    An AKARITODReader instance will, upon request, fetch a continuous sequence of data from the
    AKARI FITS TOD files, joined together across files. The data to be fetched is defined by the
    user through a formatting function.

    Arguments:

    akari_fits_dir (str): The directory in which the .fits tod files reside (should not be
        compressed)
    band (str): The band for which to load the TODs. Should be '160', '140', '065' or '090'. (N160,
        WIDE-L, N60 or WIDE-S, respectively)
    fits2output_formatter (function): A function defining the output format of the TOD. The function
        signature should be (f, f_start, f_end, band), where f is an open fitsio.FITS instance
        pointing to one of the TOD files, f_start is the starting index of the current file, f_end
        is the ending index of the current file, and band is the same type of string as above. The
        function should then fetch the desired data corresponding to the range [f_start:f_end] (in
        Python slicing format), as well as any other data desired, in a dictionary. The dictionary
        should not contain sub-dictionaries. Each entry should either be a single scalar or a numpy
        array of length (f_end - f_start).
    load_idx_file_mapping (bool): The TODReader creates an internal mapping between files and the
        global TOD index represented by that file. This option allows the user to load this mapping
        from an existing .pkl file rather than to compute them. The directory from which to load the
        file is given by the 'mapping_dir' argument.
    save_idx_file_mapping(bool): Whether to save the internal mapping to file (not necessary to use
        together with 'load_idx_file_mapping'). The 'mapping_dir' argument is where the file will be
        saved.
    mapping_dir (str): The directory in which a potential mapping file will be saved.
    """
    def __init__(self, akari_fits_dir, band, fits2output_formatter, load_idx_file_mapping=False,
                 save_idx_file_mapping=True, mapping_dir=None):
        if band in ('160', '140'):
            file_identifier = 'LW'
        elif band in ('065', '090'):
            file_identifier = 'SW'
        else:
            raise ValueError(
                f"Unrecognized band {band}. Should be one of '160, '140', '090', or '065'")
        self.fitsdir = akari_fits_dir
        self.filelist = glob.glob(self.fitsdir + f'/*/*{file_identifier}*.fits')
        self.filelist.sort()
        self.filelist = self.filelist[:60]
        self.fits2output_formatter = fits2output_formatter
        if load_idx_file_mapping and mapping_dir is None:
            raise(ValueError(f"load_idx_file_mapping is True but no mapping_dir is given"))
        if save_idx_file_mapping and mapping_dir is None:
            raise(ValueError(f"save_idx_file_mapping is True but no mapping_dir is given"))
        if load_idx_file_mapping:
            with open(f'{mapping_dir}/idx_file_mapping_{band}.pkl', 'rb') as pkl_file:
                self.idx_file_mapping = pickle.load(pkl_file)
        else:
            self.idx_file_mapping = self._calculate_idx_file_mapping()
        if save_idx_file_mapping:
            with open(f'{mapping_dir}/idx_file_mapping_{band}.pkl', 'wb') as pkl_file:
                pickle.dump(self.idx_file_mapping, pkl_file)
        self.band = band

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

    def _calculate_idx_file_mapping(self):
        """Calculates the mapping between the fits files and their global tod index.

        Returns:
            A list of ints, in which each element corresponds to the same indexed-element in
            self.filelist. Each int represents the final index of that file, so that the current
            file global indices are file_mapping[i-1]:file_mapping[i]
        """
        self.idx_file_mapping = []
        curr_idx = 0
        num_files = len(self.filelist)
        year_passed = False
        for i, file in enumerate(self.filelist):
            f = fitsio.FITS(file)
            curr_idx += f[1].get_nrows()
            f.close()
            self.idx_file_mapping.append(curr_idx)
        return self.idx_file_mapping


    def get_file_index_range(self, filename):
        """Given a .fits file, returns the global tod index range for that file.

        Arguments:
            filename (str): The filename (including the full path) of the file for which to get
                indices.

        Returns:
            start_idx, end_idx (int): the unique global index range for the file in question, so
                that [start_idx:end_idx] is the (pythonic) range of the TODs within the file.

        """
        file_idx = self.filelist.index(filename)
        if file_idx == 0:
            start_idx = 0
        else:
            start_idx = self.idx_file_mapping[file_idx-1]
        end_idx = self.idx_file_mapping[file_idx]
        return [start_idx, end_idx]


    def get_data(self, start_idx, end_idx, should_compress=False):
        """Given global indices, gets the contiguous data corresponding to that index range.
        
        The output data are formatted according to what's given by fits2output_formatter.

        Arguments:
            start_idx, end_idx (int): The start and end indices of the data desired, so that it
                represents a slicing (potentially across fits files) equal to start_idx:end_idx.

        Returns:
            A dict with the same keys as output by fits2output_formatter, where each scalar value
            output byt fits2output_formatter is still a scalar, but where each numpy array output
            has a total length of end_idx - start_idx.
        """

        data = None
        tot_tod_size = end_idx - start_idx
        file_indices = self._get_file_indices(start_idx, end_idx)
        if file_indices[0] == 0:
            start_firstfile = 0
        else:
            start_firstfile = self.idx_file_mapping[file_indices[0]-1]
        first_file_tod_idx = start_idx - start_firstfile
        if file_indices[-1] == 0:
            last_file_tod_idx = min(self.idx_file_mapping[file_indices[-1]], end_idx)
        else:
            last_file_tod_idx = min(self.idx_file_mapping[file_indices[-1]], end_idx) - self.idx_file_mapping[file_indices[-1]-1]
        
        tot_idx = 0
        for file_index in file_indices:
            with fitsio.FITS(self.filelist[file_index]) as currfile:
                if file_index == file_indices[0]:
                    file_start_idx = first_file_tod_idx
                else:
                    file_start_idx = 0
                if file_index == file_indices[-1]:
                    file_end_idx = last_file_tod_idx
                else:
                    file_end_idx = currfile[1].get_nrows()
                n_points = file_end_idx - file_start_idx
                curr_data = self.fits2output_formatter(currfile, file_start_idx,
                                                       file_end_idx, self.band,
                                                       should_compress=should_compress)
            if data is None:
                data = {}
                for field, field_data in curr_data.items():
                    try:
                        ndim = len(field_data.shape)
                    except AttributeError:
                        pass
                    else:
                        output_shape = list(field_data.shape)
                        output_shape[0] = tot_tod_size
                        data[field] = np.zeros(output_shape, dtype=field_data.dtype)
            for field, field_data in curr_data.items():
                try:
                    field_data.shape
                except AttributeError:
                    data[field] = field_data
                else:
                    data[field][tot_idx:tot_idx+n_points] = field_data
            tot_idx += n_points

        return data

    def _get_file_indices(self, start_idx, end_idx):
        """Finds the internal file indices corresponding to global TOD indices.

        Arguments:
            start_idx, end_idx (int): The global tod indices to be queried.

        Output:
            A list of all the file indices that are 'touched' by the input global TOD indices.
        """
            
        start_fileind = np.searchsorted(self.idx_file_mapping, start_idx, side='right')
        end_fileind = np.searchsorted(self.idx_file_mapping, end_idx, side='left')
        file_indices = [i for i in range(start_fileind, end_fileind+1)]

        return file_indices


    @staticmethod
    def _ring_spin_axis(vecs):
        vecs = np.array(vecs)
        outAng = np.array([0., 0.])
        nsamps = min(100, len(vecs))
        pair1 = random.sample(range(len(vecs)), nsamps)
        pair2 = random.sample(range(len(vecs)), nsamps)
        for a1, a2 in zip(pair1, pair2):
            crossP = np.cross(vecs[a1], vecs[a2])
            if(crossP[0] < 0):
                crossP *= -1
            theta1, phi1 = hp.vec2ang(crossP)
            if(not math.isnan(theta1) and not math.isnan(phi1)):
                outAng[0] += theta1[0]/nsamps
                outAng[1] += phi1[0]/nsamps
            else:
                outAng[0] += outAng[0]/nsamps
                outAng[1] += outAng[1]/nsamps
        return outAng
