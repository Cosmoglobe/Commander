import functools
from cosmoglobe.tod_tools import CommanderDataAdapter
import glob
from commander_tools.tod_tools.akari_native_tod_reader import AKARITODReader
#from akari_commander_data_adapter import AKARITODReader
from astropy.coordinates import SkyCoord
import astropy.coordinates as coords
import astropy.units as u
import numpy as np
import healpy
from astropy.time import Time, TimeDelta
import os.path
import pickle

""" 
This module's purpose is to define the AKARICommanderDataAdapter, which is the link between the
AKARI TOD reader code and the CommanderHDFWriter class. The latter class takes an instance of this
class as its input.

This is the file to change if there is to be any changes in the HDF file layout, such as flags,
detector names etc.
"""

BANDS = ['065', '090', '140', '160']

NUM_DETECTORS = {
    '065': 40, #N60
    '090': 60, #WIDE-S
    '140': 45, #WIDE-L 
    '160': 30  #N160
}

START_DET_INDS = { # These are the starting indices for each detector inside the FITS files.
    '065': 60, #N60
    '090': 0,  #WIDE-S
    '140': 0,  #WIDE-L 
    '160': 45  #N160
}

# These are the implied starting indices of each detector for the lon/lat pkl files created by
# Ryosuke. I currently assume that these are 'wrongly' indexed compared to the 'standard' starting
# indices (above) because this wrong assumption was made throughout, and I don't believe these
# lonlat files depend on the fits files any further than just getting the boresight pointing, which
# is common for all detectors. Thus, I believe that Ryosuke used the implied indexing scheme below
# when he calculated the lon/lats of each detector pixel. This we should take away when a proper
# calculation is implemented rather than using the pkl files.
# I also don't assume that the below will be the correct indexing for the pkl *flag* array, because
# those flags actually come from the .fits files and have just wrongly indexed when added to the HDF
# files, I believe.
START_DET_INDS_LONLAT_PKL_FILES = { 
    '065': 0,
    '090': 40,
    '140': 0,
    '160': 45
}

BAND_DETS = {}


# What follows are mappings between flag names and their position in the FITS files. Most of these
# can be found in the Data User Manual, exceptions are noted.


# These are the per-detector flags given in the first extension of the FITS files.
PIX_FLAG_MAP = {
    'bad': 0,
    'undef_anom': 1,
    'arith_err': 2,
    'dead': 3,
    'saturate': 4,
    'reset': 5,
    'rstanom': 6,
    'no_diff': 7,
    'no_rp_corr': 8,
    'no_dk_corr': 9,
    'no_dccal_corr': 10,
    'no_tr_corr': 11,
    'no_flat_field': 12,
    'no_gpgl': 13,
    'no_mtgl': 14,
    'no_abscal': 15,
    'gpgl_type1': 16,  # Glitch detection using Gaussian Processing method
    'gpgl_type2': 17,
    'gpgl_type3': 18,
    'gpgl_type4': 19,
    'gpgl_tail': 20,
    'mtgl_type1': 21, # Glitch detection using median transformation method
    'mtgl_type2': 22,
    'mtgl_type3': 23,
    'mtgl_type4': 24,
    'mtgl_tail': 25,
    'no_peri_corr': 26,
    'sx_peak': 27,       # For SUSSEXtractor
    'sx_source': 28,
    'sx_obj': 29,
    'flutter': 30, # This is not in the FIS data manual, but in the tsd_def.pdf
    'blank': 31
}

# These are the 'frame' flags that are true for all detectors simultaneously, given in the first
# FITS extension.
FRAME_FLAG_MAP = {
    'bad_frame': 0,
    'undef_anom_frame': 1,
    'blank': 2,
    'in_saa': 3,
    'near_moon': 4,
    'untrusted_frame': 5,
}

# These flags should hopefully not be loaded from the pkl files soon. They contain attitude
# information + moon avoidance angle information
GADS_FLAG_MAP = {
    'att_cnt_md_b1': 0, # 0 for ASS and maneuver modes, 1 for pointed observations
    'att_cnt_md_b0': 1, # 0 for ASS mode, 1 for maneuver mode
    'aa_lun2' : 2       # 1 when lunar avoidance angle < 21.5
}

# These are the 'status' flags that are true for all detectors simultaneously, given in the first
# FITS extension.
STATUS_FLAG_MAP = {
    'creon': 0,
    'shtop': 1,
    'fwposon': 2,
    'fwpos_b1': 3,
    'fwpos_b0': 4,
    'mposon': 5,
    'mpos_b1': 6,
    'mpos_b0': 7,
    'rstwidelon': 8,
    'rstwideson': 9,
    'rstn170on': 10,
    'rstn60on': 11,
    'lwbooston': 12,
    'swbooston': 13,
    'lwbiason': 14,
    'swbiason': 15,
    'calalon': 16,
    'calason': 17,
    'calbon': 18,
    'sinalon': 19,
    'sinason': 20
}

# This dictionary indicates 1) which flags to include in the final HDF, and 2) which bit the flag
# should be represented by (thus, a given bit should just be covered by one of the desired flags).
# The ground-level keys are 'frame_flags', 'pix_flags' and 'status_flags', each of which points to
# another dict that contains the relevant flags. Each flag points to the bit number corresponding to
# the flag.
DESIRED_FLAGS = {
    'frame_flags': {
        'bad_frame': 0,
        'near_moon': 1
    },
    'pix_flags': {
        'bad': 2,
        'dead': 3,
        'saturate': 4,
        'reset': 5,
        'rstanom': 6,
        'mtgl_tail': 7,
        'flutter': 8,
        'no_rp_corr': 9,
        'blank': 10,
        'mtgl_type1': 11,
        'mtgl_type2': 12,
        'mtgl_type3': 13,
        'mtgl_type4': 14,
    },
    'status_flags': {
        'calalon': 15,
        'calason': 16,
        'shtop': 17, # Shutter open flag. NB!!!!! This is inverted! When 1 in the FITS files, should
                     # be 0 in the Commander HDF files (i.e. in the Commander HDF the interpretation
                     # is 'the shutter closed flag'
    },
    'gads_flags': {  #Flags that hopefully will be generated by this code eventually rather than
                     #being loaded from the old .pkl files.
        'att_cnt_md_b1': 18, # 0 for ASS and maneuver modes, 1 for pointed observations
        'att_cnt_md_b0': 19, # 0 for ASS mode, 1 for maneuver mode
        'aa_lun2' : 20       # 1 when lunar avoidance angle < 21.5

    },
    'mode': {  #Special type of flag defined by the packet_id data.
        'cds': 21,
        'coadd': 22,
        'maneuver': 23
    }
}

# How many FITS files should be included in a single HDF file
NUM_FITS_FILES_PER_SEGMENT = 200
#NUM_FITS_FILES_PER_SEGMENT = 4

REFERENCE_TIME = Time("2000-01-01", format="isot", scale="utc")

PATH_TO_GADS_PKL_FILES = '/mn/stornext/d23/cmbco/akari/akari_TSD_pkl/gads/all/'
PATH_TO_LONLAT_PKL_FILES = '/mn/stornext/d23/cmbco/akari/akari_TSD_pkl/'

for band, ndets in NUM_DETECTORS.items():
    BAND_DETS[band] = []
    for i in range(1, ndets+1):
        detname = f'AKARI_{band}-{i:02d}'
        BAND_DETS[band].append(detname)


def load_gads_flags(fitsfilename, band, pkl_path=PATH_TO_GADS_PKL_FILES):
    """This function loads certain flags that have been extracted by Tamaki and
    Ryosuke from the GADS data (see data manual)

    Arguments:
        fitsfilename (str): The filename of the (original) TOD FITS file for which we
            want to extract the flags. Not actually opened, just used to figure
            out which section of data we're interested in.
        band (str): The band name for which we want to extract the flags.
        pkl_path (str): The path to where the GADS PKL files reside.

    Returns:
        A dict where each key (the detector names) points to a bool numpy array
            with the GADS flags set.
    """
    ff = os.path.basename(fitsfilename)
    split_ff = ff.split('_')
    time = split_ff[2]
    outflags = {}
    start_ind = START_DET_INDS[band]
    numdets = NUM_DETECTORS[band]
    detnames = BAND_DETS[band]
    for flagidx, flag in enumerate(['att_cnt_md_b1', 'att_cnt_md_b0', 'aa_lun2']):
        with open(pkl_path + f'/{time[:10]}/{"_".join(split_ff[:3])}_{flag}.pkl', 'rb') as f:
            currflags = np.array(pickle.load(f)).astype(bool)[start_ind:start_ind+numdets]
        for detidx, detname in enumerate(detnames):
            if detname not in outflags:
                outflags[detname] = np.zeros((len(currflags[0]), 3), dtype=bool)
            outflags[detname][:, flagidx] = currflags[detidx, :]
    return outflags


def load_lonlat_pkl(fitsfilename, band, pkl_path=PATH_TO_LONLAT_PKL_FILES):
    """This function loads the lon/lat files that have been created by Ryosuke
    and which provides accurate pointing for each of the detector pixels.

    Arguments:
        fitsfilename (str): The filename of the (original) TOD FITS file for which we
            want to extract the lon/lat.
        band (str): The band name for which we want to extract the lon/lat.
        pkl_path (str): The path to where the lonlat PKL files reside.

    Returns:
        Two dicts (lon, lat) where each dict key (the detector names) point to
            a numpy array containing the lon or lat, respectively, for that
            detector.
    """
    ff = os.path.basename(fitsfilename)
    split_ff = ff.split('_')
    time = split_ff[2]
    start_ind = START_DET_INDS_LONLAT_PKL_FILES[band]
    numdets = NUM_DETECTORS[band]
    detnames = BAND_DETS[band]
    outlons = {}
    outlats = {}
    lonfname = pkl_path + f'/lon/all/{time[:10]}/{"_".join(split_ff[:3])}_gb_lon_gads.pkl' 
    latfname = pkl_path + f'/lat/all/{time[:10]}/{"_".join(split_ff[:3])}_gb_lat_gads.pkl' 
    with open(lonfname, 'rb') as f_lon, open(latfname, 'rb') as f_lat:
        currlon = np.array(pickle.load(f_lon))[start_ind:start_ind+numdets]
        currlat = np.array(pickle.load(f_lat))[start_ind:start_ind+numdets]
        for detidx, detname in enumerate(detnames):
            outlons[detname] = currlon[detidx, :]
            outlats[detname] = currlat[detidx, :]
    return outlons, outlats

# This is a callback function to be used with the AKARITODReader. It formats the output data from
# the TOD reader in a format that is useful for the HDF file generation.
def fits2output_formatter(file, start_index, end_index, band,
                          start_det_inds=START_DET_INDS,
                          num_detectors=NUM_DETECTORS, band_dets=BAND_DETS,
                          reftime=REFERENCE_TIME,
                          should_compress=False):
    """Formats the data in an AKARI fits file in the way we need for the Commander HDFs.

    Arguments:
        file (fitsio.FITS): an open fitsio.FITS file instance, pointing to a TOD fits file.
        start_index, end_index (int): The indices (relative to the start of the file) to fetch.
        band (str): The band in question (is needed because each fits file contains data for two
                    bands)

    Returns:
        dictionary containing the data we want from each fits file.
    """
    out_data = {}
    # There is an implicit assumption here that we only want contiguous detectors. This seems to be
    # an okay assumption for now, but I just note it here.
    start_det_ind = start_det_inds[band]
    end_det_ind = start_det_ind + num_detectors[band]
    tot_adu = file[1]['DET'].read()[start_index:end_index, start_det_ind:end_det_ind]
    tot_flux = file[1]['FLUX'].read()[start_index:end_index, start_det_ind:end_det_ind]
    tot_pixflag = file[1]['PIX_FLAG'].read()[start_index:end_index,
                                             start_det_ind:end_det_ind, :].astype(bool)
    for local_det_idx, detname in enumerate(band_dets[band]):
        out_data[f'{detname}/tod'] = tot_flux[:, local_det_idx]
        out_data[f'{detname}/todz'] = tot_adu[:, local_det_idx]
        out_data[f'{detname}/pixel_flag'] = tot_pixflag[:, local_det_idx, :]
    out_data['status_flag'] = file[1]['STATUS'].read()[start_index:end_index].astype(bool)
    # NB! This is necessary because we interpret the shtop flag inversely internally in commander
    shtop_idx = STATUS_FLAG_MAP['shtop']
    out_data['status_flag'][:, shtop_idx] = ~out_data['status_flag'][:, shtop_idx]
    out_data['frame_flag'] = file[1]['FLAG'].read()[start_index:end_index].astype(bool)
    out_data['aftime'] = file[1]['AFTIME'].read()[start_index:end_index]
    out_data['ra'] = file[5]['RA'].read()[start_index:end_index]
    out_data['dec'] = file[5]['DEC'].read()[start_index:end_index]
    out_data['packet_id'] = file[1]['PACKETID'].read()[start_index:end_index].astype(int)
#    out_data['start_satpos'] = 0
#    out_data['end_satpos'] = 0
#    out_data['start_earthpos'] = 0
#    out_data['end_earthpos'] = 0

    return out_data


class AKARICommanderDataAdapter(CommanderDataAdapter):
    """An AKARI-specific implementation of the CommanderDataAdapter.

    All the non-underlined functions are defined as in that class.
    """
    def __init__(self, akari_fits_dir, nside, bands=BANDS,
                 num_detectors=NUM_DETECTORS, band_dets=BAND_DETS,
                 reference_time=REFERENCE_TIME,
                 extend_reset_flag=False):
        self.bands = bands
        self.num_detectors = num_detectors
        self.band_dets = band_dets
        self.nside = nside
        self.extend_reset_flag = extend_reset_flag


        self.fsamp = {
            '065': 25.28,  #N60
            '090': 25.28,  #WIDE-S
            '140': 16.86,  #WIDE-L 
            '160': 16.86   #N160
        }
        self.band2filetype = {
            '065': 'SW',
            '090': 'SW',
            '140': 'LW',
            '160': 'LW'
        }

        self.detectors = {}
        self.todreaders = {}
        self.akari_fits_dir = akari_fits_dir
        self.filelist = {}
        for band in self.bands:
            self.todreaders[band] = AKARITODReader(
                akari_fits_dir, band, fits2output_formatter, load_idx_file_mapping=True,
                save_idx_file_mapping=False,
                mapping_dir='/mn/stornext/u3/eirikgje/data/akari_analysis/')
            self.filelist[band] = self.todreaders[band].filelist
            self.detectors[band] = []
            for i in range(1, self.num_detectors[band]+1):
                self.detectors[band].append(f'AKARI_{band}-{i:02d}')

        self.chunk_file_map = self._calculate_chunk_files()

        self.nchunks = {}
        for band in bands: 
            self.nchunks[band] = len(self.chunk_file_map[band])
        self.reference_time = reference_time


    def _calculate_chunk_files(self):
        """
        Internal function that calculates a mapping between chunk index and fits files. This is
        temporary, as we might want to do something more sophisticated about how we define the
        chunks.
        """
        chunk_file_map = {}
        for band in self.bands:
            chunk_file_map[band] = {}
            curr_chunk_idx = 1
            for file in self.filelist[band]:
                chunk_file_map[band][curr_chunk_idx] =  []
                chunk_file_map[band][curr_chunk_idx].append(file)
                curr_chunk_idx += 1
        return chunk_file_map

    def get_num_segments(self, band: str):
        return int(np.ceil(self.nchunks[band] / NUM_FITS_FILES_PER_SEGMENT))

    def get_version(self):
        return 3

    def get_experiment_name(self):
        return "AKARI"

    # Number of psi bins. This is the default we seem to be using for experiments.
    def get_npsi(self):
        return 64

    def get_nside(self, band: str):
        return self.nside

    def get_polangs(self, band: str):
        ndet = self.num_detectors[band]
        return [0]*ndet

    def get_mbangs(self, band: str):
        ndet = self.num_detectors[band]
        return [0]*ndet

    def get_fsamp(self, band: str):
        fsamp = self.fsamp[band]
        return fsamp

    def get_detector_names(self, band: str):
        return self.detectors[band]

    def get_chunk_indices(self, band: str, segment:int):
        start_idx = (segment-1) * NUM_FITS_FILES_PER_SEGMENT + 1
        end_idx = segment * NUM_FITS_FILES_PER_SEGMENT + 1
        if end_idx > self.nchunks[band] + 1:
            end_idx = self.nchunks[band] + 1
        return [i for i in range(start_idx, end_idx)]

    def _process_akari_tod_reader_data(self):
        """Internal function to process the data coming from the AKARI TOD reader and getting it
        into the format needed for the HDFwriter"""
        todreader_data = {}
        for band in self.bands:
            todreader_data[band] = {}
            files = self.chunk_file_map[band][self.chunk_idx]
#            print(f'files: {files}')
            # This does not yet include any code for when the start and end happens
            # within a file rather than at the beginning of the file.
            start_idx = self.todreaders[band].get_file_index_range(files[0])[0]
            end_idx = self.todreaders[band].get_file_index_range(files[-1])[1]
            curr_data = self.todreaders[band].get_data(start_idx, end_idx)
            mode = curr_data['packet_id']
            if len(files) > 1:
                raise NotImplementedError("""We currently don't support loading gads flags for several
                                          different fits files.""")
            gads_flags = load_gads_flags(files[0], band)
            detlons, detlats = load_lonlat_pkl(files[0], band)
#            print(gads_flags)

            spin_axis_arr = np.zeros((2, len(self.band_dets[band])))
            for detidx, det in enumerate(self.band_dets[band]):
                todreader_data[band][det] = {}
                todreader_data[band][det]['tod'] = curr_data[f'{det}/tod']
                todreader_data[band][det]['todz'] = curr_data[f'{det}/todz']
                todreader_data[band][det]['flag'] = self._process_akari_flags(
                    curr_data[f'{det}/pixel_flag'],
                    curr_data['frame_flag'],
                    curr_data['status_flag'], gads_flags[det], mode,
                    extend_reset_flag=self.extend_reset_flag)
                c = SkyCoord(ra=detlons[det]*u.deg, dec=detlats[det]*u.deg, frame='icrs', distance=1*u.AU)
                todreader_data[band][det]['pix'] = healpy.ang2pix(self.nside,
                                                                  c.galactic.l.value,
                                                                  c.galactic.b.value,
                                                                  lonlat=True)

                theta = 0.5 * np.pi - c.hcrs.dec.radian
                phi = c.hcrs.ra.radian
                todreader_data[band][det]['pix_solarcentric'] = healpy.ang2pix(
                    self.nside, theta, phi)
                vecs = c.galactic.cartesian.get_xyz().value.T
                spin_axis_arr[:,detidx] = self.todreaders[band]._ring_spin_axis(vecs)

            if np.any(~np.isfinite(spin_axis_arr)):
                mu = np.nanmean(spin_axis_arr, axis=1)
                print('Got another nan', mu, start_idx,
                        self.chunk_file_map[band][self.chunk_idx], c)
            else:
                mu = spin_axis_arr.mean(axis=1)
            todreader_data[band]['spin_axis'] = mu

            ra = curr_data['ra']
            dec = curr_data['dec']
            c = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs', distance=1*u.AU)
            todreader_data[band]['pix'] = healpy.ang2pix(self.nside,
                                                         c.galactic.l.value,
                                                         c.galactic.b.value,
                                                         lonlat=True)

            theta = 0.5 * np.pi - c.hcrs.dec.radian
            phi = c.hcrs.ra.radian
            todreader_data[band]['pix_solarcentric'] = healpy.ang2pix(
                self.nside, theta, phi)



            if 'start_time' not in todreader_data:
                starttime = curr_data['aftime'][0]
                endtime = curr_data['aftime'][-1]
                todreader_data['start_time'] = (self.reference_time +
                                                TimeDelta(curr_data[f'aftime'][0], format='sec',
                                                          scale='tai'))
                todreader_data['end_time'] = (self.reference_time +
                                                TimeDelta(curr_data[f'aftime'][-1], format='sec',
                                                          scale='tai'))
                earthpos_start = coords.get_body(
                    'earth', todreader_data['start_time'],
                    ephemeris='builtin').transform_to(coords.HeliocentricMeanEcliptic).cartesian.xyz.to_value(u.AU)
                earthpos_end = coords.get_body(
                    'earth', todreader_data['end_time'],
                    ephemeris='builtin').transform_to(coords.HeliocentricMeanEcliptic).cartesian.xyz.to_value(u.AU)
                todreader_data['start_satpos'] = earthpos_start
                todreader_data['end_satpos'] = earthpos_end
                todreader_data['start_earthpos'] = earthpos_start
                todreader_data['end_earthpos'] = earthpos_end

        return todreader_data

    
    def _process_akari_flags(self, pixflag_array, frameflag_array,
                             statusflag_array, gadsflag_array,
                             mode, pixflagmap=PIX_FLAG_MAP, frameflagmap=FRAME_FLAG_MAP,
                             statusflagmap=STATUS_FLAG_MAP, gadsflagmap=GADS_FLAG_MAP,
                             desired_flags=DESIRED_FLAGS, extend_reset_flag=False):
        """
        Internal function that takes the flags coming from the fits files and turns them into the
        bitmasks needed for the HDF file"""
        outflags = np.zeros(len(pixflag_array), dtype=int)
        for flagtype, flag_arr, flagmap in zip(('pix_flags',
                                                'frame_flags', 
                                                'status_flags',
                                                'gads_flags'),
                                               (pixflag_array,
                                                frameflag_array,
                                                statusflag_array,
                                                gadsflag_array),
                                               (pixflagmap,
                                                frameflagmap,
                                                statusflagmap,
                                                gadsflagmap)):
            for flag_name, target_bit in desired_flags[flagtype].items():
                flag_idx = flagmap[flag_name]
                curr_mask = flag_arr[:, flag_idx]
                if extend_reset_flag and flagtype == 'pix_flags' and flag_name == 'reset':
                    # If this option is set, we want to extend the reset flag
                    # to cover four samples rather than one.
                    for i in range(1, 4):
                        curr_mask[i:] = curr_mask[i:] | curr_mask[:-i]
                outflags[curr_mask] += 2 ** int(target_bit)
        modeflags = desired_flags['mode']
        # These numbers are the modes of the pixflag_array field, listed in the AKARI fits data user
        # manual.
        is_cds = mode == 65
        is_coadd = mode == 66
        is_maneuver = mode == 82
        for mode, modearr in zip(('cds', 'coadd', 'maneuver'),
                                 (is_cds, is_coadd, is_maneuver)):
            target_bit = modeflags[mode]
            outflags[modearr] += 2 ** int(target_bit)

        return outflags

    def get_bands(self):
        return self.bands

    def set_chunk_index(self, chunk_idx: int):
        self.chunk_idx = chunk_idx
        self.all_chunk_data = None #Reset chunk

    def get_chunk(func):
        """Decorator function to make sure that a chunk is always preloaded"""

        @functools.wraps(func)
        def inner(self, *args):
            if self.all_chunk_data is None:
                self.all_chunk_data = self._process_akari_tod_reader_data()
            return func(self, *args)

        return inner

    @get_chunk
    def get_chunk_data(self, band: str, det:str):
        # The third entry is just zeroes - i.e. currently we're not operating with a 'psi'.
        tod = self.all_chunk_data[band][det]['tod']
        todz = self.all_chunk_data[band][det]['todz']
        psi_arr = np.zeros_like(self.all_chunk_data[band][det]['pix'], dtype=int)
        return (tod,
                todz,
                np.array([self.all_chunk_data[band][det]['pix'],
                          self.all_chunk_data[band][det]['pix_solarcentric']]),
                psi_arr,
                self.all_chunk_data[band][det]['flag'])

    @get_chunk
    def get_chunk_start_time(self):
        return self.all_chunk_data['start_time']

    @get_chunk
    def get_chunk_end_time(self):
        return self.all_chunk_data['end_time']

    @get_chunk
    def get_chunk_start_satpos(self):
        return self.all_chunk_data['start_satpos']

    @get_chunk
    def get_chunk_end_satpos(self):
        return self.all_chunk_data['end_satpos']

    @get_chunk
    def get_chunk_start_earthpos(self):
        return self.all_chunk_data['start_earthpos']

    @get_chunk
    def get_chunk_end_earthpos(self):
        return self.all_chunk_data['end_earthpos']

    @get_chunk
    def get_chunk_satvel(self):
        return [0]*3


    def get_gain(self, band: str, detector: str):
        return 1

    def get_alpha(self, band:str, detector:str):
        return -1

    def get_fknee(self, band:str, detector:str):
        return  1 / (10 * 60)

    def get_sigma0(self, band:str, detector:str):
        return 1

    def get_spinaxis(self, band:str):
        '''
        Position in radians of the spin axis of the satellite.
        Referred to as outP at various points in the code.
        '''
        return self.all_chunk_data[band]['spin_axis']
