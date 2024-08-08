#================================================================================
#
# Copyright (C) 2024 Institute of Theoretical Astrophysics, University of Oslo.
#
# This file is part of Commander3.
#
# Commander3 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Commander3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Commander3. If not, see <https://www.gnu.org/licenses/>.
#
#================================================================================

import h5py
import os
import healpy as hp
import numpy as np
import math
import argparse
from astropy.io import fits
from commander_tools.tod_tools.spider import spider
from commander_tools.tod_tools import commander_instrument as inst

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--out-dir', type=str, action='store', help='output directory', default='/mn/stornext/d16/cmbco/bp/mathew/spider_data')

    parser.add_argument('--det-lists', type=str, action='store', help='directory containing the detector list files', default='/mn/stornext/d16/cmbco/bp/mathew/spider_data')

    parser.add_argument('--bandpass-dir', type=str, action='store', help='directory with the bandpass files', default='/mn/stornext/d16/cmbco/spider/bandpass')

    args = parser.parse_args()
    outDir = args.out_dir

    version = 4

    inst_file = inst.commander_instrument(outDir, spider.instrument_filename(version), version, 'w')

    for freq in spider.freqs:

        detlist = np.genfromtxt(os.path.join(args.det_lists, spider.detlist_name(freq)), dtype=str)

        bandpass = np.loadtxt(os.path.join(args.bandpass_dir, 'bandpass_spider_'+str(freq).zfill(3)+'GHz_transpose_threshold.txt'))

        inst_file.add_bandpass('SPIDER_' + str(freq), bandpass[:,0], bandpass[:,1])

        for det in detlist:
            inst_file.add_bandpass(det, bandpass[:,0], bandpass[:,1])

            #fake alms for beam and sidelobe
            inst_file.add_alms(det, 'beam', 0, 0, [0], [0], [0])
            inst_file.add_alms(det, 'sl',   0, 0, [0], [0], [0])

            inst_file.add_field(det + '/centFreq', data=freq)

            #beam parameters spider doesn't use
            inst_file.add_field(det +'/elip', data=0)
            inst_file.add_field(det + '/mbeam_eff', data=1)
            inst_file.add_field(det + '/psi_ell', data=0)


    inst_file.finalize()
    spider.verify_instrument_file(outDir, version)

if __name__ == '__main__':
    main()
