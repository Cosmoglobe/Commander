#================================================================================
#
# Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
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

import numpy as np
import h5py as h5

# Creates instrument files for LiteBIRD, output in three different files for the three Telescopes LFT, MFT, HFT
# Uses delta bandpasses

f=h5.File("LFT_instrument.h5", "w")
f2=h5.File("MFT_instrument.h5", "w")
f3=h5.File("HFT_instrument.h5", "w")

freq_name_LFT = np.array(["040", "050", "060", "068", "078", "089", "100_1", "119_1", "140_1"])
freq_name_MFT = np.array(["100_2", "119_2", "140_2", "166", "195_1"])
freq_name_HFT = np.array(["195_2", "235", "280", "337", "402"])
freq_LFT = np.array([40., 50., 60., 68., 78., 89., 100., 119., 140.])
freq_MFT = np.array([100., 119., 140., 166., 195.])
freq_HFT = np.array([195., 235., 280., 337., 402.])


for i in range(0,np.size(freq_name_LFT)):
    grp = f.create_group(freq_name_LFT[i])
    yval = grp.create_dataset("bandpass", (1,), dtype='f')
    xval = grp.create_dataset("bandpassx", (1,), dtype='f')
    yval[0]=1.
    xval[0]=freq_LFT[i]

for i in range(0,np.size(freq_name_MFT)):
    grp = f2.create_group(freq_name_MFT[i])
    yval = grp.create_dataset("bandpass", (1,), dtype='f')
    xval = grp.create_dataset("bandpassx", (1,), dtype='f')
    yval[0]=1.
    xval[0]=freq_MFT[i]

for i in range(0,np.size(freq_name_HFT)):
    grp = f3.create_group(freq_name_HFT[i])
    yval = grp.create_dataset("bandpass", (1,), dtype='f')
    xval = grp.create_dataset("bandpassx", (1,), dtype='f')
    yval[0]=1.
    xval[0]=freq_HFT[i]


subdetector040 = np.array(["0001a","0001b","0002a","0002b","0003a","0003b"])

for i in range(0, np.size(subdetector040)):
    grp = f.create_group(subdetector040[i])
    yval = grp.create_dataset("bandpass", (1,), dtype='f')
    xval = grp.create_dataset("bandpassx", (1,), dtype='f')
    yval[0]=1.
    xval[0]=40.
    b = grp.create_group("beam")
    bT = b.create_dataset("T", (1,), dtype = "f")
    bT[:]=0.
    bB = b.create_dataset("B", (1,), dtype = "f")
    bB[:]=0.
    bE = b.create_dataset("E", (1,), dtype = "f")
    bE[:]=0.
    blmax = grp.create_dataset("beamlmax", (1,), dtype='f')
    blmax[0]=2400
    bmmax= grp.create_dataset("beammmax", (1,), dtype='f')
    bmmax[0]=100
    el = grp.create_dataset("elip", (1,), dtype='f')
    el[0]=1.
    fw = grp.create_dataset("fwhm", (1,), dtype='f')
    fw[0]=69.3
    psi = grp.create_dataset("psi_ell", (1,), dtype='f')
    psi[0] = 1.
    sl= grp.create_group("sl")
    slT = sl.create_dataset("T", (1,), dtype = "f")
    slT[:]=0.
    slB = sl.create_dataset("B", (1,), dtype = "f")
    slB[:]=0.
    slE = sl.create_dataset("E", (1,), dtype = "f")
    slE[:]=0.
    sll= grp.create_dataset("sllmax", (1,), dtype='f')
    sll[0]=512
    slm= grp.create_dataset("slmmax", (1,), dtype='f')
    slm[0]=100
    eff= grp.create_dataset("mbeam_eff", (1,), dtype='f')
    eff[0]=1.

# Also need 'bandpass', 'bandpassx', 'beam', 'beamlmax', 'beammmax', 'elip', 'fwhm', 'psi_ell', 'sl', 'sllmax', 'slmmax' for each single detector
