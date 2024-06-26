# Inputs include the engineering data and scientific data:
# FDQ_ENG=/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng.h5
# FDQ_SDF=/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf.h5
# Current files are direct numpy array dumps

import h5py


data = h5py.File('/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf.h5')
print(data.keys())

