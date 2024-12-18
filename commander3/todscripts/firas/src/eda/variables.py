import h5py
import numpy as np

import matplotlib.pyplot as plt

fdq_sdf = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5")
fdq_eng = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5")

channels = {"lh": 2, "ll": 3, "rh": 0, "rl": 1}

print("adds_per_group")

for channel, chan in channels.items():
    # allowed values for adds per group: 1, 2, 3, 8, 12
    adds_per_group = np.array(fdq_sdf[f"fdq_sdf_{channel}/sci_head/sc_head9"]).tolist()
    unique_values = sorted(set(adds_per_group))

    print(f"Channel: {channel} - sc_head9")
    for val in unique_values:
        print(val, adds_per_group.count(val))

    up_adds_per_group = np.array(fdq_eng[f"chan/up_adds_per_group"][:, chan]).tolist()
    unique_values = sorted(set(up_adds_per_group))

    print(f"Channel: {channel} - up_adds_per_group")
    for val in unique_values:
        print(val, up_adds_per_group.count(val))

print("sweeps")

for channel, chan in channels.items():
    # allowed values for sweeps: 1, 4, 16?
    sweeps = np.array(fdq_sdf[f"fdq_sdf_{channel}/sci_head/sc_head11"]).tolist()

    unique_values = sorted(set(sweeps))

    print(f"Channel: {channel} - sc_head11")
    for val in unique_values:
        print(val, sweeps.count(val))

print("gain")

for channel, chan in channels.items():
    # allowed values for gain: 1, 3, 10, 30, 100, 300, 1000, 3000
    gain = np.array(fdq_sdf[f"fdq_sdf_{channel}/sci_head/gain"]).tolist()

    unique_values = sorted(set(gain))

    print(f"Channel: {channel} - gain")
    for val in unique_values:
        print(val, gain.count(val))

    sci_gain = np.array(fdq_eng[f"chan/sci_gain"][:, chan]).tolist()
    unique_values = sorted(set(sci_gain))

    print(f"Channel: {channel} - sci_gain")
    for val in unique_values:
        print(val, sci_gain.count(val))


print("investigate sweeps")
data = h5py.File("./../../data/df_v14.h5", "r")

sweeps_ll = np.array(data["df_data/sweeps_ll"][()])
sweeps_rh = np.array(data["df_data/sweeps_rh"][()])

plt.plot(sweeps_ll, sweeps_rh)
plt.show()
