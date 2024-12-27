import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

fdq_sdf = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5")
fdq_eng = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5")
data = h5py.File("./../../data/df_v14.h5", "r")

channels = {"lh": 2, "ll": 3, "rh": 0, "rl": 1}

print("adds_per_group")

adds_per_group = {}
for channel, chan in channels.items():
    # allowed values for adds per group: 1, 2, 3, 8, 12
    adds_per_group[channel] = np.array(
        fdq_sdf[f"fdq_sdf_{channel}/sci_head/sc_head9"]
    ).tolist()
    unique_values = sorted(set(adds_per_group[channel]))

    print(f"Channel: {channel} - sc_head9")
    for val in unique_values:
        print(val, adds_per_group[channel].count(val))

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

sweeps_ll = np.array(data["df_data/sweeps_ll"][()])
sweeps_rh = np.array(data["df_data/sweeps_rh"][()])
sweeps_lh = np.array(data["df_data/sweeps_lh"][()])
sweeps_rl = np.array(data["df_data/sweeps_rl"][()])

fig, ax = plt.subplots(4, 4, figsize=(10, 10), sharex=True, sharey=True)
ax[0, 0].plot(sweeps_ll, sweeps_lh, "o")
ax[0, 1].plot(sweeps_ll, sweeps_rl, "o")
ax[0, 2].plot(sweeps_ll, sweeps_rh, "o")
ax[0, 3].plot(sweeps_lh, sweeps_rl, "o")
ax[1, 0].plot(sweeps_lh, sweeps_rh, "o")
ax[1, 1].plot(sweeps_lh, sweeps_rl, "o")
ax[1, 2].plot(sweeps_rl, sweeps_rh, "o")
ax[1, 3].plot(sweeps_rl, sweeps_ll, "o")
ax[2, 0].plot(sweeps_rh, sweeps_ll, "o")
ax[2, 1].plot(sweeps_rh, sweeps_lh, "o")
ax[2, 2].plot(sweeps_rh, sweeps_rl, "o")
ax[2, 3].plot(sweeps_rh, sweeps_rh, "o")
ax[3, 0].plot(sweeps_rl, sweeps_ll, "o")
ax[3, 1].plot(sweeps_rl, sweeps_lh, "o")
ax[3, 2].plot(sweeps_rl, sweeps_rl, "o")
ax[3, 3].plot(sweeps_rl, sweeps_rh, "o")
fig.suptitle("Sweeps")
plt.show()

gain_ll = np.array(data["df_data/gain_ll"][()])
gain_rh = np.array(data["df_data/gain_rh"][()])
gain_lh = np.array(data["df_data/gain_lh"][()])
gain_rl = np.array(data["df_data/gain_rl"][()])

fig, ax = plt.subplots(4, 4, figsize=(10, 10), sharex=True, sharey=True)
ax[0, 0].plot(gain_ll, gain_lh, "o")
ax[0, 1].plot(gain_ll, gain_rl, "o")
ax[0, 2].plot(gain_ll, gain_rh, "o")
ax[0, 3].plot(gain_lh, gain_rl, "o")
ax[1, 0].plot(gain_lh, gain_rh, "o")
ax[1, 1].plot(gain_lh, gain_rl, "o")
ax[1, 2].plot(gain_rl, gain_rh, "o")
ax[1, 3].plot(gain_rl, gain_ll, "o")
ax[2, 0].plot(gain_rh, gain_ll, "o")
ax[2, 1].plot(gain_rh, gain_lh, "o")
ax[2, 2].plot(gain_rh, gain_rl, "o")
ax[2, 3].plot(gain_rh, gain_rh, "o")
ax[3, 0].plot(gain_rl, gain_ll, "o")
ax[3, 1].plot(gain_rl, gain_lh, "o")
ax[3, 2].plot(gain_rl, gain_rl, "o")
ax[3, 3].plot(gain_rl, gain_rh, "o")
fig.suptitle("Gain")
plt.show()


# check adds_per_group for nans
for channel in channels:
    if np.isnan(adds_per_group[channel]).any():
        print("There are nans in adds_per_group (variables.py)")
        # how many nans are there?
        print(np.sum(np.isnan(adds_per_group[channel])))
    else:
        print("No nans in adds_per_group (variables.py)")

# check for nans in adds_per_group in data
adds_per_group_data = {}
for channel in channels:
    adds_per_group[channel] = data[f"df_data/adds_per_group_{channel}"][()]
    if np.isnan(adds_per_group[channel]).any():
        print("There are nans in adds_per_group (main.py)")
        # how many nans are there?
        print(np.sum(np.isnan(adds_per_group[channel])))
    else:
        print("No nans in adds_per_group (main.py)")


print("data_quality")

data_quality = {}

for channel in channels:
    total_values = {}
    data_quality[channel] = np.array(fdq_sdf[f"fdq_sdf_{channel}/dq_data/data_quality"])
    print(f"{channel}: {data_quality[channel].shape}")

    for i in range(len(data_quality[channel])):
        unique_values = sorted(set(data_quality[channel][i].tolist()))

        for val in unique_values:
            total_values[val] = total_values.get(val, 0) + data_quality[channel][
                i
            ].tolist().count(val)

    # print total values sorted by the keys
    for key in sorted(total_values.keys()):
        print(f"{key}: {total_values[key]}")
