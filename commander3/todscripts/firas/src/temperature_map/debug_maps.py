import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

# m = hp.fitsfunc.read_map(
#     "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/maps/fits/0179.fits",
# )

# hp.mollview(m, coord="G", title="179", unit="MJy/sr", min=0, max=500)
# plt.show()


sky = np.load("../../output/data/sky.npy")
data = h5py.File(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v1.h5",
    "r",
)

print(data["df_data"].keys())

spec_bad = np.empty((len(sky)))
spec_good = np.empty((len(sky)))
adds_per_group_bad = np.empty((len(sky)))
adds_per_group_good = np.empty((len(sky)))
bol_cmd_bias_bad = np.empty((len(sky)))
bol_cmd_bias_good = np.empty((len(sky)))
bol_volt_bad = np.empty((len(sky)))
bol_volt_good = np.empty((len(sky)))
dihedral_bad = np.empty((len(sky)))
dihedral_good = np.empty((len(sky)))
gain_bad = np.empty((len(sky)))
gain_good = np.empty((len(sky)))
ical_bad = np.empty((len(sky)))
ical_good = np.empty((len(sky)))
mtm_length_bad = np.empty((len(sky)))
mtm_length_good = np.empty((len(sky)))
mtm_speed_bad = np.empty((len(sky)))
mtm_speed_good = np.empty((len(sky)))
sweeps_bad = np.empty((len(sky)))
sweeps_good = np.empty((len(sky)))

spec_bad[np.max(sky[:, 1:]) > 1000] = np.mean(np.abs(sky[:, 1:]), axis=1)
spec_bad[np.max(sky[:, 1:]) <= 1000] = np.nan
spec_good[np.max(sky[:, 1:]) <= 1000] = np.mean(np.abs(sky[:, 1:]), axis=1)
spec_good[np.max(sky[:, 1:]) > 1000] = np.nan

adds_per_group_bad[np.max(sky[:, 1:]) > 1000] = data["df_data/adds_per_group_ll"][()]
adds_per_group_bad[np.max(sky[:, 1:]) <= 1000] = np.nan
adds_per_group_good[np.max(sky[:, 1:]) <= 1000] = data["df_data/adds_per_group_ll"][()]
adds_per_group_good[np.max(sky[:, 1:]) > 1000] = np.nan

bol_cmd_bias_bad[np.max(sky[:, 1:]) > 1000] = data["df_data/bol_cmd_bias_ll"][()]
bol_cmd_bias_bad[np.max(sky[:, 1:]) <= 1000] = np.nan
bol_cmd_bias_good[np.max(sky[:, 1:]) <= 1000] = data["df_data/bol_cmd_bias_ll"][()]
bol_cmd_bias_good[np.max(sky[:, 1:]) > 1000] = np.nan

bol_volt_bad[np.max(sky[:, 1:]) > 1000] = data["df_data/bol_volt_ll"][()]
bol_volt_bad[np.max(sky[:, 1:]) <= 1000] = np.nan
bol_volt_good[np.max(sky[:, 1:]) <= 1000] = data["df_data/bol_volt_ll"][()]
bol_volt_good[np.max(sky[:, 1:]) > 1000] = np.nan

dihedral_bad[np.max(sky[:, 1:]) > 1000] = data["df_data/dihedral"][()]
dihedral_bad[np.max(sky[:, 1:]) <= 1000] = np.nan
dihedral_good[np.max(sky[:, 1:]) <= 1000] = data["df_data/dihedral"][()]
dihedral_good[np.max(sky[:, 1:]) > 1000] = np.nan

gain_bad[np.max(sky[:, 1:]) > 1000] = data["df_data/gain_ll"][()]
gain_bad[np.max(sky[:, 1:]) <= 1000] = np.nan
gain_good[np.max(sky[:, 1:]) <= 1000] = data["df_data/gain_ll"][()]
gain_good[np.max(sky[:, 1:]) > 1000] = np.nan

ical_bad[np.max(sky[:, 1:]) > 1000] = data["df_data/ical"][()]
ical_bad[np.max(sky[:, 1:]) <= 1000] = np.nan
ical_good[np.max(sky[:, 1:]) <= 1000] = data["df_data/ical"][()]
ical_good[np.max(sky[:, 1:]) > 1000] = np.nan

mtm_length_bad[np.max(sky[:, 1:]) > 1000] = data["df_data/mtm_length"][()]
mtm_length_bad[np.max(sky[:, 1:]) <= 1000] = np.nan
mtm_length_good[np.max(sky[:, 1:]) <= 1000] = data["df_data/mtm_length"][()]
mtm_length_good[np.max(sky[:, 1:]) > 1000] = np.nan

mtm_speed_bad[np.max(sky[:, 1:]) > 1000] = data["df_data/mtm_speed"][()]
mtm_speed_bad[np.max(sky[:, 1:]) <= 1000] = np.nan
mtm_speed_good[np.max(sky[:, 1:]) <= 1000] = data["df_data/mtm_speed"][()]
mtm_speed_good[np.max(sky[:, 1:]) > 1000] = np.nan

sweeps_bad[np.max(sky[:, 1:]) > 1000] = data["df_data/sweeps_ll"][()]
sweeps_bad[np.max(sky[:, 1:]) <= 1000] = np.nan
sweeps_good[np.max(sky[:, 1:]) <= 1000] = data["df_data/sweeps_ll"][()]
sweeps_good[np.max(sky[:, 1:]) > 1000] = np.nan

fig, ax = plt.subplots(5, 2, figsize=(10, 10), sharex=True)
ax[0, 0].plot(spec_bad, color="red")
ax[0, 0].plot(spec_good, color="green")
ax[0, 0].set_title("spec")
ax[0, 1].plot(adds_per_group_bad, color="red")
ax[0, 1].plot(adds_per_group_good, color="green")
ax[0, 1].set_title("adds_per_group")
ax[1, 0].plot(bol_cmd_bias_bad, color="red")
ax[1, 0].plot(bol_cmd_bias_good, color="green")
ax[1, 0].set_title("bol_cmd_bias")
ax[1, 1].plot(bol_volt_bad, color="red")
ax[1, 1].plot(bol_volt_good, color="green")
ax[1, 1].set_title("bol_volt")
ax[2, 0].plot(dihedral_bad, color="red")
ax[2, 0].plot(dihedral_good, color="green")
ax[2, 0].set_title("dihedral")
ax[2, 1].plot(gain_bad, color="red")
ax[2, 1].plot(gain_good, color="green")
ax[2, 1].set_title("gain")
ax[3, 0].plot(ical_bad, color="red")
ax[3, 0].plot(ical_good, color="green")
ax[3, 0].set_title("ical")
ax[3, 1].plot(mtm_length_bad, color="red")
ax[3, 1].plot(mtm_length_good, color="green")
ax[3, 1].set_title("mtm_length")
ax[4, 0].plot(mtm_speed_bad, color="red")
ax[4, 0].plot(mtm_speed_good, color="green")
ax[4, 0].set_title("mtm_speed")
ax[4, 1].plot(sweeps_bad, color="red")
ax[4, 1].plot(sweeps_good, color="green")
ax[4, 1].set_title("sweeps")

plt.tight_layout()
plt.show()
