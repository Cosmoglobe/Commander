"""
This script compares the calculated temperature from the versions 4 and 5 of the pre-processing script.
"""

import astropy.units as u
import h5py
import matplotlib.pyplot as plt
import numpy as np

import globals as g

channel = "ll"

# Use 'r' mode and avoid loading unnecessary data
with h5py.File(f"{g.PREPROCESSED_DATA_PATH_SKY}", "r") as f:
    dataprep4 = f["df_data"]
    # Load only what you need and do it efficiently
    time4 = (dataprep4[f"midpoint_time_{channel}"][:] * (100 * u.ns)).to("s")
    # ical4 = (dataprep4["a_ical"][:] + dataprep4["b_ical"][:]) / 2
    # dihedral4 = (dataprep4["a_dihedral"][:] + dataprep4["b_dihedral"][:]) / 2
    # refhorn4 = (dataprep4["a_refhorn"][:] + dataprep4["b_refhorn"][:]) / 2
    # skyhorn4 = (dataprep4["a_skyhorn"][:] + dataprep4["b_skyhorn"][:]) / 2
    # collimator4 = (dataprep4["a_collimator"][:] + dataprep4["b_collimator"][:]) / 2
    # bolometer4 = (
    #     dataprep4[f"a_bol_assem_{channel}"][:] + dataprep4[f"b_bol_assem_{channel}"][:]
    # ) / 2
    ical4 = dataprep4["ical"][:]
    dihedral4 = dataprep4["dihedral"][:]
    refhorn4 = dataprep4["refhorn"][:]
    skyhorn4 = dataprep4["skyhorn"][:]
    collimator4 = dataprep4["collimator"][:]
    bolometer4 = dataprep4[f"bol_assem_{channel}"][:]

dataprep5 = np.load(f"{g.PREPROCESSED_DATA_PATH}/sky_{channel}.npz", allow_pickle=True)

# Extract data from npz
time5 = dataprep5["midpoint_time_s"]
ical5 = dataprep5["ical"]
dihedral5 = dataprep5["dihedral"]
refhorn5 = dataprep5["refhorn"]
skyhorn5 = dataprep5["skyhorn"]
collimator5 = dataprep5["collimator"]
bolometer5 = dataprep5[f"bolometer"]

print(f"Loaded {len(ical4)} records from dataprep4 and {len(ical5)} from dataprep5")

fig, ax = plt.subplots(3, 2, figsize=(12, 10), sharex=True)

ax.flatten()[0].plot(time4, ical4)
ax.flatten()[0].plot(time5, ical5, alpha=0.7)
ax.flatten()[0].set_ylabel("Ical (K)")
ax.flatten()[0].set_title(f"Ical Comparison for Channel {channel.upper()}")
ax.flatten()[0].grid(alpha=0.3)

ax.flatten()[1].plot(time4, dihedral4)
ax.flatten()[1].plot(time5, dihedral5, alpha=0.7)
ax.flatten()[1].set_ylabel("Dihedral (K)")
ax.flatten()[1].set_title(f"Dihedral Comparison for Channel {channel.upper()}")
ax.flatten()[1].grid(alpha=0.3)

ax.flatten()[2].plot(time4, refhorn4)
ax.flatten()[2].plot(time5, refhorn5, alpha=0.7)
ax.flatten()[2].set_ylabel("Refhorn (K)")
ax.flatten()[2].set_title(f"Refhorn Comparison for Channel {channel.upper()}")
ax.flatten()[2].grid(alpha=0.3)

ax.flatten()[3].plot(time4, skyhorn4)
ax.flatten()[3].plot(time5, skyhorn5, alpha=0.7)
ax.flatten()[3].set_ylabel("Skyhorn (K)")
ax.flatten()[3].set_title(f"Skyhorn Comparison for Channel {channel.upper()}")
ax.flatten()[3].grid(alpha=0.3)

ax.flatten()[4].plot(time4, collimator4)
ax.flatten()[4].plot(time5, collimator5, alpha=0.7)
ax.flatten()[4].set_ylabel("Collimator (K)")
ax.flatten()[4].set_title(f"Collimator Comparison for Channel {channel.upper()}")
ax.flatten()[4].grid(alpha=0.3)

ax.flatten()[5].plot(time4, bolometer4)
ax.flatten()[5].plot(time5, bolometer5, alpha=0.7)
ax.flatten()[5].set_ylabel("Bolometer (K)")
ax.flatten()[5].set_title(f"Bolometer Comparison for Channel {channel.upper()}")
ax.flatten()[5].grid(alpha=0.3)

plt.savefig(f"eda/output/compare_grts_channel_{channel}.png", dpi=300)
