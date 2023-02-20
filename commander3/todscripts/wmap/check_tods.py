import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
import h5py
import healpy as hp

from glob import glob

prefix = "/mn/stornext/d16/cmbco/bp/wmap/"

"""
I have noticed that the TODs after I process them are not behaving in the way I
expect them to. I will take some time to look at K1 to see if they have the sort
of behavior I would expect them to.
"""


files = glob(prefix + "tod/new/*.fits")
files.sort()
files = np.array(files)
file_input = files[0]


alpha = -1
fknee = 0.1

nside = 256
ntodsigma = 100
npsi = 2048
psiBins = np.linspace(0, 2 * np.pi, npsi)
fsamp = 30 / 1.536  # A single TOD record contains 30 1.536 second major science frames
chunk_size = 1875
nsamp = chunk_size * fsamp
chunk_list = np.arange(25)


gain_guesses = np.array(
    [
        -0.974,
        +0.997,
        +1.177,
        -1.122,
        +0.849,
        -0.858,
        -1.071,
        +0.985,
        +1.015,
        -0.948,
        +0.475,
        -0.518,
        -0.958,
        +0.986,
        -0.783,
        +0.760,
        +0.449,
        -0.494,
        -0.532,
        +0.532,
        -0.450,
        +0.443,
        +0.373,
        -0.346,
        +0.311,
        -0.332,
        +0.262,
        -0.239,
        -0.288,
        +0.297,
        +0.293,
        -0.293,
        -0.260,
        +0.281,
        -0.263,
        +0.258,
        +0.226,
        -0.232,
        +0.302,
        -0.286,
    ]
)


band = "K1"

labels = ["K113", "K114", "K123", "K124"]


data = fits.open(file_input)

TODs = []
for key in labels:
    TODs.append(data[2].data[key])


fig1, axes1 = plt.subplots(nrows=4, ncols=1, sharex=True, num="Raw")
fig2, axes2 = plt.subplots(nrows=4, ncols=1, sharex=True, num="Gain-corrected")
fig3, axes3 = plt.subplots(nrows=4, ncols=1, sharex=True, num="Gain-corrected - Mean")
DAs = []
for i in range(4):
    tod = np.zeros(TODs[i].size)
    for n in range(len(TODs[i][0])):
        tod[n :: len(TODs[i][0])] = TODs[i][:, n]
    axes1[i].plot(tod[:5000])
    axes2[i].plot(tod[:5000] / gain_guesses[i])
    dsub = (tod[:5000] / gain_guesses[i]) - tod[:5000].mean() / gain_guesses[i]
    axes3[i].plot(dsub)
    DAs.append(dsub)

    axes1[i].set_ylabel(labels[i])
    axes2[i].set_ylabel(labels[i])
fig1.suptitle("Raw data (du)")
fig2.suptitle("Gain-corrected data (mK)")
fig3.suptitle("Gain-corrected mean-subtracted (mK)")


fig1.savefig("raw_tods.png", bbox_inches="tight")
fig2.savefig("temp_tods.png", bbox_inches="tight")
fig3.savefig("temp_minus_mu_tods.png", bbox_inches="tight")

DAs = np.array(DAs)
d = DAs.mean(axis=0)
plt.figure("Intensity signal")
plt.plot(d[:5000])
plt.title("Intensity signal")
plt.savefig("total_intensity.png", bbox_inches="tight")

d1 = 0.5 * (DAs[0] + DAs[1])
d2 = 0.5 * (DAs[2] + DAs[3])
fig, axes = plt.subplots(nrows=2, sharex=True, sharey=False, num="DA Combinations")
axes[0].plot(d1[:5000])
axes[1].plot(d2[:5000])
axes[0].set_ylabel("K113+K114")
axes[1].set_ylabel("K123+K124")
plt.savefig("horns.png", bbox_inches="tight")

plt.show()
