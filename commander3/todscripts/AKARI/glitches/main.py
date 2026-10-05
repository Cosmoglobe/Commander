import h5py
import numpy as np
from eirik_flags import load_flags
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import detection

# open a random TOD
path = "/mn/stornext/d23/cmbco/globe/akari/tod/eirik_newhdf/v7/n8192/"
file = "AKARI_090_000048.h5"
tod_name = "009531"
detector = "AKARI_090-27"

save_path = "/mn/stornext/d23/cmbco/akari/aimartin/figures/"

nsamples = 2000
samprate = 25.28 # Hz

with h5py.File(path + file, "r") as f:
    # print the keys in the file
    tod = f[f"{tod_name}/{detector}/tod"][:]

# get flags
flagnums = np.arange(0, 21)
# remove the flags that are related to glitches
flagnums = np.delete(flagnums, [7, 11, 12, 13, 14])
flags = np.array(load_flags(path + file, flagnums=flagnums, detector=detector, chunk=tod_name))

print(flags)
plt.imshow(flags, aspect='auto', interpolation='none')
plt.xlabel("Sample")
plt.ylabel("Flag Number")
plt.title("Flags for the TOD")
plt.colorbar(label='Flag Value')
plt.savefig(save_path + "flags.png")
plt.clf()

allflags = np.any(flags, axis=0) # combine all flags into one boolean array
print(allflags)

# plot the TOD + flags + possible glitch
# t = np.arange(len(tod)) / samprate

# plt.vlines(np.array(flags[:nsamples]), ymin=np.min(tod[:nsamples]), ymax=np.max(tod[:nsamples]),
#            color='grey', alpha=0.5)
# plt.plot(t[:nsamples], tod[:nsamples])

# circle = plt.Circle((29, 1), 2 , color='red', fill=False)
# plt.gca().add_patch(circle)

# # plt.ylim(-2, 2)
# plt.xlabel("Time (s)")
# plt.ylabel("Signal")
# plt.title("Raw TOD")

# plt.savefig(save_path + "raw_tod.png")
# plt.clf()

filtered = detection.threepointfilter(tod[:nsamples])
detection.search(filtered, allflags)