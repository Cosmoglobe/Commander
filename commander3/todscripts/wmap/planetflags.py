import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import h5py
import huffman

band = "K1"
version = 11


fnames = glob(f"/mn/stornext/d16/cmbco/bp/wmap/data/wmap_{band}_*v{version}.h5")
fnames.sort()

fname = np.random.choice(fnames)
labels = [f"{band}13", f"{band}14", f"{band}23", f"{band}24"]
f = h5py.File(fname, "r")
obsid = str(list(f.keys())[0])
huffTree = f[obsid + "/common/hufftree"]
huffSymb = f[obsid + "/common/huffsymb"]
h = huffman.Huffman(tree=huffTree, symb=huffSymb)


flags = [[], [], [], []]
for num, label in enumerate(labels):
    flag = h.Decoder(np.array(f[obsid + "/" + label + "/flag"]))
    flags[num] = flags[num] + flag.tolist()

flags = np.array(flags)

for i in range(len(flags)):
    plt.plot(flags[i])
    plt.title(labels[i])

plt.figure()
mus = np.mean(flags, axis=0)
sds = np.std(flags, axis=0)
plt.plot(np.log(mus) / np.log(2), ".")

mus = np.where(mus > 0, 1, 0)
plt.figure()
plt.plot(mus, ".")
print(np.unique(flags))
print(np.sum(flags == 0) / flags.size)
plt.show()
