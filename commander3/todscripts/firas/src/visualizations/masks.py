"""
Script to check what each mask looks like.
"""
import os

import healpy as hp
import matplotlib.pyplot as plt

path = "/mn/stornext/d16/cmbco/ola/masks"

dir_list = os.listdir(path)
for file in dir_list:
    if file.endswith(".fits"):
        print(f"Loading mask: {file}")
        hp.mollview(hp.read_map(os.path.join(path, file)), title=file, unit="1", cmap="gray")
        hp.graticule()
        plt.savefig("./masks/" + file.replace(".fits", ".png"), dpi=300)
        plt.clf()