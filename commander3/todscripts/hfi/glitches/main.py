import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import Commander.commander3.todscripts.hfi.glitches.detection as detection
import Commander.commander3.todscripts.hfi.glitches.globals as g
import matplotlib.pyplot as plt
import Commander.commander3.todscripts.hfi.glitches.templates as templates
from astropy.io import fits

chains_path = "/mn/stornext/d5/data/aimartin/hfi/chains/"
current_chain = "run3_nosample_owl44/"

# open a random TOD
scan_num = 500
iter = 1
filename = f"tod_{scan_num:06d}_samp{iter:06d}"

save_path = "/mn/stornext/d5/data/aimartin/hfi/figures/"

samprate = 180.3737  # Hz

data = h5py.File(f"{chains_path}{current_chain}{filename}.h5", "r")

res = data["res"][:]

for nsamp in range(0, res.shape[1], g.WINDOWSIZE):
    for i, band in enumerate(g.BAND_LIST):
        x_axis = np.linspace(i * g.WINDOWSIZE / samprate, (i + 1) * g.WINDOWSIZE / samprate,
                             num=g.WINDOWSIZE)
        plt.plot(res[i][nsamp : nsamp + g.WINDOWSIZE])
        plt.savefig(f"{save_path}debug/res_{band}_{nsamp}.png")
        plt.close()

        filtered = detection.threepointfilter(res[i][:nsamp + g.WINDOWSIZE])
        biggest_glitch = detection.search(filtered, x_axis)

        print("Biggest glitch sample in this window:", biggest_glitch)
        templates.fit_glitch(res[i], biggest_glitch, x_axis)

        quit()