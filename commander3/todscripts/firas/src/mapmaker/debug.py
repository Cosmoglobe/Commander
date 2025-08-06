import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from sim import sim_dust

dust_map_downgraded_mjy, frequencies, signal = sim_dust()
dust_map = dust_map_downgraded_mjy[:, np.newaxis] * signal[np.newaxis, :]

d = np.load("test_output/ifgs.npz")["ifg"]
m = np.abs(np.fft.rfft(d, axis=1))

ratio = dust_map/m
print(ratio.shape)
print(ratio)
print(f"max: {np.max(ratio)}, min: {np.min(ratio)}, mean: {np.mean(ratio)}, std: {np.std(ratio)}")

hp.mollview(
    ratio[:, 40],
    title=f"Ratio at {int(frequencies[40].value):04d} GHz",
    unit="MJy/sr",
    min=0.5,
    max=1.5,
)
plt.show()