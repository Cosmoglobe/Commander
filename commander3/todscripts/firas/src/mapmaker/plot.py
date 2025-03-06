import healpy as hp
import matplotlib.pyplot as plt
from sim import sim_dust

dust_map_downgraded_mjy, frequencies, signal = sim_dust()

# plot map for each frequency
for i in range(len(frequencies)):
    print(f"Plotting dust map for frequency {i}")
    dust_map = dust_map_downgraded_mjy * signal[i]
    hp.mollview(dust_map, title=f"{int(frequencies[i]):04d} GHz", unit="MJy/sr", min=0, max=200)
    plt.savefig(f"tests/dust_maps/{int(frequencies[i]):04d}.png")
    plt.close()
    plt.clf()
