import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

import matplotlib.colors as col
import matplotlib as mpl

wmap_K = hp.read_map(
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_Q1_v5.fits",
    field=(0, 1, 2),
)

transparent = False


comm_K = hp.read_map(
    "/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_noise_test/tod_040-WMAP_Q1_map_c0001_k000010.fits",
    field=(0, 1, 2),
)

cmap = col.ListedColormap(
    np.loadtxt("/mn/stornext/u3/duncanwa/c3pp/src/planck_cmap.dat") / 255.0, "planck"
)
mpl.cm.register_cmap(name="planck", cmap=cmap)


fontsize = 20
# dpi = 300
dpi = 100
mpl.rcParams.update({"font.size": fontsize})

# Official T sky maps for K and V
plt.figure(figsize=(16, 5))
hp.mollview(
    wmap_K[0],
    min=-0.25,
    max=0.25,
    cmap="planck",
    cbar=False,
    title="K",
    sub=121,
    unit="mK",
)

plt.savefig("wmap_T.png", bbox_inches="tight", transparent=transparent, dpi=dpi)

# Commander T sky maps for K and V, with big PRELIMINARY
fig = plt.figure(figsize=(16, 5))
hp.mollview(
    comm_K[0],
    min=-3.5,
    max=3.5,
    cmap="planck",
    cbar=False,
    title="K",
    sub=121,
    unit="mK",
)

# plt.text(0.5, 0.5, "PRELIMINARY", size=30, rotation=15., ha="center",
# va="center", transform=fig.transFigure,
#         bbox=dict(boxstyle="round",
#                  ec='black',
#                  fc='white',
#                  alpha=0.5))

plt.savefig("comm_T.png", bbox_inches="tight", transparent=transparent, dpi=dpi)

comm_K[0] = hp.remove_dipole(comm_K[0], gal_cut=30)

mu_WK = wmap_K[0].mean()
mu_CK = comm_K[0].mean()

fig = plt.figure(figsize=(16, 5))
hp.mollview(
    comm_K[0] - mu_CK + mu_WK,
    min=-0.25,
    max=0.25,
    cmap="planck",
    cbar=False,
    title="K",
    sub=121,
    unit="mK",
)

# plt.text(0.5, 0.5, "PRELIMINARY", size=30, rotation=15., ha="center",
# va="center", transform=fig.transFigure,
#         bbox=dict(boxstyle="round",
#                  ec='black',
#                  fc='white',
#                  alpha=0.5))

plt.savefig("comm_T_dipsub.png", bbox_inches="tight", transparent=transparent, dpi=dpi)


# Official P sky maps for K and V

wmap_K = hp.smoothing(wmap_K, fwhm=1 * np.pi / 180)

plt.figure(figsize=(12, 8))
hp.mollview(
    wmap_K[1], cmap="planck", sub=221, title="K $Q$", min=-0.1, max=0.1, unit="mK"
)
hp.mollview(
    wmap_K[2], cmap="planck", sub=222, title="K $U$", min=-0.1, max=0.1, unit="mK"
)
plt.savefig("wmap_P.png", bbox_inches="tight", transparent=transparent, dpi=dpi)

# Commander P sky maps for K and V, with big PRELIMINARY
comm_K = hp.smoothing(comm_K, fwhm=1 * np.pi / 180)
fig = plt.figure(figsize=(12, 8))
hp.mollview(
    comm_K[1], cmap="planck", sub=221, title="K $Q$", min=-0.1, max=0.1, unit="mK"
)
hp.mollview(
    comm_K[2], cmap="planck", sub=222, title="K $U$", min=-0.1, max=0.1, unit="mK"
)

# plt.text(0.5, 0.5, "PRELIMINARY", size=30, rotation=45., ha="center",
# va="center", transform=fig.transFigure,
#         bbox=dict(boxstyle="round",
#                  ec='black',
#                  fc='white',
#                  alpha=0.5))

plt.savefig("comm_P.png", bbox_inches="tight", transparent=transparent, dpi=dpi)

fig = plt.figure(figsize=(12, 8))
hp.mollview(
    comm_K[1] - wmap_K[1],
    cmap="planck",
    sub=221,
    title="K $Q$",
    min=-0.05,
    max=0.05,
    unit="mK",
)
hp.mollview(
    comm_K[2] - wmap_K[2],
    cmap="planck",
    sub=222,
    title="K $U$",
    min=-0.05,
    max=0.05,
    unit="mK",
)
plt.savefig("diff_P.png", bbox_inches="tight", transparent=transparent, dpi=dpi)

plt.close("all")
m = comm_K[0] - mu_CK + mu_WK - wmap_K[0]
m = hp.remove_dipole(m)
mu = m.mean()
sd = m.std()
hp.mollview(
    m * 1e3,
    min=-25,
    max=25,
    cmap="planck",
    cbar=True,
    title="",
    unit=r"$\mathrm{\mu K}$",
)
plt.savefig("diff.png", bbox_inches="tight", dpi=dpi, transparent=transparent)
plt.show()

# data = hp.read_map('/mn/stornext/d16/cmbco/bp/delivery/v8.00/BP8/goodness/BP_res_061-WMAP_V_P_QU_full_n16_0arcmin_uK_v1.fits', field=(0,1))
# I = np.zeros_like(data[0])
# m = np.array([I, data[0], data[1]])
# hp.write_map('res_V.fits', m)
# map2png res_V.fits -sig 2 -sig 3 -bar -range 10

plt.show()
