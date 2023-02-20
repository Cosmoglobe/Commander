import healpy as hp
import cosmoglobe
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt


w = 5
f = {"rlabel": 6, "llabel": 6}

f = 12


def mylabel(rlabel, fontsize):
    ax = plt.gca()
    plt.text(
        0.025,
        1.05,
        rlabel,
        ha="left",
        va="center",
        fontsize=fontsize,
        transform=ax.transAxes,
    )


try:
    d_100_s = hp.read_map("npipe_100_2deg.fits", field=(0, 1, 2))
except FileNotFoundError:
    d_100 = hp.read_map(
        "/mn/stornext/d16/cmbco/ola/npipe/freqmaps/npipe6v20_100_map_K.fits",
        field=(0, 1, 2),
    )
    d_100_s = 1e6 * hp.ud_grade(hp.smoothing(d_100, fwhm=2 * np.pi / 180), 512)
    hp.write_map("npipe_100_2deg.fits", d_100_s)


DIR = "/mn/stornext/d16/www_cmb/dwatts/v0"
d_W1 = hp.read_map(f"{DIR}/BP_090-WMAP_W1_IQU_n0512_v0.fits", field=(0, 1, 2, 6, 7, 8))
d_W2 = hp.read_map(f"{DIR}/BP_090-WMAP_W2_IQU_n0512_v0.fits", field=(0, 1, 2, 6, 7, 8))
d_W3 = hp.read_map(f"{DIR}/BP_090-WMAP_W3_IQU_n0512_v0.fits", field=(0, 1, 2, 6, 7, 8))
d_W4 = hp.read_map(f"{DIR}/BP_090-WMAP_W4_IQU_n0512_v0.fits", field=(0, 1, 2, 6, 7, 8))
d_W = (
    d_W1[:3] / d_W1[3:] ** 2
    + d_W2[:3] / d_W2[3:] ** 2
    + d_W3[:3] / d_W3[3:] ** 2
    + d_W4[:3] / d_W4[3:] ** 2
) / (1 / d_W1[3:] ** 2 + 1 / d_W2[3:] ** 2 + 1 / d_W3[3:] ** 2 + 1 / d_W4[3:] ** 2)


d_W_s = hp.smoothing(d_W, fwhm=2 * np.pi / 180)
d_W_s = 1e3 * d_W_s

d_W_orig = hp.read_map(
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_band_iqusmap_r9_9yr_W_v5.fits",
    field=(0, 1, 2),
)
d_W_o_s = 1e3 * hp.smoothing(d_W_orig, fwhm=2 * np.pi / 180)

cosmoglobe.standalone_colorbar(
    "planck", ticks=[-10, 0, 10], unit=r"$\mathrm{\mu K}$", extend="both"
)
plt.savefig("cbar.png", bbox_inches="tight", dpi=300)

cosmoglobe.plot(-d_100_s + d_W_o_s, sig=1, min=-10, max=10, cbar=False, width=w)
mylabel("$\mathit{WMAP9}-\mathit{Planck}\ 100\,\mathrm{GHz}$", f)
plt.savefig("diff_orig.png", bbox_inches="tight", dpi=300)
cosmoglobe.plot(-d_100_s + d_W_s, sig=1, min=-10, max=10, width=w, cbar=False)
mylabel(r"$\mathrm{Watts\ et\ al.}-\mathit{Planck}\ 100\,\mathrm{GHz}$", f)
plt.savefig("diff_cg.png", bbox_inches="tight", dpi=300)
cosmoglobe.plot(
    d_W_s, sig=1, min=-10, max=10, unit=r"$\mathrm{\mu K}$", cbar=False, width=w
)
mylabel(r"$\mathrm{Watts\ et\ al.}$", f)
plt.savefig("W_cg.png", bbox_inches="tight", dpi=300)
cosmoglobe.plot(
    d_W_o_s,
    sig=1,
    min=-10,
    max=10,
    unit=r"$\mathrm{\mu K}$",
    cbar=False,
    width=w,
    llabel=r"\mathit{WMAP9}",
)
plt.savefig("W_WMAP_cg.png", bbox_inches="tight", dpi=300)
cosmoglobe.plot(
    d_100_s,
    sig=1,
    min=-10,
    max=10,
    unit=r"$\mathrm{\mu K}$",
    cbar=True,
    width=w,
    llabel=r"\mathit{Planck}\ 100\,\mathrm{GHz}",
)
plt.savefig("npipe_100.png", bbox_inches="tight", dpi=300)
