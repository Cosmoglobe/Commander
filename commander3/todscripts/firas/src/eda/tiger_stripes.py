import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
import globals as g
import my_utils as mu

data = np.load(g.PROCESSED_DATA_PATH, allow_pickle=True)

channels = ["ll", "rh"]
max_amp = {"ll": 25, "rh": 200}

ical = data['ical_ss'] 
pix_gal = data[f"pix_gal_ss"]
gmt = data["gmt_ss"]

# filter for ical bigger than tcmb
# test temps ical vs tcmb
ical_mask = ical > g.T_CMB
print(f"Filtering {np.sum(ical_mask)} pixels out of {len(ical)} where the temperature of the ICAL is bigger than T_CMB")
pix_gal_ical = pix_gal[ical_mask]

# filter for tcmb bigger than ical
xcal_mask = g.T_CMB > ical
print(f"Filtering {np.sum(xcal_mask)} pixels out of {len(ical)} where the temperature of the ICAL is smaller than T_CMB")
pix_gal_xcal = pix_gal[xcal_mask]


start = "89-326-1130"
start_dt = datetime.strptime(start, "%y-%j-%H%M")
mission_periods = np.array(
    [
        "89-327-2359",
        "89-343-0151",
        "90-019-0204",
        "90-080-0114",
        "90-128-2359",
        "90-139-1534",
        "90-193-1849",
        "90-207-1103",
        "90-208-1119",
        "90-220-0459",
        "90-264-0936",
    ]
)
mission_periods_dt = [datetime.strptime(period, "%y-%j-%H%M") for period in mission_periods]

period_filter = np.empty((len(mission_periods), len(gmt)), dtype=bool)
period_filter[0] = (gmt < mission_periods_dt[0]) & (
    gmt > start_dt
)

for i in range(len(mission_periods_dt) - 1):
    period_filter[i + 1] = (gmt < mission_periods_dt[i + 1]) & (
        gmt > mission_periods_dt[i]
    )

for channel in channels:
    sky = data[f"{channel}_ss"]
    f_ghz = mu.generate_frequencies(channel, "ss")

    sky_ical = sky[ical_mask]
    sky_xcal = sky[xcal_mask]

    m_ical = np.zeros((g.NPIX, len(f_ghz)))
    dd_ical = np.zeros(g.NPIX)

    for i, pix in enumerate(pix_gal_ical):
        m_ical[pix] += np.abs(sky_ical[i])
        dd_ical[pix] += 1
    mask_ical = dd_ical == 0
    m_ical[mask_ical] = np.nan
    monopole = mu.planck(f_ghz, np.array(g.T_CMB))
    m_ical[~mask_ical] = m_ical[~mask_ical] / dd_ical[~mask_ical, np.newaxis] - monopole

    for freqi in range(len(f_ghz)):
        hp.mollview(m_ical[:, freqi], title=f"{int(f_ghz[freqi]):04d} GHz", min=1, max=max_amp[channel], unit="MJy/sr", norm='log')
        plt.savefig(f"{g.SAVE_PATH}eda/{channel}/ical_{int(f_ghz[freqi]):04d}.png")
        plt.close()

    m_xcal = np.zeros((g.NPIX, len(f_ghz)))
    dd_xcal = np.zeros(g.NPIX)

    for i, pix in enumerate(pix_gal_xcal):
        m_xcal[pix] += np.abs(sky_xcal[i])
        dd_xcal[pix] += 1
    mask_xcal = dd_xcal == 0
    m_xcal[mask_xcal] = np.nan
    monopole = mu.planck(f_ghz, np.array(g.T_CMB))
    m_xcal[~mask_xcal] = m_xcal[~mask_xcal] / dd_xcal[~mask_xcal, np.newaxis] - monopole

    for freqi in range(len(f_ghz)):
        hp.mollview(m_xcal[:, freqi], title=f"{int(f_ghz[freqi]):04d} GHz", min=1, max=max_amp[channel], unit="MJy/sr", norm='log')
        plt.savefig(f"{g.SAVE_PATH}eda/{channel}/xcal_{int(f_ghz[freqi]):04d}.png")
        plt.close()

    # split up maps by their mission splits
    for j, period in enumerate(period_filter):
        m = np.zeros((g.NPIX, len(f_ghz)))
        dd = np.zeros(g.NPIX)

        sky_j = sky[period]
        for i, pix in enumerate(pix_gal[period]):
            m[pix] += np.abs(sky_j[i])
            dd[pix] += 1
        mask = dd == 0
        m[mask] = np.nan
        monopole = mu.planck(f_ghz, np.array(g.T_CMB))
        m[~mask] = m[~mask] / dd[~mask, np.newaxis] - monopole

        save_path = f"{g.SAVE_PATH}eda/{channel}/period{j}_"
        mu.save_mollview(m, f_ghz, save_path, max_amp[channel], norm='log')
