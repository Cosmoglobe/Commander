import numpy as np
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric, get_body, get_moon
import h5py

n_days = 1000
n = n_days * 24
mjd_start = 58400
mjd_arr = np.linspace(mjd_start, mjd_start + n_days, n + 1)

outname = "ephem_data.h5"
f2 = h5py.File(outname, "w")

bodies = ["jupiter", "moon", "mars", "sun", "venus", "saturn"]
n_bodies = len(bodies)
ra = np.zeros((n + 1, n_bodies))
dec = np.zeros_like(ra)
dis = np.zeros_like(ra)
with solar_system_ephemeris.set("builtin"):
    loc = coord.EarthLocation(lon=-118.283 * u.deg, lat=37.2313 * u.deg)
    for i in range(n + 1):
        t = Time(mjd_arr[i], format="mjd")
        for j, body in enumerate(bodies):
            jup = get_body(body, t, loc)
            ra[i, j] = jup.ra.deg
            dec[i, j] = jup.dec.deg
            dis[i, j] = jup.distance.AU
            if i % 24 == 0:
                print(i // 24, body)

for j, body in enumerate(bodies):
    data = np.transpose(np.array((mjd_arr, ra[:, j], dec[:, j], dis[:, j])))
    f2.create_dataset(body, data=data)

f2.close()
