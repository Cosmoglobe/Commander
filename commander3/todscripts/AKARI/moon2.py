import numpy as np
import healpy as hp
from healpy.rotator import angdist

'''
Adapted from lines 1070--1145 of LongCorrect.cc

There are three basic components to the model;

1. An inner radius within which the Moon contamination is too bright, 21.5 degrees
2. Three rings with varying radii, centered roughly 120 degrees away from each other, maximum 34 degrees away from the moon.
3. A phase-dependent cut of scattered light, between 23 and 27 degrees away, and between 27 and 63 degrees in azimuth.

From figure 3 of Doi et al. 2015, it appears that the inner radius could be reduced to around 7 degrees without introducing too much contamination, as well as reducing the outer radius slightly. However this should be checked.
'''

DEG2RAD = np.pi / 180.0
RAD2DEG = 180.0 / np.pi

def wrap_pi(x):
    return (x + np.pi) % (2.0 * np.pi) - np.pi

def moon_flag_like_longcorrect_angdist(
    nside,
    moon_lon_deg, moon_lat_deg,
    # This essentially is a flag for whether or not the telescope is scanning
    # towards the NEP versus away from it. In the original code, it was
    # calculated by taking the difference between 
    to_nep=True,
    outer_cut = 35 * DEG2RAD,
    inner_cut = 21.5*DEG2RAD,
    theta1 = 24.5*DEG2RAD,
    phi1 = 59.5*DEG2RAD,
    theta2 = 29.4*DEG2RAD,
    phi2 = 179.4*DEG2RAD,
    theta3 = 24.5*DEG2RAD,
    phi3 = -59*DEG2RAD,
    d1l = 23.5*DEG2RAD, 
    d1u = 25.7*DEG2RAD,
    d2l = 27.9*DEG2RAD, 
    d2u = 29.2*DEG2RAD,
    d3l = 23.5*DEG2RAD, 
    d3u = 25.5*DEG2RAD,
):
    """
    Returns:
      moon_in: bool array [npix], True where sample/pixel should be masked by moon logic
      moon_avoid: angular distance from moon [rad]
    """
    npix = hp.nside2npix(nside)
    ipix = np.arange(npix)

    # Pixel angles in radians (theta=co-lat, phi=lon)
    theta, phi = hp.pix2ang(nside, ipix)

    # Convert to ecliptic-style lon/lat
    lon = phi
    lat = 0.5 * np.pi - theta

    # Moon position (rad)
    mlon = moon_lon_deg * DEG2RAD
    mlat = moon_lat_deg * DEG2RAD
    mtheta = 0.5 * np.pi - mlat
    mphi = mlon
    moon_vec = np.array(hp.ang2vec(mtheta, mphi))

    # Explicit healpy angular distance from moon
    moon_avoid = angdist([mtheta, mphi], [theta, phi])  # radians, shape (npix,)

    # Constants from LongCorrect for specific asymmetric feature
    moon23 = 23.0 * DEG2RAD
    moon27 = 26.5 * DEG2RAD
    phase27 = 27.0 * DEG2RAD
    phase63 = 63.0 * DEG2RAD


    moon_in = np.zeros(npix, dtype=int)

    in_calc = moon_avoid < outer_cut

    #moon_in |= in_calc & (moon_avoid < inner_cut)
    moon_in += (in_calc & (moon_avoid < inner_cut)).astype(int)
    print(moon_in)

    mid = in_calc & (moon_avoid < outer_cut) 
    if np.any(mid):
        sample_vecs = np.vstack(hp.ang2vec(theta[mid], phi[mid]))

        ref = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(ref, moon_vec)) > 0.99:
            ref = np.array([1.0, 0.0, 0.0])

        e1 = np.cross(ref, moon_vec)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(moon_vec, e1)

        x = sample_vecs.dot(e1)
        y = sample_vecs.dot(e2)
        phase = np.arctan2(y, x)
        if not to_nep:
            phase = wrap_pi(phase + np.pi)

        m = moon_avoid[mid]

        c1 = np.array([theta1, phi1])
        c2 = np.array([theta2, phi2])
        c3 = np.array([theta3, phi3])

        dirs = np.vstack((m, phase))
        d1 = hp.rotator.angdist(dirs, c1)
        d2 = hp.rotator.angdist(dirs, c2)
        d3 = hp.rotator.angdist(dirs, c3)

        lobe_hit = ((d1 > d1l) & (d1 < d1u)) | ((d2 > d2l) & (d2 < d2u)) | ((d3 > d3l) & (d3 < d3u))

        # extra phase gate in 23--27 deg radial range
        phase_gate = (m > moon23) & (m < moon27) & (phase > phase27) & (phase < phase63)

        #moon_in[mid] = lobe_hit | phase_gate
        # The numerical factors make it easier to visualize the asymmetry.
        moon_in[mid] += (phase_gate).astype(int)
        moon_in[mid] += ((d1 > d1l) & (d1 < d1u)).astype(int)
        moon_in[mid] += ((d2 > d2l) & (d2 < d2u)).astype(int)*2
        moon_in[mid] += ((d3 > d3l) & (d3 < d3u)).astype(int)*3

    return moon_in, moon_avoid


# Example usage
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    nside = 512
    moon_lon_deg = 0.0
    moon_lat_deg = 90.0
    to_nep = True

    # Note that to_nep is shorthand for the roll angle.

    mask, moon_avoid = moon_flag_like_longcorrect_angdist(
        nside=nside,
        moon_lon_deg=moon_lon_deg,
        moon_lat_deg=moon_lat_deg,
        to_nep=to_nep
    )

    hp.gnomview(mask.astype(float), rot=(moon_lon_deg, moon_lat_deg), reso=5,
                xsize=1000, cmap='turbo', title='To NEP')
    hp.graticule()
    plt.savefig('to_nep.png')
    plt.close()

    hp.gnomview(mask.astype(bool), rot=(moon_lon_deg, moon_lat_deg), reso=5,
                xsize=1000, cmap='binary', title='Moon flag')
    hp.graticule()
    plt.savefig('moon_flag.png')
    plt.close()

    mask, moon_avoid = moon_flag_like_longcorrect_angdist(
        nside=nside,
        moon_lon_deg=moon_lon_deg,
        moon_lat_deg=moon_lat_deg,
        to_nep=False
    )

    hp.gnomview(mask.astype(float), rot=(moon_lon_deg, moon_lat_deg), reso=5,
                xsize=1000,
                cmap='turbo',
                title='From NEP')
    hp.graticule()
    plt.savefig('from_nep.png')
    plt.close()
