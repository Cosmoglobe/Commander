import numpy as np 
import matplotlib.pyplot as plt
from astropy.io import fits
import tools 
import healpy as hp

data = fits.open('lambda_chipass_healpix_r10.fits')

print(data[0].header)

print(data[1].header)

# nside = 512  # 256
# hpxmap = np.load('maskmap_%i.npy' % nside) 
# #hpxmap = np.load('map_%i.npy' % nside) 
# # np.save('nhit_%i' % nside, nhit)
# maxval = 2                                                                                                                        
# hp.mollview(np.log10(hpxmap + 2), min=0, max=maxval) 
# #plt.savefig('chipass_binned.pdf')
# plt.show()

# t = 123 * 60 + 14.54

# h = t // 3600
# m = (t - h * 3600) // 60
# s = (t - h * 3600 - m * 60) 
# print(h, m, s)

# data_raw = fits.open('1999-07-26_1427_344p042229_17a_l1.sdfits')
# data_cal = fits.open('1999-07-26_1427_344p042229_17a_l2.sdfits')


# def median_filter(tod, n_med=50):
#     n = tod.shape[0]
#     out = np.zeros_like(tod)
#     print(out.shape)
#     n_half = n_med // 2
#     for i in range(n):
#         low = max(i-n_half, 0)
#         hi = min(i+n_half, n)
#         out[i] = np.nanmedian(tod[low:hi], axis=0)
#     return out
# date = np.array(data_raw[1].data['date-obs']).reshape(-1,13)[:, 0]
# t = data_cal[1].data['time'].reshape(-1,13)[:, 0]
# print(t / 3600)
# h = t // 3600
# m = (t - h * 3600) // 60
# s = (t - h * 3600 - m * 60) 
# print(s)
# timestr = []
# for i in range(len(h)):
#     timestr.append(date[0] + 'T' + '%02i:%02i:%06.3f' % (h[i], m[i], s[i]))
# print(timestr)

# az = data_cal[1].data['AZIMUTH'].reshape(-1,13)[:, 0]
# el = data_cal[1].data['ELEVATIO'].reshape(-1,13)[:, 0]
# ra = data_cal[1].data['obj-ra'].reshape(-1,13)[:, 0]
# dec = data_cal[1].data['obj-dec'].reshape(-1,13)[:, 0]
# # print(az.shape)
# # sys.exit()
# bla = np.array(data_cal[1].data['Data'][:,0,0,:,:])
# cal_cube = np.reshape(bla, (-1, 13, 2, 1024))

# bla = np.array(data_raw[1].data['Data'][:,0,0,:,:])
# raw_cube = np.reshape(bla, (-1, 13, 2, 1024))

# tsys = np.array(data_raw[1].data['Tsys'])
# tsys_raw = np.reshape(tsys, (-1, 2, 13)).transpose((0, 2, 1))
# # tsys_raw = np.reshape(tsys, (-1, 13, 2))  # this should be the correct one

# # Lat/Lon of Parkes telescope -32.99778, 148.26292
# loc = coord.EarthLocation(lon=148.26292 * u.deg, lat=-32.99778 * u.deg)
# time = Time(timestr, format='isot', scale='utc')
# aa = AltAz(location=loc, obstime=time)
# newAltAzcoordiantes = SkyCoord(alt=el * u.deg, az=az * u.deg, obstime = time, frame = 'altaz', location = loc)
# print(newAltAzcoordiantes.galactic, el[0], az[0], ra[0], dec[0])
# print(tsys_raw.shape)
# print(raw_cube.shape)
# ratio = raw_cube[:, :, :, :] / tsys_raw[:, :, :, None]
# bst = median_filter(ratio)
# bs = median_filter(raw_cube[:, :, :, :])
# bt = median_filter(tsys_raw)
# s_prime = raw_cube / bst - bt[:, :, :, None]
# # for i in range(20, 40, 5):
# #     plt.plot(s_prime[i, 0, 0, :], 'r')
# #     plt.plot(cal_cube[i, 0, 0, :], 'b')
# #     plt.ylim(-1.5, 1.5)
# #     print(i)
# #     plt.show()
# # plt.plot(raw_cube[:, 0, 0, 50])

# plt.plot(cal_cube[:, 5, 0, 50])
# plt.plot(s_prime[:, 5, 0, 50])
# x = np.linspace(-5, 5, 100)
# plt.figure()
# plt.scatter(s_prime.flatten(), cal_cube.flatten(), s=1, rasterized=True)
# plt.plot(x, x, '--k', lw=3)
# plt.xlabel('our_calibration')
# plt.ylabel('their')
# plt.xlim(-15, 25)
# plt.ylim(-15, 25)
# plt.savefig('wrong_tsys.pdf')
# plt.show()

# nu = 1.3945e9 # Hz
# area = 2.0 * np.pi * (14.0 / 60 / 180 * np.pi / np.sqrt(8 * np.log(2))) ** 2

# dnu = 64e6 / 1024
# jy2k = tools.jy2kcmb(1, nu, area)
# print(jy2k)
# data_cal = fits.open('1999-07-26_1427_344p042229_17a_l2.sdfits')
# t = data_cal[1].data['time'].reshape(-1,13)[:, 0]
# t = t - t[0]
# dt = t[1] - t[0]
# ex = data_cal[1].data['object'].reshape(-1,13)[:, 0]
# print(ex)
# print(dt)
# az = data_cal[1].data['AZIMUTH'].reshape(-1,13)[:, 0]
# el = data_cal[1].data['ELEVATIO'].reshape(-1,13)[:, 0]
# object = data_cal[1].data['object'].reshape(-1,13)[:, 2]
# print(object)
# ra = data_cal[1].data['obj-ra'].reshape(-1,13)[:, 2]
# dec = data_cal[1].data['obj-dec'].reshape(-1,13)[:, 0]
# print(ra.shape, az.shape)
# print(data_cal[1].columns)
# bla = data_cal[1].data['Data'][:,0,0,0,:]
# cal_cube = np.reshape(bla, (-1, 13, 1024))
# tsys = data_cal[1].data['Tsys']
# # print("tsys", tsys.shape)
# tsys_cube = np.reshape(tsys, (-1, 13, 2))
# # print()
# print(cal_cube.shape)
# print(tsys_cube.shape)
# # print(ra)
# # plt.scatter(ra, dec)

# plt.figure()
# plt.plot(t, (el - el.mean()), label='elevation - <el>')
# plt.plot(t, (az - az.mean()), label='azimuth - <az>')
# plt.plot(t, (ra - ra.mean()), label='ra - <ra>')
# plt.plot(t, (dec - dec.mean()), label='dec - <dec>')
# plt.plot(t, cal_cube[:, 0, 400] * jy2k, label='single channel')
# plt.plot(t, cal_cube[:, 0, 250:].mean(1) * jy2k, label='band average')
# plt.ylabel(r'Data [K${}_\mathrm{cmb}$]')
# plt.xlabel('Time [s]')
# plt.legend()
# plt.savefig('chipass_timestream.pdf', bbox_inches='tight')

# plt.figure()
# plt.plot(tsys_cube[:, :, 0])
# plt.show()

# nu = 1.4e9 # Hz

# flux = 1 # Jy

# area = 2.0 * np.pi * (14.0 / 60 / 180 * np.pi / np.sqrt(8 * np.log(2))) ** 2

# print(tools.jy2kcmb(flux, nu, area))