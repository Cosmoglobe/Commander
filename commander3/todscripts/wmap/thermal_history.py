import numpy as np
import matplotlib.pyplot as plt

from astropy.time import Time

from scipy.interpolate import interp1d

t_i = 2452131.50000
t_f = 2455418.50000


t_out = np.linspace(t_i, t_f, 10000)

t, T_rxb = np.loadtxt('gain_params/extracted_rxb.txt', delimiter=',').T
# starts on day 222 of 2001
t = Time(t + 222/365, format='jyear')

#t = Time(t, format='jd')
plt.plot(t.datetime, T_rxb, '.', ms=1)
plt.gcf().autofmt_xdate()  # orient date labels at a slant  
plt.figure()

plt.plot(t.jd, T_rxb, '.', ms=1)
plt.gcf().autofmt_xdate()  # orient date labels at a slant  


f = interp1d(t.jd, T_rxb)
plt.plot(t_out, f(t_out))


t, T_fpa = np.loadtxt('gain_params/extracted_fpa.txt', delimiter=',').T
t = Time(t + 222/365, format='jyear')
plt.figure()
plt.plot(t.jd, T_fpa, '.', ms=1)
plt.gcf().autofmt_xdate()  # orient date labels at a slant  
f = interp1d(t.jd, T_fpa, fill_value='extrapolate')
plt.plot(t_out, f(t_out))

plt.show()
