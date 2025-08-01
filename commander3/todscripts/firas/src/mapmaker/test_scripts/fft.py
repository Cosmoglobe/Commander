import astropy.units as u
import funcs as f
import globals as g
import matplotlib.pyplot as plt
import numpy as np

z=[[1,2,34,45],[1,2,5,6],[7,8,9,10]]
x1=np.fft.irfft(z)
x2=np.fft.rfft(x1)
print(x2)
print(z)

x = np.linspace(0, 10, 100)
y = np.sin(x)
x1 = np.fft.irfft(y)
x2 = np.fft.rfft(x1)
# print(x2)
# print(y)

difference = np.abs(x2 - y)
plt.plot(x, difference)
plt.show()

dnu = 13.604162
frequencies = np.linspace(1e-5, dnu * g.SPEC_SIZE, g.SPEC_SIZE)
y = f.planck(frequencies, np.array(2.73))

x1 = np.fft.irfft(y)
x2 = np.fft.rfft(x1)

difference = np.abs(x2 - y)
plt.plot(frequencies, difference)
plt.show()

nu0_dust = 545 * u.GHz # Planck 2015
A_d = 163 * u.uK
T_d = 21 * u.K
beta_d = 1.53
y= f.dust(frequencies * u.GHz, A_d, nu0_dust, beta_d, T_d).value

x1 = np.fft.irfft(y)
x2 = np.fft.rfft(x1)
difference = np.abs(x2 - y)
plt.plot(frequencies, difference)
plt.show()