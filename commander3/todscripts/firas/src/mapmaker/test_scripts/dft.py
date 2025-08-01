"""
Script to do some experiments on how to use the DFT matrix.
"""

import matplotlib.pyplot as plt
import numpy as np

import globals as g

W = np.zeros((g.IFG_SIZE, g.IFG_SIZE), dtype=complex)
W[0, :] = 1
W[:, 0] = 1
omega = np.exp(-2j * np.pi / g.IFG_SIZE)
for xi in range(1, g.IFG_SIZE):
    for nui in range(1, g.IFG_SIZE):
        W[nui, xi] = omega ** ((xi * nui) % g.IFG_SIZE) # the mod operator just avoids calculating high exponents
W = W #/ np.sqrt(g.IFG_SIZE)

plt.imshow(W.real, cmap='gray')
plt.show()

Wnsq = np.zeros((g.SPEC_SIZE, g.IFG_SIZE), dtype=complex)
Wnsq[0, :] = 1
Wnsq[:, 0] = 1
omega = np.exp(-2j * np.pi / g.IFG_SIZE)
for xi in range(1, g.IFG_SIZE):
    for nui in range(1, g.SPEC_SIZE):
        Wnsq[nui, xi] = omega ** ((xi * nui) % g.IFG_SIZE) # the mod operator just avoids calculating high exponents
Wnsq = Wnsq #/ g.IFG_SIZE

plt.imshow(Wnsq.real, cmap='gray')
plt.show()

x = np.linspace(0, 10, g.IFG_SIZE)
y = np.sin(x)

# compare W matrix with numpy fft and rfft
fft = np.fft.rfft(y)
freq = np.fft.fftfreq(g.IFG_SIZE, d=g.IFG_SIZE / g.SPEC_SIZE)
rfft = np.fft.rfft(y, n=g.IFG_SIZE)
rfreq = np.fft.rfftfreq(g.IFG_SIZE, d=g.IFG_SIZE / g.SPEC_SIZE)
Wfft = np.dot(W, y)
Wnsqfft = np.dot(Wnsq, y)

print(Wfft.real)
# print shapes
print(f"fft: {fft.shape}")
print(f"freq: {freq.shape}")
print(f"rfft: {rfft.shape}")
print(f"rfreq: {rfreq.shape}")
print(f"Wfft: {Wfft.shape}")

# plot
fig, ax = plt.subplots(2, 1)
ax[0].plot(freq, Wfft.real, label='Wfft', alpha=0.5)
# ax[0].plot(freq, fft.real, label='fft', alpha=0.5)
ax[1].plot(rfreq, Wnsqfft.real, label='Wnsq', alpha=0.5)
ax[1].plot(rfreq, rfft.real, label='rfft', alpha=0.5)
fig.legend()
plt.show()

assert np.allclose(Wnsqfft.real, rfft.real)

# now compare the inverse operations
IW = np.zeros((g.IFG_SIZE, g.IFG_SIZE), dtype=complex)
IW[0, :] = 1
IW[:, 0] = 1
omega = np.exp(2j * np.pi / g.IFG_SIZE)
for xi in range(1, g.IFG_SIZE):
    for nui in range(1, g.IFG_SIZE):
        IW[xi, nui] = omega ** ((xi * nui) % g.IFG_SIZE) # the mod operator just avoids calculating high exponents
IW = IW / g.IFG_SIZE

plt.imshow(IW.real, cmap='gray')
plt.colorbar()
plt.title('IW')
plt.show()

# check if the square dtf matrices are identity when multiplied
plt.imshow(np.dot(IW, W).real, cmap='gray')
plt.colorbar()
plt.show()

assert np.allclose(np.dot(IW, W), np.identity(g.IFG_SIZE))

IWnsq = np.zeros((g.IFG_SIZE, g.SPEC_SIZE), dtype=complex)
IWnsq[0, :] = 1
# IWnsq[:, 0] = 1
omega = np.exp((2j * np.pi) / g.IFG_SIZE)
for xi in range(1, g.IFG_SIZE):
    for nui in range(0, g.SPEC_SIZE):
        IWnsq[xi, nui] = omega ** ((xi * nui))#% g.IFG_SIZE) # the mod operator just avoids calculating high exponents
IWnsq = (IWnsq-1) / (g.SPEC_SIZE)

plt.imshow(IWnsq.real, cmap='gray')
plt.title('IWnsq')
plt.colorbar()
plt.show()

# check size
print(f"Wnsq: {Wnsq.shape}")
IWnsq_H = np.asarray((np.matrix(Wnsq).H -1) / g.SPEC_SIZE)
# check size
print(f"IWnsq_H: {IWnsq_H.shape}")

ifft = np.fft.irfft(fft)

# check if ifft and y are the same
difference = np.abs(ifft) - y
ratio = np.abs(ifft) / y
print(f"Difference: {difference}")
print(f"Ratio: {ratio}")
plt.plot(x, difference)
plt.title('Difference')
plt.show()
plt.plot(x, ratio)
plt.title('Ratio')
plt.show()



Wifft = np.dot(IW, Wfft)
irfft = np.fft.irfft(rfft, n=g.IFG_SIZE)
Wifftnsq = np.dot(IWnsq, Wnsqfft)
# check sizes
print(f"Wnsqfft: {Wnsqfft.shape}")
print(f"IWnsq_H: {IWnsq_H.shape}")
Wifftnsq_H = np.squeeze(np.dot(IWnsq_H, Wnsqfft))
# check size
print(f"Wifftnsq_H: {Wifftnsq_H.shape}")

fig, ax = plt.subplots(2, 1)
ax[0].plot(x, Wifft.real, label='Wifft')
ax[0].plot(x, ifft.real, label='ifft')
ax[1].plot(x, Wifftnsq.real, label='Wifftnsq')
ax[1].plot(x, irfft.real, label='irfft')
ax[1].plot(x, Wifftnsq_H.real, label='Wifftnsq_H')
fig.legend()
plt.show()

fig, ax = plt.subplots(2, 1)
ax[0].plot(x, Wifftnsq.real - irfft.real)
ax[0].set_title('Wifftnsq - irfft')
ax[1].plot(x, Wifftnsq_H.real - irfft.real)
ax[1].set_title('Wifftnsq_H - irfft')
plt.show()

# check if it's identity
print(np.allclose(np.dot(Wnsq, IWnsq), np.identity(g.SPEC_SIZE)))
print(np.dot(Wnsq, IWnsq))
# print properly
matrix = np.dot(IWnsq, Wnsq)
plt.imshow(matrix.real, cmap='gray')
plt.colorbar()
plt.show()

matrix2 = np.dot(IWnsq_H, Wnsq)
plt.imshow(matrix2.real, cmap='gray')
plt.colorbar()
plt.show()

# assert np.allclose(Wifftnsq.real[1:], y[1:])
assert np.allclose(Wifftnsq_H.real[1:], ifft.real[1:])

# assert np.allclose(Wifftnsq.real[1:], irfft.real[1:])

assert np.allclose(np.dot(Wifftnsq_H, Wnsq)[1:, 1:], np.identity(g.IFG_SIZE)[1:, 1:])