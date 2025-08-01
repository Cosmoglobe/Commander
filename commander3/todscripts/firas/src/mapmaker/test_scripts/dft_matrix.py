import matplotlib.pyplot as plt
import numpy as np

import globals as g

N = g.IFG_SIZE

W = np.zeros((N, N), dtype=complex)
W[0, :] = 1
W[:, 0] = 1
omega = np.exp(2j * np.pi / N)
for i in range(1, N):
    for j in range(1, N):
        W[i, j] = omega ** ((i * j))

IW = np.zeros((N, N), dtype=complex)
IW[0, :] = 1
IW[:, 0] = 1
for i in range(1, N):
    for j in range(1, N):
        IW[i, j] = omega ** (- ((i * j)))
IW = IW / N

print("Sine wave")

x = np.linspace(0, 10, N)
y = np.sin(x)

# fourier transform
yhat = np.dot(W, y)
yhat_fft = np.fft.fft(y)

# compare
print(f"Are the FFT and the DTF matrix performing the same operation? Real part: {np.allclose(yhat.real, yhat_fft.real)}, imaginary part: {np.allclose(yhat.imag, yhat_fft.imag)}")

fig, ax = plt.subplots(2, 1)
ax[0].plot(x, yhat.real, label='DFT Matrix')
ax[0].plot(x, yhat_fft.real, label='FFT')
ax[0].legend()
ax[0].set_title('Real Part')
ax[1].plot(x, yhat.imag, label='DFT Matrix')
ax[1].plot(x, yhat_fft.imag, label='FFT')
ax[1].legend()
ax[1].set_title('Imaginary Part')
# plt.plot(x, yhat, label = 'DFT Matrix')
# plt.plot(x, yhat_fft, label = 'FFT')
plt.legend()
plt.show()

# inverse fourier transform
y_new = np.dot(IW, yhat)
y_new_fft = np.fft.ifft(yhat_fft)

print(f"Are the inverse FFT and the inverse DTF matrix performing the same operation? {np.allclose(y_new.real, y_new_fft.real)}")
print(f"Is the original signal recovered? {np.allclose(y, y_new.real)}")

plt.plot(x, y_new, label = 'Inverse DFT Matrix')
plt.plot(x, y_new_fft, label = 'Inverse FFT')
plt.plot(x, y, label = 'Original Signal')
plt.legend()
plt.show()