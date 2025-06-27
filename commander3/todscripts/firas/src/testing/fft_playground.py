import matplotlib.pyplot as plt
import numpy as np

# generate white noise
white_noise =  np.random.normal(0, 10, 512)

plt.plot(white_noise)
plt.title("White Noise")
plt.xlabel("Sample")
plt.ylabel("Amplitude")
plt.grid()
plt.show()

white_noise_fft = np.fft.fft(white_noise)
plt.plot(np.abs(white_noise_fft))
plt.title("FFT of White Noise")
plt.xlabel("Frequency Bin")
plt.ylabel("Magnitude")
plt.grid()
plt.show()

white_noise_fft = np.random.normal(30, 15, 257)
plt.plot(white_noise_fft)
plt.title("White Noise FFT")
plt.xlabel("Frequency Bin")
plt.ylabel("Magnitude")
plt.grid()
plt.show() 

white_noise_ifg = np.fft.ifft(white_noise_fft)
plt.plot(np.abs(white_noise_ifg))
plt.title("Inverse FFT of White Noise FFT")
plt.xlabel("Sample")
plt.ylabel("Amplitude")
plt.grid()
plt.show()