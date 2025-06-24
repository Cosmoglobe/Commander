import os
import sys

import matplotlib.pyplot as plt
import numpy as np

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
import my_utils as mu

frequency = mu.generate_frequencies("ll", "ss", 257)
spectrum = mu.planck(frequency, np.array([2.76]))[0]
spectrum43 = spectrum.copy()
cutoff = 5
spectrum43[:cutoff] = 0
spectrum43[cutoff+43:] = 0
spectrumcut = spectrum43[spectrum43>0]

plt.plot(frequency, spectrum, label="2.76 K")
plt.plot(frequency, spectrum43, label="2.76 K with 43 removed")
plt.plot(frequency[cutoff:cutoff+43], spectrumcut, label="2.76 K with 43 removed and cut off")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Brightness (MJy/sr)")
plt.title("Brightness Spectrum at 2.76 K")
plt.legend()
plt.grid()
plt.show()

# get irfft of each to see the difference
spectrum_ifft = np.fft.irfft(spectrum)
spectrum43_ifft = np.fft.irfft(spectrum43)
spectrumcut_ifft = np.fft.irfft(spectrumcut)
plt.plot(spectrum_ifft, label="2.76 K")
plt.plot(spectrum43_ifft, label="2.76 K with 43 removed")
plt.plot(spectrumcut_ifft, label="2.76 K with 43 removed and cut off")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Brightness (MJy/sr)")
plt.title("Brightness Spectrum at 2.76 K")
plt.legend()
plt.grid()
plt.show()
