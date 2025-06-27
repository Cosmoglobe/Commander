# compare outputs from simulations and residuals
import matplotlib.pyplot as plt
import numpy as np

sim = np.load("/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/src/simulations/ifgsim.npy")
data = np.load("/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/src/calibration/ifgdata.npy")

print(data)

plt.plot(sim[0], label="Simulation IFG")
plt.plot(data, label="Data IFG")
plt.title("Comparison of Simulation and Data IFG")
plt.xlabel("Sample")
plt.ylabel("Amplitude")
plt.legend()
plt.grid()
plt.show()