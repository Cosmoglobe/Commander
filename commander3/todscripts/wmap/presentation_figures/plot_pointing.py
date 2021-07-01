import matplotlib.pyplot as plt
import numpy as np

t = np.linspace(0.1, 8*np.pi-0.1, 100)
p1 = ((1+np.cos(t))*10).astype('int')
p2 = ((1+np.cos(t+np.pi/3))*10).astype('int')

t = np.linspace(0.1, 4*np.pi-0.1, 110)
p1 = ((np.cos(t)**2)*20).astype('int')
p2 = ((np.cos(t+np.pi/3)**2)*20).astype('int')
P = np.zeros((len(t), 20))



plt.figure(figsize=(3, 5))

for i in range(len(t)):
  P[i,p1[i]] += 1
plt.imshow(P, vmin=-1, vmax=1, cmap='coolwarm')
plt.title(r'$\mathsf{P}$')
plt.xlabel(r'$N_\mathrm{pix}$')
plt.ylabel(r'$N_\mathrm{TOD}$')
plt.xticks([])
plt.yticks([])
plt.savefig('sing_point.png', transparent=True, dpi=300)

plt.figure(figsize=(3,3))
plt.imshow(P.T.dot(P), cmap='coolwarm', vmin=-5, vmax=5)
plt.title(r'$\mathsf{P}^T\mathsf{N}^{-1}\mathsf{P}$')
plt.xlabel(r'$N_\mathrm{pix}$')
plt.ylabel(r'$N_\mathrm{pix}$')
plt.xticks([])
plt.yticks([])
plt.savefig('sing_mat.png', transparent=True, dpi=300)

plt.figure(figsize=(3, 5))
for i in range(len(t)):
  P[i,p2[i]] -= 1
plt.imshow(P, vmin=-1, vmax=1, cmap='coolwarm')
plt.title(r'$\mathsf{P}$')
plt.xlabel(r'$N_\mathrm{pix}$')
plt.ylabel(r'$N_\mathrm{TOD}$')
plt.xticks([])
plt.yticks([])
plt.savefig('diff_point.png', transparent=True, dpi=300)

plt.figure(figsize=(3,3))
plt.imshow(P.T.dot(P), cmap='coolwarm', vmin=-5, vmax=5)
plt.title(r'$\mathsf{P}^T\mathsf{N}^{-1}\mathsf{P}$')
plt.xlabel(r'$N_\mathrm{pix}$')
plt.ylabel(r'$N_\mathrm{pix}$')
plt.xticks([])
plt.yticks([])
plt.savefig('diff_mat.png', transparent=True, dpi=300)
plt.show()
