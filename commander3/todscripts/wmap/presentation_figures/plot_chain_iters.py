import matplotlib.pyplot as plt
import numpy as np

import h5py

data = h5py.File('/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_WMAP_all/chain_c0001.h5', 'r')


bands=['023-WMAP_K', 
       '060-WMAP_V2']
label_list = [
         ['K113', 'K114', 'K123', 'K124'],
         ['V213', 'V214', 'V223', 'V224']]

#bands = [f'090-WMAP_W{i}' for i in range(1,5)]
#
#data = h5py.File('chains_WMAP_W_exp/chain_c0001.h5', 'r')
#label_list = [ [f'W{i}13', f'W{i}14', f'W{i}23', f'W{i}24'] for i in range(1,5)]
for band, labels in zip(bands, label_list):
    print('\n')
    print(band, labels)
    for j in range(4):
        print('\n')
        maxs =-np.ones(5)*np.inf
        mins = np.ones(5)*np.inf

        inds = (data[str(0).zfill(6)+f'/tod/{band}/accept'][j] == 1)
        g = data[str(0).zfill(6)+f'/tod/{band}/gain'][j][inds]
        mins[0] = min([g[0], g[-1]])
        maxs[0] = max([g[0], g[-1]])
        for i in range(30,len(data.keys())-1):
            inds = (data[str(i).zfill(6)+f'/tod/{band}/accept'][j] == 1) & \
            (abs(data[str(i).zfill(6)+f'/tod/{band}/chisq'][j]) < 100) & \
            (data[str(i).zfill(6)+f'/tod/{band}/gain'][j] > 0)
            inds = np.array(inds)
            g = data[str(i).zfill(6)+f'/tod/{band}/gain'][j][inds]
            if min(g) < mins[0]:
              mins[0] = min(g)
            if max(g) > maxs[0]:
              maxs[0] = max(g)

            if min(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][1][j][inds]) < mins[1]:
              mins[1] = min(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][1][j][inds])
            if max(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][1][j][inds]) > maxs[1]:
              maxs[1] = max(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][1][j][inds])
            if min(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][2][j][inds]) < mins[2]:
              mins[2] = min(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][2][j][inds])
            if max(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][2][j][inds]) > maxs[2]:
              maxs[2] = max(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][2][j][inds])
            if min(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][0][j][inds]) < mins[3]:
              mins[3] = min(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][0][j][inds])
            if max(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][0][j][inds]) > maxs[3]:
              maxs[3] = max(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][0][j][inds])
            if min(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][0][j][inds]/g) < mins[4]:
              mins[4] = min(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][0][j][inds]/g)
            #if max(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][0][j][inds]/g) > maxs[4]:
            #  maxs[4] = max(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][0][j][inds]/g)

            #b = data[str(i).zfill(6)+f'/tod/{band}/baseline'][j][inds]
            #if min((b-b.mean())/g) < mins[5]:
            #  mins[5] = min((b-b.mean())/g)
            #if max((b-b.mean())/g) > maxs[5]:
            #  maxs[5] = max((b-b.mean())/g)
            if min(data[str(i).zfill(6)+f'/tod/{band}/chisq'][j][inds]) < mins[4]:
              mins[4] = min(data[str(i).zfill(6)+f'/tod/{band}/chisq'][j][inds]) 
            if max(data[str(i).zfill(6)+f'/tod/{band}/chisq'][j][inds]) > maxs[4]:
              maxs[4] = max(data[str(i).zfill(6)+f'/tod/{band}/chisq'][j][inds]) 

        for i in range(5,len(data.keys())-1):
            fig, axes = plt.subplots(figsize=(8, 10), sharex=True, nrows=5)
            c = 'k'
            inds = (data[str(i).zfill(6)+f'/tod/{band}/accept'][j] == 1) & \
            (abs(data[str(i).zfill(6)+f'/tod/{band}/chisq'][j]) < 100) & \
            (data[str(i).zfill(6)+f'/tod/{band}/gain'][j] > 0)
            inds = np.array(inds)
            t = np.arange(len(inds))

            g = data[str(1).zfill(6)+f'/tod/{band}/gain'][j][inds]
            axes[0].plot(t[inds], g, '.', color='r', ms=1)

            g = data[str(i).zfill(6)+f'/tod/{band}/gain'][j][inds]
            b = data[str(i).zfill(6)+f'/tod/{band}/baseline'][j][inds]
            axes[0].plot(t[inds], g, '.', color=c, ms=1)
            axes[1].plot(t[inds], data[str(i).zfill(6)+f'/tod/{band}/xi_n'][1][j][inds], '.', color=c, ms=1)
            axes[2].plot(t[inds], data[str(i).zfill(6)+f'/tod/{band}/xi_n'][2][j][inds], '.', color=c, ms=1)
            axes[3].plot(t[inds], data[str(i).zfill(6)+f'/tod/{band}/xi_n'][0][j][inds], '.', color=c, ms=1)
            #axes[4].plot(t[inds], data[str(i).zfill(6)+f'/tod/{band}/xi_n'][0][j][inds]/data[str(i).zfill(6)+f'/tod/{band}/gain'][j][inds],
            #    '.', color=c, ms=1)
            #axes[5].plot(t[inds], (b-b.mean())/g, '.', color=c, ms=1)
            axes[4].plot(t[inds], data[str(i).zfill(6)+f'/tod/{band}/chisq'][j][inds], color=c, ms=1)
    
            axes[4].axhline(0, color='r', linestyle=':', lw=0.5)
    
            print(np.mean(data[str(i).zfill(6)+f'/tod/{band}/chisq'][j][inds]))

            for num in range(5):
              axes[num].set_ylim(mins[num], maxs[num])
    
            #axes[1].set_ylim([0.0002, 0.002])
            #axes[2].set_ylim([-0.9, -0.5])
            #axes[6].set_ylim([-20, 0])
            #if j < 2:
            #  axes[0].set_ylim([0.85, 1.35])
            #  axes[3].set_ylim([2.5, 3.5])
            #  axes[4].set_ylim([2.65, 2.85])
            #  #axes[1].hlines(0.40e-3, xmin=0, xmax=468)
            #  #axes[1].hlines(6.13e-3, xmin=0, xmax=468)
            #else:
            #  axes[0].set_ylim([1.10, 1.35])
            #  axes[3].set_ylim([3.2, 4.4])
            #  axes[4].set_ylim([2.9, 3.2])
            #  #axes[1].hlines(0.51e-3, xmin=0, xmax=468)
            #  #axes[1].hlines(5.37e-3, xmin=0, xmax=468)
            #if j == 2:
            #  axes[0].set_ylim([1.1, 1.4])
            #  #axes[3].set_ylim([3.5, 4.5])
            #  #axes[4].set_ylim([3.0, 3.3])
            #  axes[5].set_ylim([-5, 5])
            #if j == 3:
            #  axes[0].set_ylim([1.05, 1.35])
            #  axes[3].set_ylim([3.25, 4.25])
            #  axes[4].set_ylim([2.9, 3.2])
            #  axes[5].set_ylim([-5, 5])
            axes[0].set_ylabel(r'$g$ [du/mK]')
            axes[1].set_ylabel(r'$f_\mathrm{k}$ [Hz]')
            axes[1].set_yscale('log')
            axes[2].set_ylabel(r'$\alpha$')
            axes[3].set_ylabel(r'$\sigma_0$ [du]')
            #axes[4].set_ylabel(r'$\sigma_0$ [mK]')
            #axes[5].set_ylabel(r'$\Delta b/g$ [mK]')
            axes[4].set_ylabel(r'$(\chi^2-n_\mathrm{tod})/\sqrt{2n_\mathrm{tod}}$')
            axes[4].set_xlabel('Scan number')
    
    
            plt.suptitle(i, size=16)
            if i > 0:
              plt.savefig(f'plots/{labels[j]}_iter_{str(i).zfill(3)}.png',
                  bbox_inches='tight')
              plt.close()
    
