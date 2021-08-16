import matplotlib.pyplot as plt
import numpy as np

import h5py

data = h5py.File('/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_WMAP_all/chain_c0001.h5', 'r')

data = h5py.File('/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_WMAP_beamtest/chain_c0001.h5', 'r')


gain = {}
gain['K113'] = 0.97
gain['K114'] = 0.9938
gain['K123'] = 1.1745
gain['K124'] = 1.12

gain['Ka113'] = 0.8668
gain['Ka114'] = 0.8753
gain['Ka123'] = 1.0914
gain['Ka124'] = 1.0033

gain['Q113'] = 1.053
gain['Q114'] = 0.9834
gain['Q123'] = 0.4914
gain['Q124'] = 0.5365

gain['Q213'] = 0.9882
gain['Q214'] = 1.0173
gain['Q223'] = 0.8135
gain['Q224'] = 0.7896

gain['V113'] = 0.4896
gain['V114'] = 0.5380
gain['V123'] = 0.5840
gain['V124'] = 0.5840

gain['V213'] = 0.4948
gain['V214'] = 0.4872
gain['V223'] = 0.4096
gain['V224'] = 0.3802

gain['W113'] = 0.3888
gain['W114'] = 0.4139
gain['W123'] = 0.3290
gain['W124'] = 0.3003

gain['W213'] = 0.3587
gain['W214'] = 0.3701
gain['W223'] = 0.3655
gain['W224'] = 0.3666

gain['W313'] = 0.3255
gain['W314'] = 0.3517
gain['W323'] = 0.3291
gain['W324'] = 0.3225

gain['W413'] = 0.2841
gain['W414'] = 0.2918
gain['W423'] = 0.3796
gain['W424'] = 0.3591

bands=['023-WMAP_K', 
       '060-WMAP_V2']
label_list = [
         ['K113', 'K114', 'K123', 'K124'],
         ['V213', 'V214', 'V223', 'V224']]

#data = h5py.File('/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_WMAP_opt/chain_c0001.h5', 'r')
#bands=['023-WMAP_K', 
#       '030-WMAP_Ka',
#       '040-WMAP_Q1',
#       '040-WMAP_Q2',
#       '060-WMAP_V1',
#       '060-WMAP_V2',
#       '090-WMAP_W1',
#       '090-WMAP_W2',
#       '090-WMAP_W3',
#       '090-WMAP_W4']
#label_list = [
#         ['K113', 'K114', 'K123', 'K124'],
#         ['Ka113', 'Ka114', 'Ka123', 'Ka124'],
#         ['Q113', 'Q114', 'Q123', 'Q124'],
#         ['Q213', 'Q214', 'Q223', 'Q224'],
#         ['V113', 'V114', 'V123', 'V124'],
#         ['V213', 'V214', 'V223', 'V224'],
#         ['W113', 'W114', 'W123', 'W124'],
#         ['W213', 'W214', 'W223', 'W224'],
#         ['W313', 'W314', 'W323', 'W324'],
#         ['W413', 'W414', 'W423', 'W424']]
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
        #for i in range(len(data.keys())-2, len(data.keys())-10, -1):
        #for i in range(len(data.keys())-20, len(data.keys())-1):
        for i in range(1, len(data.keys())-1):
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

            if min(data[str(i).zfill(6)+f'/tod/{band}/chisq'][j][inds]) < mins[4]:
              mins[4] = min(data[str(i).zfill(6)+f'/tod/{band}/chisq'][j][inds]) 
            if max(data[str(i).zfill(6)+f'/tod/{band}/chisq'][j][inds]) > maxs[4]:
              maxs[4] = max(data[str(i).zfill(6)+f'/tod/{band}/chisq'][j][inds]) 

        for i in range(1,len(data.keys())-1):
            fig, axes = plt.subplots(figsize=(8, 10), sharex=True, nrows=5)
            c = 'k'
            inds = (data[str(i).zfill(6)+f'/tod/{band}/accept'][j] == 1) & \
            (abs(data[str(i).zfill(6)+f'/tod/{band}/chisq'][j]) < 100) & \
            (data[str(i).zfill(6)+f'/tod/{band}/gain'][j] > 0)
            inds = np.array(inds)
            t = np.arange(len(inds))

            g = data[str(1).zfill(6)+f'/tod/{band}/gain'][j][inds]
            axes[0].plot(t[inds], g, '.', color='r', ms=1)
            #axes[0].axhline(gain[labels[j]], color='r')

            g = data[str(i).zfill(6)+f'/tod/{band}/gain'][j][inds]
            b = data[str(i).zfill(6)+f'/tod/{band}/baseline'][j][inds]
            axes[0].plot(t[inds], g, '.', color=c, ms=1)
            axes[1].plot(t[inds], data[str(i).zfill(6)+f'/tod/{band}/xi_n'][1][j][inds], '.', color=c, ms=1)
            axes[2].plot(t[inds], data[str(i).zfill(6)+f'/tod/{band}/xi_n'][2][j][inds], '.', color=c, ms=1)
            axes[3].plot(t[inds], data[str(i).zfill(6)+f'/tod/{band}/xi_n'][0][j][inds], '.', color=c, ms=1)
            axes[4].plot(t[inds], data[str(i).zfill(6)+f'/tod/{band}/chisq'][j][inds], color=c, ms=1)
    
            axes[4].axhline(0, color='r', linestyle=':', lw=0.5)
    
            #print(np.mean(data[str(i).zfill(6)+f'/tod/{band}/chisq'][j][inds]))
            #print(np.mean(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][2][j][inds]))
            print(np.mean(data[str(i).zfill(6)+f'/tod/{band}/xi_n'][1][j][inds]))

            for num in range(5):
              axes[num].set_ylim(mins[num], maxs[num])
    
            axes[0].set_ylabel(r'$g$ [du/mK]')
            axes[1].set_ylabel(r'$f_\mathrm{k}$ [Hz]')
            axes[1].set_yscale('log')
            axes[2].set_ylabel(r'$\alpha$')
            axes[3].set_ylabel(r'$\sigma_0$ [du]')
            #axes[4].set_ylabel(r'$\sigma_0$ [mK]')
            #axes[5].set_ylabel(r'$\Delta b/g$ [mK]')
            axes[4].set_ylabel(r'$(\chi^2-n_\mathrm{tod})/\sqrt{2n_\mathrm{tod}}$')
            axes[4].set_xlabel('Scan number')
    
    
            plt.suptitle(i, size=16, ha='right')
            if i > 0:
              plt.savefig(f'plots/{labels[j]}_iter_{str(i).zfill(3)}.png',
                  bbox_inches='tight')
            plt.close()
    
