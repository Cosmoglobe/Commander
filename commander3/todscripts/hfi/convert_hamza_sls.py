import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

beams_dir = '/mn/stornext/d23/cmbco/globe/orig/planck/aux/beams'

freq = '857'
dets = ['1', '2', '3', '4']
versions = ['1', '2', '3', '4', '5']
nsides = {'1':'16', '2':'32', '3':'32', '4':'32', '5':'64'}

for det in dets:

    map_data = np.zeros(hp.nside2npix(64))

    for version in versions:

        sub_map = hp.read_map(beams_dir + '/FSLbeam_final_v' + version + '_' + freq + '-' + det + '_nside=' + nsides[version] + '.fits')

        if(nsides[version] != '64'):
            sub_map = hp.ud_grade(sub_map, 64)

        map_data += sub_map


    plt.figure()
    hp.projview(map_data)
    plt.savefig(beams_dir + '/FSL_' + freq + '-' + det + '_map_hamza_summed.png')

    alms = hp.map2alm(map_data, pol=False)

    hp.write_alm(beams_dir + '/FSL_' + freq + '-' + det + '_alms_hamza_summed.fits', alms, overwrite=True)


