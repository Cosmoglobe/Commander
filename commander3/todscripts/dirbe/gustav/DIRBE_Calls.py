from DIRBE_FITS2HDF5 import *
from DIRBE_Functions import *
from DIRBE_Plotting_and_Figures import *

if __name__=="__main__":

    # The DIRBE reformatting project function calls
    parse_directory()

    start = timer()
    binmap_all_bands()
    end = timer()
    print('Binmapping all bands and detectors took %.2f minutes' % ((end - start)/60))

    binmap_pretty_plot_all(save=True)

    diff_map_all()

    pretty_Y_plot_all(save=True)

    TT_plot_all()
