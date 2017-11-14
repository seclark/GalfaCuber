from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import itertools
import re

import sys 
sys.path.insert(0, '../../FITSHandling/code')
import cutouts

import galfa_vel_helpers
import galfa_cuber

class NewCube():
    """
    New data cube from RA, DEC min, max inputs
    specify coordinates
    methods fill with data
    """
    def __init__(self, RA_min=None, RA_max=None, DEC_min=None, DEC_max=None, 
                       verbose=False):
        
        if any(coord is None for coord in [RA_min, RA_max, DEC_min, DEC_max]):
            print("Please specify all coordinates.")
            
        self.RA_min = np.float(RA_min)
        self.RA_max = np.float(RA_max)
        self.DEC_min = np.float(DEC_min)
        self.DEC_max = np.float(DEC_max)
        
        # create a flag for 180-wrapped RA. Deal with this later.
        if self.RA_max < self.RA_min:
            self.wrapRA = True
        else:
            self.wrapRA = False
        
        # Define dimensions of new cube. 
        
        # All the existing data cubes
        all_center_DECs = [2.35, 10.35, 18.35, 26.35, 34.35]
        all_center_RAs = [ra for ra in np.arange(4, 360, 8)]
        cube_halflength = 4.275
        all_min_DECs = [dec - cube_halflength for dec in all_center_DECs]
        all_max_DECs = [dec + cube_halflength for dec in all_center_DECs]
        all_min_RAs = np.asarray([ra - cube_halflength for ra in all_center_RAs])
        all_max_RAs = np.asarray([ra + cube_halflength for ra in all_center_RAs])

        # All the center RAs and DECs of the constituent cubes needed to make new cube
        my_center_RAs = [ra for i, ra in enumerate(all_center_RAs) if all_min_RAs[i] > self.RA_min and all_max_RAs[i] < self.RA_max]
        my_center_DECs = [dec for i, dec in enumerate(all_center_DECs) if all_min_DECs[i] > self.DEC_min and all_max_DECs[i] < self.DEC_max]

        [str(a)+"."+str(b) for (a, b) in itertools.product([5, 6, 3], [3,4])]
        self.all_RADEC_strs = ["RA+DEC_{:06.2f}+{:05.2f}".format(ra, dec) for (ra, dec) in itertools.product(my_center_RAs, my_center_DECs)]
