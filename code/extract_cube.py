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
        print(self.DEC_min, self.DEC_max)
        
        # create a flag for wrapped RA. Deal with this later.
        if self.RA_max < self.RA_min:
            self.wrapRA = True
        else:
            self.wrapRA = False
        
        # Define dimensions of new cube. 
        
        # All the existing data cubes
        all_center_DECs = [2.35, 10.35, 18.35, 26.35, 34.35]
        all_center_RAs = [ra for ra in np.arange(4, 360, 8)]
        cube_halflength = 4.275 # in degrees
        all_min_DECs = [dec - cube_halflength for dec in all_center_DECs]
        all_max_DECs = [dec + cube_halflength for dec in all_center_DECs]
        all_min_RAs = np.asarray([ra - cube_halflength for ra in all_center_RAs])
        all_max_RAs = np.asarray([ra + cube_halflength for ra in all_center_RAs])

        # All the center RAs and DECs of the constituent cubes needed to make new cube 
        self.my_center_RAs = [ra for i, ra in enumerate(all_center_RAs) if all_max_RAs[i] >= self.RA_min and all_min_RAs[i] <= self.RA_max]
        self.my_center_DECs = [dec for i, dec in enumerate(all_center_DECs) if all_max_DECs[i] >= self.DEC_min and all_min_DECs[i] <= self.DEC_max]

        # RA+DEC strings for all constituent cubes
        [str(a)+"."+str(b) for (a, b) in itertools.product([5, 6, 3], [3,4])]
        self.all_RADEC_strs = ["RA+DEC_{:06.2f}+{:05.2f}".format(ra, dec) for (ra, dec) in itertools.product(self.my_center_RAs, self.my_center_DECs)]
        self.all_RADEC_str_pairs = [("{:06.2f}".format(ra), "{:05.2f}".format(dec)) for (ra, dec) in itertools.product(self.my_center_RAs, self.my_center_DECs)]
        
        print(self.all_RADEC_str_pairs)
        
        self.n_cubes_dec = len(self.my_center_DECs)
        self.n_cubes_ra = len(self.my_center_RAs)
        print(self.n_cubes_dec, self.n_cubes_ra)
        print(self.my_center_DECs)
        print(self.my_center_RAs)

        # corner RA DEC strings
        self.blc = "RA+DEC_{:06.2f}+{:05.2f}".format(np.nanmin(self.my_center_RAs), np.nanmin(self.my_center_DECs))
        self.ulc = "RA+DEC_{:06.2f}+{:05.2f}".format(np.nanmin(self.my_center_RAs), np.nanmax(self.my_center_DECs))
        self.brc = "RA+DEC_{:06.2f}+{:05.2f}".format(np.nanmax(self.my_center_RAs), np.nanmin(self.my_center_DECs))
        self.urc = "RA+DEC_{:06.2f}+{:05.2f}".format(np.nanmax(self.my_center_RAs), np.nanmax(self.my_center_DECs))
    
    
        #for radec_str in self.all_RADEC_strs:
        #    cube = galfa_cuber.Cube(RA=ra, DEC=dec
        
        # the lazy way -- make big cube, then refine
        #bigcube = 
        

cc = NewCube(RA_min=50., RA_max=54, DEC_min=11, DEC_max=13)
print(cc.all_RADEC_strs)


