from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import itertools
import re
import copy

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
        
        # number of overlap segments -- 32 pixels wide each -- is 1 less than number of cubes in each dimension. Subtract 16 = half overlap region
        self.naxis1 = 512
        self.naxis2 = 512
        self.max_len_RA = (self.naxis1 * self.n_cubes_ra) - ((self.n_cubes_ra - 1) * 16)
        self.max_len_DEC = (self.naxis2 * self.n_cubes_dec) - ((self.n_cubes_dec - 1) * 16)
        print(self.max_len_RA)
    
        # RA incr to left? yes- should round min up, max down to recover max area 
        # must step first through dec, then ra -- can't do them at the same time because will double-count
        for r in self.my_center_RAs:
            ra = "{:06.2f}".format(r)
            ppv_cube = galfa_cuber.Cube(RA=ra, DEC="{:05.2f}".format(self.my_center_DECs[0]))
            
            ppv_cube_wcs = cutouts.make_wcs(ppv_cube.ppv_cube_fn)
            print(ra, dec)
            
            xmin, ymin = cutouts.radec_to_xy(self.RA_min, self.DEC_min, ppv_cube_wcs)
            xmax, ymax = cutouts.radec_to_xy(self.RA_max, self.DEC_max, ppv_cube_wcs)
            
            xmin = min(np.ceil(xmin), self.naxis1)
            xmax = max(np.floor(xmax), 0)
            
            self.max_len_RA -= (self.naxis1 - (xmin - xmax))    
            
            print("xmin, xmax, xmin-xmax = ", xmin, xmax, xmin-xmax)
                
            print("xmin, ymin", xmin, ymin)
            print("xmax, ymax", xmax, ymax)
            print(self.max_len_RA)
            print(self.max_len_DEC)
            
        for d in self.my_center_DECs:
            dec = "{:05.2f}".format(d)
            ppv_cube = galfa_cuber.Cube(RA="{:06.2f}".format(self.my_center_RAs[0]), DEC=dec)
            
            ppv_cube_wcs = cutouts.make_wcs(ppv_cube.ppv_cube_fn)
            print(ra, dec)
            
            xmin, ymin = cutouts.radec_to_xy(self.RA_min, self.DEC_min, ppv_cube_wcs)
            xmax, ymax = cutouts.radec_to_xy(self.RA_max, self.DEC_max, ppv_cube_wcs)
            
            ymin = max(np.floor(ymin), 0)
            ymax = min(np.ceil(ymax), self.naxis2)
            
            self.max_len_DEC -= (self.naxis2 - (ymax - ymin))
            
            print("xmin, xmax, xmin-xmax = ", xmin, xmax, xmin-xmax)
                
            print("xmin, ymin", xmin, ymin)
            print("xmax, ymax", xmax, ymax)
            print(self.max_len_RA)
            print(self.max_len_DEC)
            
        self.max_len_RA = np.int(self.max_len_RA)
        self.max_len_DEC = np.int(self.max_len_DEC)

    def make_RHT_XYT_cube(self, rht_velstart="0974", rht_velstop="0978"):
        """
        create new Rtheta cube
        """
        
        self.nthets = 165
        print("DEC, RA", self.max_len_DEC, self.max_len_RA)
        self.RHT_XYT_cube = np.zeros((self.nthets, self.max_len_DEC, self.max_len_RA), np.float_)
        
        for ra, dec in self.all_RADEC_str_pairs:
            # load small xyt cube
            ppv_cube = galfa_cuber.Cube(RA=ra, DEC=dec)
            ppv_cube.load_RHT_XYT_cube(rht_velstart=rht_velstart, rht_velstop=rht_velstop)
            rht_xyt_smallcube = ppv_cube.get_RHT_XYT_cube(ashdulist=False)
            
            xyt_cube_wcs = cutouts.make_wcs(ppv_cube.xyt_fn)
            print(ra, dec)
            
            # find start and end points in small cube
            xmin, ymin = cutouts.radec_to_xy(self.RA_min, self.DEC_min, xyt_cube_wcs)
            xmax, ymax = cutouts.radec_to_xy(self.RA_max, self.DEC_max, xyt_cube_wcs)
            
            xmin = min(np.ceil(xmin), self.naxis1)
            xmax = max(np.floor(xmax), 0)
            ymin = max(np.floor(ymin), 0)
            ymax = min(np.ceil(ymax), self.naxis2)
            
            # find start and end points in large new cube
            large_cube_wcs = copy.copy(xyt_cube_wcs)
            large_cube_wcs.wcs.naxis1 = self.max_len_RA
            large_cube_wcs.wcs.naxis2 = self.max_len_DEC
            new_xmin, new_ymin = cutouts.radec_to_xy(self.RA_min, self.DEC_min, large_cube_wcs)
            new_xmax, new_ymax = cutouts.radec_to_xy(self.RA_max, self.DEC_max, large_cube_wcs)
            new_xmin = np.ceil(new_xmin)
            new_xmax = np.floor(new_xmax)
            new_ymin = np.floor(new_ymin)
            new_ymax = np.ceil(new_ymax)
            
            print("new x y coords: ", new_xmin, new_ymin, new_xmax, new_ymax)
            self.RHT_XYT_cube[:, new_ymin:new_ymax, new_xmax:new_xmin] = rht_xyt_smallcube[:, ymin:ymax, xmax:xmin]

cc = NewCube(RA_min=50., RA_max=154, DEC_min=11, DEC_max=14)
cc.make_RHT_XYT_cube()
print(cc.all_RADEC_strs)

gg=galfa_cuber.Cube(RA="004.00", DEC="02.35")
gg_cube_wcs = cutouts.make_wcs(gg.ppv_cube_fn)
ra, dec = cutouts.xy_to_radec(512, 512, gg_cube_wcs)
print(ra, dec)

gg2=galfa_cuber.Cube(RA="012.00", DEC="10.35")
gg2_cube_wcs = cutouts.make_wcs(gg2.ppv_cube_fn)
ra2, dec2 = cutouts.xy_to_radec(0, 0, gg2_cube_wcs)
print(ra2, dec2)

xx, yy = cutouts.radec_to_xy(ra, dec, gg2_cube_wcs)
print(xx, yy)
xx, yy = cutouts.radec_to_xy(ra2, dec2, gg_cube_wcs)
print(xx, yy)
