from __future__ import division, print_function
import numpy as np
from astropy import wcs
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
            
        self.RA_min_orig = np.float(RA_min)
        self.RA_max_orig = np.float(RA_max)
        self.DEC_min_orig = np.float(DEC_min)
        self.DEC_max_orig = np.float(DEC_max)
        
        # create a flag for wrapped RA. Deal with this later.
        if self.RA_max_orig < self.RA_min_orig:
            self.wrapRA = True
        else:
            self.wrapRA = False
        
        # Define dimensions of new cube. 
        
        # get allsky RHT data header
        self.path_to_rht_thetaslices = "/disks/jansky/a/users/goldston/susan/Wide_maps/single_theta_maps/"
        self.galfa_allsky_hdr = fits.getheader(self.path_to_rht_thetaslices+"S0974_0978/intrht_S0974_0978.fits")
        
        #startpix_ras = [(10800.5 - pix)*0.0166667 + 180.0 for pix in np.arange(0, 21600)]
                
        #startpix_ras = np.arange(1/180., 360, 0.0166667)[::-1] #np.linspace(0, 360, 21599)
        
        self.allsky_crval1 = self.galfa_allsky_hdr['CRVAL1']
        self.allsky_crval2 = self.galfa_allsky_hdr['CRVAL2']
        self.allsky_cdelt1 = self.galfa_allsky_hdr['CDELT1']
        self.allsky_cdelt2 = self.galfa_allsky_hdr['CDELT2']
        self.allsky_naxis1 = self.galfa_allsky_hdr['NAXIS1']
        self.allsky_naxis2 = self.galfa_allsky_hdr['NAXIS2']
        self.allsky_crpix1 = self.galfa_allsky_hdr['CRPIX1']
        self.allsky_crpix2 = self.galfa_allsky_hdr['CRPIX2']
        
        startpix_ras =  self.allsky_crval1 +  self.allsky_cdelt1 * (np.arange(self.allsky_naxis1) + 1 -  self.allsky_crpix1)
        startpix_decs =  self.allsky_crval2 +  self.allsky_cdelt2 * (np.arange(self.allsky_naxis2) + 1 -  self.allsky_crpix2)
        print("first startpixra", startpix_ras[0])
        print("first startpixdec", startpix_decs[0])
        
        # redefine RA, DEC min and max s.t. they are integer pixels.
        self.RA_min = startpix_ras[np.argmin(np.abs(startpix_ras - self.RA_min_orig))]
        self.RA_max = startpix_ras[np.argmin(np.abs(startpix_ras - self.RA_max_orig))]
        print("RA min orig = {}, rounded = {}".format(self.RA_min_orig, self.RA_min))
        self.DEC_min = startpix_decs[np.argmin(np.abs(startpix_decs - self.DEC_min_orig))]
        self.DEC_max = startpix_decs[np.argmin(np.abs(startpix_decs - self.DEC_max_orig))]
        print("DEC min orig = {}, rounded = {}".format(self.DEC_min_orig, self.DEC_min))
        
        
        print("DEC max - min = {} in pixels = {}".format(self.DEC_max - self.DEC_min, (self.DEC_max - self.DEC_min)/self.allsky_cdelt2))
        
        # translate cube corners to allsky x, y
        allsky_w = cutouts.make_wcs(self.galfa_allsky_hdr)
        allsky_x1, allsky_y1 = cutouts.radec_to_xy(self.RA_max, self.DEC_min, allsky_w) # LLH corner
        allsky_x2, allsky_y2 = cutouts.radec_to_xy(self.RA_min, self.DEC_max, allsky_w) # URH corner

        # ROUND to integer x, y. Necessary because allsky header does not quite match up. Valid per Yong's check.
        self.allsky_xstart = np.int(np.floor(allsky_x1))
        self.allsky_xstop = np.int(np.ceil(allsky_x2))
        self.allsky_ystart = np.int(np.floor(allsky_y1))
        self.allsky_ystop = np.int(np.ceil(allsky_y2))
        
        self.newcube_xlen = (self.allsky_xstop - self.allsky_xstart) +1
        self.newcube_ylen = (self.allsky_ystop - self.allsky_ystart) + 1
        print("new cube xlen, ylen = {}, {}".format(self.newcube_xlen, self.newcube_ylen))
        self.newcube_centerRA, self.newcube_centerDEC = cutouts.xy_to_radec(self.allsky_xstart + self.newcube_xlen/2.0, self.allsky_ystart + self.newcube_ylen/2.0 + 1, allsky_w)
        
        # define new 2D wcs object
        self.new_cube_flat_wcs = wcs.WCS(naxis=2)
        self.new_cube_flat_wcs.wcs.crpix = [self.newcube_xlen/2.0 + 1, self.newcube_ylen/2.0 + 1]
        self.new_cube_flat_wcs.wcs.cdelt = np.array([-0.0166667, 0.0166667]) # assume format is cdelt1, cdelt2 (ra, dec)
        self.new_cube_flat_wcs.wcs.crval = [self.newcube_centerRA, self.newcube_centerDEC]
        self.new_cube_flat_wcs.wcs.ctype = ["RA      ", "DEC     "]
        
        print("old ra max = {}, ra min = {}, dec max = {}, dec min = {}".format(self.RA_max, self.RA_min, self.DEC_max, self.DEC_min))
        self.zRA_max, self.zDEC_min = cutouts.xy_to_radec(0, 0, self.new_cube_flat_wcs)
        self.zRA_min, self.zDEC_max = cutouts.xy_to_radec(self.newcube_xlen, self.newcube_ylen, self.new_cube_flat_wcs)
        print("new ra max = {}, ra min = {}, dec max = {}, dec min = {}".format(self.RA_max, self.RA_min, self.DEC_max, self.DEC_min))
        
        testx, testy = cutouts.radec_to_xy(52.0, 2.35, self.new_cube_flat_wcs)
        print("COORD TEST: 52, 2.35 are at x={}, y={}".format(testx, testy))
        
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
        print("n cubes des = {}, n cubes ra = {}".format(self.n_cubes_dec, self.n_cubes_ra))
        #print(self.my_center_DECs)
        #print(self.my_center_RAs)
        
        # number of overlap segments -- 32 pixels wide each -- is 1 less than number of cubes in each dimension. Subtract 16 = half overlap region
        self.naxis1 = 512
        self.naxis2 = 512
        #self.max_len_RA = (self.naxis1 * self.n_cubes_ra) - ((self.n_cubes_ra - 1) * 16)
        #self.max_len_DEC = (self.naxis2 * self.n_cubes_dec) - ((self.n_cubes_dec - 1) * 16)
        #print(self.max_len_RA)
    

    def make_RHT_XYT_cube(self, rht_velstart="0974", rht_velstop="0978"):
        """
        create new Rtheta cube
        """
        
        self.nthets = 165
        self.RHT_XYT_cube = np.zeros((self.nthets, self.newcube_ylen, self.newcube_xlen), np.float_)
        
        # strep through all necessary constituent cubes
        for ra, dec in self.all_RADEC_str_pairs:
            # load small xyt cube
            ppv_cube = galfa_cuber.Cube(RA=ra, DEC=dec)
            ppv_cube.load_RHT_XYT_cube(rht_velstart=rht_velstart, rht_velstop=rht_velstop)
            rht_xyt_smallcube = ppv_cube.get_RHT_XYT_cube(ashdulist=False)
            
            xyt_cube_wcs = cutouts.make_wcs(ppv_cube.xyt_fn)
            print(ra, dec)
            
            testx, testy = cutouts.radec_to_xy(52.0, 2.35, xyt_cube_wcs)
            print("COORD TEST: 52, 2.35 are at x={}, y={}".format(testx, testy))
            
            # find start and end pixels in small cube
            xmin, ymin = cutouts.radec_to_xy(self.RA_min, self.DEC_min, xyt_cube_wcs)
            xmax, ymax = cutouts.radec_to_xy(self.RA_max, self.DEC_max, xyt_cube_wcs)
            
            xmin = np.int(min(np.ceil(xmin), self.naxis1))
            xmax = np.int(max(np.floor(xmax), 0))
            ymin = np.int(max(np.floor(ymin), 0))
            ymax = np.int(min(np.ceil(ymax), self.naxis2))
            
            # find start and end points in large new cube
            #large_cube_wcs = copy.copy(xyt_cube_wcs)
            #large_cube_wcs.wcs.naxis1 = self.max_len_RA
            #large_cube_wcs.wcs.naxis2 = self.max_len_DEC
            
            # find start and end pixels in new cube
            ramin, decmin = cutouts.xy_to_radec(xmin, ymin, xyt_cube_wcs)
            ramax, decmax = cutouts.xy_to_radec(xmax, ymax, xyt_cube_wcs)
            new_xmin, new_ymin = cutouts.radec_to_xy(ramin, decmin, self.new_cube_flat_wcs)
            new_xmax, new_ymax = cutouts.radec_to_xy(ramax, decmax, self.new_cube_flat_wcs)
            new_xmin = np.int(min(np.ceil(new_xmin), self.newcube_xlen))
            new_xmax = np.int(max(np.floor(new_xmax), 0))
            new_ymin = np.int(max(np.floor(new_ymin), 0))
            new_ymax = np.int(min(np.ceil(new_ymax), self.newcube_ylen))
            
            print("new x y coords: ", new_xmin, new_ymin, new_xmax, new_ymax)
            print("insert x y coords: ", xmin, ymin, xmax, ymax)
            print("new x len:", new_xmin-new_xmax + 1)
            print("insert x len:", xmin-xmax + 1)
            print("new y len:", new_ymax-new_ymin + 1)
            print("insert y len:", ymax-ymin + 1)
            print(rht_xyt_smallcube.shape)
            self.RHT_XYT_cube[:, new_ymin:new_ymax+1, new_xmax:new_xmin+1] = rht_xyt_smallcube[:, ymin:ymax+1, xmax:xmin+1]

cc = NewCube(RA_min=50., RA_max=55., DEC_min=2.35-1./60, DEC_max=2.35+1./60)
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
