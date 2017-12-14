from __future__ import division
import numpy as np
from astropy.io import fits
import re

import sys 
sys.path.insert(0, '../../FITSHandling/code')
import cutouts

import galfa_vel_helpers
import galfa_cuber


# test
#RA = "156.00"
#DEC = "26.35"
#rht_velstr="S0974_0978"
#cube = Cube(RA=RA, DEC=DEC)
#cube.get_cube_coordinates_in_allsky()
#cube.make_RHT_XYT_cube(rht_velstr=rht_velstr)
#hdulist = cube.get_RHT_XYT_cube(ashdulist = True)
#hdulist.writeto("../testdata/testrht_velcube_"+rht_velstr+"RA+DEC_"+RA+"+"+DEC+".fits")

if __name__ == "__main__":
    all_DECs = ["02.35", "10.35", "18.35", "26.35", "34.35"]
    #all_RAs = ["{0:0=3d}.00".format(ra) for ra in np.arange(12, 350, 8)]
    #all_RAs = ["{0:0=3d}.00".format(ra) for ra in np.arange(180, 360, 8)]
    all_RAs = ["{0:0=3d}.00".format(ra) for ra in [4, 356]]

    for ra in all_RAs:
        for dec in all_DECs:
            galfa_cuber.make_single_cube_rtheta(RA=ra, DEC=dec, rht_velstart="1039", rht_velstop="1043", verbose=True)
            #galfa_cuber.make_single_cube_IQU(RA=ra, DEC=dec, verbose=True)

    #RA = "156.00"
    #DEC = "26.35"
    #make_single_cube_IQU(RA=RA, DEC=DEC, verbose=True)
    