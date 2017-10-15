from __future__ import division
import numpy as np
from astropy.io import fits

import sys 
sys.path.insert(0, '../../FITSHandling/code')
import cutouts

path_to_galfa_cubes = "/disks/jansky/a/users/goldston/DR2W_RC5/Wide/"
path_to_rht_thetaslices = "/disks/jansky/a/users/goldston/susan/Wide_maps/single_theta_maps/"

# eventually we will step through these cubes -- start w/ 1 test
galfa_cube_name = "GALFA_HI_RA+DEC_356.00+34.35_W"
galfa_cube_fn = path_to_galfa_cubes + galfa_cube_name + ".fits"

galfa_cube_hdr = fits.getheader(galfa_cube_fn)
galfa_allsky_hdr = fits.getheader(path_to_rht_thetaslices+"S0974_0978/intrht_S0974_0978.fits")

# R(theta) cube dimensions
nthets = 165
rht_data_cube = np.zeros((nthets, galfa_cube_hdr['NAXIS2'], galfa_cube_hdr['NAXIS1']), np.float_)


# Obtain x, y values in *big* GALFA data. These should be integers.
cube_crval1 = galfa_cube_hdr['CRVAL1']
cube_crval2 = galfa_cube_hdr['CRVAL2']
cube_crpix1 = galfa_cube_hdr['CRPIX1']

print(cube_crval1, cube_crval2)
cube_w = cutouts.make_wcs(galfa_cube_fn)
cube_edgex1, cube_edgey1 = cutouts.radec_to_xy(cube_crval1, cube_crval2, cube_w)
print(cube_edgex1, cube_edgey1)
cube_edgera1, cube_edgedec1 = cutouts.xy_to_radec(cube_crval1, cube_crval2, cube_w)
print(cube_edgera1, cube_edgedec1)

allsky_w = cutouts.make_wcs(galfa_allsky_hdr)
cutout_x1, cutout_y1 = cutouts.radec_to_xy(cube_edgera1, cube_edgedec1, allsky_w)
print(cutout_x1, cutout_y1)

xycut_hdr, xycut_data = cutouts.xycutout_data(big_data, big_hdr, xstart = 0, xstop = None, ystart = 0, ystop = None)


# construct a 2D header from galfa cube to project each theta slice to
new_header = fits.getheader(galfa_cube_fn)




new_header.remove('CRPIX3') # remove all 3rd axis keywords from fits header
new_header.remove('CTYPE3')
new_header.remove('CRVAL3')
new_header.remove('CDELT3')
new_header.remove('NAXIS3')
new_header.remove('CROTA3')
new_header['NAXIS'] = 2

print(new_header)