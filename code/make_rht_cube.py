from __future__ import division
import numpy as np
from astropy.io import fits

import sys 
sys.path.insert(0, '../../FITSHandling/code')
import cutouts

path_to_galfa_cubes = "/disks/jansky/a/users/goldston/DR2W_RC5/Wide/"
path_to_rht_thetaslices = "/disks/jansky/a/users/goldston/susan/Wide_maps/"

# eventually we will step through these cubes -- start w/ 1 test
galfa_cube_name = "GALFA_HI_RA+DEC_356.00+34.35_W"
galfa_cube_fn = path_to_galfa_cubes + galfa_cube_name + ".fits"

galfa_cube_hdr = fits.getheader(galfa_cube_fn)

# R(theta) cube dimensions
nthets = 165
rht_data_cube = np.zeros((nthets, galfa_cube_hdr['NAXIS2'], galfa_cube_hdr['NAXIS1']), np.float_)

# construct a 2D header from galfa cube to project each theta slice to
new_header = fits.getheader(galfa_cube_fn)

print(new_header)


new_header.remove('CRPIX3') # remove all 3rd axis keywords from fits header
new_header.remove('CTYPE3')
new_header.remove('CRVAL3')
new_header.remove('CDELT3')
new_header.remove('NAXIS3')
new_header.remove('CROTA3')
new_header['NAXIS'] = 2

print(new_header)