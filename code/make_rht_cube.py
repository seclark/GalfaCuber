from __future__ import division
import numpy as np
from astropy.io import fits

import sys 
sys.path.insert(0, '../../FITSHandling/code')
import cutouts

path_to_galfa_cubes = "/disks/jansky/a/users/goldston/DR2W_RC5/Wide/"
path_to_rht_thetaslices = "/disks/jansky/a/users/goldston/susan/Wide_maps/single_theta_maps/"
path_to_nhi_allskymaps = "/disks/jansky/a/users/goldston/zheng/151019_NHImaps_SRcorr/data/GNHImaps_SRCORR_final/NHImaps/"

# eventually we will step through these cubes -- start w/ 1 test
#galfa_cube_name = "GALFA_HI_RA+DEC_356.00+34.35_W"
galfa_cube_name = "GALFA_HI_RA+DEC_044.00+02.35_W"
galfa_cube_fn = path_to_galfa_cubes + galfa_cube_name + ".fits"

galfa_cube_hdr = fits.getheader(galfa_cube_fn)
galfa_allsky_hdr = fits.getheader(path_to_rht_thetaslices+"S0974_0978/intrht_S0974_0978.fits")
#galfa_allsky_hdr = fits.getheader(path_to_nhi_allskymaps+"GALFA-HI_NHI_VLSR-90+90kms.fits")

# R(theta) cube dimensions
nthets = 165
rht_velstr = "S0974_0978"
rht_data_cube = np.zeros((nthets, galfa_cube_hdr['NAXIS2'], galfa_cube_hdr['NAXIS1']), np.float_)


# Obtain x, y values in *big* GALFA data. These should be integers.
cube_crval1 = galfa_cube_hdr['CRVAL1']
cube_crval2 = galfa_cube_hdr['CRVAL2']
cube_crpix1 = galfa_cube_hdr['CRPIX1']

#test
print(cube_crval1, cube_crval2)
cube_w = cutouts.make_wcs(galfa_cube_fn)
cube_edgex1, cube_edgey1 = cutouts.radec_to_xy(cube_crval1, cube_crval2, cube_w)
print(cube_edgex1, cube_edgey1)
cube_edgera1, cube_edgedec1 = cutouts.xy_to_radec(cube_edgex1, cube_edgey1, cube_w)
print(cube_edgera1, cube_edgedec1)

print(galfa_allsky_hdr['CTYPE1'])
galfa_allsky_hdr['CTYPE1'] = 'RA      '
galfa_allsky_hdr['CTYPE2'] = 'DEC     '
allsky_w = cutouts.make_wcs(galfa_allsky_hdr)
cutout_x1, cutout_y1 = cutouts.radec_to_xy(cube_crval1, cube_crval2, allsky_w)
print(cutout_x1, cutout_y1)

# Actually grab cube edges
cube_w = cutouts.make_wcs(galfa_cube_fn)
cube_edgera1, cube_edgedec1 = cutouts.xy_to_radec(0, 0, cube_w)
cube_edgera2, cube_edgedec2 = cutouts.xy_to_radec(galfa_cube_hdr["NAXIS1"], galfa_cube_hdr["NAXIS2"], cube_w)

# translate cube corners to allsky x, y
allsky_w = cutouts.make_wcs(galfa_allsky_hdr)
cutout_x1, cutout_y1 = cutouts.radec_to_xy(cube_edgera1, cube_edgedec1, allsky_w)
cutout_x2, cutout_y2 = cutouts.radec_to_xy(cube_edgera2, cube_edgedec2, allsky_w)

# ROUND to integer x, y. Necessary because allsky header does not quite match up. Should check.
print(cutout_x1, cutout_x2, cutout_y1, cutout_y2)
print(np.int(np.round(cutout_x1)), np.int(np.round(cutout_x2)), np.int(np.round(cutout_y1)), np.int(np.round(cutout_y2)))
cutout_xstart = np.int(cutout_x1)

#xycut_hdr, xycut_data = cutouts.xycutout_data(big_data, big_hdr, xstart = 0, xstop = None, ystart = 0, ystop = None)
"""
# Grab new data
for thet_i in xrange(nthets):
    allsky_fn = path_to_rht_thetaslices + "GALFA_HI_allsky_-10_10_w75_s15_t70_thetabin_"+str(thet_i)+".fits"
    allsky_thetaslice_data = fits.getdata(allsky_fn)
    allsky_thetaslice_hdr = fits.getheader(allsky_fn)
    
    print('theta_i', thet_i)
    
    # Reproject each theta slice into appropriate theta bin in cube
    output, footprint = reproject_interp((allsky_thetaslice_data, allsky_thetaslice_hdr), new_header)
    rht_data_cube[thet_i, :, :] = output

new_hdr = copy.copy(galfa_cube_hdr)    
new_hdr['NAXIS3'] = nthets
new_hdr['CTYPE3'] = 'THETARHT'
new_hdr['CRVAL3'] = 0.000000
new_hdr['CRPIX3'] = 0.000000
new_hdr['CDELT3'] = np.pi/nthets


# Create new HDU object
hdu = fits.PrimaryHDU(rht_data_cube)
hdulist = fits.HDUList([hdu])
priheader = hdulist[0].header

print(priheader)
priheader.set('NAXIS', 2)

new_header.remove('CRPIX3') # remove all 3rd axis keywords from fits header
new_header.remove('CTYPE3')
new_header.remove('CRVAL3')
new_header.remove('CDELT3')
new_header.remove('NAXIS3')
new_header.remove('CROTA3')
new_header['NAXIS'] = 2
"""
