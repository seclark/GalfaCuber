from __future__ import division
import numpy as np
from astropy.io import fits
import re

import sys 
sys.path.insert(0, '../../FITSHandling/code')
import cutouts

import galfa_vel_helpers

class Cube():
    """
    Single data cube defined by coordinates of PPV cube
    """
    def __init__(self, RA="180.00", DEC="02.35"):
        
        # string center cube coordinates
        self.centerRA = RA
        self.centerDEC = DEC
        
        # define paths
        self.path_to_galfa_cubes = "/disks/jansky/a/users/goldston/DR2W_RC5/Wide/"
        self.ppv_cube_name = "GALFA_HI_RA+DEC_"+self.centerRA+"+"+self.centerDEC+"_W"
        self.ppv_cube_fn = self.path_to_galfa_cubes + self.ppv_cube_name + ".fits"
        
        # get PPV cube header 
        self.ppv_cube_hdr = fits.getheader(self.ppv_cube_fn)
        
        # cube dimensions
        self.naxis1 = self.ppv_cube_hdr['NAXIS1']
        self.naxis2 = self.ppv_cube_hdr['NAXIS2']
        
        # get allsky RHT data header
        self.path_to_rht_thetaslices = "/disks/jansky/a/users/goldston/susan/Wide_maps/single_theta_maps/"
        self.galfa_allsky_hdr = fits.getheader(self.path_to_rht_thetaslices+"S0974_0978/intrht_S0974_0978.fits")

    def get_cube_coordinates_in_allsky(self):
        """
        x, y start, stop in coordinates of allsky data
        """
    
        # Actually grab cube edges
        cube_w = cutouts.make_wcs(self.ppv_cube_fn)
        cube_edgera1, cube_edgedec1 = cutouts.xy_to_radec(0, 0, cube_w)
        cube_edgera2, cube_edgedec2 = cutouts.xy_to_radec(self.naxis1, self.naxis2, cube_w)

        # translate cube corners to allsky x, y
        allsky_w = cutouts.make_wcs(self.galfa_allsky_hdr)
        cutout_x1, cutout_y1 = cutouts.radec_to_xy(cube_edgera1, cube_edgedec1, allsky_w)
        cutout_x2, cutout_y2 = cutouts.radec_to_xy(cube_edgera2, cube_edgedec2, allsky_w)

        # ROUND to integer x, y. Necessary because allsky header does not quite match up. Valid per Yong's check.
        self.cutout_xstart = np.int(np.round(cutout_x1))
        self.cutout_xstop = np.int(np.round(cutout_x2))
        self.cutout_ystart = np.int(np.round(cutout_y1))
        self.cutout_ystop = np.int(np.round(cutout_y2))
    
    def make_RHT_XYT_cube(self, rht_velstr="S0974_0978", verbose=False):
        """
        make a cube of dimensions (nx, ny, ntheta)
        inputs: rht_velstr :: velocity range string for RHT channel map
        """
        
        self.rht_velstr = rht_velstr
        
        # There are 165 theta bins in RHT data
        self.nthets = 165
        
        # Empty cube dimensions
        self.rht_data_cube = np.zeros((self.nthets, self.naxis2, self.naxis1), np.float_)
        
        # Grab new data
        for thet_i in xrange(self.nthets):
            allsky_fn = self.path_to_rht_thetaslices + self.rht_velstr + "/GALFA_HI_W_"+rht_velstr+"_newhdr_SRcorr_w75_s15_t70_theta_"+str(thet_i)+".fits"
            allsky_thetaslice_data = fits.getdata(allsky_fn)
            allsky_thetaslice_hdr = fits.getheader(allsky_fn)
    
            xycut_hdr, xycut_data = cutouts.xycutout_data(allsky_thetaslice_data, allsky_thetaslice_hdr, xstart=self.cutout_xstart, xstop=self.cutout_xstop, ystart=self.cutout_ystart, ystop=self.cutout_ystop)
            if verbose:
                print(self.cutout_xstart, self.cutout_xstop, self.cutout_ystart, self.cutout_ystop)
            
            self.rht_data_cube[thet_i, :, :] = xycut_data

        # Create new HDU object
        hdu = fits.PrimaryHDU(self.rht_data_cube)
        hdulist = fits.HDUList([hdu])
        priheader = hdulist[0].header

        priheader.set('NAXIS', 3)
        priheader.set('CTYPE1', 'RA      ')
        priheader.set('CRPIX1', self.ppv_cube_hdr['CRPIX1'])
        priheader.set('CRVAL1', self.ppv_cube_hdr['CRVAL1'])
        priheader.set('CDELT1', self.ppv_cube_hdr['CDELT1'])

        priheader.set('CTYPE2', 'DEC     ')
        priheader.set('CRPIX2', self.ppv_cube_hdr['CRPIX2'])
        priheader.set('CRVAL2', self.ppv_cube_hdr['CRVAL2'])
        priheader.set('CDELT2', self.ppv_cube_hdr['CDELT2'])

        priheader.set('NAXIS3', self.nthets)
        priheader.set('CDELT3', np.pi/self.nthets)
        priheader.set('CTYPE3', 'THETARHT')
        priheader.set('CRVAL3', 0.000000)
        priheader.set('CRPIX3', 1.000000)

        with open('../text/newhistory.txt') as histtext:
            allhistory = histtext.readlines()

        # strip /n characters, 'HISTORY'
        allhistory = [x.strip() for x in allhistory] 
        allhistory = [x.replace('HISTORY ', '') for x in allhistory] 

        for line in allhistory:
            priheader.set('HISTORY', line)
            
        self.hdulist = hdulist

    def get_RHT_XYT_cube(self, ashdulist=False):
        """
        retrieve (x, y, theta) cube.
        """
        if ashdulist: 
            return self.hdulist
        else: 
            return self.rht_data_cube
            
    def make_RHT_IQU_cubes(self, verbose=False):
        """
        make a cube of dimensions (nx, ny, ntheta)
        inputs: rht_velstr :: velocity range string for RHT channel map
        """
        
        # Empty cube dimensions - there are 21 velocity channels
        self.nchannels = 21
        self.rht_I_cube = np.zeros((self.nchannels, self.naxis2, self.naxis1), np.float_)
        self.rht_Q_cube = np.zeros((self.nchannels, self.naxis2, self.naxis1), np.float_)
        self.rht_U_cube = np.zeros((self.nchannels, self.naxis2, self.naxis1), np.float_)
        
        # Grab new data
        for vel_i, rht_velstr in enumerate(galfa_vel_helpers.all_rht_velstrs):
            allsky_I_fn = self.path_to_rht_thetaslices + rht_velstr + "/intrht_"+rht_velstr+".fits"
            allsky_Q_fn = self.path_to_rht_thetaslices + rht_velstr + "/QRHT_"+rht_velstr+".fits"
            allsky_U_fn = self.path_to_rht_thetaslices + rht_velstr + "/URHT_"+rht_velstr+".fits"
            allsky_I_data = fits.getdata(allsky_I_fn)
            allsky_Q_data = fits.getdata(allsky_Q_fn)
            allsky_U_data = fits.getdata(allsky_U_fn)
            allsky_thetaslice_hdr = fits.getheader(allsky_I_fn)
    
            xycut_hdr, xycut_I = cutouts.xycutout_data(allsky_I_data, allsky_thetaslice_hdr, xstart=self.cutout_xstart, xstop=self.cutout_xstop, ystart=self.cutout_ystart, ystop=self.cutout_ystop)
            xycut_hdr, xycut_Q = cutouts.xycutout_data(allsky_Q_data, allsky_thetaslice_hdr, xstart=self.cutout_xstart, xstop=self.cutout_xstop, ystart=self.cutout_ystart, ystop=self.cutout_ystop)
            xycut_hdr, xycut_U = cutouts.xycutout_data(allsky_U_data, allsky_thetaslice_hdr, xstart=self.cutout_xstart, xstop=self.cutout_xstop, ystart=self.cutout_ystart, ystop=self.cutout_ystop)
            #if verbose:
            #    print(self.cutout_xstart, self.cutout_xstop, self.cutout_ystart, self.cutout_ystop)
            
            self.rht_I_cube[vel_i, :, :] = xycut_I
            self.rht_Q_cube[vel_i, :, :] = xycut_Q
            self.rht_U_cube[vel_i, :, :] = xycut_U
            
        # Create new HDU object
        hdu_I = fits.PrimaryHDU(self.rht_I_cube)
        hdu_Q = fits.PrimaryHDU(self.rht_Q_cube)
        hdu_U = fits.PrimaryHDU(self.rht_U_cube)
        hdulist_I = fits.HDUList([hdu_I])
        hdulist_Q = fits.HDUList([hdu_Q])
        hdulist_U = fits.HDUList([hdu_U])
        
        priheader_I = hdulist_I[0].header
        priheader_Q = hdulist_Q[0].header
        priheader_U = hdulist_U[0].header
        
        # start and end velocities in meters
        startvel = np.float(galfa_vel_helpers.get_galfa_W_truevel(974))*1000
        stopvel = np.float(galfa_vel_helpers.get_galfa_W_truevel(1078))*1000
        velstep = (stopvel - startvel)/self.nchannels
        
        for priheader in [priheader_I, priheader_Q, priheader_U]:

            priheader.set('NAXIS', 3)
            priheader.set('CTYPE1', 'RA      ')
            priheader.set('CRPIX1', self.ppv_cube_hdr['CRPIX1'])
            priheader.set('CRVAL1', self.ppv_cube_hdr['CRVAL1'])
            priheader.set('CDELT1', self.ppv_cube_hdr['CDELT1'])

            priheader.set('CTYPE2', 'DEC     ')
            priheader.set('CRPIX2', self.ppv_cube_hdr['CRPIX2'])
            priheader.set('CRVAL2', self.ppv_cube_hdr['CRVAL2'])
            priheader.set('CDELT2', self.ppv_cube_hdr['CDELT2'])

            priheader.set('NAXIS3', self.nchannels)
            priheader.set('CDELT3', velstep)
            priheader.set('CTYPE3', 'VELO-LSR')
            priheader.set('CRVAL3', startvel)
            priheader.set('CRPIX3', 1.000000)
            
            priheader.set('EQUINOX', 2000.00)
            
            with open('../text/newhistory.txt') as histtext:
                allhistory = histtext.readlines()

                # strip /n characters, 'HISTORY'
                allhistory = [x.strip() for x in allhistory] 
                allhistory = [x.replace('HISTORY ', '') for x in allhistory] 

                for line in allhistory:
                    priheader.set('HISTORY', line)
            
        self.hdulist_I = hdulist_I
        self.hdulist_Q = hdulist_Q
        self.hdulist_U = hdulist_U
            
            
    def get_RHT_IQU_cubes(self, ashdulist=False):
        """
        retrieve (x, y, I), (x, y, Q), and (x, y, U) cubes.
        """
        if ashdulist: 
            return self.hdulist_I, self.hdulist_Q, self.hdulist_U
        else: 
            return self.rht_I_cube, self.rht_Q_cube, self.rht_U_cube
            

def make_single_cube_rtheta(RA="180.00", DEC="02.35", rht_velstart="0974", rht_velstop="0978", verbose=False):
    """
    create and save R(theta) cube for single RHT velocity range.
    """
    velstart = galfa_vel_helpers.galfa_name_dict[rht_velstart]
    velstop = galfa_vel_helpers.galfa_name_dict[rht_velstop]
    rht_velstr = "S"+rht_velstart+"_"+rht_velstop
    
    if verbose:
        print("velocities {} to {}".format(velstart, velstop))
    
    cube = Cube(RA=RA, DEC=DEC)
    cube.get_cube_coordinates_in_allsky()
    cube.make_RHT_XYT_cube(rht_velstr=rht_velstr, verbose=False)
    hdulist = cube.get_RHT_XYT_cube(ashdulist = True)
    
    outroot = "/disks/jansky/a/users/goldston/susan/RHT_RC1/Rtheta_cubes/"
    outfn = outroot + rht_velstr + "/GALFA-HI_RHT_spect_v"+velstart+"_"+velstop+"_RA+DEC_"+RA+"+"+DEC+".fits"
    
    if verbose:
        print("Saving to {}".format(outfn))
    
    hdulist.writeto(outfn)

def make_single_cube_IQU(RA="180.00", DEC="02.35", verbose=False):
    """
    create and save I, Q, U cube
    """
    cube = Cube(RA=RA, DEC=DEC)
    cube.get_cube_coordinates_in_allsky()
    cube.make_RHT_IQU_cubes(verbose=verbose)
    hdulist_I, hdulist_Q, hdulist_U = cube.get_RHT_IQU_cubes(ashdulist=True)
    
    outroot = "/disks/jansky/a/users/goldston/susan/RHT_RC1/IQU_cubes/"
    outfn_I = outroot + "GALFA-HI_RHT_I_RA+DEC_"+RA+"+"+DEC+".fits"
    outfn_Q = outroot + "GALFA-HI_RHT_Q_RA+DEC_"+RA+"+"+DEC+".fits"
    outfn_U = outroot + "GALFA-HI_RHT_U_RA+DEC_"+RA+"+"+DEC+".fits"
    
    hdulist_I.writeto(outfn_I)
    hdulist_Q.writeto(outfn_Q)
    hdulist_U.writeto(outfn_U)
    

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
    all_RAs = ["{0:0=3d}.00".format(ra) for ra in np.arange(12, 350, 8)]
    #all_RAs = ["{0:0=3d}.00".format(ra) for ra in np.arange(180, 360, 8)]

    for ra in all_RAs:
        for dec in all_DECs:
            make_single_cube_rtheta(RA=ra, DEC=dec, rht_velstart="1029", rht_velstop="1033", verbose=True)
            #make_single_cube_IQU(RA=ra, DEC=dec, verbose=True)

    #RA = "156.00"
    #DEC = "26.35"
    #make_single_cube_IQU(RA=RA, DEC=DEC, verbose=True)
    
    
    
