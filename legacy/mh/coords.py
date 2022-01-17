#Library imports
import numpy as np
import matplotlib.pyplot as plt
import time
from astropy.io import fits
from astropy.table import Table,vstack
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astroquery.gaia import Gaia
from scipy.ndimage import gaussian_filter, gaussian_gradient_magnitude
from scipy import optimize
import healpy as hp
from astropy.table import Table, hstack, vstack


def in_box(xxs, yys, pa_deg, inst):
    # tells you if an xy coordinate (in the focal plane) is on the guider chip
    # based on the instrument specs (above)
    # xy coordinates follow the definitions laid out in LVM-0040_LVM_Coordinates.pdf
    #
    #input:
    #xxs and yys:  np array with x,y positions (in mm)
    #pa:   position angle (east of north) in degrees
    #
    #returns: np array set to 1 (in box) or 0 (not in box)
    #        x&y positions on the chip
    #        note! These do not account for the 6th mirror, which flips the handedness

    #convert position angle to radians
    pa = pa_deg / 180. * np.pi
    
    #find some vertices, A and B, of the box (see Photo_on13.07.20at13.58.jpg)
    Ar=inst.r_outer
    Atheta=np.arcsin(inst.chip_size_mm[0] / 2. / inst.r_outer)
    
    phi=(np.pi-Atheta) / 2.
    h1 = inst.chip_size_mm[0] / (2. * np.tan(phi))
    h2 = inst.r_outer-inst.chip_size_mm[1] - h1
    
    chi = np.arctan(h2 / (inst.chip_size_mm[0] / 2.))
    
    Br = np.sqrt(inst.chip_size_mm[0]*  inst.chip_size_mm[0] / 2. / 2. + h2 * h2)
    Btheta = np.pi / 2. - chi
            
    
    #convert from polar to cartesian
    Ay = Ar * np.cos(Atheta)
    Ax = Ar * np.sin(Atheta)
    
    By = Br * np.cos(Btheta)
    Bx = Br * np.sin(Btheta)
    
    #print(Ax,Ay,Bx,By)
    
    #are the positions to test within the guider chip?
    #derotate in pa
    rrs = np.sqrt(xxs * xxs + yys * yys)
    thetas = np.arctan2(yys, xxs)
    
    derot_thetas=thetas-pa
    
    derot_xxs = rrs * np.cos(derot_thetas)
    derot_yys = rrs * np.sin(derot_thetas)
        
        
    #compare with box edges
    flagg = np.array(len(xxs) * [False])
    
    ii=((derot_xxs < Ax) & (derot_xxs > -1. * Ax) & (derot_yys < Ay) & (derot_yys > By))
    flagg[ii] = True
        
    #return flags testing whether object is on the chip
    #also return x,y position on chip (in mm), in chip coordinates
    #origin is arbitrarily at (Bx,By) the lower left corner
    #note! this does not account for the 6th mirror, which flips the handedness
    return flagg, derot_xxs + Bx,derot_yys - By


def ad2xy(cats2, c_icrs, inst):
    # converts ra/dec positions to angular offsets from field center (c_icrs)
    #
    # inputs:
    # cats2: table of SkyCoords with ra & dec positions 
    # c_icrs: SkyCoord of the IFU field center
    # returns: np array of x and y positions in focal plane, in mm
    
    #convert ra & dec of guide stars to xy position in relation to field center (in arcsec)

    
    #Without the slow loop (Max)
    dd_y = sphdist(cats2['ra'],c_icrs.dec.deg,cats2['ra'], cats2['dec']) * 3600. #in arcsec 
    dd_y *= (2 * (cats2['dec'] > c_icrs.dec.deg).astype(float) - 1)
    
    dd_x = sphdist(c_icrs.ra.deg,cats2['dec'],cats2['ra'],cats2['dec']) * 3600. #in arcsec
    dd_x *= (2 * (cats2['ra'] < c_icrs.ra.deg).astype(float) - 1)  
    
    #convert to mm
    dd_x_mm = dd_x * inst.image_scale / 1e3
    dd_y_mm = dd_y * inst.image_scale / 1e3
    
    return dd_x_mm, dd_y_mm
    
def sphdist(ra1, dec1, ra2, dec2):
# measures the spherical distance in degrees
# The input has to be in degrees too
    dec1_r = np.deg2rad(dec1)
    dec2_r = np.deg2rad(dec2)
    ra1_r = np.deg2rad(ra1)
    ra2_r = np.deg2rad(ra2)
    return 2 * np.rad2deg(np.arcsin(np.sqrt((np.sin((dec1_r - dec2_r) / 2)) ** 2 + np.cos(dec1_r) * np.cos(dec2_r) * (np.sin((np.deg2rad(ra1 - ra2)) / 2)) ** 2)))
#    return np.rad2deg(np.arccos(np.sin(dec1_r)*np.sin(dec2_r)+np.cos(dec1_r)*np.cos(dec2_r)*np.cos(np.abs(ra1_r-ra2_r))))
    
