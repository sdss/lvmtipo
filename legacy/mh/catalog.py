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
#Instrument specs


def get_cat_using_healpix(c:SkyCoord, inst, plotflag=False):
    vec = hp.ang2vec(np.deg2rad(c.dec.value+90),np.deg2rad(c.ra.value))

    ipix_disc = hp.query_disc(nside=64, vec=vec, radius=np.deg2rad(inst.outer_search_radius),inclusive=True)
    print(ipix_disc)

    if plotflag:
        fig,ax = plt.subplots(figsize=(12,12))
        ax.set_aspect("equal")
        ax.axhline(c.dec.value)
        ax.axvline(c.ra.value)
        
    counter=0
    for ipix in ipix_disc:
        filename = inst.catalog_path + "/Gaia_Healpix_64/{:06d}.fits".format(ipix)

        hdul = fits.open(filename)
        data= Table(hdul[1].data)
        print(filename,len(data))
        #data = data.filled()
        if plotflag:
            ax.plot(data["ra"],data["dec"],".")
        if counter==0:
            data_combined = data
            counter+=1
        else:
            data_combined = vstack([data_combined, data])
    return data_combined


def get_cat_using_healpix2(c:SkyCoord, inst, plotflag=False, verbose=False):
    vec = hp.ang2vec(np.deg2rad(-c.dec.value + 90), np.deg2rad(c.ra.value))

    ipix_disc = hp.query_disc(nside=64, vec=vec, radius=np.deg2rad(inst.outer_search_radius),inclusive=True,nest=True)
    if verbose: 
        print(ipix_disc)

    if plotflag:
        fig,ax = plt.subplots(figsize=(12,12))
        ax.set_aspect("equal")
        ax.axhline(c.dec.value)
        ax.axvline(c.ra.value)
        
    counter=0
    for ipix in ipix_disc:
        filename = inst.catalog_path + "/Gaia_Healpix_6/lvl6_{:06d}.npy".format(ipix)
        data = np.load(filename)
        #hdul = fits.open(filename)
        #data= Table(hdul[1].data)
        if verbose: print(filename,len(data))
        #data = data.filled()
        if plotflag:
            ax.plot(data["ra"],data["dec"],".")
        if counter==0:
            data_combined = data
            counter+=1
        else:
            #data_combined = vstack([data_combined, data])
            data_combined = np.concatenate([data_combined, data])
    return data_combined


def calc_sn(gmag, inst, n_pix=7*7, sky_flux=10, exp_time=5):
        gaia_flux = 10**(-(gmag+inst.zp)/2.5)
        background = (sky_flux+inst.dark_current)*exp_time
        background_noise = np.sqrt(background+inst.readout_noise**2)
        signal = gaia_flux*exp_time
        noise = np.sqrt(inst.readout_noise**2+signal+n_pix*background)
        sn = signal/noise
        return sn
