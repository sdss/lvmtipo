from datetime import datetime
from math import nan, cos
import numpy as np

from astropy.coordinates import EarthLocation, SkyCoord, Angle
from astropy.time import Time
from astropy import units as u

from .ambient import Ambient
from .site import Site
from .siderostat import Siderostat


def radec2azel(raD, decD, lstD, site: Site):
    '''
        See LVM0040_LVM Coordinates for the transformations below

        Returns the Azimuth and Elevation for the supplied RA, Dec, and sidereal time
            Inputs: ra,dec   - Right Ascension and Declination in decimal degrees
                    lst      - Local Sideral Time in decimal degrees
            Returns:  az,el  - Azimuth and Elevation (Altitude) in decimal degrees
    '''
    lat_r = np.radians(site.lat)
    
    ra,dec,lst = np.radians(raD),np.radians(decD),np.radians(lstD)   # Convert to radians
    ha = lst - ra
    el = np.arcsin( np.sin(dec)*np.sin(lat_r) + np.cos(dec)*np.cos(lat_r)*np.cos(ha) )
    rat = ( np.sin(dec) - np.sin(el)*np.sin(lat_r) ) / ( np.cos(el)*np.cos(lat_r) )   # Ratio - need to pin [-1,1]
    if rat<-1.0:          # Goes wonky if roundoff puts it outside [1,1]
      rat = -1.0
    if rat > 1.0:
      rat = 1.0
    if np.sin(ha) < 0.0:
        az = np.arccos( rat )
    else:
        az = 2.0*np.pi - np.arccos( rat )
    return np.degrees(az), np.degrees(el)


def azel2sazsel(azD, elD):
    '''
        Returns the siderostat coordinates (saz, Sel) for the supplied Az-El
            Inputs:   az,el       - Azimuth and Elevation (Altitude) in decimal degrees
            Returns:  sazD,selD   - Siderostat angles in degrees
    '''
    r90 = np.radians(90.0)        # 90 deg in radians
    az, el = np.radians(azD),np.radians(elD)          # Convert to radians
    SEl = np.arccos( np.cos(el) * np.cos(az) ) - r90   # SEl in radians
    rat = np.sin(el) / np.cos(SEl)                     # Ratio
    if azD < 180.0:
        SAz = r90 - np.arcsin(rat)      # saz in radians
    else:
        SAz = np.arcsin(rat) - r90

    return np.degrees(SAz),np.degrees(SEl)           # Return values in degrees


def delta_radec2mot_axis(ref_midpoint, new_midpoint):
    sid = Siderostat(azang=0.0)
    site = Site(name='LCO')
    observing_location = EarthLocation(lat=site.lat*u.deg, lon=site.long*u.deg)

    observing_time = Time(datetime.utcnow(), scale='utc', location=observing_location)
    lst = observing_time.sidereal_time('mean')

    ref_az_d, ref_el_d = radec2azel(ref_midpoint.ra.deg, ref_midpoint.dec.deg, lst.deg, site)
#    print(f"new {ref_midpoint.ra.deg} {ref_midpoint.dec.deg} {ref_az_d} {ref_el_d}")
    
    new_az_d, new_el_d = radec2azel(new_midpoint.ra.deg, new_midpoint.dec.deg, lst.deg, site)
#    print(f"new {new_midpoint.ra.deg} {new_midpoint.dec.deg} {new_az_d} {new_el_d}")
    
    ref_saz_d, ref_sel_d = azel2sazsel(ref_az_d, ref_el_d)
    new_saz_d, new_sel_d = azel2sazsel(new_az_d, new_el_d)
    
    saz_diff_d = ref_saz_d - new_saz_d
    sel_diff_d = ref_sel_d - new_sel_d
    
    saz_diff_d *= -3600.
    sel_diff_d *= -3600.

    return Angle(saz_diff_d * u.deg), Angle(sel_diff_d * u.deg)
