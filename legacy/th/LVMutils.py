#
# LVMutils.py - Utility Subroutines for LVM
#
# This file contains useful routines for LVM calculations. In particular, it provides
# coordinate transforms:
#
#  Pointing (look direction) Transformations
#
#   RAdec2AzEl   - Returns the Azimuth and Elevation for the supplied RA and Dec
#   AzEl2RAdec   - Returns the RA/Dec for the supplied Az-El
#
#   AzEl2SazSel  - Returns the siderostat coordinates (Saz, Sel) for the supplied Az-El
#   SazSel2AzEl  - Returns the Az-El for the supplied siderostat coordinates (Saz, Sel)
#
#   calcParallactic - Calculates the parallactic angle, given RA, Dec, Latitued, and LST
#
#  Orientation Transformations (for example, which way is north? which way is perpendicular to the horizon)
#
#
# Note that most of these transformations come from the document LVM0040_LVM Coordinates.
#
# Note that we (finally) embrace Gaia-like rationality and use DECIMAL DEGREES only! Nevertheless,
#      there are handy utility routines to go back to the old way:
#
#   dcm2hms - Returns strings in the form "HH MM SS +DD MM SS" for the supplied decimal RA and Dec
#   hms2dcm - Returns decimal degrees (float) for the supplied string "HH MM SS +DD MM SS"


import numpy as np

#-----------------------------------------------------------------------------#    

# Observatory Constants

LCO_latD = -29.0182          # Latitude
LCO_lat = np.radians(LCO_latD)

lngD = -70.6915              # Longitude
lng = np.radians(lngD)

hgt = 2380.0                 # Altitude (m)

#-----------------------------------------------------------------------------#    
#                         COORDINATE TRANSFORMS
#-----------------------------------------------------------------------------#    

def RAdec2AzEl(raD,decD,lstD):

    '''
        Returns the Azimuth and Elevation for the supplied RA, Dec, and sidereal time

            Inputs: ra,dec   - Right Ascension and Declination in decimal degrees
                    lst      - Local Sideral Time in decimal degrees

            Returns:  az,el  - Azimuth and Elevation (Altitude) in decimal degrees
    '''

    ra,dec,lst = np.radians(raD),np.radians(decD),np.radians(lstD)   # Convert to radians

    # See LVM0040_LVM Coordinates for the transformations below

    ha = lst - ra

    el = np.arcsin( np.sin(dec)*np.sin(LCO_lat) + np.cos(dec)*np.cos(LCO_lat)*np.cos(ha) )

    rat = ( np.sin(dec) - np.sin(el)*np.sin(LCO_lat) ) / ( np.cos(el)*np.cos(LCO_lat) )   # Ratio - need to pin [-1,1]

    if rat<-1.0:          # Goes wonky if roundoff puts it outside [1,1]
      rat = -1.0
    if rat > 1.0:
      rat = 1.0

    if np.sin(ha) < 0.0:

        az = np.arccos( rat )

    else:

        az = 2.0*np.pi - np.arccos( rat )


    return np.degrees(az),np.degrees(el)

#-----------------------------------------------------------------------------#    

def AzEl2RAdec(azD,elD,lstD):

    '''
        Returns the RA/Dec for the supplied Az-El and sidereal time

            Inputs:   az,el    - Azimuth and Elevation (Altitude) in decimal degrees
                      lst      - Local Sideral Time in decimal degrees

            Returns:  ra,dec   - Right Ascension and Declination in decimal degrees
    '''

    az,el,lst = np.radians(azD),np.radians(elD),np.radians(lstD)   # Convert to radians

    # See LVM0040_LVM Coordinates for the transformations below

    dec = np.arcsin( np.sin(el)*np.sin(LCO_lat) + np.cos(el)*np.cos(LCO_lat)*np.cos(az) )

    rat = ( np.sin(el) - np.sin(dec)*np.sin(LCO_lat) ) / ( np.cos(dec)*np.cos(LCO_lat) )  # Ratio

    # if rat<-1.0:          # Goes wonky if roundoff puts it outside [1,1]
    #   rat = -1.0
    # if rat > 1.0:
    #   rat = 1.0

    if np.sin(az) < 0.0:

        ha = np.arccos( rat )

    else:

        ha = 2.0*np.pi - np.arccos( rat )

    ra = lst - ha

    raD,decD = np.degrees(ra),np.degrees(dec)

    if raD<0.0:
        raD+=360.0     # Normalize to 0-360

    return raD,decD

#-----------------------------------------------------------------------------#    

def AzEl2SazSel(azD,elD):

    '''
        Returns the siderostat coordinates (Saz, Sel) for the supplied Az-El

            Inputs:   az,el       - Azimuth and Elevation (Altitude) in decimal degrees
            Returns:  sazD,selD   - Siderostat angles in degrees
    '''

    r90 = np.radians(90.0)        # 90 deg in radians

    az,el = np.radians(azD),np.radians(elD)          # Convert to radians

    SEl = np.arccos( np.cos(el) * np.cos(az) ) - r90   # SEl in radians

    rat = np.sin(el) / np.cos(SEl)                     # Ratio

    if azD < 180.0:
        SAz = r90 - np.arcsin(rat)      # SAz in radians
    else:
        SAz = np.arcsin(rat) - r90


    return np.degrees(SAz),np.degrees(SEl)           # Return values in degrees

#-----------------------------------------------------------------------------#    

def SazSel2AzEl(sazD,selD):

    '''
        Returns the Az-El for the supplied siderostat coordinates (Saz, Sel)

            Inputs:   sazD,selD  - Siderostat angles in degrees 
            Returns:  azD,elD    - Azimuth and Elevation (Altitude) in decimal degrees
    '''

    saz,sel = np.radians(sazD),np.radians(selD)     # Convert to radians

    el = np.arcsin( np.cos(sel) * np.cos(saz) )     # Elevation in radians

    rat = -1.0 * np.sin(sel) / np.cos(el)           # Ratio in Xform

    if saz<0.0:
        az = np.arcos(rat)        # Azimuth in radians
    else:
        az = np.radians(360.0) - np.arcos(rat)


    return np.degrees(az),np.degrees(el)            # Return values in degrees

#-----------------------------------------------------------------------------#    

def calcParallactic(raD,decD,lstD):

    '''
        Returns the parallactic angle for the supplied RA, dec, and LST

            Inputs:   raD,decD,lstD  - RA, dec, and LST 
            Returns:  parD           - Parallactic angle in degrees
    '''

    ra,dec,lst = np.radians(raD),np.radians(decD),np.radians(lstD)   # Convert to radians

    ha = lst - ra       # Calculate hour angle (radians)

    par = np.arctan2(np.sin(ha),(np.cos(dec)*np.tan(LCO_lat)-np.sin(dec)*np.cos(ha)))

    return np.degrees(par)

#-----------------------------------------------------------------------------#    
#                         MATRIX OPTICS STUFF
#-----------------------------------------------------------------------------#    

def MirMat(nx,ny,nz):

    '''
        Returns the mirror matrix for the supplied normal vector.
    '''

    M = np.array( [ [ (1.0-2*nx*nx) , (-2.0*nx*ny)  , (-2.0*nx*nz) ],  \
                    [ (-2.0*nx*ny)  , (1.0-2*ny*ny) , (-2.0*ny*nz) ],  \
                    [ (-2.0*nx*nz)  , (-2.0*ny*nz)  , (1.0-2*nz*nz)] ] )

    return M

#----------------

def rotXmat(alpha):

    '''
        Returns the rotation matrix for angle alpha (in degrees) around the X-axis
    '''

    th = np.deg2rad(alpha)    # Convert to radians

    rX = np.array( [ [ 1.0 , 0.0  , 0.0 ],  \
                     [ 0.0 , np.cos(th) , -1.0*np.sin(th) ],  \
                     [ 0.0 , np.sin(th)  , np.cos(th) ] ] )

    return (rX)               # Send it back

#----------------   

def rotYmat(beta):

    '''
        Returns the rotation matrix for angle beta (in degrees) around the Y-axis
    '''

    th = np.deg2rad(beta)    # Convert to radians

    rY = np.array( [ [ np.cos(th) , 0.0 , np.sin(th) ], \
                     [ 0.0 , 1.0 , 0.0] , \
                     [ -1.0*np.sin(th) , 0.0 , np.cos(th) ] ] )

    return (rY)               # Send it back

#----------------   

def rotZmat(gamma):

    '''
        Returns the rotation matrix for angle gamma (in degrees) around the Z-axis
    '''

    th = np.deg2rad(gamma)    # Convert to radians

    rZ = np.array( [ [ np.cos(th) , -1.0*np.sin(th), 0.0 ], \
                     [ np.sin(th) , np.cos(th) , 0.0 ], \
                     [ 0.0 , 0.0 , 1.0 ] ] )

    return (rZ)               # Send it back
