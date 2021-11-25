#
#  plotRates.py - Visualize the LVM tracking rates
#
#  The goal here is to derive the required tracking rates for various LVM telescope
#  mechanisms, such as the siderostat axes, the ADC, and the K-mirror. We would also
#  like to confirm that things don't diverge at zenith.
#
#  The calculation strategy is pretty simple. To determine the track rate in degrees
#  per second, we calculate the angle (azimuth, for example) at a certain time and then
#  again one second later. The difference in angle is then the number of degrees/sec
#  that the mechanism needs to track.
#
#  Actually, to avoid crossing +/- 180 degrees and other anomalies, we make a step of
#  0.1 seconds and multiply by 10. This small step is sufficient (by experiment) to
#  not cause this wrapping issue with the current sampling. Obviously, if we are plotting
#  a non-differential quantity (like just the azimuth), we don't do this scaling.
#
#  This code generated the grayscale plots in Section 6


import numpy as np
import matplotlib.pyplot as plt      # For producing contour plots
from matplotlib import ticker, cm    # For plot enhancements

import LVMutils as LVMu              # Contains utility routines

lst = 0.0                            # Assume 0h LST

dT = 360.0 / 24.0 / 36000.0          # This is 1/10 second of arc in degrees

iV = np.array( [0.0, 0.0, 1.0])      # Input vector perpendicular to horizon (z unit vector)

#---- Create numpy arrays to hold calculations

X,Y = np.zeros((180)),np.zeros((180))
ZA,ElA,decA = np.zeros((180,180)),np.zeros((180,180)),np.zeros((180,180))

minmax = []   # List to hold relevant values for calculating min-max above 30 elevation

#---- Generate the base (zero angle) mirror matrices for the SAz and SEl mirrors - See Doc

nx,ny,nz = -1.0/np.sqrt(2.0) , +1.0/np.sqrt(2.0) , 0.0  # Normal unit vector to SAz mirror
mSAz = LVMu.MirMat(nx,ny,nz)                            # Calculate mirror matrix

nx,ny,nz = 1.0/np.sqrt(2.0) , 0.0 , 1.0/np.sqrt(2.0)   # Normal unit vector to SEl mirror
mSEl = LVMu.MirMat(nx,ny,nz)                           # Calculate mirror matrix

# We are going to plot the track rates as a function of Az and El. We use a rectangular
# grid to work with matplotlib
#
# Note that we offset the center by 0.5 degree in both directions to avoid math anomalies.

for i in range(180):         # Step from horizon to horizon in degrees 

    X[i] = 90.5 - i          # How far from center EW

    for j in range(180):     # Perpendicular direction

        Y[j] = 90.5 - j      # How far from center NS

        el1 = 90.0 - np.sqrt(X[i]*X[i] + Y[j]*Y[j])    # This is the elevation
        az1 = np.degrees(np.arctan2(Y[j],X[i]))        #   and azimuth

        if az1 < 0.0:        # Gotta wrap it properly - This produced the correct Fig.11
            az1+=360.0

        ra,dec = LVMu.AzEl2RAdec(az1,el1,lst)          # Calculate RA and dec of this spot
        az2,el2 = LVMu.RAdec2AzEl(ra,dec,lst+dT)       # Now look 1/10 second later

        par1 = LVMu.calcParallactic(ra,dec,lst)        # Parallactic angle for this spot and time
        par2 = LVMu.calcParallactic(ra,dec,lst+dT)     #  and 1/10 second later

        #---- Calculate the SEl and SAz siderostat angles (Section 4.3 of doc)

        sazD1,selD1 = LVMu.AzEl2SazSel(az1,el1)        # Calculate siderostat angles
        sazD2,selD2 = LVMu.AzEl2SazSel(az2,el2)        #    and one second later

        #---- Evaluate the rotated and effective mirror matrices at the two time steps

        rY1 = LVMu.rotYmat(sazD1)                      # Returns the Y rotation matrix for SAz (degrees)
        rY2 = LVMu.rotYmat(sazD2)

        rX1 = LVMu.rotXmat(selD1)                      # Returns the X rotation matrix for SEl (degrees)
        rX2 = LVMu.rotXmat(selD2)

        mpSAz1 = np.matmul(rY1,np.matmul(mSAz,rY1.T))  # Modified SAz mirror matrix - timestep 1
        mpSAz2 = np.matmul(rY2,np.matmul(mSAz,rY2.T))  #   and at timestep 2

        mpSEl1 = np.matmul(rY1,np.matmul(rX1,np.matmul(mSEl,np.matmul(rX1.T,rY1.T))))    # Modified SEl mirror matrix
        mpSEl2 = np.matmul(rY2,np.matmul(rX2,np.matmul(mSEl,np.matmul(rX2.T,rY2.T))))    #   and timestep 2

        mEff1 = np.matmul(mpSAz1,mpSEl1)               # Effective mirror matrix for system - see doc
        mEff2 = np.matmul(mpSAz2,mpSEl2)               #   and at timestep 2

        #---- And now, rotate the Z unit-vector and calculate the angle to put it up

        rV1 = np.matmul(mEff1,iV)     # These are the rotated z-unit vectors
        rV2 = np.matmul(mEff2,iV)

        zu1 = np.degrees(np.arctan2(rV1[0],rV1[2]))     # This is the angle to straight up
        zu2 = np.degrees(np.arctan2(rV2[0],rV2[2]))     #    and 1/10 second later

        nu1 = zu1 - par1         # Subtract parallactic to get angle to put North up
        nu2 = zu2 - par2         #   and 1/10 sec later

        #---- Calculate change in angle - Note that we multiply by 10 to get 1 sec (5 for K-mirror)

        dpar = (par2 - par1)*10.0         # Change in parallactic angle
        dzu  = (zu2 - zu1)*10.0           # Change in zenith-up angle
        dnu  = (nu2 - nu1)*5.0            #    and in north-up angle - Note the K-mirror doubling!

        dsaz = (sazD1 - sazD2)*10.0       # How far the SAz axis moved in one timestep
        dsel = (selD1 - selD2)*10.0       #    and same for SEl

        if el1>0.0:                       # Don't go below the horizon!
            ElA[i,j] = el1                # This is to plot reference contours for elevation
            decA[i,j] = dec               # Declination for reference lines of constant dec

            ZA[i,j] = dnu                 # This is what we are plotting

        else:                             # Below horizon - set values to zero
            ElA[i,j] = 0.0
            decA[i,j] = 0.0
            ZA[i,j] = 0.0

        #--- Accumulate stats if above 30 degrees elevation

        if el1>30.0:                      # Accumulate values for statistics
            minmax.append(ZA[i,j])

#--- Done looping through az/alt...present results

print ""
print "Range of values between El=30-90 deg have min ",min([(number) for number in minmax]),"  max absolute = ", max(([(number) for number in minmax]))
print "Absolute values between El=30-90 deg have min ",min([abs(number) for number in minmax]),"  max absolute = ", max(([abs(number) for number in minmax]))
print "Median absolute value is                      ",np.median([abs(number) for number in minmax])
print "5"

#---- Time to generate the visuals...quantity ZA in grayscale, dec in solid lines, and 30,60 degree elevation as dashed

fig = plt.figure()
print "0"

ax = fig.add_subplot(111)
print "1"
cpf = ax.contourf(X,Y,ZA,np.arange(-0.01, 0.01, 0.0005),cmap=cm.Spectral) #Greys_r)   # Grayscale contour map of ZA

# Set the colours of the contours and labels so they're white where the
# contour fill is dark (ZA < 0) and black where it's light (ZA >= 0)
print "2"

colours = ['w' if level<0 else 'k' for level in cpf.levels]

dp = ax.contour(X, Y, decA, 16, colors=colours)                        # Declination lines
cp = ax.contour(X, Y, ElA, 2, colors=colours, linestyles='dashed')     # Plot 30 and 60 degrees elevations
print "3"

ax.clabel(cp, fontsize=12, colors=colours)         # Label axes
ax.clabel(dp, fontsize=12, colors=colours)

cbar=fig.colorbar(cpf)                             # Add a colourbar for ZA
print "4"

plt.savefig('dnu_minus.pdf', bbox_inches='tight')  # Make a pdf file of the result
print "5"

plt.show()

print "6"
