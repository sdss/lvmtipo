# Calculates the LMC example in LVM_Coordinates.doc
# 

import numpy as np
import LVMutils as LVMu        # Contains utility routines

np.set_printoptions(precision=3,suppress=True)     # Minimize on-screen digits

#---------------- MAIN ROUTINE STARTS HERE -------------------------

r2d = 180.0 / np.pi       # Converts radians to degrees
r360 = 2 * np.pi          # 360 deg in radians (2pi)
r180 = np.pi              # 180 deg in radians
r90 = np.pi / 2           # 90 deg in radians

latD = LVMu.latD          # Latitude of Las Campanas in degrees

raD = (5.0 + 23.0/60 + 34.5/3600) / 24.0 * 360.0    # RA of LMC  (05h 23m 34.5s) in degrees
decD = -1.0*(69.0 + 45.0/60 + 22.0/3600)            # dec of LMC (-69 45 22) in degrees
haD = -1.0 / 24 * 360                               # HA of LMC in degrees one hour before transit

print ("lat = ",latD)
print ("RA = ",raD)
print ("Dec = ",decD)
print ("HA = ",haD)
print ("")

#---- Calculate Azimuth and Elevation (Section 4.1 of doc)

azD,elD = LVMu.RAdec2AzEl(raD,decD,raD+haD)          # Calculate Az and El (degrees) of target 1 hr before transit

print ("Az = ",azD)
print ("El = ",elD)
print ("")

#---- Calculate the SEl and SAz siderostat angles (Section 4.3 of doc)

sazD,selD = LVMu.AzEl2SazSel(azD,elD)         # Calculate siderostat angles

print ("SAz = ",sazD)
print ("SEl = ",selD)
print ("")

#---- Generate the base mirror matrices for the SAz and SEl mirrors - See Doc

nx,ny,nz = -1.0/np.sqrt(2.0) , +1.0/np.sqrt(2.0) , 0.0  # Normal unit vector to SAz mirror
mSAz = LVMu.MirMat(nx,ny,nz)                            # Calculate mirror matrix

nx,ny,nz = 1.0/np.sqrt(2.0) , 0.0 , 1.0/np.sqrt(2.0)   # Normal unit vector to SEl mirror
mSEl = LVMu.MirMat(nx,ny,nz)                           # Calculate mirror matrix

print ("mSAz ")
print (mSAz)
print ("")

print ("mSEl ")
print (mSEl)
print ("")

#---- Evaluate the rotated and effective mirror matrices

rY = LVMu.rotYmat(sazD)                       # Returns the Y rotation matrix for SAz (degrees)
mpSAz = np.matmul(rY,np.matmul(mSAz,rY.T))    # Modified SAz mirror matrix

rX = LVMu.rotXmat(selD)                       # Returns the X rotation matrix for SEl (degrees)

rXmat = np.matmul(rX,np.matmul(mSEl,rX.T))    # Check pure rotation matrix

print ("Pure X / El rotation ")
print (rXmat)
print ("")

mpSEl = np.matmul(rY,np.matmul(rX,np.matmul(mSEl,np.matmul(rX.T,rY.T))))    # See LVM Coordinates doc

mEff = np.matmul(mpSAz,mpSEl)                 # Effective mirror matrix for system - see doc

print ("rY ")
print (rY)
print ("")

print ("rX ")
print (rX)
print ("")

print ("MSAz mod ")
print (mpSAz)
print ("")

print ("MSEl mod ")
print (mpSEl)
print ("")

#--- Calculate effect on unit vector(s)

iV = np.array( [0.0, 0.0 , 1.0])     # Input vector (z unit vector)

print ("Meff")
print (mEff)
print ("")
print ("Output Z unit-vector")
print (np.matmul(mEff,iV))
