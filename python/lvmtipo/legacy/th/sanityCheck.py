'''
    sanityCheck.py is based on plotRates.py. See that script and the document
    for more information.
      
    Instead of calculating things at all azimuths and elevations, this program
    propagates the N and E vector through the system just prior to and then
    just after transit.

'''

import numpy as np
import matplotlib.pyplot as plt      # For producing contour plots

import LVMutils as LVMu              # Contains utility routines

np.set_printoptions(precision=3)     # Prett(ier) print

nV = np.array([0.0, 1.0, 0.0])       # On-sky North-pointing (y) unit vector for sanity test
eV = np.array([1.0, 0.0, 0.0])       # On-sky East-pointing (x) unit vector for sanity test

az1,el1 = 90.0,89.99                 # A spot just before transit near zenith
az2,el2 = 270.0,89.99                #  and now just after transit on the western side

#---- Generate the base (zero angle) mirror matrices for the SAz and SEl mirrors - See Doc

nx,ny,nz = -1.0/np.sqrt(2.0) , +1.0/np.sqrt(2.0) , 0.0  # Normal unit vector to SAz mirror
mSAz = LVMu.MirMat(nx,ny,nz)                            # Calculate mirror matrix

nx,ny,nz = 1.0/np.sqrt(2.0) , 0.0 , 1.0/np.sqrt(2.0)   # Normal unit vector to SEl mirror
mSEl = LVMu.MirMat(nx,ny,nz)                           # Calculate mirror matrix

#---- Calculate the SEl and SAz siderostat angles (Section 4.3 of doc)

sazD1,selD1 = LVMu.AzEl2SazSel(az1,el1)        # Calculate siderostat angles before
sazD2,selD2 = LVMu.AzEl2SazSel(az2,el2)        #    and after transit

print "\nSAz,SEl before transit ",sazD1,selD1
print "SAz,SEl after  transit ",sazD2,selD2,"\n"

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

#---- And now, propagate the unit-vectors

nV1 = np.matmul(mEff1,nV)     # These are the on-sky North-pointing vectors in LVM coordinates before
nV2 = np.matmul(mEff2,nV)     #    and after transit

eV1 = np.matmul(mEff1,eV)     # These are the on-sky East-pointing vectors in LVM coordinates before
eV2 = np.matmul(mEff2,eV)     #    and after transit


print "Projected in XZ plane\n"
print "  N before transit        ",nV1[0],nV1[2]
print "  E before transit        ",eV1[0],eV1[2],"\n"
print "  N after transit         ",nV2[0],nV2[2]
print "  E after transit         ",eV2[0],eV2[2],"\n"

