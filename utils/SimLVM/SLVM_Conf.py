
'''
    Contains parameters and a few classes and functions for SatPlot

    Classes and Functions

       MyQGS          - (class) QGraphicsScene override with mouse click capability
       GetFITSthread  - (class) Asynchronous download thread for FITS files
       ParseNewFld    - (function) Parses new field info and returns standard values
       suggestFldName - (function) Generates a field name based on supplied HH MM SS.S +DD MM SS
       overlayPolys   - (function) Returns 3 QPolygons corresponding to IFU, AG fields
       genPixMap      - (function) Generates a PixMap from the supplie image (w/ autoscale, invert)
       LoadListMD     - (class) Modal dialog to select a target from a list

'''

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.uic import loadUi    # To load .ui files

import os, datetime
import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord              # For dealing with -00 24 43.4 etc.

#-----------------------------------------------------------------------------#    

# --- Parameters ---

verbose = True   # How wordy to be

wSize = 800      # Size of QGraphicsView onscreen (will change, based on .ui file)
fSize = 950      # Default size of FITS file to download (pixels)
fDsize = 1.5 #/2.0 # Defaults size of field to download (degrees)

fldName = ""     # Current Field name
fldType = ""     # Field type - "RADEC" or "NAMED"
fldTxt  = ""     # Text of field location or name

#--- Vizier search parameters...valid servers are:
#---    'vizier.u-strasbg.fr','vizier.nao.ac.jp', 'vizier.hia.nrc.ca',
#---    'vizier.ast.cam.ac.uk','vizier.cfa.harvard.edu', 'www.ukirt.jach.hawaii.edu',
#---    'vizier.iucaa.ernet.in','vizier.china-vo.org'


vizServer = 'vizier.ast.cam.ac.uk'     # Vizier server. See above for options
survey = "Gaia-DR2"                    #  was "USNO-B1"...Which VIZIER survey to use for reference stars
refLim = 10                            # Magnitude limit for returned Vizier reference stars

#workDir = "/Users/herbst/Dropbox/LVM/SimLVM/FITS"
workDir = "FITS"

#--- Image display parameters

pAuto = True         # Pixmap autoscale
pB = 0.0             # Pixmap black value
pW = 10000.0         # Pixmap white value
pInvert = True       # Whether to invert Pixmap

pAuto_AG = False     # Autoscale AG cams
pB_AG = 0.0          # Pixmap black for AG cams
pW_AG = 150.0        #    and white
pInvert_AG = True    # Invert AG Cam grayscale

useOnAxis = True     # Whether the on-axis camera is installed

#--- GUI parameters

refRad = 10      # Radius/half-width of circle/box to draw around reference stars

#-----------------------------------------------------------------------------#    

class MyQGS (QtWidgets.QGraphicsScene):
    """ Override of QGraphicsScene to capture mouse down events
        (originally developed for All-Sky-Cam ProcViewer.py)

        Updated for Py3Qt5 - 1 May 2020
    """

    from PyQt5.QtCore import pyqtSignal

    sendMouse = pyqtSignal(str)             # Need to define this in pyqt5

    def __init__(self,parent = None):
        super(MyQGS,self).__init__(parent)

    def mousePressEvent(self, event):
        super(MyQGS, self).mousePressEvent(event)

        if event.button() == QtCore.Qt.LeftButton:
            thePos = QtCore.QPoint(event.scenePos().x(), event.scenePos().y())
            print (str(thePos.x()),str(thePos.y()))
            self.sendMouse.emit(str(thePos.x()) + " " + str(thePos.y()))

#-----------------------------------------------------------------------------# 

class GetFITSthread(QtCore.QThread):

    '''
        Download a FITS file for the specified location. Implemented as a QThread
        so that the GUI doesn't lock up.

        Inputs  loc     - Location, either RA,dec as string or an object name
                type    - Either "radec" or "name"

        Outputs fNam    - FITS file. Does not include workDir full path!

        Possible surveys are:
            R band:     'SDSSr' 'DSS1 Red' 'DSS2 Red'
            NIR         '2MASS-K' '2MASS-H' '2MASS-J'

            see http://astroquery.readthedocs.io/en/latest/skyview/skyview.html?highlight=dss for more.

        Updated for Py3Qt5 - 1 May 2020

    '''

    msg = QtCore.pyqtSignal(str)          # To send messages back
    finished = QtCore.pyqtSignal(str)     # Sent at completion

    def __init__(self,namP,locP,typeP):
        QtCore.QThread.__init__(self)

        self.nam=namP      # Field name (and FITS file root)
        self.loc=locP      # Grab location
        self.type=typeP    #  and type

        if verbose:
            self.msg.emit("GetFITSthread launched with "+self.nam+" "+self.loc+" "+self.type)

    #----------------

    def __del__(self):     # Presumably kills thread?
        self.wait()

    #----------------

    def run(self):                       # This gets called at thread.start()

        from astroquery.skyview import SkyView

        self.startT=datetime.datetime.now()                        # Set initial "last" picture time

        if verbose:
            self.msg.emit("Initiating astroquery for "+self.loc+" at "+str(self.startT))

        theFITS = SkyView.get_images(position=SkyCoord(self.loc,unit=(u.hourangle,u.deg)), survey=['DSS2 Red'], radius=fDsize*u.deg, pixels=fSize)    # NEW - Sep19

        if len(theFITS):    # Non-zero result returned - Some places on-sky fail, like 00 07 3.082 +65 38 35.50

            fNam = self.nam+".fits"                             # Assemble FITS file name

            theFITS[0].writeto(workDir+"/"+fNam,overwrite="True")   # Write FITS file to working directory

        else:

            fNam = "NoData"   # Flags the problem

            self.msg.emit("Zero-length data returned...bad location?\a\a")


        self.finished.emit(fNam)     # Send back the FITS file name
        self.exit()                  # All done...

#-----------------------------------------------------------------------------#    

def ParseNewFld():
    '''
        Based on ParseNewCat() from ChartMaker

        User has entered the parameters for a new field. Work out central
        location and other info. If named, call Simbad to get coordinates.

        Specifically, we have:
            fldType - Either "RADEC" or "NAMED"
            fldTxt  - The string of RA, Dec or the object name

        This routine uses astroquery.simbad

        Note also that Simbad occasionally returns coordinates like this:

           18 01 06 +02 54.1   instead of the usual:
           18 01 06 +02 54 06

        In this instance (i.e. Dec MM.M), we have to reassemble the standard format.

    '''

    from astroquery.simbad import Simbad

    if fldType == "RADEC":                                       # Simple central coordinates
        ctrStr = fldTxt                                          # Coordinates are already there
        c = SkyCoord(ctrStr,unit=(u.hourangle,u.deg))                # Convert RA Dec to SkyCoord

    else:                                                            # Object name...go to Simbad for coordinates
        cS = Simbad()                                                # Create custom Simbad
        cS.add_votable_fields('pmra','pmdec','flux(R)','flux(B)')    # Add proper motion, filter mags (NOT IMPLEMENTED YET)
        r=cS.query_object(fldTxt)                                # Send the query to Simbad

        print("SIMBAD", r)

        simCtrStr = str(r[0]['RA'])+" "+str(r[0]['DEC'])             # This extracts RA and dec as HH MM SS.SS +DD MM SS.S strings

        c = SkyCoord(simCtrStr,unit=(u.hourangle,u.deg))             # Convert RA Dec to SkyCoord
        v  = c.to_string('hmsdms',sep=" ",precision=3).split()       # Split up to build standard format 

        ctrStr = v[0]+" "+v[1]+" "+("%6.3f" % float(v[2]))+"  "+v[3]+" "+v[4]+" "+("%5.2f" % float(v[5]))    # Our standard HH MM SS.SSS +DD MM SS.SS

    fldName = suggestFldName(ctrStr)                      # Generate field name based on center loc

    ctrPos=np.array([c.ra.degree,c.dec.degree])           # And into decimal degrees (RA,dec)

    return fldName,ctrPos,ctrStr          # Return values

    #-----------------------------------------------------------------------------#    

#-----------------------------------------------------------------------------#       

def suggestFldName(loc):

    '''
        Based on suggestCatName of ChartMaker, but with only location strings,
        not named objects. Also, no ".cat" extension

        Suggest a catalog name, based on the supplied location. This is
        a (HH MM SS.SS +DD MM SS.S) location.

        The current (14Mar18) standard is HHMM.m+DDMM.cat (or "-" for 
        negative decs...This means a unique catalog name for each 1.5
        arcmin (RA) at the equator by 1.0 arcmin (dec). The names will
        be unique for correspondingly smaller sky areas toward the poles.

        Note that the science IFU is ~30 arcmin across.
    '''

    locS = loc.split()                      # Grab and split HH MM SS DD MM SS


    RAstr = locS[0]+locS[1]+("%4.1f" % (float(locS[2])/60))[-2:]  # Gives HHMM.m

    decStr = locS[3] + locS[4]       # Gives DDMM

    sugNam = RAstr + decStr          # Put it all together

    return sugNam

#-----------------------------------------------------------------------------#

def overlayPolys(rAn):

    '''
        Calculates QPolygons corresponding to the IFU location and the two
        Acquisition and Guider chips.

          Input:
            rAn     - Rotation angle of the focal plane
          Output:
            ifuPoly - QPolygon for Hexagon corresponding to IFU
            ag1Poly - QPolygon for AG sensor 1
            ag2Poly - QPolygon for AG sensor 2

        The QPolygons are calculated parametrically, based on the sizes and
        separations of the various fields. They should therefore remain ok
        if we change the number of image pixels, etc. Note also that we have
        to work with integers when using QPolygons.

    '''

    #--- Center point, pixels per arcsec, etc.

    


    # points = QPolygon([
    #             QPoint(10,10),
    #             QPoint(10,100),
    #             QPoint(100,10),
    #             QPoint(100,100)       

#-----------------------------------------------------------------------------#

def genPixMap(ary,wid,z1,z2,auto,invert):
    """
       Generates and returns a PixMap suitable for onscreen display, that is, of size wid
       and scaled from 0-255.

       Parameters:  ary    - numpy array of input data (any size, but must be SQUARE!)
                    wid    - width of PixMap to return
                    z1,z2  - grayscale cuts
                    auto   - whether to autoscale
                    invert - whether to invert colour map

       Returns:     pMp    - pixmap suitable for onscreen display
                    z1,z2  - grayscale cuts (changes if auto=True)
                    ok     - flag indicating success
    """

    import scipy.ndimage.interpolation     # For easy scaling
    from PyQt5.QtGui import QImage,qRgb

    ok=False       # Assume the worst

    ary = np.flipud(ary)                   # Need to y-flip Pixmap to match PCam GUI

    pMp=QtGui.QPixmap("sadFace.png")       # Sad face image

    if auto:       # Autoscale
        pMin = np.amin(ary)       # Min value
        pMax = np.amax(ary)       # Max value
    else:
        pMin = float(z1)          # User supplied 
        pMax = float(z2)

    sclFactor = 1.0*wid/ary.shape[0]                                          # Need to zoom by this factor (plus buffer)
    pS = scipy.ndimage.interpolation.zoom(ary,sclFactor,order=3)              # Zoom to fit

    # pSG = ((pS-pMin)*255.0 / (pMax-pMin)).astype(np.uint8)
    # pSG[pSG < 0] = 0
    # pSG[pSG > 255] = 255

    # https://stackoverflow.com/questions/46689428/convert-np-array-of-type-float64-to-type-uint8-scaling-values

    # print (np.finfo(pS.dtype))

    # max_value = np.finfo(pS.dtype).max       # Get the information of the incoming image type
    # pS = pS.astype(np.float64) / max_value   # normalize the data to 0 - 1
    # pS = 255 * pS # Now scale by 255

    pSS = (pS - pMin) / (pMax - pMin) * 255.0  # Scaled to 0-255
    pSS[pSS < 0] = 0
    pSS[pSS > 255] = 255

    pSG = pSS.astype(np.uint8)

    # pSG = np.clip( ((pS-pMin)*255.0 / (pMax-pMin)).astype(np.uint8), 0, 255)       #   and scale 0-255
    # pSG = np.clip( ((pS-pMin)*255.0 / (pMax-pMin)), 0, 255)       #   and scale 0-255

    import matplotlib.pyplot as plt
    plt.imsave('pSG.png', pSG,vmin=0.0,vmax=255.0)
    # print ("pMin/pMax: ",pMin,pMax)
    # print ("pSS min/max: ",np.amin(pSS),np.amax(pSS))
    # print ("pSG min/max: ",np.amin(pSG),np.amax(pSG))

    if invert:    # Reverse graymap
        gray_color_table = [qRgb(255-i, 255-i, 255-i) for i in range(256)]
    else:
        gray_color_table = [qRgb(i, i, i) for i in range(256)]

    qImg = QImage(pSG.data, pSG.shape[1], pSG.shape[0], pSG.strides[0], QImage.Format_Indexed8)
    qImg.setColorTable(gray_color_table)

    pMp = QtGui.QPixmap.fromImage(qImg)                               #   and to pixmap
    ok=True                                                           # Success!

    return pSS,pMp,pMin,pMax,ok    # Send back scaled numpy array, pixmap, clip limits, ok flag

#-----------------------------------------------------------------------------#
#
# https://stackoverflow.com/questions/45020672/convert-pyqt5-qpixmap-to-numpy-ndarray

def QPixmapToArray(pixmap):

    size = pixmap.size()      ## Get the size of the current pixmap
    h = size.width()
    w = size.height()

    
    qimg = pixmap.toImage()  ## Get the QImage Item and convert it to a byte string
    byte_str = qimg.bits().tobytes()

    ## Using the np.frombuffer function to convert the byte string into an np array
    img = np.frombuffer(byte_str, dtype=np.uint8).reshape((w,h,4))

    return img

#-----------------------------------------------------------------------------#

def NaNscale255(ary,wid,z1,z2,auto,invert,doQP):

    '''
        NaNscale255 - Scale input numpy array 0-255 integer, remove NaN's

            A modified version of genPixMap above...

            This routine accepts an input numpy array and replaces all instances
            of NaN with zero before scaling the result 0-255 (uint8) and returning.
            It also returns a QPixmap, if doQP is True.

                Inputs:  ary     - numpy array of data (any size, but must be SQUARE!)
                         doQP    - whether or not to create and return a QPixmap
                         wid     - width of PixMap to return
                         z1,z2   - grayscale cuts
                         auto    - whether to autoscale
                         invert  - whether to invert colour map

                Returns: PSG     - uint8 array of scaled input data
                         pMp     - pixmap suitable for onscreen display
                         pMn,pMx - min / max values used for scaling
                         ok      - flag indicating success

                Note: If auto is False, (pMn,pMx) are equal to (z1,z2), otherwise, the
                      min/max values of the input array are returned.

                NOTE: If using to just create numpy arrays, scaled 0-255 (to Assemble
                      composite DSS / AG pixmaps, for example), set doQP to False. This
                      will also suppress the scipy interpolation !!!
    '''

    import scipy.ndimage.interpolation     # For easy scaling
    from PyQt5.QtGui import QImage,qRgb
    import matplotlib.pyplot as plt        # To create debug PNG

    ok=False       # Assume the worst

    ary = np.flipud(ary)                   # Need to y-flip Pixmap to match

    if auto:                  # Autoscale
        pMin = np.amin(ary)       # Min value of input array
        pMax = np.amax(ary)       # Max value
    else:
        pMin = float(z1)          # User supplied cuts
        pMax = float(z2)

    #--- Replace any NaN values with 0.0 and scale to desired wid (if generating a pixmap)

    ary[np.isnan(ary)] = 0.0

    if doQP:       # Do scaling for final pixmap

        sclFactor = 1.0*wid/ary.shape[0]                              # Need to zoom by this factor (plus buffer)
        pS = scipy.ndimage.interpolation.zoom(ary,sclFactor,order=3)  # Zoom to fit

    else:

        pS = ary   # No scaling, since we are just returning a 0-255 numpy array of the input

    #--- Now scale the result to 0-255 (clipping at limits) and convert to uint8

    pSS = (pS - pMin) / (pMax - pMin) * 255.0    # Scaled to 0-255
    pSS[pSS < 0] = 0                             # Ensure no negative values
    pSS[pSS > 255] = 255                         #    or values > 255

    pSG = pSS.astype(np.uint8)                      # Convert to uint8

    plt.imsave('pSG.png', pSG,vmin=0.0,vmax=255.0)  # Save an image for debugging

    #--- And create QPixmap if desired

    if doQP:

        if invert:    # Reverse graymap
            gray_color_table = [qRgb(255-i, 255-i, 255-i) for i in range(256)]
        else:
            gray_color_table = [qRgb(i, i, i) for i in range(256)]

        qImg = QImage(pSG.data, pSG.shape[1], pSG.shape[0], pSG.strides[0], QImage.Format_Indexed8)
        qImg.setColorTable(gray_color_table)

        pMp = QtGui.QPixmap.fromImage(qImg)                   #   and to pixmap

    else:

        pMp = None      # Set QPixmap to None

    ok=True                                               # Success!

    # print ("pMin/pMax: ",pMin,pMax)
    # print ("pSS min/max: ",np.amin(pSS),np.amax(pSS))
    # print ("pSG min/max: ",np.amin(pSG),np.amax(pSG))

    return pSG,pMp,pMin,pMax,ok    # Send back scaled numpy array, pixmap, clip limits, ok flag


#-----------------------------------------------------------------------------#    

class LoadListMD(QtWidgets.QDialog):

    ''' 
       User has selected "Load List". Pop a dialog for the list file name
       and then present the targets.

       The user selects one (or cancels).
 

       Note that this code was built on Tom's template:

            self_mdDialog.py and self_mdApp.py in ~/Dropbox/Python/ModalDialog/Python_3.7_PyqQt5

    '''

    def __init__(self):    #,parent=None

        QtWidgets.QDialog.__init__(self)
        
        self.initVar()          # Initialize important variables
        self.initUI()           # Initialize this window
        self.modalUI.show()     # Show the user interface

        self.doGetFile()        # Pop the dialog for list file
                
    def initVar(self):
        
        self.gotIt = False      # Goes True when all is Ok.

    def initUI(self):               
        """ Initializes the GUI, especially attaching widgets to routines """

        self.modalUI = loadUi("LoadList_Modal.ui",self)             # Load GUI XML file into ui

        self.modalUI.Ok_Cancel_bb.accepted.connect(self.doOk)            # User pressed Ok...proceed
        self.modalUI.Ok_Cancel_bb.rejected.connect(self.doCancel)        # User cancelled

    def doGetFile(self):

        # loadName = str(QtWidgets.QFileDialog.getOpenFileName(self, 'Select List File', '','Text files (*.txt)',options=QtWidgets.QFileDialog.DontUseNativeDialog))
        # loadName, _ = str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', '/home/herbst/Dropbox/','Text files (*.txt)'))
        loadName, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select Catalog File", "", "Catalog files (*.txt)",options=QtWidgets.QFileDialog.DontUseNativeDialog)

        if (loadName==""):         # User Cancelled
            print ("Cancelled")
            self.doCancel()
        else:                      # Proceed with file read and populate list
            print ("Opening ",loadName)
            self.inList = open(loadName, "r").read().splitlines()     # Now a list

            for i in range(1,len(self.inList)):                       # Step through (ignore header)
                self.modalUI.targets_lw.addItem(self.inList[i])       #   adding targets to list widget


    def doOk(self):

        self.theTarget = self.inList[self.modalUI.targets_lw.currentRow() + 1]     # We picked this one

        self.gotIt = True    # We have what we need
        self.accept()        # Ok...tell calling program

    def doCancel(self):

        self.theTarget = ""  # We picked Nothing
        self.reject()        # Cancel...tell calling program


