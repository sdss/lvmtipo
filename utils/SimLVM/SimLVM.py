#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
    SimLVM - Simulate LVM Focal Plane

    This program, based on LN ChartMaker, presents the LVM focal plane,
    including the Science IFU and the two Acquisition and Guiding (AG)
    sensors.

    It goes to Vizier to identify and mark potential guide stars.

    A Word about Image Sizes

    The LVM focal plane is 44.5 mm in diameter, which is 4989 arcsec (1.39 deg)
    at the nominal image scale of 8.92 microns/arcsec. This is the
    diameter of the escribed circle that just touches the most distant
    corners of the AG chips. We choose to work with somewhat larger images
    that are 1.5 degrees or 5400 arcsec square. To keep image sizes manageable
    (naturally at the expense of resolution), we download 900x900 pixel FITS
    files. These will then be 6 arcsec per pixel.

    The Sony IMX432 sensors have 1608 x 1104 pixels, each of which is 9 microns
    square. This corresponds to 14.472 x 9.936 mm and 1622 x 1114 arcsec
    (27.04 x 18.57 arcmin). For the SimLVM images (at 6 arcsec/pixel) the guide
    chips are 270 x 186 onscreen pixels. We will see if this is good enough...

    Note that the file "overlay.png" contains the overlay exposing the IFU and
    AG fields. It is also 900 x 900.

"""

# ####  Import Needed Packages

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.uic import loadUi

import sys,os,subprocess
import numpy as np

from astropy import wcs                     # To deal with WCS
from astropy import units as u
from astropy.coordinates import SkyCoord    # For dealing with -00 24 43.4 etc.

import SLVM_Conf as SLC
from SimLvmScraper import QtLvmScraper, CommandStatus, flatten, unpack

class AppWin(QtWidgets.QMainWindow):

    """ 
        Class for Main application window. 
    """
    
    def __init__(self):                         # Initialiaze self
        super(AppWin, self).__init__()          # Initialize inherited class
#        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)   # Set to delete on close

        self.initVar()                          # Initialize important variables
        self.initUI()                           # Initialize this window

        self.ui.hide()        # Added Feb23 (fixed missing cursor at startup)
        self.ui.show()        # Show the user interface
    
    #-----------------------------------------------------------------------------#    
    
    def initVar(self):
        """ Initializes important variables (i.e. current stage positions, limits)"""
        
        self.webT = None         # No web search thread yet!
        self.gotViz = False      # No vizier references yet!
        self.showViz = True      # Assume user wants to see Vizier stars (when there)
        self.whatImg = "DSS"     # Initially assume DSS only display
        self.gotAG = False       # Whether we have acquired AG images
        self.doBlinkDSS = False  # Tick-Tock flag to blink between DSS and AG

        self.zoomInFactor = 1.25                    # Scale factor for zooming in
        self.zoomOutFactor = 1/self.zoomInFactor    #   and out
        self.curIzoom = 0.8                         # Need to keep track of zoom factor

        self.ovPM = QtGui.QPixmap("overlay.png")    # Semi-transparent overlay

        self.WCSkeys = None      # To prevent Annotate errors before FITS download complete

    #-----------------------------------------------------------------------------#    

    def initUI(self):               
        """ Initializes the GUI, especially attaching widgets to routines """

        self.ui = loadUi("SimLVM.ui",self)     # Load GUI XML file into ui

        #--- Push Buttons and Check Boxes ---

        self.ui.dummy_AG_pb.clicked.connect(self.dummy_AG)          # CHANGE! - dummy AG Frame trigger


        self.ui.quit_pb.clicked.connect(self.shutdown)              # Quit button

        self.ui.loadList_pb.clicked.connect(self.loadList)          # Load a list of target coordinates
        self.ui.fetch_pb.clicked.connect(self.loadField)            # Fetch new field button
        self.ui.launchCDS_pb.clicked.connect(self.launchCDS)        # Open CDS with current target

        self.ui.searchViz_pb.clicked.connect(self.searchViz)        # Search VIZIER for Ref Stars
        self.ui.searchViz_pb.setEnabled(False)                      # Disabled without a field (yet)

        self.ui.showViz_cb.toggled.connect(self.togShowVizier)      # Show Vizier on chart check box
        self.ui.showViz_cb.setChecked(True)                         # Visibility toggled on

        self.ui.Zin_pb.clicked.connect(self.imgZin)                 # Zoom buttons
        self.ui.Zout_pb.clicked.connect(self.imgZout)

        self.ui.Autoscale_cb.toggled.connect(self.togAutoscale)         # Autoscale check box (DSS)
        self.ui.Invert_cb.toggled.connect(self.togInvert)               #  and Invert check box

        self.ui.Autoscale_AG_cb.toggled.connect(self.togAutoscale_AG)   # Autoscale check box (AGC)
        self.ui.Invert_AG_cb.toggled.connect(self.togInvert_AG)         #  and Invert check box

        self.ui.Blink_pb.clicked.connect(self.blinkImg)             # Blink between DSS and Both mode

        self.ui.DSS_rb.clicked.connect(self.doShowRadioButton)
        self.ui.AG_rb.clicked.connect(self.doShowRadioButton)
        self.ui.Both_rb.clicked.connect(self.doShowRadioButton)
        self.ui.None_rb.clicked.connect(self.doShowRadioButton)

        #--- Line Edits ---

        self.ui.locNam_le.returnPressed.connect(self.loadField)     # White background...ready to go

        self.ui.vizLim_le.setValidator(QtGui.QDoubleValidator())    # Only floats for Vizier limit
        self.ui.vizLim_le.returnPressed.connect(self.setVizLim)     # Connect

        self.ui.Bk_le.setValidator(QtGui.QDoubleValidator())        # DSS Cuts Only floats for cuts
        self.ui.Bk_le.returnPressed.connect(self.setBk)             # Connect
        self.ui.Wh_le.setValidator(QtGui.QDoubleValidator())        # Only floats for cuts
        self.ui.Wh_le.returnPressed.connect(self.setWh)             # Connect

        self.ui.Bk_AG_le.setValidator(QtGui.QDoubleValidator())     # AG Cuts Only floats for cuts
        self.ui.Bk_AG_le.returnPressed.connect(self.setBk_AG)       # Connect
        self.ui.Wh_AG_le.setValidator(QtGui.QDoubleValidator())     # Only floats for cuts
        self.ui.Wh_AG_le.returnPressed.connect(self.setWh_AG)       # Connect

        #--- List Widgets ---

        self.ui.vizier_lw.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)   # Select max 1 item
        self.ui.vizier_lw.clicked.connect(self.vizLWclick)                            # Process user click

        self.imgScn=SLC.MyQGS(self)                              # Set up Qgraphics scene

        self.imgScn.sendMouse.connect(self.MousePress)           # Connect mouse press

        SLC.wSize = self.ui.img_qgv.width()                      # Size of QGraphicsView

        imgMap=QtGui.QPixmap("imgSplash.png")                    # Eye candy to start
        self.imgScn.addPixmap(imgMap) #.scaled(SPC.wSize*0.9,SPC.wSize*0.9,QtCore.Qt.KeepAspectRatio))

        self.ui.img_qgv.setScene(self.imgScn)                       # Set scene to QGraphicsView

        self.ui.statusBar().showMessage('Ready')           # Show status bar

        self.UpdateGUI()      # Stuff initial parameters

    #-----------------------------------------------------------------------------#    

    def loadList(self):

        '''
            User has pressed "Load List". Present dialog to select file and then allow user
            to choose a target from the list.

            Written in the field (LVM telescope enclosure, 13 Feb 2023)
        '''

        self.loadList_modal = SLC.LoadListMD()    # Pops modal for target seletion

        self.loadList_modal.exec_()               # Run it

        if self.loadList_modal.gotIt:         # User did not press Ok with valid input

            target = self.loadList_modal.theTarget         # Text line containing field
            self.fieldName = target.split()[0]             # Target name
            RD = target.split()
            RA_Dec_Txt = RD[1]+" "+RD[2]+" "+RD[3]+" "+RD[4]+" "+RD[5]+" "+RD[6]   # Assemble name
            self.ui.fieldName_lbl.setText(self.fieldName)
            self.ui.locNam_le.setText(RA_Dec_Txt)

    #-----------------------------------------------------------------------------#    

    def loadField(self):

        """ User pressed the Load button. We assume that there is a valid target
            field defined!

        """

        if (self.ui.locNam_le.text()=="") :              # No entry - warn user and return
            reply = QtWidgets.QMessageBox.critical(self, 'Warning',"Please enter both a catalog name and a location/object name",QtWidgets.QMessageBox.Ok)
    
            return

        locS = str(self.ui.locNam_le.text()).split()   # Grab and split location or object name

        if (len(locS)==6):                             # We have an RA Dec target

            rah,ram,ras,ded,dem,des = locS             # And we have already split them up

            cc = SkyCoord(rah+" "+ram+" "+ras,ded+" "+dem+" "+des,unit=(u.hourangle,u.deg))
            v  = cc.to_string('hmsdms',sep=" ",precision=3).split()

            theTxt = v[0]+" "+v[1]+" "+("%6.3f" % float(v[2]))+"  "+v[3]+" "+v[4]+" "+("%5.2f" % float(v[5]))

            SLC.fldType = "RADEC"           # User entered RA and Dec...Format nicely
            SLC.fldTxt = theTxt             #   and add to "globals"

        else:                                                # Named target
            SLC.fldType = "NAMED"                            # Record type of target
            SLC.fldTxt = str(self.ui.locNam_le.text())       #   and its name

        self.fieldName,self.ctrPos,self.ctrStr = SLC.ParseNewFld()      # Get field name & center (numbers, string)

            #--- Fetch the FITS file. 

        self.fetchFITS()       # Check if already on disk and if not, launch thread to get FITS

        self.searchViz()       # And as a convenience, request Vizier sources


#-----------------------------------------------------------------------------#    

    def fetchFITS(self):

        '''
            Checks if already downloading a FITS file.
            If not, launch GetFITSthread to go online and grab the FITS file. The
            thread itself launches "GotFITS" if successful.

        '''

        import sys

        if (os.path.isfile('FITS/'+self.fieldName+'.fits')):   # FITS file already there

            self.ui.statusBar().showMessage('Found local FITS file')        # Inform user
            self.GotFITS(self.fieldName+'.fits')                            # Load up and display

        self.ui.statusBar().showMessage('Fetching FITS file')           # Inform user
        self.setStyleSheet("QMainWindow{background:rgb(255,225,225)}")  # Set window pink

        if self.webT is None:                                           # First time for thread or has been killed

            self.webT=SLC.GetFITSthread(self.fieldName,self.ctrStr,'radec')            # Launch thread to get FITS file

            self.webT.msg.connect(self.thrdMsg)                         # (Intermediate) Message
            self.webT.finished.connect(self.GotFITS)                    # Connect to finished routine
            
            self.webT.start()                                           # Start web search


#-----------------------------------------------------------------------------#

    def GotFITS(self,fName):

        '''
            A FITS file is in place. Load up and display it. Note
            that this routine can be called directly if the FITS file
            exists on disk. It can also be called by the web thread
            which goes online to get the file (and places it on disk).

            Note that some locations return zero data (no survey there?).
            We catch this with the fName="NoData".
        '''

        from astropy.io import fits       # To read and write FITS files

        if (fName=="NoData"):             # Ooops! a problem

            if SLC.verbose:
                self.logTxt("Zero data returned...bad location?",True)         # Log the problem

            self.setStyleSheet("QMainWindow{background:rgb(255,255,255)}")     # Back to white

        else:

            if SLC.verbose:
                self.logTxt("GotFITS " + SLC.workDir+fName,True)       # Log that we have it

            self.ui.statusBar().showMessage("Got FITS file " + fName)
            self.setStyleSheet("QMainWindow{background:rgb(255,255,255)}")  # Back to white

            self.fitsName = str(SLC.workDir+"/"+fName)                      # Full path

            hdulist = fits.open(self.fitsName)   # Read in FITS file
            self.image = hdulist[0].data         # Extract image data

            print ("FRESH self.image max, median :",np.median(self.image),np.median(self.image))

            self.WCSkeys = wcs.WCS(hdulist[0].header)    # Pull out WCS keywords

            self.NAXIS1,self.NAXIS2 = hdulist[0].header['NAXIS1'],hdulist[0].header['NAXIS2']    # Number of pixels
            self.CRVAL1,self.CRVAL2 = hdulist[0].header['CRVAL1'],hdulist[0].header['CRVAL2']    # Central pixel long and lat
            self.CDELT1,self.CDELT2 = hdulist[0].header['CDELT1'],hdulist[0].header['CDELT2']    # Pixel scale (degrees)
            self.CRPIX1,self.CRPIX2 = hdulist[0].header['CRPIX1'],hdulist[0].header['CRPIX2']    # Reference pixel

            if (abs(self.CDELT1)!=abs(self.CDELT2)):       # Oooops!
                self.logTxt('WARNING - Image pixels not square',False)
                self.logTxt('   CDELT1 = '+str(self.CDELT1)+"  CDELT2 = "+str(self.CDELT2),True)

            self.pxScale = abs(self.CDELT1) * 3600.0  # Pixel scale (arcsec/px) NOTE: CDELT1 can be negative!

            hdulist.close()            # Close input file

            self.GenPixMap("DSS")      # Create PixMap of DSS data

            self.DispImg()             # And display it

            self.ui.searchViz_pb.setEnabled(True)    # User can now search VIZIER

        if self.webT != None:                                      # Thread still around somehow
            # self.logTxt("Thread already exists, killing",True)     # Log that we will kill it
            self.webT = None                                       # Zorch zombie!

#-----------------------------------------------------------------------------#    

    def DispImg(self):

        '''
            More or less cloned from ChartMaker.
            Generates a Pixmap from self.image (current chart image data)
            and loads it into the QGraphicsView.

            13 Feb 23 - Moved much of the functionality to GenPixMap and RedrawChart
        '''

        # try:           # Don't do anything if no data!
        ##### UNCOMMENT THE TRY EXCEPT BLOCK AND INDENT WHEN WORKING 

        self.ui.imgScn.clear()           # Presumably, this will remove all previous items?

        # Fwid = SLC.fSize                 # Size of FITS file
        # Wwid = SLC.wSize                 #   and of QGraphicsView
        # zerOff = -1 * (Fwid - Wwid)/2.0  # Need to offset pixmaps this much to center

        #--- Generate PixMap of DSS image and / or AG Image ---
        #---   using NaNscale255(ary,wid,z1,z2,auto,invert,doQP)

        # ok = True        # Assume the best for now

        # print ("")

        # if self.whatImg=="DSS":          # DSS only
        #     print ("DSS only")

        #     self.pMap = self.GenPixMap("DSS")

            # pImg = np.copy(self.image)       # Set pixel map image to DSS

            # if (SLC.pAuto==True):            # Autoscale
            #     gDSS,self.pMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,0.0,0.0,True,SLC.pInvert,True)    # Generate Autoscaled onscreen pixmap
            #     SLC.pB,SLC.pW = pMin,pMax                                                           # Store new values
            # else:
            #     gDSS,self.pMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,SLC.pB,SLC.pW,False,SLC.pInvert,True)   # Or fixed scaling


        # elif self.whatImg=="AG":            # AG only
        #     print("AG only")
        #     pImg = np.copy(self.image_AG)       # Set pixel map image to AG

        #     if (SLC.pAuto_AG==True):            # Autoscale
        #         gAG,self.pMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,0.0,0.0,True,SLC.pInvert_AG,True)                # Generate Autoscaled onscreen pixmap
        #         SLC.pB_AG,SLC.pW_AG = pMin,pMax                                                                         # Store new values
        #     else:
        #         gAG,self.pMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,SLC.pB_AG,SLC.pW_AG,False,SLC.pInvert_AG,True)   # Or fixed scaling

        # elif self.whatImg=="Both":
        #     print ("Both")

        #     ###--- Create DSS numpy array with holes, scaled 0-255
            
        #     pImg = np.copy(self.image) * self.mask_AG       # Set to DSS and zero AG pixel locations

        #     if (SLC.pAuto==True):            # Autoscale
        #         compDSS,self.pMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,0.0,0.0,True,SLC.pInvert,False)          # Generate Autoscaled onscreen pixmap
        #         SLC.pB,SLC.pW = pMin,pMax                                                                         # Store new values
        #     else:
        #         compDSS,self.pMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,SLC.pB,SLC.pW,False,SLC.pInvert,False)   # Fixed scaling

        #     print ("HERE ")

        #     ###-- Now create AG numpy array with zeroes elsewhere

        #     pImg = np.copy(self.image_AG)

        #     if (SLC.pAuto_AG==True):            # Autoscale
        #         compAG,self.pMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,0.0,0.0,True,SLC.pInvert_AG,False)                # Generate Autoscaled onscreen pixmap
        #         SLC.pB_AG,SLC.pW_AG = pMin,pMax                                                                           # Store new values
        #     else:
        #         compAG,self.pMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,SLC.pB_AG,SLC.pW_AG,False,SLC.pInvert_AG,False)   # Or fixed scaling

        #     ###-- And finally add them together to generate the final pixmap (note scaling 0-255)

        #     pImg = np.flipud(compDSS + compAG)
        #     # print ("composite max / median  ",np.amax(pImg),np.median(pImg))

        #     comp,self.pMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,0.0,255.0,False,SLC.pInvert,True)    # Generate Autoscaled onscreen pixmap

        # else:
        #     print ("None")

        # if ok:      # Successful (or no) Pixmap conversions...update things

        self.RedrawChart()      # Puts up PixMap, Overlay, etc.

        self.ui.fieldName_lbl.setText(self.fieldName)

        self.AnnotateChart()  # Annotate the chart

        self.UpdateGUI()      # To update any changes        

#-----------------------------------------------------------------------------#

    def GenPixMap(self,which):

        '''
            This routine was added on 13 Feb 23 to speed up blinking etc. It
            basically creates the following entities:

                self.DSSpMap  - Pixmap of current DSS image
                self.AGpMap   - Pixmap of current AG image(s)
                self.BothPMap - Pixmap of both DSS and AG overlaid

            The pixmaps depend on the current Autoscale / Cuts, and hence this
            routine should be called when new data appears or when the user changes
            settings.

            Inputs:     which  - Flag which is either "DSS", "AG", or "Both"

            Returns:    Nothing (it creates self.DSSpMap, etc.)

            Note that much of this code was in DispImg before...

        '''

        Fwid = SLC.fSize                 # Size of FITS file
        Wwid = SLC.wSize                 #   and of QGraphicsView
        zerOff = -1 * (Fwid - Wwid)/2.0  # Need to offset pixmaps this much to center

        if which=="DSS":                 # Make DSS Pixmap

            pImg = np.copy(self.image)       # Set pixel map image to DSS

            if (SLC.pAuto==True):            # Autoscale
                gDSS,self.DSSpMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,0.0,0.0,True,SLC.pInvert,True)    # Generate Autoscaled onscreen pixmap
                SLC.pB,SLC.pW = pMin,pMax                                                           # Store new values
            else:
                gDSS,self.DSSpMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,SLC.pB,SLC.pW,False,SLC.pInvert,True)   # Or fixed scaling

        elif which=="AG":                # AG only

            pImg = np.copy(self.image_AG)       # Set pixel map image to AG

            if (SLC.pAuto_AG==True):            # Autoscale
                gAG,self.AGpMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,0.0,0.0,True,SLC.pInvert_AG,True)                # Generate Autoscaled onscreen pixmap
                SLC.pB_AG,SLC.pW_AG = pMin,pMax                                                                         # Store new values
            else:
                gAG,self.AGpMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,SLC.pB_AG,SLC.pW_AG,False,SLC.pInvert_AG,True)   # Or fixed scaling

        elif which=="Both":             # PixMap with both DSS and AG

            ###--- Create DSS numpy array with holes, scaled 0-255
            
            pImg = np.copy(self.image) * self.mask_AG       # Set to DSS and zero AG pixel locations

            if (SLC.pAuto==True):            # Autoscale
                compDSS,thePMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,0.0,0.0,True,SLC.pInvert,False)          # Generate Autoscaled onscreen pixmap
                SLC.pB,SLC.pW = pMin,pMax                                                                         # Store new values
            else:
                compDSS,thePMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,SLC.pB,SLC.pW,False,SLC.pInvert,False)   # Fixed scaling

            ###-- Now create AG numpy array with zeroes elsewhere

            pImg = np.copy(self.image_AG)

            if (SLC.pAuto_AG==True):            # Autoscale
                compAG,thePMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,0.0,0.0,True,SLC.pInvert_AG,False)                # Generate Autoscaled onscreen pixmap
                SLC.pB_AG,SLC.pW_AG = pMin,pMax                                                                           # Store new values
            else:
                compAG,thePMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,SLC.pB_AG,SLC.pW_AG,False,SLC.pInvert_AG,False)   # Or fixed scaling

            ###-- And finally add them together to generate the final pixmap (note scaling 0-255)

            pImg = np.flipud(compDSS + compAG)
            # print ("composite max / median  ",np.amax(pImg),np.median(pImg))

            comp,self.BothPMap,pMin,pMax,ok = SLC.NaNscale255(pImg,Fwid,0.0,255.0,False,SLC.pInvert,True)    # Generate Autoscaled onscreen pixmap

        else:
            print ("Invalid Option")
            self.logTxt("Invalid Pixmap Selection",True)

#-----------------------------------------------------------------------------#

    def blinkImg(self):

        '''
            blinkImg - Blink between DSS and Both modes

            This utility allows the user to compare the AG and DSS frames.
            It simply toggles back and forth between the two while the user
            clicks the button over and over (klunky but better than nothing).

        '''

        if self.gotAG:              # Only start if AG frames present

            if self.doBlinkDSS:        # Initiate blinking

                print ("blink DSS")
                self.whatImg="DSS"
                self.DispImg()
                self.doBlinkDSS=False

            else:                   # doBlink was off...turn it on

                print ("blink Both")
                self.whatImg="Both"
                self.DispImg()
                self.doBlinkDSS=True


#-----------------------------------------------------------------------------#

    def searchViz(self):
        '''
            Again based on the ChartMaker original code.

            A new field has been fetched, or the user has pressed "Search VIZIER",
            or the mag limit has changed. Initiate a web search using the astroquery
            package for reference stars. 

            Note that there are multiple VIZIER servers around. We pick SLC.vizServer.
            The options are:

                'vizier.u-strasbg.fr','vizier.nao.ac.jp', 'vizier.hia.nrc.ca',
                'vizier.ast.cam.ac.uk','vizier.cfa.harvard.edu', 'www.ukirt.jach.hawaii.edu',
                'vizier.iucaa.ernet.in','vizier.china-vo.org'

            NOTE: The default row limit for Vizier searches is 50 items. For fainter magnitudes
                  (say above 15 or so, depending on the sky location), this can mean that you
                  don't get all objects, and they are all on one side. We increase the limit
                  to avoid this.

            NOTE: I have eliminated the colour field for now (13Mar18)

        '''

        from astroquery.vizier import Vizier      # To use Vizier server

        # Vizier.ROW_LIMIT = -1 #250              # To make sure we get enough

        if self.fieldName=="":     # Makes no sense with no central coordinate
            return  

        if not self.showViz:       #    or with flag off
            return  

        Vizier.VIZIER_SERVER = SLC.vizServer      # Select VIZIER server

        #--- To restrict the columns that we download, we need to instantiate and modify the Vizier class

        if SLC.survey == "USNO-B1":

            v = Vizier(columns=['_RAJ2000', '_DEJ2000','e_RAJ2000','e_DEJ2000','B2mag', 'R2mag'], column_filters={"R2mag":"<"+str(SLC.refLim)})
            v.ROW_LIMIT = -1
            result = v.query_region(SkyCoord(ra=self.ctrPos[0],dec=self.ctrPos[1], unit=(u.deg, u.deg)),width="1.5d",catalog=["USNO-B1"] )

        elif SLC.survey == "Gaia-DR2":

            v = Vizier(columns=['RA_ICRS', 'DE_ICRS', 'pmRA' , 'pmDE', 'phot_bp_mean_mag', 'phot_rp_mean_mag'], column_filters={"phot_rp_mean_mag":"<"+str(SLC.refLim)})
            v.ROW_LIMIT = -1
            result = v.query_region(SkyCoord(ra=self.ctrPos[0],dec=self.ctrPos[1], unit=(u.deg, u.deg)),width="1.5d",catalog=["I/345/gaia2"] )

        else:

            reply = QtWidgets.QMessageBox.critical(self, 'Error',"No survey defined...Check SLVM_Conf.py",QtWidgets.QMessageBox.Ok)
            return

        if len(result)==0:         # Ooops...something went wrong
            reply = QtWidgets.QMessageBox.critical(self, 'Warning',"Zero-length result. This usually means no reference stars found...Try again.",QtWidgets.QMessageBox.Ok)

            return

        #--- Now reset the Vizier lists

        self.vNam = []          # List of Vizier "names"
        self.vStr = []          #   and their sky coordinates in hh mm ss +dd mm ss format
        self.vPos = []          # These will be a list of numpy decimal degrees RA/Dec
        self.vPM = []           #   and proper motion
        self.vMag = []          #   and vMag
        self.vCol = []          #   and colour

        #--- And stuff arrays

        start = True                    # Need to grab first

        self.logTxt("Retrieved "+str(len(result[0]))+" reference stars",False)

        for i in range(len(result[0])):

            if SLC.survey == "USNO-B1":
                ra,de,pr,pd,bm,rm = result[0][i][0],result[0][i][1],result[0][i][2],result[0][i][3],result[0][i][4],result[0][i][5]    # Blue and Red mags
            elif SLC.survey == "Gaia-DR2":
                ra,de,pr,pd,bm,rm = result[0][i][0],result[0][i][1],result[0][i][2],result[0][i][3],result[0][i][4],result[0][i][5]    # Now also for Gaia!

            self.vNam.append("V"+str(i).zfill(3))          # Creates star name V001, etc.   

            vLS = SkyCoord(ra=ra*u.degree, dec=de*u.degree).to_string('hmsdms',sep=" ",precision=3).split()                         # Loc string,split            
            self.vStr.append(vLS[0]+" "+vLS[1]+" "+("%6.3f" % float(vLS[2]))+"  "+vLS[3]+" "+vLS[4]+" "+("%5.2f" % float(vLS[5])))  # .cat format

            self.vPos.append(np.array([ra,de]))
            self.vPM.append(np.array([pr,pd]))
            self.vMag.append(np.array([rm]))
            self.vCol.append(np.array([bm-rm]))                    

            # if SLC.verbose:                               # Visual feedback
            #     print self.vNam[i]+" "+self.vStr[i]+" "+str(self.vPos[i])+" "+str(self.vPM[i])+" "+str(self.vMag[i])+" "+str(self.vCol[i])

        #--- And stuff GUI widget

        self.ui.vizier_lw.clear()                           # Empty previous Vizier sources
        self.ui.vizier_lw.addItem("Targ  "+self.ctrStr)     # Add field center

        for i in range(len(self.vNam)):
            self.ui.vizier_lw.addItem(self.vNam[i]+"  "+self.vStr[i]+"  "+("%7.2f" % self.vPM[i][0])+" "+("%7.2f" % self.vPM[i][1])+"  "+("%5.2f" % self.vMag[i]) +"  "+("%6.2f" % self.vCol[i]))

        #--- Finally, we have Vizier sources! Annotate them

        self.gotViz = True        # Set the flag

        self.RedrawChart()        # To clear old if going to brighter stars
        self.AnnotateChart()      # Annotate the chart

#-----------------------------------------------------------------------------#    

    def RedrawChart(self):

        '''
            The field has changed, or the user has toggled visibility
            of the Vizier targets or DSS image. The pixmap is still good,
            so just clear the imgScn and redraw the chart.
        '''

        self.ui.imgScn.clear()           # Presumably, this will remove all previous items?

        Fwid = SLC.fSize                 # Size of FITS file
        Wwid = SLC.wSize                 #   and of QGraphicsView

        ok = True        # Assume the best for now

        if self.DSSpMap!=None:           # Avoid crash on first field
            if self.whatImg=="DSS":
                thePMap = self.DSSpMap
            elif self.whatImg=="AG":
                thePMap = self.AGpMap
            elif self.whatImg=="Both":
                thePMap = self.BothPMap
            else:
                print ("Invalid (or None) PMAp in RedrawChart")
                self.logTxt("Invalid (or None) PMAp in RedrawChart",True)

        try:
            zerOff = -1 * (Fwid - Wwid)/2.0  # Need to offset this much to center
            pixItem = self.ui.imgScn.addPixmap(thePMap.scaled(Fwid,Fwid,QtCore.Qt.KeepAspectRatio))     # Add to scene to display and get ptr
            pixItem.setPos(zerOff,zerOff)    # And center it

        except:
            print("Unable to add PixMap in RedrawChart")
            self.logTxt("Unable to add PixMap in RedrawChart",True)

            ok = False

        self.AnnotateChart()             # And annotate it

#-----------------------------------------------------------------------------#    

    def AnnotateChart(self):
        '''
            Place annotations on the chart, based on current status.
            Specifically, the code puts the following onscreen:

               - Boxes red for Vizier objects
        '''

        #--- Mark Science, GWS, and HWS fields

        # try:

        if self.WCSkeys != None:     # Only annotate if we have WCS info

            rX,rY = SLC.wSize/2.0,SLC.wSize/2.0    # Central reference pixel (QGraphicsView)
            rR = SLC.refRad                        # Half-width of ref star squares

            if SLC.pInvert:                        # Want black markings
                col = QtGui.QColor(0,0,0)
            else:                                  # Want white markings
                col = QtGui.QColor(255,255,255)

            #--- Mark Survey stars (if present)

            if self.gotViz and self.showViz:           # If we have them and want them onscreen      

                col = QtGui.QColor(192,0,0)            # Red for Vizier

                nViz = len(self.vPos)                  # We have this many sources

                for i in range(nViz):              # Step through Vizier targets

                    rLoc = self.WCSkeys.wcs_world2pix(np.array([self.vPos[i]]),1)          # Location of object (pix)
                    cX,cY = rLoc[0][0] - self.CRPIX1 + rX, self.CRPIX1 - rLoc[0][1] + rY   # Center of object (QGV coord)
                    self.imgScn.addRect(cX-rR,cY-rR,2*rR,2*rR,QtGui.QPen(col))             # Draw Vizier box
                    lbl = self.imgScn.addText(self.vNam[i], QtGui.QFont('Arial', 11, QtGui.QFont.Light)) # Label for star number
                    lbl.setPos(cX-1.75*rR,cY+0.5*rR)                                       # Position label
                    lbl.setDefaultTextColor(col)                                           # And set colour

        # except:
        #     print("AnnotateChart bypassed")

#-----------------------------------------------------------------------------#    

    def MousePress(self,theText):
        """ User clicked mouse """

        self.xClk, self.yClk = float(theText.split()[0]),float(theText.split()[1])   # Extract position

        print ("User clicked at ",self.xClk,self.yClk)

#-----------------------------------------------------------------------------#    

    def vizLWclick(self):

        ''' User has clicked on an item in the Vizier list. 
            Do nothing for now...
        '''

        row = self.ui.vizier_lw.currentRow()                       # This is the selected row
        sTx = str(self.ui.vizier_lw.selectedItems()[0].text())     # And its text

        t = sTx.split()                                               # Grab individual text fields
        radecStr = t[1]+" "+t[2]+" "+t[3]+" "+t[4]+" "+t[5]+" "+t[6]  # Assemble Ra Dec string

        self.logTxt("Clicked on "+radecStr,False)

#-----------------------------------------------------------------------------#    

    def launchCDS(self):       # Launch CDS website on source

        '''
            Launch a web browser to CDS with the central target (how cool is that)

            For positive declination 09 51 58.56 +06 28 37.80 the location is:
            http://cdsportal.u-strasbg.fr/?target=09%2051%2058.56%20%2B06%2028%2037.80

            and for negative dec     09 51 58.56 -06 28 37.80 it is:
            http://cdsportal.u-strasbg.fr/?target=09%2051%2058.56%20-06%2028%2037.80
            Note the subtle difference in address.
        '''

        import webbrowser      # This makes it all possible

        rh,rm,rs,dd,dm,ds = self.ctrStr.split()[0],self.ctrStr.split()[1],self.ctrStr.split()[2],  \
                            self.ctrStr.split()[3],self.ctrStr.split()[4],self.ctrStr.split()[5]

        ddf = self.ctrStr.split()[3]   # Grab declination to check if negative

        if "-" in ddf :                # Need to formulate negative decs correctly
            lStr = rh+"%20"+rm+"%20"+rs+"%20-"+dd[-2:]+"%20"+dm+"%20"+ds
        else:
            lStr = rh+"%20"+rm+"%20"+rs+"%20%2B"+dd[-2:]+"%20"+dm+"%20"+ds


        webbrowser.open('http://cdsportal.u-strasbg.fr/?target='+lStr)     # Pretty simple!#-----------------------------------------------------------------------------#    

#-----------------------------------------------------------------------------#    

    def UpdateGUI(self):
        '''
            Update the values and controls in the GUI
        '''

        self.ui.VizierSurvey_lbl.setText(SLC.survey)  # Which VIZIER survey to search for ref stars

        if SLC.survey == "USNO-B1":        # We have Rmag and colour
            self.ui.VizListHead_lbl.setText("Object  RA(2000)    Dec(2000)        PM        Rmag   B-R")        
        if SLC.survey == "Gaia-DR2":       # We have only gmag
            self.ui.VizListHead_lbl.setText("Object  RA(2000)    Dec(2000)        PM        rpMag  B-R")        

        self.ui.vizLim_le.setText(str(SLC.refLim))    # Mag limit for returned VIZIER stars

        self.ui.Bk_le.setText(str(SLC.pB))            # Grayscale cuts for DSS
        self.ui.Wh_le.setText(str(SLC.pW))
        self.ui.Bk_AG_le.setText(str(SLC.pB_AG))      #   and for AG Cams
        self.ui.Wh_AG_le.setText(str(SLC.pW_AG))

        if (SLC.pAuto):                               # Autoscale check box (DSS)
            self.ui.Autoscale_cb.setChecked(True)
        else:
            self.ui.Autoscale_cb.setChecked(False)

        if (SLC.pInvert):                             # Autoscale check box (DSS)
            self.ui.Invert_cb.setChecked(True)
        else:
            self.ui.Invert_cb.setChecked(False)

        if (SLC.pAuto_AG):                            # Autoscale check box (AG)
            self.ui.Autoscale_AG_cb.setChecked(True)
        else:
            self.ui.Autoscale_AG_cb.setChecked(False)

        if (SLC.pInvert_AG):                          # Autoscale check box (AG)
            self.ui.Invert_AG_cb.setChecked(True)
        else:
            self.ui.Invert_AG_cb.setChecked(False)

        # print ("Update_GUI self.whatImg = ",self.whatImg)

#-----------------------------------------------------------------------------#    

    def imgZin(self):                             # Zoom in image
        self.ui.img_qgv.scale(self.zoomInFactor,self.zoomInFactor)
        self.curIzoom = self.curIzoom * self.zoomInFactor

        self.DispImg()      # Added Feb23

#-----------------------------------------------------------------------------#    

    def imgZout(self):                            # Zoom out image
        self.ui.img_qgv.scale(self.zoomOutFactor,self.zoomOutFactor)
        self.curIzoom = self.curIzoom * self.zoomOutFactor

        self.DispImg()      # Added Feb23

#-----------------------------------------------------------------------------#    

    def togAutoscale(self):                           # Toggle autoscale (DSS)

        if (self.ui.Autoscale_cb.isChecked()==True):  # If checked
            SLC.pAuto = True                          # Toggle flag
        else:
            SLC.pAuto = False                         # Toggle flag

        try:                       # Avoids errors at startup
            self.GenPixMap("DSS")  # Create new PixMap
            self.DispImg()         # And redisplay

        except:
            pass

#-----------------------------------------------------------------------------#    

    def togInvert(self):                          # Toggle Invert grayscale 

        if (self.ui.Invert_cb.isChecked()==True): # If checked
            SLC.pInvert = True                    # Toggle flag
        else:
            SLC.pInvert = False                   # Toggle flag

        try:                       # Avoids errors at startup
            self.GenPixMap("DSS")  # Create new PixMap
            self.DispImg()         # And redisplay
        except:
            pass

#-----------------------------------------------------------------------------#    

    def togAutoscale_AG(self):                    # Toggle autoscale (AG Cameras)

        if (self.ui.Autoscale_AG_cb.isChecked()==True):   # If checked
            SLC.pAuto_AG = True                           # Toggle flag
        else:
            SLC.pAuto_AG = False                          # Toggle flag

        try:                       # Avoids errors at startup
            self.GenPixMap("AG")   # Create new PixMap
            self.DispImg()         # And redisplay

        except:
            pass

#-----------------------------------------------------------------------------#    

    def togInvert_AG(self):                   # Toggle Invert grayscale (AG Cameras)

        if (self.ui.Invert_AG_cb.isChecked()==True):  # If checked
            SLC.pInvert_AG = True                     # Toggle flag
        else:
            SLC.pInvert_AG = False                    # Toggle flag

        try:                       # Avoids errors at startup
            self.GenPixMap("AG")   # Create new PixMap
            self.DispImg()         # And redisplay
        except:
            pass

#-----------------------------------------------------------------------------#    

    def togShowVizier(self):                          # Toggle visibility of Vizier 

        if (self.ui.showViz_cb.isChecked()==True):    # If checked
            self.showViz = True                       # Toggle flag
        else:
            self.showViz = False                      # Toggle flag

        try:                       # Avoids errors at startup
            self.RedrawChart()     # redisplay chart with or without Vizier stars
        except:
            pass

#-----------------------------------------------------------------------------#    

    def doShowRadioButton(self):
        ''' User clicked one of the Show Image: radio buttons.
            Info: https://www.pythontutorial.net/pyqt/pyqt-qradiobutton/

            This sets the variable self.whatImg to "DSS", "AG", "Both", or "None"
        '''        
        rb = self.sender()
        if rb.isChecked():
            if rb.text()=="DSS":
                self.whatImg = "DSS"
            if rb.text()=="AG Cams":
                if self.gotAG:                  # Can't do this if no AG data
                    self.whatImg = "AG"
                else:
                    self.ui.DSS_rb.setChecked = True         # Go back to DSS
                    self.logTxt("No AG frames yet !", True)  #   and inform user
            if rb.text()=="Both":
                if self.gotAG:                  # Can't do this if no AG data
                    self.whatImg = "Both"
                else:
                    self.ui.DSS_rb.setChecked = True         # Go back to DSS
                    self.logTxt("No AG frames yet !", True)  #   and inform user

            if rb.text()=="None":
                self.whatImg = "None"

        # try:                       # Avoids errors at startup

        self.DispImg()         # redisplay chart with new wishes
        # except:
            # pass

#-----------------------------------------------------------------------------#    

    def yelVizLim(self):   # User has not pressed <CR>...not accepted yet!

        self.ui.vizLim_le.setStyleSheet("QLineEdit{background:rgb(255,255,177)}")  # Yellow

#-----------------------------------------------------------------------------#    

    def setVizLim(self):   # Set Vizier reference limit and do a search

        SLC.refLim = float(self.ui.vizLim_le.text())
        self.ui.vizLim_le.setStyleSheet("QLineEdit{background:rgb(255,255,255)}")  # Back to white

        self.ui.searchViz_pb.click()    # Issue a click by default

#-----------------------------------------------------------------------------#    

    def setBk(self):   # Set DSS black cut

        SLC.pB = float(self.ui.Bk_le.text())
        self.GenPixMap("DSS")   # Create new PixMap
        self.DispImg()          # And redisplay

#-----------------------------------------------------------------------------#    

    def setWh(self):   # Set DSS white cut

        SLC.pW = float(self.ui.Wh_le.text())
        self.GenPixMap("DSS")   # Create new PixMap
        self.DispImg()          # And redisplay

#-----------------------------------------------------------------------------#    

    def setBk_AG(self):   # Set AG black cut

        SLC.pB_AG = float(self.ui.Bk_AG_le.text())
        self.GenPixMap("AG")   # Create new PixMap
        self.DispImg()         # And redisplay

#-----------------------------------------------------------------------------#    

    def setWh_AG(self):   # Set AG white cut

        SLC.pW_AG = float(self.ui.Wh_AG_le.text())
        self.GenPixMap("AG")   # Create new PixMap
        self.DispImg()         # And redisplay

#-----------------------------------------------------------------------------#    

    def logTxt(self,msgT,nLin):      # Log something to the Main message window, optional blank line

        self.ui.msg_lw.addItem(msgT)      # Add to onscreen log
        if nLin:                          #  and newline if wanted
          self.ui.msg_lw.addItem("")          

#-----------------------------------------------------------------------------#    

    def thrdMsg(self,msgT):             # Log something to the  message window
        self.ui.msg_lw.addItem(msgT)      # Add to onscreen log

#-----------------------------------------------------------------------------#    

    def shutdown(self):          # Run shutdown routines

        self.close()             # Close App window

#-----------------------------------------------------------------------------#    

    def dummy_AG(self):          # Create dummy AG Camera frames

        '''
            This routine opens the files "FITS/E.fits" etc. and reprojects
            them into the WCS of the DSS image ()
        '''

        from astropy.io import fits
        from astropy.utils.data import get_pkg_data_filename
        import reproject
        from astropy.wcs import WCS

        hduD = fits.open('FITS/'+self.fieldName+'.fits')[0]      # Grab the DSS image HDU

        #--- Fetch dummy Pleiades HDU's, WCS, data:

        self.eFITS = "DummyFITS/E.fits"                       # Dummy files - FLORIAN
        self.wFITS = "DummyFITS/W.fits"

        hduE = fits.open(get_pkg_data_filename(self.eFITS))[0]    # Extract HDU's
        hduW = fits.open(get_pkg_data_filename(self.wFITS))[0]

        Erep,Efoot = reproject.reproject_interp(hduE, hduD.header)      # Reproject into DSS WCS
        Wrep,Wfoot = reproject.reproject_interp(hduW, hduD.header)

        if SLC.useOnAxis:   # On-axis camera in place

            self.cFITS = "DummyFITS/C.fits"                             # Dummy FITS file - FLORIAN
            hduC = fits.open(get_pkg_data_filename(self.cFITS))[0]      # Get HDU
            Crep,Cfoot = reproject.reproject_interp(hduC, hduD.header)  # Reproject into DSS WCS
            stack = np.array([Crep,Erep,Wrep])                          #   and stack
            mask  = np.array([Cfoot,Efoot,Wfoot])

        else:
            stack = np.array([Erep,Wrep])                               # Stack only E-W frames
            mask  = np.array([Efoot,Wfoot])

        self.image_AG = np.nanmedian(stack,axis=0)     # Create stacked frame
        self.mask_AG = 1.0 - np.sum(mask,axis=0)       #   and 1's / 0's indicating coverage

        self.gotAG = True                              # We have AG frames!

        self.GenPixMap("AG")      # Create PixMap of AG data
        self.GenPixMap("Both")    #   and composite of DSS/AG data

        hduM = fits.PrimaryHDU(self.image_AG)          # Create HDU of new data
        hduM.writeto('mosaic.fits',overwrite=True)     #  and write it out!

        import matplotlib.pyplot as plt

        plt.imsave('mosaic.png', self.image_AG,vmin=0.0,vmax=1000.0)


    def handle_data(self, data):
        # gets executed on scraper_event
#        print(f" {data.timestamp} {data.sender} {data.command_status} {data}")

        if "pwi" in data.sender:
            print(f" {data.timestamp} ra, dec {data.unpack('ra_j2000_hours', 'dec_j2000_degs')}")
        elif "agcam" in data.sender and data.command_status == CommandStatus.DONE:
            print(f" {data.timestamp} {data.flatten().unpack('*.filenames')}")


######################## MAIN APPLICATION LOOP #########################

if __name__ == '__main__':

    app = QtWidgets.QApplication(sys.argv)  # Create an application object

    aWin = AppWin()                         # Instantiate appWin

    scraper = QtLvmScraper()
    scraper.start(app, aWin.handle_data)

    print("done")
#    sys.exit(app.exec_())

