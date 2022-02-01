# A model for visualizing the light output from a star + planet or an eclipsing binary system of two stars. 

from argparse import ZERO_OR_MORE
from re import S
import sys
import math
from tkinter.constants import ANCHOR
import numpy as np
import tkinter
from tkinter import CENTER, NW, StringVar, messagebox
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg)
from matplotlib.ticker import (MultipleLocator)
from systemModel import SystemModel
from os.path import exists
import time

# planet radius RJupiter to RSun
RJupiterToRSun = 0.10047

# create the system model
pm = SystemModel(
    1.0, # star radius in units of Solar radius
    2.0, # planet orbital radius in units of R_Star
    1.0 * RJupiterToRSun, # planet radius in units of Jupiter radius
    0.8, # planet luminosity(0 = planet, > 1 eclipsing binary star)
    15.0, # inclination view angle (degrees) 
    0.5 # star limb-darkening factor (0 = none)
    )


class App_Window(tkinter.Tk):
    def __init__(self, parent):
        tkinter.Tk.__init__(self,parent)
        self.initialize()
    def appQuit(self):
        sys.exit()
    def loadInputs(self):
        # Load the model inputs from a local text file.
        # First ensure that the modelInputs file exists
        if exists("modelInputs.txt"):
            inFile = open("modelInputs.txt","r")
            for line in inFile:
                # ignore comment lines
                if line.index('#') == 0:
                    continue
                split1 = line.split('#')
                split2 = split1[0].split(':')
                if line.startswith('diskOuterRadius'):
                    pm.diskOuterRadius = float(split2[1])
                elif line.startswith('diskInnerRadius'):
                    pm.diskInnerRadius = float(split2[1])
                elif line.startswith('diskBright'):
                    pm.diskBright = float(split2[1])
                elif line.startswith('diskInclination'):
                    pm.diskInclination = float(split2[1])
                elif line.startswith('diskOrbitAngle'):
                    pm.diskOrbitAngle = float(split2[1])
                elif line.startswith('planRadius'):
                    pr = float(split2[1])
                    pm.planRadius = pr*RJupiterToRSun
                    self.planRadEntry.config(textvariable=tkinter.StringVar(value=pr))
                elif line.startswith('planOrbitRadius'):
                    pm.planOrbitRadius = float(split2[1])
                    self.planOrbitRadiusEntry.config(textvariable = tkinter.StringVar(value=pm.planOrbitRadius))
                elif line.startswith('planBright'):
                    pm.planBright = float(split2[1])
                    self.planBrightEntry.config(textvariable = tkinter.StringVar(value=pm.planBright))
                elif line.startswith('starRadius'):
                    pm.starRadius = float(split2[1])
                    self.starRadEntry.config(textvariable=tkinter.StringVar(value=pm.starRadius))
                elif line.startswith('LDF'):
                    pm.LDF = float(split2[1])
                    self.ldfEntry.config(textvariable=tkinter.StringVar(value=pm.LDF))
                else:
                    print("unrecognized line",line)
            inFile.close()
            print("Loaded system model inputs from 'modelInputs.txt'")
            self.plotLDF()
            self.resetCanvas2()
    def initialize(self):
        self.title("ExoTransit_Visualization V2.0")
        self.geometry("1100x800")
        # scaling from brightness to pixel value = uint16/2
        self.pixelScale = 65535
        # define arrays for the orbital angle and the total luminosity
        pm.systemLuminosity = pm.starLuminosity + pm.planetLuminosity
        self.luminosities = np.zeros(pm.nAngleBins, float)
        # the starting angle is chosen so that the planet is at the minimum X position, at -orbitalRadius
        step = 360./pm.nAngleBins
        self.angles = np.arange(0.,360.,step,float)
       # place the entry fields
        ypo = 30
        xpo1 = 0
        xpo2 = 200
        xpo3 = 300
        self.lbl0 = tkinter.Label(self, text = "Star radius (R_star)")
        self.lbl01 = tkinter.Label(self, text = "in units of R_Sun")
        self.starRadEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=pm.starRadius), width = 8)
        self.starRadEntry.bind("<Return>",self.validateStarRadius)
        self.lbl0.place(x=xpo1,y=ypo)
        self.starRadEntry.place(x=xpo2, y=ypo)
        self.lbl01.place(x=xpo3,y=ypo)
        ypo += 25
        self.lbl1 = tkinter.Label(self, text = "Limb-darkening factor (LDF)")
        self.lbl11 = tkinter.Label(self, text = "0 (none) < LDF < 1")
        self.ldfEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=pm.LDF), width = 8)
        self.ldfEntry.bind("<Return>",self.validateLDF)
        self.lbl1.place(x=xpo1, y=ypo)
        self.ldfEntry.place(x=xpo2, y=ypo)
        self.lbl11.place(x=xpo3, y=ypo)
        #
        ypo += 25
        self.lbl3 = tkinter.Label(self, text = "Planet radius (R_planet)")
        self.lbl31 = tkinter.Label(self, text = "R_Jupiter. (R_Jupiter= 0.1 R_Sun)")
        self.planRadEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=pm.planRadius/RJupiterToRSun), width = 8)
        self.planRadEntry.bind("<Return>",self.validatePlanRadius)
        self.lbl3.place(x=xpo1, y=ypo)
        self.planRadEntry.place(x=xpo2, y=ypo)
        self.lbl31.place(x=xpo3, y=ypo)
        #
        ypo += 25
        self.lbl2 = tkinter.Label(self, text = "Planet orbital radius")
        self.lbl21 = tkinter.Label(self, text = "R_Star. Note: R_Mercury = 83 R_Sun")
        self.planOrbitRadiusEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=pm.planOrbitRadius), width = 8)
        self.planOrbitRadiusEntry.bind("<Return>",self.validatePlanOrbitRadius)
        self.lbl2.place(x=xpo1, y=ypo)
        self.planOrbitRadiusEntry.place(x=xpo2, y=ypo)
        self.lbl21.place(x=xpo3, y=ypo)
        #
        ypo += 25
        self.lbl4 = tkinter.Label(self, text = "Planet central brightness")
        self.lbl41 = tkinter.Label(self, text = "0 for planet or >0 for a binary companion star\u2020")
        self.planBrightEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=pm.planBright), width = 8)
        self.planBrightEntry.bind("<Return>",self.validatePlanBright)
        self.lbl4.place(x=xpo1, y=ypo)
        self.planBrightEntry.place(x=xpo2, y=ypo)
        self.lbl41.place(x=xpo3, y=ypo)
        #
        ypo += 25
        self.lbl5 = tkinter.Label(self, text = "Planet orbit inclination angle")
        self.lbl51 = tkinter.Label(self, text = "degrees")
        self.planInclEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=pm.inclination), width = 8)
        self.planInclEntry.bind("<Return>",self.validateInclination)
        self.lbl5.place(x=xpo1, y=ypo)
        self.planInclEntry.place(x=xpo2, y=ypo)
        self.lbl51.place(x=xpo3, y=ypo)
        # noise
        ypo += 25
        self.lbl6 = tkinter.Label(self, text = "Noise rms")
        self.lbl61 = tkinter.Label(self, text = "(about 0.0005 for real NGTS data)")
        self.noiseEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=pm.noise), width = 8)
        self.noiseEntry.bind("<Return>",self.validateNoise)
        self.lbl6.place(x=xpo1, y=ypo)
        self.noiseEntry.place(x=xpo2, y=ypo)
        self.lbl61.place(x=xpo3, y=ypo)
        # run button
        ypo += 30
        self.runBut=tkinter.Button(self, text='Run', command=self.run)
        self.runBut.place(x=xpo1, y=ypo)
        # transit zoom
        self.zoomText = tkinter.StringVar()
        self.zoomText.set("Zoom")
        self.zoombutton=tkinter.Button(self, textvariable=self.zoomText, command=self.zoomTransit)
        self.zoombutton.place(x=xpo1+50, y=ypo)
        # add features like dust disks and stellar luminosity asymmetry
        self.featureText = tkinter.StringVar()
        self.featureText.set("Add feature")
        self.featurebutton=tkinter.Button(self, textvariable=self.featureText, command=self.addFeature)
        self.featurebutton.place(x=xpo1+200, y=ypo)
        # quit button at the top
        self.b3=tkinter.Button(self, text="Quit", command=self.appQuit)
        self.b3.pack(anchor="nw")
        self.b4=tkinter.Button(self, text="Load inputs", command=self.loadInputs)
        self.b4.place(x=190,y=0)
        # read in system inputs
        # add some informational text fields
        ypo += 30
        daggerTxt = tkinter.Label(self, text = "\u2020 Note: Star central brightness, B(r=0) = 0, by definition.  (Luminosity = \u222B B(r) LDF(r) dA)")
        daggerTxt.place(x=xpo1, y=ypo)
        ypo += 20
        self.starLumTxt = tkinter.Label(self, text = "Star luminosity = UNDEFINED")
        self.starLumTxt.place(x=xpo1, y=ypo)
        self.planLumEntryTxt = tkinter.Label(self, text = "Planet luminosity = UNDEFINED")
        self.planLumEntryTxt.place(x=xpo1+200, y=ypo)
        ypo += 20
        self.systLumEntryTxt = tkinter.Label(self, text = "System luminosity(0 deg) = UNDEFINED")
        self.systLumEntryTxt.place(x=xpo1, y=ypo)
        ypo += 20
        self.ang90LumTxt = tkinter.Label(self, text = "System luminosity(90 deg) = UNDEFINED")
        self.ang90LumTxt.place(x=xpo1, y=ypo)
        self.dipTxt = tkinter.Label(self, text = "Normalized luminosity dip = UNDEFINED")
        self.dipTxt.place(x=xpo1+250, y=ypo)
        # add the Brightness vs Radius figure
        Fig = Figure(figsize=(5.5,3.5),dpi=100)
        FigSubPlot = Fig.add_subplot(111)
        FigSubPlot.set_xlabel("Radius Fraction")
        FigSubPlot.set_ylabel("LDF Brightness")
        FigSubPlot.set_xlim(0, 1.1)
        FigSubPlot.set_ylim(0, 1.1)
        Fig.tight_layout()
        ra = np.linspace(0, 1, num=50, endpoint=True)
        ya = np.zeros(50)
        indx = 0
        for r in ra:
            r2 = r * r
            ya[indx] = self.limbDarkFactor(r2)
            indx += 1
        self.line1, = FigSubPlot.plot(ra,ya,"r-")
        self.canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(Fig, master=self)
        self.canvas.get_tk_widget().pack(side=tkinter.TOP,anchor="ne", expand=True)
        # create the system cartoon and transit figure in another canvas
        Fig2 = Figure(figsize=(10,5))
        self.canvas2 = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(Fig2, master=self)
        self.F2Cartoon = Fig2.add_subplot(211)
        self.updateStarImage(pm)
        self.updatePlanetImage(pm)
        self.updateSystemImage(pm,0)
        self.ims = self.F2Cartoon.imshow(pm.systemImage, cmap="gray", vmin=0, vmax=65535)
        self.F2Cartoon.axis("off")
        self.F2Transit = Fig2.add_subplot(212)
        self.F2Transit.autoscale(True,axis='both',tight=True)
        angla = np.linspace(0, 360, num=60, endpoint=True)
        intena = np.zeros(60)
        self.line2, = self.F2Transit.plot(angla,intena,"o")
        self.canvas2.get_tk_widget().pack(side=tkinter.BOTTOM, expand=True)
        self.canvas2.draw()
        self.UpdateRunInfo()
        # finish up
        self.resizable(True, True)
#        self.update()
    def UpdateRunInfo(self):
        # Updates the detailed results below the Run button
        self.starLumTxt.configure(text="Star luminosity = {:.4f}".format(pm.starLuminosity))
        self.planLumEntryTxt.configure(text="Planet luminosity = {:.4f}".format(pm.planetLuminosity))
        lumMax = np.double(pm.starLuminosity + pm.planetLuminosity)
        self.systLumEntryTxt.configure(text="System luminosity(0 deg) = {:.4f}".format(lumMax))
        self.updateSystemImage(pm,90.)
        pm.luminosity90 = pm.systemLuminosity
        self.ang90LumTxt.configure(text="System luminosity(90 deg) = {:.4f}".format(pm.luminosity90))
        dip = np.double((lumMax - pm.luminosity90) / lumMax)
        self.dipTxt.configure(text="Normalized luminosity dip = {:.4f}".format(dip))
        # move the planet back to 0
        self.updateSystemImage(pm,0.)
    def updateStarImage(self, pm):
        # there is no need to initialize the star image to 0's because the radius is always 1
        centerRow = int(pm.starImage.shape[0]/2)
        centerCol = int(pm.starImage.shape[1]/2)
        maxBright = 1.
        if pm.planBright > maxBright: maxBright = pm.planBright
        if pm.diskBright > maxBright: maxBright = pm.diskBright
        self.pixelScale = 65535./maxBright
        # determine the step size per pixel for a circle of radius = 1
        step = pm.stepSize
        xs = np.double(-1)
        while xs < 1:
            ymax = np.sqrt(1 - xs*xs)
            ys = -ymax
            xPixel = int(centerRow + xs/step)
            while ys < ymax:
                # distance^2 from star center
                rs2 = xs*xs + ys*ys
                yPixel = centerCol + int(ys/step)
                pixelValue = int(self.pixelScale*self.limbDarkFactor(rs2))
                if pixelValue <= 0: pixelValue = 1
                elif pixelValue > 65535: pixelValue = 65535
                pm.starImage[xPixel][yPixel] = np.uint16(pixelValue)
                ys += step
            xs += step
        # define the luminosity normalization if it hasn't been done
        if pm.luminosityNorm < 0:
            pm.luminosityNorm = np.double(1/(self.pixelScale*np.double(np.count_nonzero(pm.starImage))))
        pm.starLuminosity = np.double(np.sum(pm.starImage))*pm.luminosityNorm
    def updatePlanetImage(self,pm):
        # reset to 0
        pm.planetImage.fill(0)
        maxBright = 1.
        if pm.planBright > maxBright: maxBright = pm.planBright
        if pm.diskBright > maxBright: maxBright = pm.diskBright
        self.pixelScale = 65535./maxBright
        centerRow = pm.planetImage.shape[0]/2
        centerCol = pm.planetImage.shape[1]/2
        step = pm.stepSize
        # scale the planet radius to the star radius
        planRad = pm.planRadius / pm.starRadius
        pr2 = planRad * planRad
        xp = -planRad
        isStar = (pm.planBright > 0)
        while xp < planRad:
            xPixel = int(centerRow + xp/step)
            ymax = np.sqrt(pr2 - xp*xp)
            yp = -ymax
            while yp < ymax:
                yPixel = int(centerCol + yp/step)
                # make the planet not-totally-dark to distinguish it from a deep-space background
                pixelValue = 1
                if isStar: 
                    # apply limb-darkening to the 2ndry star
                    rp2 = xp*xp + yp*yp
                    pixelValue = int(self.pixelScale*self.limbDarkFactor(rp2/pr2)*pm.planBright)
                    if pixelValue <= 0: pixelValue = 1
                    elif pixelValue > 65535: pixelValue = 65535
                pm.planetImage[xPixel][yPixel] = np.uint16(pixelValue)
                yp += step
            xp += step
        # Add a dustdisk?
        if pm.featureState == 1 and pm.diskOuterRadius > 0:
            diskRad = pm.diskOuterRadius * pm.starRadius
            dr2 = diskRad*diskRad
            diskIn = pm.diskInnerRadius * pm.starRadius
            di2 = diskIn*diskIn
            # define the rotation cosines
            angle1 = pm.diskOrbitAngle * np.pi / 180.
            cs1 = np.cos(angle1)
            sn1 = np.sin(angle1)
            angle2 = pm.diskInclination * np.pi / 180.
            cs2 = np.cos(angle2)
            sn2 = np.sin(angle2)
            xd = -diskRad
            pixelValue = int(self.pixelScale*pm.diskBright)
            if pixelValue <= 0: pixelValue = 1
            elif pixelValue > 65535: pixelValue = 65535
            while xd < diskRad:
                ymax = np.sqrt(dr2 - xd*xd)
                yd = -ymax
                while yd < ymax:
                    rd2 = yd*yd + xd*xd
                    if rd2 < di2:
                        yd += step
                        continue
                    xr1 =  cs1*xd + sn1*yd
                    yr1 = -sn1*xd + cs1*yd
                    xr2 =  cs2*xr1 + sn2*yr1
                    yr2 = -sn2*xr1 + cs1*yr1
                    xPixel = int(centerRow + xr2/step)
                    yPixel = int(centerCol + yr2/step)
                    pm.planetImage[xPixel][yPixel]=np.uint16(pixelValue)
                    yd += step
                xd += step
            print("Added a dust disk to the planet image")
        # dust disk done
        sum = np.sum(pm.planetImage)
        pixCnt = np.count_nonzero(pm.planetImage)
        # subtract the pixel value = 1 dum = pixCnt and normalize to the number of pixels in the star
        pm.planetLuminosity = np.double(sum - pixCnt)*pm.luminosityNorm
        # normalize to the star radius
        pm.planetLuminosity /= pm.starRadius
    def setAngleRange(self):
        # Finds one or two transits without updating the image arrays
        # Default for zoomState == 0
        if pm.zoomState == 0 or (pm.zoomState == 2 and pm.planetLuminosity == 0):
            self.angles = np.arange(0., 360., 360./pm.nAngleBins)
        else:
            start = -1
            nt = 0
            minAngle = 0.
            maxAngle = 360.
            if pm.planetLuminosity == 0:
                maxAngle = 180
            # adjust the step size so that there will be about 20 occluded points in the transit curve
            step = pm.starRadius/pm.planOrbitRadius
            step /= 20.
            step *= 180. / np.pi
            for angleDeg in np.arange(0.,maxAngle,step):
                angle = angleDeg * np.pi/180.
                xcp = -pm.planOrbitRadius * np.cos(angle)
                ycp = pm.planOrbitRadius * np.sin(angle) * np.tan(pm.inclination*np.pi/180)
                rp2 = xcp*xcp + ycp*ycp
                radCut2 = (1.0+pm.planRadius)*(1.0+pm.planRadius)
                if rp2 < radCut2 and start < 0:
                    # the planet center position is within the star disk
                    start = angleDeg
                elif rp2 > radCut2 and start > 0:
                    nt += 1
                    if nt == pm.zoomState:
                        # pad the angle range
                        pad = 0.3 * (angleDeg - start)
                        print("transit",start,angleDeg,pad)
                        start = round(start-pad)
                        end = round(angleDeg+pad)
                        if start == angleDeg:
                            print("Failed to find an expected dip")
                            return
                        step = (end - start)/pm.nAngleBins
                        for bin in np.arange(0,pm.nAngleBins,1):
                            self.angles[bin] = start + bin * step
                        break
                    start = -1
        # update the transit curve
#        print("angle range",self.angles[0],self.angles[pm.nAngleBins-1])
        self.line2.set_xdata(self.angles)
        self.F2Transit.set_xlim([self.angles[0],self.angles[pm.nAngleBins-1]])
        self.resetTransitCurve()
        self.canvas2.draw()
        self.canvas2.flush_events()
    def updateSystemImage(self, pm, angleDeg):
        # updates the image of the star/planet system and sets the system luminosity
        # The first step is to see if the planet image is within the cartoon canvas.
        step = pm.stepSize
        planetInFront = bool(angleDeg < 180)
        angle = angleDeg * np.pi/180.
        # x, y position of the planet center viewed at the inclination angle in
        # coordinate system where RStar = 1.
        xcp = -pm.planOrbitRadius * np.cos(angle)
        # see if the planet center is outside the cartoon window
        planetEdge = abs(xcp - pm.planRadius)
        if planetEdge > 2.5:
            pm.systemImage = pm.starImage.copy()
            pm.systemLuminosity = pm.starLuminosity + pm.planetLuminosity
            return
        ycp = pm.planOrbitRadius * np.sin(angle) * np.tan(pm.inclination*np.pi/180)
        # find the number of pixels to shift the planet position into the system image
        iyShift = int(xcp/step) + int(pm.systemImage.shape[1]/2 - pm.planetImage.shape[1]/2)
        ixShift = int(ycp/step)
        # get a tuple array of non-zero planet pixel indices
        nzPixels = np.nonzero(pm.planetImage)
        maxX = pm.systemImage.shape[0] 
        maxY = pm.systemImage.shape[1] 
        if planetInFront:
            # copy the star image into the system image
            pm.systemImage = pm.starImage.copy()
            # shift and paste non-zero planet pixels
            for ii in range(nzPixels[0].shape[0]):
                ix = nzPixels[0][ii]
                iy = nzPixels[1][ii]
                ixs = ix+ixShift
                iys = iy+iyShift
                if ixs < 0 or ixs >= maxX: continue
                if iys < 0 or iys >= maxY: continue
                pm.systemImage[ixs][iys] = pm.planetImage[ix][iy]
        else:
            # The planet is behind the star. Clear the system image
            pm.systemImage.fill(0)
            # paste the shifted planet image pixels
            for ii in range(nzPixels[0].shape[0]):
                ix = nzPixels[0][ii]
                iy = nzPixels[1][ii]
                ixs = ix+ixShift
                iys = iy+iyShift
                if ixs < 0 or ixs >= maxX: continue
                if iys < 0 or iys >= maxY: continue
                pm.systemImage[ixs][iys]= pm.planetImage[ix][iy]
            # Get a list of the non-zero star pixels
            nzPixels = np.nonzero(pm.starImage)
            # paste them into the system image
            for ii in range(nzPixels[0].shape[0]):
                ix = nzPixels[0][ii]
                iy = nzPixels[1][ii]
                pm.systemImage[ix][iy]= pm.starImage[ix][iy]
        pm.systemLuminosity = np.double(np.sum(pm.systemImage))*pm.luminosityNorm
    def run(self):
        # Run the model and update the screen with the results.
        if pm.zoomState > 1 and pm.planetLuminosity == 0:
            print("No secondary transit is expected")
            return
        self.runBut.configure(text= '...')
        self.resetTransitCurve()
        start = time.time()
        sysCnt = np.count_nonzero(pm.starImage)+np.count_nonzero(pm.planetImage)
        isPlanet = bool(pm.planetLuminosity == 0)
        self.luminosities.fill(pm.systemLuminosity)
        # guess at the transit curve Y limits
        ymax = pm.systemLuminosity + 0.001
        self.updateSystemImage(pm,90.)
        ymin = pm.systemLuminosity
        if pm.planetLuminosity > 0:
            # find low value for a bright companion
            self.updateSystemImage(pm,270.)
            if pm.systemLuminosity < ymin:
                ymin = pm.systemLuminosity
        ymin -= 0.001
        self.F2Transit.set_ylim(ymin, ymax)
        self.canvas2.draw()
        self.canvas2.flush_events()
        cnt = 0
        for angle in self.angles:
            if isPlanet and angle > 180:
                break
            self.updateSystemImage(pm,angle)
            self.ims.set_data(pm.systemImage)
            self.luminosities[cnt] = pm.systemLuminosity
            self.line2.set_ydata(self.luminosities)
            self.canvas2.draw()
            self.canvas2.flush_events()
            cnt += 1
        end = time.time()
        print("cpu time {:.2f} s".format(end-start))
        self.UpdateRunInfo()
        self.runBut.configure(text='Run')
        self.updateTransitCurve()
    def resetCanvas2(self):
        self.updateStarImage(pm)
        self.updatePlanetImage(pm)
        self.updateSystemImage(pm, 0.)
        self.ims.set_data(pm.systemImage)
        self.UpdateRunInfo()
        self.resetTransitCurve()
        self.setAngleRange()
    def resetTransitCurve(self):
        self.luminosities.fill(0.)
        self.updateTransitCurve()
    def updateTransitCurve(self):
        # update the transit curve
        self.line2.set_xdata(self.angles)
        transitY = self.luminosities
        if pm.noise > 0:
            transitY = self.luminosities + np.random.normal(0, pm.noise, pm.nAngleBins)
        self.line2.set_ydata(transitY)
        maxLum = np.amax(transitY)
        minLum = np.amin(transitY)
        if (maxLum == minLum): 
            maxLum = minLum + 0.01
        yRange = maxLum - minLum
        loY = minLum - 0.1*yRange
        hiY = maxLum + 0.1*yRange
        self.F2Transit.set_ylim(loY, hiY)
        self.canvas2.draw()
        self.canvas2.flush_events()
    def plotLDF(self):
        # Get the entries in the text fields and validate
        ra = np.linspace(0, 1, 50)
        ya = np.zeros(50)
        indx = 0
        for r in ra:
            r2 = r * r
            ya[indx] = self.limbDarkFactor(r2)
            indx += 1
        self.line1.set_ydata(ya)
        self.canvas.draw()
        self.canvas.flush_events()
    def addFeature(self):
        pm.featureState = (pm.featureState+1) % 3
        if pm.featureState == 0:
            self.featureText.set("No feature")
        elif pm.featureState == 1:
            if pm.diskOuterRadius == 0.:
                print("Dust disk model inputs haven't been loaded.")
            else:
                self.featureText.set("Dust disk")
        else:
            self.featureText.set("Future feature")
#        self.updatePlanetImage(pm)
#        self.updateSystemImage(pm, 0.)
        self.resetCanvas2()
    def zoomTransit(self):
        pm.zoomState = (pm.zoomState+1) % 3
        if pm.zoomState == 0:
            self.zoomText.set("Zoom")
        elif pm.zoomState == 1:
            self.zoomText.set("Zoom Primary transit")
        elif pm.zoomState == 2:
            self.zoomText.set("Zoom 2ndry transit")
        else:
            self.zoomText.set("???")
        # testing
        # Set the angle range for the current zoom state
        self.setAngleRange()
    def validateInclination(self, event):
        # inclination angle
        try:
            inclin = float(self.planInclEntry.get())
        except:
            self.planInclEntry.config(bg="red")
            return False
        if (inclin < 0 or inclin > 90):
            print("The inclination angle must be between 0 and 90 degrees")
            self.planInclEntry.config(bg="red")
            return False
        maxInclination = 180*math.atan((pm.starRadius+pm.planRadius)/pm.planOrbitRadius)/np.pi
        if (inclin > maxInclination):
            mess = "The inclination angle should be less than {:.0f} degrees to see a transit".format(maxInclination)
            print(mess)
            self.planInclEntry.config(bg="red")
            return False
        self.planInclEntry.config(bg="white")
        pm.inclination = inclin
        self.planInclEntry.config(textvariable=tkinter.StringVar(value=pm.inclination))
        self.resetTransitCurve()
        return True
    def validateNoise(self, event):
        # noise
        try:
            noise = float(self.noiseEntry.get())
        except:
            self.noiseEntry.config(bg="red")
            return False
        if (noise < 0 or noise > 0.005):
            print("The noise rms must be between 0 and 0.005")
            self.noiseEntry.config(bg="red")
            return False
        self.noiseEntry.config(bg="white")
        pm.noise = noise
        self.noiseEntry.config(textvariable=tkinter.StringVar(value=pm.noise))
        self.updateTransitCurve()
        return True
    def validatePlanBright(self, event):
        # planet brightness
        try:
            planBright = float(self.planBrightEntry.get())
        except:
            self.planBrightEntry.config(bg="red")
            return False
        if (planBright < 0):
            print("The planet brightness can't be < 0")
            self.planBrightEntry.config(bg="red")
            return False
        pm.planBright = planBright
        self.planBrightEntry.config(bg="white")
        self.planBrightEntry.config(textvariable=tkinter.StringVar(value=pm.planBright))
        self.resetCanvas2()
        return True
    def validatePlanOrbitRadius(self, event):
        # orbital radius
        try:
            orbRad = float(self.planOrbitRadiusEntry.get())
        except:
            self.planOrbitRadiusEntry.config(bg="red")
            return False
        # ensure that the planet isn't inside the star
        if orbRad < pm.starRadius + pm.planRadius:
            print("The planet orbital radius must be > R_Star + R_Planet")
            self.planOrbitRadiusEntry.config(bg="red")
            return False
        pm.planOrbitRadius = orbRad
        self.planOrbitRadiusEntry.config(bg="white")
        self.planOrbitRadiusEntry.config(textvariable=tkinter.StringVar(value=orbRad))
        self.resetCanvas2()
        return True
    def validatePlanRadius(self, event):
        # planet radius
        try:
            planRad = float(self.planRadEntry.get())
        except:
            self.planRadEntry.config(bg="red")
            return False
        planRad *= RJupiterToRSun
        if planRad < 0 or planRad > pm.starRadius:
            print("The planet radius must in the range 0 to R_star")
            self.planRadEntry.config(bg="red")
            return False
        self.planRadEntry.config(bg="white")
        pm.planRadius = planRad
        self.planRadEntry.config(textvariable=tkinter.StringVar(value=pm.planRadius/RJupiterToRSun))
        self.resetCanvas2()
        return True
    def validateStarRadius(self, event):
        try:
            srad = float(self.starRadEntry.get())
        except:
            self.starRadEntry.config(bg="red")
            return False
        # next ensure that the radius value is valid
        if (srad < 0.1): 
            print("The star radius should be greater than 0.1 R_Sun")
            self.starRadEntry.config(bg="red")
            return False
        else:
            self.starRadEntry.config(bg="white")
        pm.starRadius = srad
        self.starRadEntry.config(textvariable=tkinter.StringVar(value=srad))
        self.resetCanvas2()
        return True
    def validateLDF(self, event):
        try:
            ldf = float(self.ldfEntry.get())
        except:
            self.ldfEntry.config(bg="red")
            return False
        if (ldf < 0 or ldf > 1): 
            print("The limb-darkening factor must be between 0 and 1")
            self.ldfEntry.config(bg="red")
            return False
        else:
            self.ldfEntry.config(bg="white")
        print("LDF old",pm.LDF,"new",ldf)
        pm.LDF = ldf
        self.ldfEntry.config(textvariable=tkinter.StringVar(value=ldf))
        self.resetCanvas2()
        self.plotLDF()
        return True
    def limbDarkFactor(self, radialDistance2):
        # calculate the limb-darkening factor at given radial distance^2 with a given limbDark factor.
        if (radialDistance2 < 0 or radialDistance2 > 1): return 0
        arg = 1 - radialDistance2
        arg = 1 - pm.LDF * (1 - np.sqrt(arg))
        return arg

if __name__ == "__main__":
    MainWindow = App_Window(None)
    '''
    str = "Hit the return key after changing a variable to update the model.\n\n"
    str += "Hit the Run button after you are done to start the animation.\n\n"
    str += "Hitting the return key in the 'noise rms' field randomizes the noise."
    tkinter.messagebox.showinfo(message=str)
    '''
    MainWindow.mainloop()

