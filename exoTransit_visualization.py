# A model for visualizing the light output from a star + planet or an eclipsing binary system of two stars. 

from argparse import ZERO_OR_MORE
import sys
import math
from tkinter.constants import ANCHOR
import numpy as np
import tkinter
from tkinter import CENTER, NW, messagebox
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg)
from matplotlib.ticker import (MultipleLocator)

from transitCode import *
from systemModel import SystemModel

# normalization factors
# planet radius RJupiter to RSun
RJupiterToRSun = 0.10047
# planet orbital radius in AU to RSun
AU_ToRSun = 1/0.0046524726

# default values for the planet model
starRad = 1.0 # RSun
orbRad = 0.01  # AU -> RSun
planRad = 1.0 # RJupiter -> RSun
planBright = 1.
inclin = 0
ldf = 0.
noise = 0.0
# create the system model
pm = SystemModel(
    starRad, # star radius in units of Solar radius
    orbRad * AU_ToRSun, # planet orbital radius in units of AU
    planRad * RJupiterToRSun, # planet radius in units of Jupiter radius
    planBright, # planet luminosity(0 = planet, > 1 eclipsing binary star)
    inclin, # inclination view angle (degrees) 
    ldf # star limb-darkening factor (0 = none)
    )

class App_Window(tkinter.Tk):
    def __init__(self, parent):
        tkinter.Tk.__init__(self,parent)
        self.initialize()
    def appQuit(self):
        sys.exit()
    def initialize(self):
        self.title("ExoTransit_Visualization V1.1")
        self.geometry("1000x800")
        # define arrays for the orbital angle and the total luminosity
        nBins = 60
        self.angles = np.zeros(nBins)
        # the starting angle is chosen so that the planet is at the minimum X position, at -orbitalRadius
        step = 360 / nBins
        indx = 0
        for iang in range(nBins):
            self.angles[indx]=iang * step
            indx += 1
        self.luminosities = np.zeros(nBins)
        pm.needsUpdate = True
       # place the entry fields
        ypo = 30
        xpo1 = 0
        xpo2 = 200
        xpo3 = 300
        self.lbl0 = tkinter.Label(self, text = "Star radius (R_star)")
        self.lbl01 = tkinter.Label(self, text = "in units of R_Sun")
        self.starRadEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=starRad), width = 8)
        self.lbl0.place(x=xpo1,y=ypo)
        self.starRadEntry.place(x=xpo2, y=ypo)
        self.lbl01.place(x=xpo3,y=ypo)
        ypo += 25
        self.lbl1 = tkinter.Label(self, text = "Limb-darkening factor (LDF)")
        self.lbl11 = tkinter.Label(self, text = "0 (none) < LDF < 1")
        self.ldfEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=ldf), width = 8)
        self.lbl1.place(x=xpo1, y=ypo)
        self.ldfEntry.place(x=xpo2, y=ypo)
        self.lbl11.place(x=xpo3, y=ypo)
        #
        ypo += 25
        self.lbl3 = tkinter.Label(self, text = "Planet radius (R_planet)")
        self.lbl31 = tkinter.Label(self, text = "R_Jupiter. (R_Jupiter= 0.1 R_Sun)")
        self.planRadEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=planRad), width = 8)
        self.lbl3.place(x=xpo1, y=ypo)
        self.planRadEntry.place(x=xpo2, y=ypo)
        self.lbl31.place(x=xpo3, y=ypo)
        #
        ypo += 25
        self.lbl2 = tkinter.Label(self, text = "Planet orbital radius")
        self.lbl21 = tkinter.Label(self, text = "AU (R_Earth = 1 AU, R_Mercury = 0.387 AU)")
        self.orbRadEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=orbRad), width = 8)
        self.lbl2.place(x=xpo1, y=ypo)
        self.orbRadEntry.place(x=xpo2, y=ypo)
        self.lbl21.place(x=xpo3, y=ypo)
        #
        ypo += 25
        self.lbl4 = tkinter.Label(self, text = "Planet central brightness")
        self.lbl41 = tkinter.Label(self, text = "0 for planet or >0 for a binary companion \u2020")
        self.planBrightEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=planBright), width = 8)
        self.lbl4.place(x=xpo1, y=ypo)
        self.planBrightEntry.place(x=xpo2, y=ypo)
        self.lbl41.place(x=xpo3, y=ypo)
        #
        ypo += 25
        self.lbl5 = tkinter.Label(self, text = "Planet orbit inclination angle")
        self.lbl51 = tkinter.Label(self, text = "degrees")
        self.planInclEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=inclin), width = 8)
        self.lbl5.place(x=xpo1, y=ypo)
        self.planInclEntry.place(x=xpo2, y=ypo)
        self.lbl51.place(x=xpo3, y=ypo)
        # noise
        ypo += 25
        self.lbl6 = tkinter.Label(self, text = "Noise rms")
        self.lbl61 = tkinter.Label(self, text = "(about 0.0005 for real NGTS data)")
        self.noiseEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=noise), width = 8)
        self.lbl6.place(x=xpo1, y=ypo)
        self.noiseEntry.place(x=xpo2, y=ypo)
        self.lbl61.place(x=xpo3, y=ypo)
        # run button
        ypo += 30
        self.b1=tkinter.Button(self, text='Run', command=self.run)
        self.b1.place(x=xpo1, y=ypo)
        self.plotbutton=tkinter.Button(self, text="Update LDF Plot", command=self.plotLDF)
        self.plotbutton.place(x=xpo1+50, y=ypo)
        self.zoomText = tkinter.StringVar()
        self.zoomText.set("Zoom None")
        self.zoombutton=tkinter.Button(self, textvariable=self.zoomText, command=self.zoomTransit)
        self.zoombutton.place(x=xpo1+180, y=ypo)
        # quit button at the top
        self.b3=tkinter.Button(self, text='Quit', command=self.appQuit)
        self.b3.pack(anchor="nw")
        # add some informational text fields
        ypo += 30
        self.starLumTxt = tkinter.Label(self, text = "Star luminosity = UNDEFINED")
        self.starLumTxt.place(x=xpo1, y=ypo)
        ypo += 20
        self.planBrightEntryTxt = tkinter.Label(self, text = "Planet luminosity = UNDEFINED")
        self.planBrightEntryTxt.place(x=xpo1, y=ypo)
        ''''
        ypo += 20
        self.ang0LumTxt = tkinter.Label(self, text = "System luminosity(0 deg) = UNDEFINED")
        self.ang0LumTxt.place(x=xpo1, y=ypo)
        ypo += 20
        self.ang90LumTxt = tkinter.Label(self, text = "System luminosity(90 deg) = UNDEFINED")
        self.ang90LumTxt.place(x=xpo1, y=ypo)
        ypo += 20
        self.ang270LumTxt = tkinter.Label(self, text = "System luminosity(270 deg) = UNDEFINED")
        self.ang270LumTxt.place(x=xpo1, y=ypo)
        ypo += 20
        self.planRadEntryEstTxt = tkinter.Label(self, text = "Estimated planet radius = UNDEFINED")
        self.planRadEntryEstTxt.place(x=xpo1, y=ypo)
        ypo += 20
        self.dagger = tkinter.Label(self, text = "\u2020 The 'planet' is actually a binary companion star if the central brightness is set > 0. Star central brightness = 1.")
        self.dagger.place(x=xpo1, y=ypo)
        ypo += 20
        self.dagger2 = tkinter.Label(self, text = "  Luminosity is the weighted average brightness summed over the star or planet disk area, including the effects of limb-darkening.")
        self.dagger2.place(x=xpo1, y=ypo)
        ypo += 20
        self.dagger1 = tkinter.Label(self, text = "  The binary companion luminosity is calculated assuming it has the same limb-darkening factor as the primary star.")
        self.dagger1.place(x=xpo1, y=ypo)
        '''''
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
            ya[indx] = limbDarkFactor(pm, r2)
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
        self.updateSystemImage(pm, 0)
        self.F2Cartoon.imshow(pm.systemImage, cmap="gray", vmin=0, vmax=255)
        self.F2Cartoon.axis("off")
        F2Transit = Fig2.add_subplot(212)
        angla = np.linspace(0, 360, num=60, endpoint=True)
        intena = np.zeros(60)
        self.line2, = F2Transit.plot(angla,intena,"o")
        self.canvas2.get_tk_widget().pack(side=tkinter.BOTTOM, expand=True)
        self.canvas2.draw()
        self.UpdateRunInfo()
        # finish up
        self.resizable(True, True)
        self.refreshLDF(ra,ya)
        self.update()
    def refreshLDF(self,x,y):
        # Update the light-darkening factor plot
        self.line1.set_xdata(x)
        self.line1.set_ydata(y)
        self.canvas.draw()
        self.canvas.flush_events()
    def drawCartoon(self):
        # This updates the cartoon in the middle of the screen.
        # Note that here we display light levels using brightness
        yp = self.CFigScale * pm.orbitalRadius * np.tan(pm.inclination*np.pi/180)
        self.planet90.center = 0, yp
        self.planet90.radius = 0.1 * self.CFigScale * pm.planetRadius / pm.starRadius
        # define the color (actually the gray scale)
        clr = pm.planetBrightness * self.starColorScale
        if (clr > 1): clr = 1
        strClr = "{:.1f}".format(clr)
        self.planet90.set_facecolor(strClr)
        self.canvas.draw()
        self.canvas.flush_events()
    def UpdateRunInfo(self):
        # Updates the detailed results below the Run button
        self.starLumTxt.configure(text="Star luminosity = {:.3f}".format(pm.starLuminosity))
        self.planBrightEntryTxt.configure(text="Planet luminosity = {:.3f}".format(pm.planetLuminosity))
    def updateStarImage(self, pm):
        # there is no need to initialize the star image to 0's because the radius is always 1
        if pm.starImageOK:
            return
        xc = pm.starImage.shape[0]/2
        yc = pm.starImage.shape[1]/2
        # determine the step size per pixel for a circle of radius = 1
        step = pm.stepSize
        xs = -1
        while xs < 1:
            ymax = np.sqrt(1 - xs*xs)
            ys = -ymax
            xPixel = int(xc + xs/step)
            while ys < ymax:
                # distance^2 from star center
                rs2 = xs*xs + ys*ys
                yPixel = int(yc + ys/step)
                pm.starImage[xPixel][yPixel] = int(200*limbDarkFactor(pm,rs2))
                ys += step
            xs += step
        pm.starLuminosity = np.sum(pm.starImage) / (200 * np.count_nonzero(pm.starImage))
        pm.starImageOK = True
#        print("starLuminosity {:.3f}".format(pm.starLuminosity))
    def updatePlanetImage(self,pm):
        if pm.planetImageOK:
            return
        # reset to 0
        pm.planetImage.fill(0)
        xc = pm.planetImage.shape[0]/2
        yc = pm.planetImage.shape[1]/2
        step = pm.stepSize
        # scale the planet radius to the star radius
        planRad = pm.planetRadius / pm.starRadius
        pr2 = planRad * planRad
        xp = -planRad
        isStar = (pm.planetBrightness > 0)
        cnt = 0
        while xp < planRad:
            xPixel = int(xc + xp/step)
            ymax = np.sqrt(pr2 - xp*xp)
            yp = -ymax
            while yp < ymax:
                yPixel = int(yc + yp/step)
                # make the planet not-totally-dark to distinguish it from a black background
                pixelValue = 1
                if isStar: 
                    rp2 = xp*xp + yp*yp
                    pixelValue = int(200*limbDarkFactor(pm,rp2/pr2)*pm.planetBrightness)
                    if pixelValue == 0: pixelValue = 1
                pm.planetImage[xPixel][yPixel] = pixelValue
                yp += step
            xp += step
        sum = np.sum(pm.planetImage)
        pixCnt = np.count_nonzero(pm.planetImage)
#        print("chk", sum, pixCnt, np.count_nonzero(pm.starImage),"planetRadius", pm.planetRadius)
        # subtract the pixel value = 1 dum = pixCnt and normalize to the number of pixels in the star
        pm.planetLuminosity = float(sum - pixCnt) / float(200*np.count_nonzero(pm.starImage))
        # normalize to the star radius
        pm.planetLuminosity /= pm.starRadius
        pm.planetImageOK = True
#        print("planetLuminosity {:.3f}".format(pm.planetLuminosity), pm.starRadius)
    def updateSystemImage(self, pm, angleDeg):
        # Shift the planet by angleDeg to update pm.systemImage
        # NOTE ROTATED IMAGE shape[0] = y, shape[1] = x
        # star image center in pixel coordinates
        ycs = pm.systemImage.shape[0]/2
        xcs = pm.systemImage.shape[1]/2
        step = 2 / (0.75 * ycs)
#        isPlanet = bool(pm.planetBrightness == 0)
        planetInFront = bool(angleDeg < 180)
        angle = angleDeg * np.pi/180.
        # x, y position of the planet center viewed at the inclination angle in
        # coordinate system where RStar = 1.
        xcp = -pm.orbitalRadius * np.cos(angle)
        ycp = pm.orbitalRadius * np.sin(angle) * np.tan(pm.inclination*np.pi/180)
        # transform to pixel coordinates in the system image
        ixShift = int(ycp/step - pm.planetImage.shape[0]/2 + pm.starImage.shape[0]/2)
        iyShift = int(xcp/step - pm.planetImage.shape[1]/2 + pm.starImage.shape[1]/2)
        # over-write the old system image with the star image
        pm.systemImage = pm.starImage
        # add the planet
        ixpRng = np.arange(0,pm.planetImage.shape[0],dtype=int)
        iypRng = np.arange(0,pm.planetImage.shape[1],dtype=int)
        cnt = 0
        for ixp in ixpRng:
            ixs = int(ixShift+ixp)
            if ixs < 0 or ixs >= pm.systemImage.shape[0]:
                continue
            for iyp in iypRng:
                iys = int(iyShift+iyp)
                if iys < 0 or iys >= pm.systemImage.shape[1]:
                    continue
                if pm.planetImage[ixp][iyp] > 0:
                    # this assumes planetInFront
                    if pm.systemImage[ixs][iys] == 0:
                        pm.systemImage[ixs][iys]=pm.planetImage[ixp][iyp]
                    elif planetInFront:
                        pm.systemImage[ixs][iys]=pm.planetImage[ixp][iyp]
                    cnt += 1
    def run(self):
        # Validate the inputs, run the model and update the screen with the results.
        # Validate the user inputs and put them in the system model
        if not self.validateInputs(): return
        # update the LDF plot
        self.plotLDF()
        self.updateStarImage(pm)
        self.updatePlanetImage(pm)
        self.updateSystemImage(pm, 0)
        self.F2Cartoon.imshow(pm.systemImage, cmap="gray", vmin=0, vmax=255)
        self.UpdateRunInfo()
        self.canvas2.draw()
        self.canvas.flush_events()
    def updateTransitCurve(self, values):
        # update the transit curve
        self.line2.set_xdata(self.angles)
        self.line2.set_ydata(values)
        minLum = 1000.
        maxLum = 0.
        for value in values:
            if (value < minLum): minLum = value
            if (value > maxLum): maxLum = value
        if (maxLum == minLum): maxLum = minLum + 0.01
        yRange = maxLum - minLum
        loY = minLum - 0.1*yRange
        hiY = maxLum + 0.1*yRange
        self.TFigSubPlot.set_ylim(loY, hiY)
        self.TFigSubPlot.xaxis.set_major_locator(MultipleLocator(30))
        # I don't know how to update text and text position so just make all old
        # text objects invisible
        for txt in self.TFigSubPlot.texts: txt.set_visible(False)
        # update the luminosity drop at theta = 90 degrees
        if (maxLum > minLum):
            frac =  yRange / maxLum
            text = "dip = {:.4f} ".format(frac)
            textY = loY + 0.1 * (hiY - loY)
            self.TFigSubPlot.text(130,textY,text,fontsize=12)
    def plotLDF(self):
        # Get the entries in the text fields and validate
        self.validateInputs()
        ldf = float(self.ldfEntry.get())
        ra = np.linspace(0, 1, 50)
        ya = np.zeros(50)
        indx = 0
        for r in ra:
            r2 = r * r
            ya[indx] = limbDarkFactor(pm, r2)
            indx += 1
        self.refreshLDF(ra,ya)
    def zoomTransit(self):
        pm.zoomState = (pm.zoomState+1) % 3
        if pm.zoomState == 0:
            self.zoomText.set("Zoom None")
        elif pm.zoomState ==1:
            self.zoomText.set("Zoom Primary")
        elif pm.zoomState ==2:
            self.zoomText.set("Zoom 2ndry")
        else:
            self.zoomText.set("unknown")
    def validateInputs(self):
        # returns true if the user inputs are valid.
        # Assume that the entries that were changed don't require re-rerunning the transit model
        valueChanged = False
        # *****************************
        # Star radius. First ensure that the field is float-like
        try:
            srad = float(self.starRadEntry.get())
        except:
            self.starRadEntry.config(bg="red")
            return False
        # next ensure that the radius value is valid
        if (srad < 0): 
            self.starRadEntry.config(bg="red")
            return False
        else:
            self.starRadEntry.config(bg="white")
        if srad != pm.starRadiusOld: 
            pm.starImageOK = False
            pm.starRadiusOld = srad
        pm.starRadius = srad
        # *****************************
        # limb-darkening factor
        try:
            ldf = float(self.ldfEntry.get())
        except:
            self.ldfEntry.config(bg="red")
            return False
        # next ensure that the LDF value is valid
        if (ldf < 0 or ldf > 1): 
            self.ldfEntry.config(bg="red")
            return False
        else:
            self.ldfEntry.config(bg="white")
        if ldf != pm.LDF:
            pm.LDFOld = ldf
            pm.starImageOK = False
            pm.planetImageOK = False
            pm.systemImageOK = False
        pm.LDF = ldf
        # *****************************
        # planet radius
        try:
            planRad = float(self.planRadEntry.get())
        except:
            self.planRadEntry.config(bg="red")
            return False
        # orbital radius
        try:
            orbRad = float(self.orbRadEntry.get())
        except:
            self.orbRadEntry.config(bg="red")
            return False
        self.orbRadEntry.config(bg="white")
        # scale from AU to solar radius.
        orbRad *= AU_ToRSun
        if orbRad != pm.orbitalRadiusOld:
            pm.orbitalRadiusOld = orbRad
            pm.systemImageOK = False
        pm.orbitalRadius = orbRad
        self.planRadEntry.config(bg="white")
        planRad *= RJupiterToRSun
        if planRad != pm.planetRadius: valueChanged = True
        pm.planetRadius = planRad
        # *****************************
        # planet brightness
        try:
            planBright = float(self.planBrightEntry.get())
        except:
            self.planBrightEntry.config(bg="red")
            return False
        if (planBright < 0):
            self.planBrightEntry.config(bg="red")
            return False
        self.planBrightEntry.config(bg="white")
        if planBright != pm.planetBrightness: valueChanged = True
        pm.planetBrightness = planBright
        # *****************************
        # inclination angle
        try:
            inclin = float(self.planInclEntry.get())
        except:
            self.planInclEntry.config(bg="red")
            return False
        if (inclin < 0 or inclin > 90):
            self.planInclEntry.config(bg="red")
            return False
        maxInclination = 180*math.atan((1+pm.planetRadius)/pm.orbitalRadius)/np.pi
        if (inclin > maxInclination):
            mess = "The inclination angle should be less than {:.0f} degrees to see a transit".format(maxInclination)
            self.planInclEntry.config(bg="red")
            return False
        self.planInclEntry.config(bg="white")
        if inclin != pm.inclination: valueChanged = True
        pm.inclination = inclin
        # *****************************
        # noise
        try:
            noise = float(self.noiseEntry.get())
        except:
            self.noiseEntry.config(bg="red")
            return False
        if (noise < 0 or noise > 0.005):
            self.noiseEntry.config(bg="red")
            return False
        self.noiseEntry.config(bg="white")
        pm.noise = noise
        if (pm.starLuminosity < 0): 
            pm.needsUpdate = True
        else:
            pm.needsUpdate = valueChanged
        return True

if __name__ == "__main__":
    MainWindow = App_Window(None)
    MainWindow.mainloop()

