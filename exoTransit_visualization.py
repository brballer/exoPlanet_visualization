# A model for visualizing the light output from a star + planet or an eclipsing binary system of two stars. 

import sys
import math
from tkinter.constants import ANCHOR
import numpy as np
import tkinter
from tkinter import messagebox
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg)
from matplotlib.ticker import (MultipleLocator)
import matplotlib.artist as martist

from transitCode import *
from systemModel import SystemModel

# default values for the planet model
starRad = 1.0 # solar radii
orbRad = 0.01 # AU
planRad = 1.0 # Jupiter radii
# this should be 0 or the starting cartoon grayscale will be off
planLum = 0.
inclin = 10
ldf = 0.4
noise = 0.0
# create the system model
pm = SystemModel(
    starRad, # star radius in units of Solar radius
    orbRad, # planet orbital radius in units of AU
    planRad, # planet radius in units of Jupiter radius
    planLum, # planet luminosity(0 = planet, > 1 eclipsing binary star)
    inclin, # inclination view angle (degrees) 
    ldf # star limb-darkening factor (0 = none)
    )
# normalize
pm.orbitalRadius *= pm.starRadius / 0.00465

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
        self.planLumEntry = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=planLum), width = 8)
        self.lbl4.place(x=xpo1, y=ypo)
        self.planLumEntry.place(x=xpo2, y=ypo)
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
        self.plotbutton=tkinter.Button(self, text="Update LDF Brightness Plot", command=self.plotLDF)
        self.plotbutton.place(x=xpo1+50, y=ypo)
        self.b3=tkinter.Button(self, text='Quit', command=self.appQuit)
        self.b3.pack(anchor="nw")
        # add some informational text fields
        ypo += 30
        self.starLumTxt = tkinter.Label(self, text = "Star luminosity = UNDEFINED")
        self.starLumTxt.place(x=xpo1, y=ypo)
        ypo += 20
        self.planLumEntryTxt = tkinter.Label(self, text = "Planet luminosity = UNDEFINED")
        self.planLumEntryTxt.place(x=xpo1, y=ypo)
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
        # define a blank transit figure object
        # now make the transit figure
        self.TransitFig = plt.Figure(figsize=(30,30),dpi=100)
        # add one subplot for the cartoon
        CFigSubPlot = self.TransitFig.add_subplot(211, facecolor="black")
        CFigSubPlot.set_aspect(1)
        CFigSubPlot.set_xlim(-35, 35)
        CFigSubPlot.set_ylim(-35, 35)
        CFigSubPlot.xaxis.set_visible(False)
        CFigSubPlot.yaxis.set_visible(False)
        CFigSubPlot.set_title("Occultation at Orbit Angle = 90 degrees")
        # scaling for the sun and planet90
        self.CFigScale = 15
        # first draw the star disk. Set the color less than white to allow distinguising
        # between a brighter eclipsing binary star
        self.starColorScale = 0.7
        self.star = plt.Circle((0,0), self.CFigScale, color="{:.2f}".format(self.starColorScale))
        CFigSubPlot.add_patch(self.star)
        # then the planet disk
        yp = self.CFigScale * pm.orbitalRadius * np.tan(pm.inclination*np.pi/180)
        # initialize with a real planet (facecolor=0) but give it a bright edge for visibility.
        # the star is scaled to the radius of the sun (= 1), so scale the cartoon planet
        plRad = self.CFigScale * 0.1 * pm.planetRadius / pm.starRadius
        self.planet90 = plt.Circle((0,yp),plRad, facecolor="0", edgecolor="1")
        CFigSubPlot.add_patch(self.planet90)
        # add another subplot for the transit curve
        self.TFigSubPlot = self.TransitFig.add_subplot(212)
        self.TFigSubPlot.xaxis.set_major_locator(MultipleLocator(30))
        self.TFigSubPlot.set_xlabel("Orbit Angle (degrees)")
        self.TFigSubPlot.set_ylabel("Flux")
        angla = np.linspace(0, 360, num=60, endpoint=True)
        intena = np.zeros(60)
        self.line2, = self.TFigSubPlot.plot(angla,intena,"o")
        self.canvas2 = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.TransitFig, master=self)
        self.canvas2.get_tk_widget().pack(side=tkinter.TOP, expand=True)
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
        self.planLumEntryTxt.configure(text="Planet luminosity = {:.3f}".format(pm.planetLuminosity))
        self.ang0LumTxt.configure(text="System luminosity(0 deg) = {:.4f}".format(self.luminosities[0]))
        self.ang90LumTxt.configure(text="System luminosity(90 deg) = {:.4f}".format(self.luminosities[14]))
        self.ang270LumTxt.configure(text="System luminosity(270 deg) = {:.4f}".format(self.luminosities[44]))
        # the estimated planet radius only makes sense if the planet disk was inside the star disk
        # at max occultation
        yp = pm.orbitalRadius * np.tan(pm.inclination*np.pi/180) + 0.1 * pm.planetRadius
        pre = 9.73 * pm.starRadius * np.sqrt(self.luminosities[0] - self.luminosities[14])
        pre *= 1.06 # correction for grid size error
        if (pm.planetLuminosity > 0):
            self.planRadEntryEstTxt.configure(text="Can't estimate star companion radius")
        elif (yp < pm.starRadius):
            self.planRadEntryEstTxt.configure(text="Estimated planet radius = {:.1f} * R_Jupiter using R_planet = 9.73 * R_star * \u221adip".format(pre))
        else:
            self.planRadEntryEstTxt.configure(text="Planet radius cannot be estimated because it is not fully within the star disk")
    def run(self):
        # Validate the inputs, run the model and update the screen with the results.
        # Validate the user inputs and put them in the system model
        if not self.validateInputs(False): return
        # update the LDF plot
        self.plotLDF()
        # Find the luminosities of the unoccluded star and unoccluded planet
        setLuminosities(pm)
        if (pm.starLuminosity < 0 or pm.planetLuminosity < 0):
            messagebox.showerror("ooops starLuminosity {:.2f}".format(pm.starLuminosity), \
                " planetLuminosity {:.2f}".format(pm.planetLuminosity))
            return
        # draw the star-planet cartoon before running doOrbit
        self.drawCartoon()
        # Move the planet around the star and fill the angle array and luminosities array
        if pm.needsUpdate:
            doOrbit(pm, self.angles, self.luminosities)
            # normalize
            norm = self.luminosities[0]
            pm.planetLuminosity /= norm
            pm.needsUpdate = False
        # update the run info text
        self.UpdateRunInfo()
        angla = np.linspace(0, 360, num=60, endpoint=True)
        intena = np.zeros(60)
        if(pm.noise > 0) :
            # create a noise array
            noisyY = self.luminosities + np.random.normal(0, pm.noise, 60)
            self.updateTransitCurve(noisyY)
        else:
            self.updateTransitCurve(self.luminosities)
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
        self.validateInputs(True)
        ldf = float(self.ldfEntry.get())
        ra = np.linspace(0, 1, 50)
        ya = np.zeros(50)
        indx = 0
        for r in ra:
            r2 = r * r
            ya[indx] = limbDarkFactor(pm, r2)
            indx += 1
        self.refreshLDF(ra,ya)
    def validateInputs(self, justLDF):
        # returns true if the user inputs are valid.
        # Assume that the entries that were changed don't require re-rerunning the transit model
        valueChanged = False
        # limb-darkening factor. First ensure that the field is float-like
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
        if srad != pm.starRadius: valueChanged = True
        pm.starRadius = srad
        # limb-darkening factor. First ensure that the field is float-like
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
        if ldf != pm.LDF: valueChanged = True
        pm.LDF = ldf
        # only check the LDF Entry?
        if justLDF: return True
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
        # scale from AU to star radius. RSun = 0.00465 AU
        orbRad *= pm.starRadius/0.00465
        if (orbRad < starRad + 0.1*planRad):
            self.orbRadEntry.config(bg="red")
            return False
        if orbRad != pm.orbitalRadius: valueChanged = True
        pm.orbitalRadius = orbRad
        # Require the planet radius be greater than 0.2 Jupiter radius
        if (planRad < 0.2):
            self.planRadEntry.config(bg="red")
            return False
        self.planRadEntry.config(bg="white")
        if planRad != pm.planetRadius: valueChanged = True
        pm.planetRadius = planRad
        # planet luminosity
        try:
            planLum = float(self.planLumEntry.get())
        except:
            self.planLumEntry.config(bg="red")
            return False
        if (planLum < 0):
            self.planLumEntry.config(bg="red")
            return False
        self.planLumEntry.config(bg="white")
        if planLum != pm.planetBrightness: valueChanged = True
        pm.planetBrightness = planLum
        # inclination angle
        try:
            inclin = float(self.planInclEntry.get())
        except:
            self.planInclEntry.config(bg="red")
            return False
        if (inclin < 0 or inclin > 90):
            self.planInclEntry.config(bg="red")
            return False
        maxInclination = 180*math.atan((1+0.1*planRad)/orbRad)/np.pi
        if (inclin > maxInclination):
            mess = "The inclination angle should be less than {:.0f} degrees to see a transit".format(maxInclination)
            self.planInclEntry.config(bg="red")
            return False
        self.planInclEntry.config(bg="white")
        if inclin != pm.inclination: valueChanged = True
        pm.inclination = inclin
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

