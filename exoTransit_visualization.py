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
orbRad = 2.0
planRad = 0.5
# this should be 0 or the starting cartoon grayscale will be off
planLum = 0.
inclin = 20
ldf = 0.5
# create the system model
pm = SystemModel(
    orbRad, # planet orbital radius
    planRad, # planet radius relative to star radius
    planLum, # planet luminosity(0 = planet, > 1 eclipsing binary star)
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
        self.title("ExoTransit_Visualization V1.0")
        self.geometry("1000x800")
        # define arrays for the orbital angle and the total luminosity
        self.angles = np.zeros(60)
        self.luminosities = np.zeros(60)
        self.lbl1 = tkinter.Label(self, text = "Limb-darkening factor (LDF)")
        self.lbl11 = tkinter.Label(self, text = "0 < LDF < 1")
        ypo = 30
        xpo1 = 0
        xpo2 = 200
        xpo3 = 300
        self.ldf = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=ldf), width = 8)
#        self.ldf.bind("<Enter>", self.rePlotLDF)
        self.lbl1.place(x=xpo1, y=ypo)
        self.ldf.place(x=xpo2, y=ypo)
        self.lbl11.place(x=xpo3, y=ypo)
        #
        ypo += 25
        self.lbl2 = tkinter.Label(self, text = "Planet orbital radius")
        self.lbl21 = tkinter.Label(self, text = "(> 1")
        self.orbrad = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=orbRad), width = 8)
        self.lbl2.place(x=xpo1, y=ypo)
        self.orbrad.place(x=xpo2, y=ypo)
        self.lbl21.place(x=xpo3, y=ypo)
        #
        ypo += 25
        self.lbl3 = tkinter.Label(self, text = "Planet radius")
        self.lbl31 = tkinter.Label(self, text = "Range 0.02 to 1")
        self.planrad = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=planRad), width = 8)
        self.lbl3.place(x=xpo1, y=ypo)
        self.planrad.place(x=xpo2, y=ypo)
        self.lbl31.place(x=xpo3, y=ypo)
        #
        ypo += 25
        self.lbl4 = tkinter.Label(self, text = "Planet brightness")
        self.lbl41 = tkinter.Label(self, text = "0 for planet, >0 for secondary star")
        self.planlum = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=planLum), width = 8)
        self.lbl4.place(x=xpo1, y=ypo)
        self.planlum.place(x=xpo2, y=ypo)
        self.lbl41.place(x=xpo3, y=ypo)
        #
        ypo += 25
        self.lbl5 = tkinter.Label(self, text = "Planet inclination angle")
        self.lbl51 = tkinter.Label(self, text = "degrees")
        self.planIncl = tkinter.Entry(bd=3, textvariable = tkinter.StringVar(value=inclin), width = 8)
        self.lbl5.place(x=xpo1, y=ypo)
        self.planIncl.place(x=xpo2, y=ypo)
        self.lbl51.place(x=xpo3, y=ypo)
        # Create the Run button and LDF plot button
        ypo += 30
        self.b1=tkinter.Button(self, text='Run', command=self.run)
        self.b1.place(x=xpo1, y=ypo)
        self.plotbutton=tkinter.Button(self, text="Update LDF Brightness Plot", command=self.plotLDF)
        self.plotbutton.place(x=xpo1+50, y=ypo)
        self.b3=tkinter.Button(self, text='Quit', command=self.appQuit)
        self.b3.pack(anchor="nw")
        # add some informational text fields
        ypo += 30
        self.starLumTxt = tkinter.Label(self, text = "Star Luminosity = UNDEFINED")
        self.starLumTxt.place(x=xpo1, y=ypo)
        ypo += 20
        self.planLumTxt = tkinter.Label(self, text = "Planet Luminosity = UNDEFINED")
        self.planLumTxt.place(x=xpo1, y=ypo)
        ypo += 20
        self.ang0LumTxt = tkinter.Label(self, text = "System Luminosity(0 deg) = UNDEFINED")
        self.ang0LumTxt.place(x=xpo1, y=ypo)
        ypo += 20
        self.ang90LumTxt = tkinter.Label(self, text = "System Luminosity(90 deg) = UNDEFINED")
        self.ang90LumTxt.place(x=xpo1, y=ypo)
        ypo += 20
        self.ang270LumTxt = tkinter.Label(self, text = "System Luminosity(270 deg) = UNDEFINED")
        self.ang270LumTxt.place(x=xpo1, y=ypo)
        # add the Brightness vs Radius figure
        Fig = Figure(figsize=(5.5,3.5),dpi=100)
        FigSubPlot = Fig.add_subplot(111)
        FigSubPlot.set_xlabel("Radius")
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
        # initialize with a real planet (facecolor=0) but give it a bright edge for visibility
        self.planet90 = plt.Circle((0,yp),self.CFigScale*pm.planetRadius, facecolor="0", edgecolor="1")
        CFigSubPlot.add_patch(self.planet90)
        # add another subplot for the transit curve
        self.TFigSubPlot = self.TransitFig.add_subplot(212)
        self.TFigSubPlot.xaxis.set_major_locator(MultipleLocator(30))
        self.TFigSubPlot.set_xlabel("Orbit Angle (degrees)")
        self.TFigSubPlot.set_ylabel("Luminosity")
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
        self.planet90.radius = self.CFigScale*pm.planetRadius
        # define the color (actually the gray scale)
        clr = pm.planetBrightness * self.starColorScale
        if (clr > 1): clr = 1
        strClr = "{:.1f}".format(clr)
        self.planet90.set_facecolor(strClr)
        self.canvas.draw()
        self.canvas.flush_events()
    def UpdateRunInfo(self):
        # Updates the somewhat detailed results below the Run button
        self.starLumTxt.configure(text="Star Luminosity = {:.3f}".format(pm.starLuminosity))
        self.planLumTxt.configure(text="Planet Luminosity = {:.3f}".format(pm.planetLuminosity))
        self.ang0LumTxt.configure(text="System Luminosity(0 deg) = {:.3f}".format(self.luminosities[0]))
        self.ang90LumTxt.configure(text="System Luminosity(90 deg) = {:.3f}".format(self.luminosities[14]))
        self.ang270LumTxt.configure(text="System Luminosity(270 deg) = {:.3f}".format(self.luminosities[44]))
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
        doOrbit(pm, self.angles, self.luminosities)
        # update the run info text
        self.UpdateRunInfo()
        # Determine the plot limits for the transit curve
        minLum = 1000
        maxLum = 0
        for luminosity in self.luminosities:
            if (luminosity < minLum): minLum = luminosity
            if (luminosity > maxLum): maxLum = luminosity
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
            frac =  (minLum - maxLum) / maxLum
            text = "Luminosity change {:.3f} ".format(frac)
            textY = loY + 0.1 * (hiY - loY)
            self.TFigSubPlot.text(130,textY,text,fontsize=12)
        angla = np.linspace(0, 360, num=60, endpoint=True)
        intena = np.zeros(60)
        # update the transit curve
        self.line2.set_xdata(self.angles)
        self.line2.set_ydata(self.luminosities)
        self.canvas2.draw()
        self.canvas.flush_events()
    def plotLDF(self):
        # Get the entry in the ldf text field
        self.validateInputs(True)
        ldf = float(self.ldf.get())
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
        # limb-darkening factor. First ensure that the field is float-like
        try:
            ldf = float(self.ldf.get())
        except:
            self.ldf.config(bg="red")
            return False
        # next ensure that the LDF value is valid
        if (ldf < 0 or ldf > 1): 
            self.ldf.config(bg="red")
            return False
        else:
            self.ldf.config(bg="white")
        pm.LDF = ldf
        # only check the LDF Entry?
        if justLDF: return True
        # orbital radius
        try:
            orbRad = float(self.orbrad.get())
        except:
            self.orbrad.config(bg="red")
            return False
        self.orbrad.config(bg="white")
        pm.orbitalRadius = orbRad
        # planet radius
        try:
            planRad = float(self.planrad.get())
        except:
            self.planrad.config(bg="red")
            return False
        if (planRad < 0.02 or planRad > 1):
            self.planrad.config(bg="red")
            return False
        self.planrad.config(bg="white")
        pm.planetRadius = planRad
        # planet luminosity
        try:
            planLum = float(self.planlum.get())
        except:
            self.planlum.config(bg="red")
            return False
        if (planLum < 0):
            self.planlum.config(bg="red")
            return False
        self.planlum.config(bg="white")
        pm.planetBrightness = planLum
        # inclination angle
        try:
            inclin = float(self.planIncl.get())
        except:
            self.planIncl.config(bg="red")
            return False
        if (inclin < 0 or inclin > 90):
            self.planIncl.config(bg="red")
            return False
        minInclination = 180*math.atan((1+planRad)/orbRad)/np.pi
        if (inclin > minInclination):
            mess = "The inclination angle should be less than {:.0f} degrees to see a transit".format(minInclination)
            messagebox.showinfo("Hmmm",mess)
            self.planIncl.config(bg="red")
            return False
        self.planIncl.config(bg="white")
        pm.inclination = inclin
        return True

if __name__ == "__main__":
    MainWindow = App_Window(None)
    MainWindow.mainloop()

