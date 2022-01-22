
import math
import numpy as np

class SystemModel():
    def __init__(self, starRadIn, orbRadIn, planRadIn, planBrightIn, inclinIn, ldfIn):
        # class constructor of a star and single planet system.
        # The star has a radius of 1 and a central brightness of 1. The integrated star luminosity is reduced if a
        # limb-darkening factor LDF is > 0. A "planet" in the code is a true planet if the planetBrightness = 0.
        # Eclipsing binary stars are modeled by setting planetBrightness > 0. 
        self.systemImage = np.zeros((512,1024), np.uint8)
        self.starImage = np.zeros((512,1024), np.uint8)
        self.planetImage = np.zeros((512,512), np.uint8)
        self.starRadius = starRadIn # in units of RSun
        self.starRadiusOld = -1.
        self.planetRadius = planRadIn # in units of RSun
        self.planetRadiusOld = -1.
        self.orbitalRadius = orbRadIn # in units of RSun
        self.orbitalRadiusOld = -1.
        # brightness is the amount of light emitted per unit area
        self.planetBrightness = planBrightIn
        self.planetBrightnessOld = -1
        # limb-darkening factor
        self.LDF = ldfIn
        self.LDFOld = -1
        self.inclination = inclinIn # in degrees
        self.inclinationOld = -1
        self.stepSize = 0.00488281 # for mapping cartesian (x,y) to pixel
        self.zoomState = 0 # 0 = none, 1 = primary transity, 2 = secondary transit
        # Luminosity is the total amount of light emitted. These can't be determined until the
        # system model is defined by the user.
        self.starLuminosity = -1
        self.planetLuminosity = -1
        self.noise = 0.
        self.starImageOK = False # set true if the image needs updating
        self.planetImageOK = False
        self.systemImageOK = False
        self.F2CartoonOK = False

