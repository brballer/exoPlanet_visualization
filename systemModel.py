
import math
import numpy as np

class SystemModel():
    def __init__(self, starRadIn, orbRadIn, planRadIn, planBrightIn, inclinIn, ldfIn):
        # class constructor of a star and single planet system.
        # The star has a radius of 1 and a central brightness of 1. The integrated star luminosity is reduced if a
        # limb-darkening factor LDF is > 0. A "planet" in the code is a true planet if the planetBrightness = 0.
        # Eclipsing binary stars are modeled by setting planetBrightness > 0. 
        self.systemImage = np.zeros((256,512), np.uint16)
        self.starImage = np.zeros((256,512), np.uint16)
        self.planetImage = np.zeros((256,256), np.uint16)
        self.starRadius = np.double(starRadIn) # in units of RSun
        self.planRadius = np.double(planRadIn) # in units of RSun
        self.planOrbitRadius = np.double(orbRadIn) # in units of RSun
        self.planBright = planBrightIn
        self.LDF = np.float(ldfIn)
        self.inclination = np.float(inclinIn) # in degrees
        self.diskOuterRadius = 0
        self.diskInnerRadius = 0
        self.diskBright = 0
        self.diskInclination = 0 # degrees
        self.diskOrbitAngle = 0
        self.stepSize = np.double(5.0/512) # The cartoon is 5 R_Sun wide with 512 pixels
        self.luminosityNorm= np.double(-1.)
        self.luminosity90 = np.double(-1.)
        self.nAngleBins = 60
        self.zoomState = 0 # 0 = none, 1 = primary transit, 2 = secondary transit
        self.featureState = 0 # 0 = none, 1 = circumsecondary disk, etc
        # Luminosity is the total amount of light emitted. These can't be determined until the
        # system model is defined by the user.
        self.starLuminosity = np.double(-1.)
        self.planetLuminosity = np.double(-1.)
        self.systemLuminosity = np.double(-1.) # this is the total luminosity at a specifc orbit angle
        self.noise = 0.

