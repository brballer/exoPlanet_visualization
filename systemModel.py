
import math

class SystemModel():
    def __init__(self, orbitalRadius, planetRadius, planetBrightness, inclination, ldf):
        # class constructor of a star and single planet system.
        # The star has a radius of 1 and a central brightness of 1. The integrated star luminosity is reduced if a
        # limb-darkening factor LDF is > 0. A "planet" in the code is a true planet if the planetBrightness = 0.
        # Eclipsing binary stars are modeled by setting planetBrightness > 0. 
        self.isOK = True
        self.planetRadius = planetRadius
        self.orbitalRadius = orbitalRadius
        # brightness is the amount of light emitted per unit area
        self.planetBrightness = planetBrightness
        # limb-darkening factor
        self.LDF = ldf
        self.inclination = inclination
        # Luminosity is the total amount of light emitted. These can't be determined until the
        # system model is defined by the user.
        self.starLuminosity = -1
        self.planetLuminosity = -1
        # mask off the star light to measure the planet luminosity
        self.maskStar = False
