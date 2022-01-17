#
# All of the transit code lives here
#
import numpy as np

from systemModel import SystemModel as pm

def limbDarkFactor(pm, radialDistance2):
    # calculate the limb-darkening factor at given radial distance^2 with a given limbDark factor.
    if (radialDistance2 < 0 or radialDistance2 > 1): return 0
    arg = 1 - radialDistance2
    arg = 1 - pm.LDF * (1 - np.sqrt(arg))
    return arg

def setLuminosities(pm):
    # fudge the planet position and brightness to determine the star
    # luminsotity
    pm.starLuminosity = -1
    pm.planetLuminosity = -1

    # first the star
    pm.maskStar = False
    planLumSav = pm.planetBrightness
    pm.planetBrightness = 0
    # put the planet at 0 degrees in it's orbit to ensure there is no transit
    pm.starLuminosity = transitLuminosity(pm, 0)
    pm.planetBrightness = planLumSav
    
    # now do the planet
    # mask off the star light 
    pm.maskStar = True
    # Find the luminosity at 0 degrees
    pm.planetLuminosity = transitLuminosity(pm, 0)
    pm.maskStar = False
    
def transitLuminosity(pm, angleDeg):
    # Sum the light of the star + planet (or secondary star) when the planet center is
    # at the specified angle (degrees). 

    # First section: The light from the star is summed over a grid of points within the star disk, 
    # weighted by limb-darkening, and excluding points that are occluded by the "planet" disk. If the
    # "planet" is a true planet (planetBrightness = 0), then all the light from the transit has been
    # summed.
    #
    # Seccond section: This section is only needed if the "planet" is actually a secondary star in a
    # binary star system (planetBrightness == 0). When the "planet" is in front of the star,
    # the light from the "planet" is summed. When the "planet" is not in front, the light from any
    # point on the "planet" disk may come from the "planet" or the star
    #
    # This method uses starLuminosity and planetLuminosity but it is also used to define those
    # variables, so set a boolean that these variables are defined
    isPlanet = bool(pm.planetBrightness == 0)
    planetInFront = bool(angleDeg <= 180)
    angle = angleDeg * np.pi/180.
    # x, y position of the planet center viewed at the inclination angle
    xp = pm.orbitalRadius * np.cos(angle)
    yp = pm.orbitalRadius * np.sin(angle) * np.tan(pm.inclination*np.pi/180)
    # distance from the star as viewed
    rp2 = xp*xp + yp*yp
    # scale the planet radius, which is in units of Jupiter radius, to solar radius where
    # R star = 1
    planRad = 0.1 * pm.planetRadius 
    # and then normalize to the star radius
    planRad /= pm.starRadius
    # planet radius^2
    pr2 = planRad * planRad
    if (pm.starLuminosity > 0 and pm.planetLuminosity >= 0):
        transitRadius2 = (1 + planRad) * (1 + planRad) 
        # return with the star + planet Luminosity if this is NOT a transit
        if (rp2 > transitRadius2): return pm.starLuminosity + pm.planetLuminosity
    step = 0.01
    if(planRad < 0.08): step = 0.005
    # iterate over the full star surface in a rectangular grid of points with size step. The
    # variable sdl (star disk luminosity) is the sum of all star or planet Luminosity points
    #  within the star disk
    sdl = 0
    cnts = 0
    # xs is the x position of a grid point within the star disk (-1 < xs < 1). Only those points that
    # are inside the star disk and outside the planet disk are considered in this section
    if (not pm.maskStar):
        xs = -1
        while (xs < 1):
            ymax = np.sqrt(1 - xs*xs)
            # ys is the y position of a grid point in the star disk at xs (-ymax < ys < +ymax)
            ys = -ymax
            while (ys < ymax):
                # find distance^2 of this point to the planet center. 
                dx = xs - xp
                dy = ys - yp
                # distance^2 from planet center
                rp2 = dx*dx + dy*dy
                # distance^2 from star center
                rs2 = xs*xs + ys*ys
                cnts += 1
                # See if this point is outside the planet disk
                if (rp2 > pr2):
                    # This point is outside the planet disk.
                    sdl += limbDarkFactor(pm, rs2)
                else:
                    # This point is inside the planet disk
                    if (planetInFront):
                        # normalize the planet radius^2 to 1 and planet brightness
                        sdl += limbDarkFactor(pm, rp2/pr2) * pm.planetBrightness
                    else:
                        sdl += limbDarkFactor(pm, rs2)
                ys += step
            xs += step
        # We are done if the planet is indeed dark. Now add the light from a non-dark "planet", aka
        # a secondary star
        # normalize to the number of grid points
        sdl /= cnts
        if (isPlanet):
            return sdl
    # sum the planet disk Luminosity for points that are outside the star disk
    pdl = 0
    cntp= 0
    # Iterate over the grid of points in the "planet" disk
    x = -planRad
    while (x < planRad):
        ymax = np.sqrt(pr2 - x*x)
        y = -ymax
        while (y < ymax):
            # Find the distance^2 of this (x,y) point to the star center
            dx = x + xp
            dy = y + yp
            rs2 = dx*dx + dy*dy
            cntp += 1
            # see if this point is outside the star disk. 
            if (rs2 > 1): 
                # find the distance^2 to the planet center and apply limb-darkening
                rp2 = x*x + y*y
                # This point is outside the star disk so sum the "planet" light
                # normalize the planet radius^2 to 1
                pdl += limbDarkFactor(pm, rp2/pr2) * pm.planetBrightness
            y += step
        x += step
    # normalize by planet disk area / star disk area
    if (cntp > 0): pdl *= pr2 /cntp
    return (sdl + pdl)

def doOrbit(pm, angles, intensities):
    # The star is located at (0,0,0). The planet moves in a circular orbit in the Z = 0 plane and
    # is viewed by us by a 2D coordinate system which is a projection of the 3D system rotated by
    # the inclination angle around the Y axis. 
    indx = 0
    for angle in angles:
        intensities[indx]=transitLuminosity(pm, angle)
        indx += 1
