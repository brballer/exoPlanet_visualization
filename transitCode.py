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
