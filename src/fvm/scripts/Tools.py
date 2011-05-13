from numpy import *

def findMax(array):
    result = 0.0
    min = array.min(axis=0)
    max = array.max(axis=0)    
    if array.ndim == 1:
        if fabs(min) > fabs(max):
            result = fabs(min)
        else:
            result = fabs(max)
    elif array.ndim == 2:
        if fabs(min[2]) > fabs(max[2]):
            result = fabs(min[2])
        else:
            result = fabs(max[2])

    return result


def calculateVoltage(vol, epsilon0, thickness0, epsilon1, thickness1):
    #epsilon0 is air permittivity
    #thickness0 is air gap
    #epsilon1 is dielectric permittivity
    #thickness1 is dielectric thickness
    if thickness1 > 1e-9:
        alpha0 = epsilon0/thickness0
        alpha1 = epsilon1/thickness1
        bias = vol * alpha1 / (alpha0+alpha1)
        return bias
    else:
        return vol
