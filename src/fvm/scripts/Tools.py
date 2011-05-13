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
