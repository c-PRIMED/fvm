from math import *
from numpy import *

########################################################################
def computeTimeStepPrep(gap, A0, AN):
    beta = 1 - (A0-AN)/gap
    NP = int(log(AN/A0)/log(beta))
    A = zeros(NP)
    R = zeros(NP)
    sum = 0.0
    A[0] = A0
    sum = A0
    R[0] = sum
    for i in range(1, NP):
        A[i] = A[i-1] * beta
        sum = sum + A[i]
        R[i] = sum
    return R

########################################################################
def computeTravelDistance(distance, gap, rmin=0.1e-9, rmax=50e-9):
    #default: mininum displacement is 0.1 nm
    #default: maximum displacement is 50 nm
    R = computeTimeStepPrep(gap, rmin, rmax)
    NP = R.size
    # find the current position in R array
    found = False
    if distance > R[NP-1]:
        Rinit = gap
        Rend = R[NP-1]
        found = True
    elif distance < R[0]:
        Rinit = R[0]
        Rend = 0.01e-9
        found = True
    else:
        left = 0
        right = NP-1        
        while found == False:
            mid = int((left+right)*0.5)
            if distance > R[mid]:
                left = mid
            if distance <= R[mid]:
                right = mid
            if right-left == 1:
                Rinit = R[right]
                Rend = R[left]
                found = True
    
    deltaR = fabs(Rinit-Rend)
    return deltaR
########################################################################    
def computeTimeStep(dr, vel, acc, rmin=0.1e-9, rmax=50e-9):
    timeStep = 0
    if fabs(acc) > 1e-10:
    	#if vel<0:
    	#     dr = -dr
    	discr = vel*vel + 2*acc*dr
    	if discr > 0:
    	     discr = sqrt(discr)
    	     dt1 = (-vel + discr) / acc
    	     dt2 = (-vel - discr) / acc
    	     	
    	     if dt1>0 and dt2>0:
    	     	timeStep = min(dt1, dt2)
    	     elif dt1 > 0:
    	     	timeStep = dt1
    	     elif dt2 > 0:
    	     	timeStep = dt2
    	     else:
    	     	timeStep = -1
    	        	     
        else:
             timeStep = fabs(2*vel/acc)
    else:
        timeStep = -1    
    return timeStep

########################################################################

