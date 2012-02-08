import random
from math import *

H = 0.215e-18
c1 = 1.2e3
c2 = 0.5e9
c3 = 1.0e3
dc = 110e-9
tR = 8.5e-9

H_min = 0.18e-18
H_max = 0.25e-18
dc_min = 50e-9
dc_max = 200e-9
tR_min = 5e-9
tR_max = 10e-9
c1_min = 0.6e3
c1_max = 1.5e3
c2_min = 0.3e9
c2_max = 15e9
c3_min = 0.5e3
c3_max = 2e3

ds = 30e-9
d_bar = 1e-9

def func(x):
    dx = 0
    if x > 0:
        dx = 0
    else:
        dx = -x
    return dx

def generateRandomNumber(N):
    rn = []
    for i in range(0,N):
        rn.append(random.random())
    return rn

def computeContactForceOnDemand(x,y,lx,ly,d, i, rn):
    #if distance < 1e-10:
    #    distance = 1e-10
    H = (H_max-H_min)*rn[i] + H_min
    dc = (dc_max-dc_min)*rn[i] + dc_min
    tR = (tR_max-tR_min)*rn[i] + tR_min
    c1 = (c1_max-c1_min)*rn[i] + c1_min
    c2 = (c2_max-c2_min)*rn[i] + c2_min
    c3 = (c3_max-c3_min)*rn[i] + c3_min

    attractive_force = -H/(6*pi*(pow(d,3)+pow(ds,3))) -c3*func(d-(dc+tR))/d_bar
    repulsive_force = c1*pow((func(d-dc)/d_bar),1.5) + c2*exp(-d/tR)
    contact_force = attractive_force + repulsive_force

    return contact_force
