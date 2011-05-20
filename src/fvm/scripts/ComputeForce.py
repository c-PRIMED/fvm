from numpy import *
from math import *
import pdb

### constants ###
epsilon0 = 8.854187817e-12
H = 0.23e-20
B = 3529e3
alpha = 0.1127
gamma = 22.69e9
alpha01 = 1.6e-9
alpha02 = 1.99e-9

########################################################################
def computeElectrostaticForce(voltage, distance):
    # assume the distance is not smaller than 1nm
    if distance < 1e-9:
        distance = 1e-9
    ef = voltage/distance
    elecForce = 0.5 * ef * ef * epsilon0
    return elecForce

########################################################################
def computeContactForce(distance):
    #if distance < 1e-10:
    #    distance = 1e-10
    repulsive_force =  B*exp(-(distance-alpha02)*gamma) 
    attractive_force = -H/(6*pi) * ((1-alpha)/pow(distance,3) + \
                         alpha/pow(distance-alpha01,3))
    contactForce = repulsive_force #+ attractive_force 
    return contactForce

########################################################################
def computeDampForce(velocity, coeff):
    dampForce = -coeff * velocity
    return dampForce

########################################################################
def computeDampCoeff(gap, width, thickness, pressure=1.01325e5):
    A = 10.39
    B = 1.374
    c = 3.100
    d = 1.825
    e = 0.966

    temperature = 300
    k = 1.38e-23 
    diameter = 4e-10
    area = pi * diameter * diameter
    gama = k*temperature/area
    mfp = gama / pressure

    x1 = width / gap
    x2 = mfp / width
    Cf = A*pow(x1,c) * thickness / (1+B*pow(x1,d)*pow(x2,e))

    return Cf
