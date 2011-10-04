from numpy import *
from math import *
import pdb

### constants ###
epsilon0 = 8.854187817e-12
"""
old parameters
H = 0.23e-20
B = 3529e3
alpha = 0.1127
gamma = 22.69e9
alpha01 = 1.6e-9
alpha02 = 1.99e-9
"""
H = 0.215e-18
B = 0.738362e9
xx = 0.1
gamma = 6.35225
x01 = 60e-9
x02 = 42.1714e-9
alpha = 2e-9

def func(x):
    dx = 0
    if x > 0:
        dx = pow((x*x*x+alpha*alpha*alpha), 1./3.)
    else:
        dx = alpha
    return dx
########################################################################
def computeElectrostaticForce(voltage, distance, stop=100e-9):
    # assume the distance is not smaller than 1nm
    if distance < stop:
        distance = stop
    ef = voltage/distance
    elecForce = 0.5 * ef * ef * epsilon0
    return elecForce

########################################################################
def computeContactForce(x):
    #if distance < 1e-10:
    #    distance = 1e-10
    #repulsive_force =  B*exp(-(distance-alpha02)*gamma) 
    #attractive_force = -H/(6*pi) * ((1-alpha)/pow(distance,3) + \
    #                     alpha/pow(distance-alpha01,3))
    #contactForce = repulsive_force #+ attractive_force 
    attractive_force = -H/(6*pi) * ((1-xx)/pow(func(x),3) + xx/pow(func(x-x01),3))
    repulsive_force = B*exp(-func(x-x02)/x02*gamma)
    contact_force = attractive_force + repulsive_force
    return contact_force

########################################################################
def computeDampForce(velocity, coeff):
    dampForce = -coeff * velocity
    return dampForce

########################################################################
def computeDampCoeff(gap, width, thickness, pressure=1.01325e5, stop=100e-9):
    A = 10.39
    B = 1.374
    c = 3.100
    d = 1.825
    e = 0.966
    if gap < stop:
    	gap = stop
    temperature = 300
    k = 1.3806503e-23 
    diameter = 4.17e-10
    area = pi * diameter * diameter
    gama = k*temperature/area
    mfp = gama / pressure
    Tref = 273.15
    omega = 0.81

    Kn = k * temperature / (gap*sqrt(2)*area*pressure*pow(Tref/temperature, omega-0.5))

    x1 = width/gap
    x2 = Kn/x1
    Cf = A*pow(x1,c) * thickness / (1+B*pow(x1,d)*pow(x2,e))

    return Cf
