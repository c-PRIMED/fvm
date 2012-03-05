from numpy import *
from math import *
import pdb

### constants ###
epsilon0 = 8.854187817e-12
H = 0.215e-18
c1 = 1.2e3
c2 = 0.5e9
c3 = 1.0e3
dc = 110e-9
ds = 30e-9
d_bar = 1e-9
tR = 8.5e-9
ta = 2e-9

def func(x):
    dx = 0
    if x > 0:
        dx = 0
    else:
        dx = -x
    return dx
########################################################################
def computeElectrostaticForce(voltage, distance, stop=200e-9):
    # assume the distance is not smaller than 1nm
    if distance < stop:
        distance = stop
    ef = voltage/distance
    elecForce = 0.5 * ef * ef * epsilon0
    return elecForce
########################################################################
# voltage: the voltage difference between beam and substrate (including dielectric)
# dc0:  air dielectric constant, i.e. dc0 = 1
# thickness0: the gap between beam bottom and dielectric top surface
# dc1: dielectric constant, i.e. dc1 = 7.9
# thickness1: dielectric thickness, i.e. thickness1 = 200e-9
# src: charge density * QE 
def computeElectrostaticForceWithCharge(voltage, dc0=1.0, thickness0=0., 
                                        dc1=7.9, thickness1=200e-9, src=0., stop=100e-9):
    if thickness0 < stop:
        thickness0 = stop
    alpha0 = thickness0/dc0
    alpha1 = thickness1/dc1
    source = 0.5 * src * thickness1 * thickness1 / (dc1*epsilon0)
    ef = (voltage+source)/(alpha0+alpha1)
    elecForce = 0.5 * ef * ef * epsilon0
    
    return elecForce

########################################################################
def computeContactForce(x):
    #if distance < 1e-10:
    #    distance = 1e-10
    attractive_force = -H/(6*pi*(pow(x,3)+pow(ds,3))) -c3*func(x-(dc+ta))/d_bar
    repulsive_force = c1*pow((func(x-dc)/d_bar),1.5) + c2*exp(-x/tR)
    contact_force = attractive_force + repulsive_force
    return contact_force

########################################################################
def computeDampForce(velocity, coeff):
    dampForce = -coeff * velocity
    return dampForce

########################################################################
def computeDampCoeff_Old(gap, width, thickness, pressure=1.01325e5, stop=200e-9):
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
#############################################################################
def computeDampCoeff(gap, width, thickness, pressure=1.01325e5):
    temperature = 300
    k = 1.3806503e-23 
    dref = 4.17e-10
    Tref = 273.15
    omega = 0.74
    mfp = k*temperature /(sqrt(2.0)*pi*dref*dref*pressure*pow(Tref/temperature, omega-0.5))
    Kn = mfp/gap

    x_cut =  0.015
    mu    =  1.98e-5
    c1    =  5.5964
    c2    =  (48.62*Kn)/(1+0.3405*Kn)
    c3    =  (1+210.79*Kn)/(67.59*Kn)

    x1    =  width/gap
    x1_sqr = x1*x1
    x1_cube = x1*x1*x1

    if (1.0/x1 >= x_cut):
        Cf = mu * x1_cube / (1 + c1*6.0*Kn) * (1+6.0*c2/x1+12.0*c3/x1_sqr)
    else:
        A1 = mu / (1+c1*6.0*Kn) * (-3.0/pow(x_cut, 4)-12.0*c2/pow(x_cut, 3)-12.0*c3/pow(x_cut, 2))
        A2 = mu / (1+c1*6.0*Kn) * (1/pow(x_cut,3) + 6.0*c2/pow(x_cut,2) + 12.0*c3/x_cut )-x_cut*A1
        Cf = A1/x1 + A2
    return Cf
