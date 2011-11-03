#!/usr/bin/env python

################################
#Importing the relevant modules#
################################

import sys
import fvm 
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.exporters_atyped_double as exporters
import time
from numpy import *
import fvm.phonon_atyped_double as pa
import fvm.phononbaseExt as pext
from FluentCase import FluentCase
fvm.set_atype('double')

print " "


def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat if it is not specified."
    sys.exit(1)

############################################################
#Read the mesh, and extract the information into the MEMOSA#
#data structures                                           # 
############################################################
print "1"
reader = FluentCase("/home/brad/memosa/src/fvm/src/modules/phononbase/testmesh/Toy2.cas")
reader.read();
print "2"
meshes = reader.getMeshList()
print "3"
geomFields =  fvm.models.GeomFields('geom')
metricsCalculator = fvm.models.MeshMetricsCalculatorA(geomFields,meshes)
metricsCalculator.init()

###################################################################
#Specify the group velocity and the scattering rate.  Since       #
#the only parameter that matters is the Knudsen number (=vg*tau/L)#
#they needn't be physical values, so long as the desired Knudsen  #
#is obtained.  It is also assumed that the Brillioun zone is      #
#isotropic and the group velocity is the same direction as the    #
#K-space vector.  Because of this, to discretize k-space, you     #
#simply specify ntheta (number of discretizations in the theta    #
#direction) and nphi (number of discretizations in the phi        #
#direction).  For more explanation, see Chai, Lee, and Patankar   #
#JTHT.                                                            #
###################################################################

vg=.01
tau=100.
ntheta=4
nphi=8

                  #1 is for the
                  #number of
                  #polarizations
K_space=pa.KspaceA(1,tau,vg,4e12,ntheta,nphi)
                          #4e12 is the phonon
                          #frequency.  Don't
                          #worry about it.

pmacro=pext.PhononMacro("e_dubprime")

#create the model
pmodel=pa.PhononModelA(meshes,geomFields,K_space,pmacro)
#get the options
popts=pmodel.getOptions()
#get the boundary map
bcMap=pmodel.getBCs()

T_1=300
T_2=310

#initial temperature guess
popts["initialTemperature"]=305.
#set tolerances
popts.absTolerance=1e-7
popts.relTolerance=1e-5
#how often you want the residual to be printed
popts.showResidual=10

##################################################################
#This section sets the variables for specific meshes that I made.#
#L=1 unless the mesh is one of the 4 thin meshes.                #
##################################################################



######################################################################
#This section sets the boundary conditions.  Boundaries can either be#
#"temperature" or "reflecting."  For "reflecting," 1 means specular  #
#and 0 means diffuse.                                                #
######################################################################

rightbc=bcMap[6]
rightbc.bcType="temperature"
rightbc["specifiedTemperature"]=T_2
leftbc=bcMap[9]
leftbc.bcType="temperature"
leftbc["specifiedTemperature"]=T_1
topbc1=bcMap[8]
topbc1.bcType="reflecting"
topbc1["specifiedTemperature"]=T_2
topbc1["specifiedReflection"]=1.
topbc2=bcMap[5]
topbc2.bcType="reflecting"
topbc2["specifiedTemperature"]=T_2
topbc2["specifiedReflection"]=1.
botbc1=bcMap[10]
botbc1.bcType="reflecting"
botbc1["specifiedTemperature"]=T_1
botbc1["specifiedReflection"]=1.
botbc2=bcMap[7]
botbc2.bcType="reflecting"
botbc2["specifiedTemperature"]=T_1
botbc2["specifiedReflection"]=1.

L=1
L_side = 1
ac_th=L/(vg*tau)
Kn=vg*tau/L

#calculate the specific heat of the material (just for reference)
Cp=K_space.calcSpecificHeat((T_1+T_2)/2)
Cp1=K_space.calcSpecificHeat(T_1)
Cp2=K_space.calcSpecificHeat(T_2)

print "Specific Heat at Average Boundary Temperature: ",Cp
print "Specific Heat at Left Wall:",Cp1
print "Specific Heat at Right Wall:",Cp2
print "Acoustic Thickness:",ac_th
print "Knudsen Number:",Kn
print " "

#initialize
pmodel.init()
#initialize the boundary conditions
pmodel.callBoundaryConditions()
#iterate (argument is the max number of iterations)
pmodel.advance(100)

###################################################
#The simulation is finished, the rest is just post#
#processing                                       #
###################################################

writer = exporters.VTKWriterA(geomFields,meshes,
                              "temperatureProfile.vtk",
                              "Phonon Transport",
                              False,0)
writer.init()
writer.writeScalarField(pmacro.temperature,"Temperature")
writer.writeVectorField(pmacro.heatFlux,"HeatFlux")
writer.finish()

topflux=pmodel.HeatFluxIntegral(meshes[0],5)+pmodel.HeatFluxIntegral(meshes[1],8)
botflux=pmodel.HeatFluxIntegral(meshes[0],7)+pmodel.HeatFluxIntegral(meshes[1],10)
leftflux=pmodel.HeatFluxIntegral(meshes[1],9)
rightflux=pmodel.HeatFluxIntegral(meshes[0],6)
ballistic=Cp*vg*(T_2-T_1)/4.*L_side

print " "
print "=================================="
print "Acoustic Thickness: ",ac_th
print "Knudsen Number: ",Kn
if ac_th>=10.:
    psi=(4./3.)/(1.42089+ac_th)
    print "psi: ",psi
print "=================================="
print " "
print "Top Heat Flux: ", topflux
print "Bottom Heat Flux: ", botflux
print "Left Heat Flux: ", leftflux
print "Right Heat Flux: ", rightflux
print "Heat Balance: ",topflux+botflux+leftflux+rightflux
print "Ballistic Heat Flux: ",ballistic
print "Dimensionless (left): ",leftflux/ballistic
print " "
