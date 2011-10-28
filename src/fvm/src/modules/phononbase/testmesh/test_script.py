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

#Specify the name of the file (it could be a .msh or a .cas)
filename="test80x80"
extension=".cas"

fileBase="/home/james/memosa/src/fvm/src/modules/phononbase/testmesh/"+filename

def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat if it is not specified."
    sys.exit(1)

############################################################
#Read the mesh, and extract the information into the MEMOSA#
#data structures                                           # 
############################################################

reader = FluentCase(fileBase+extension)
reader.read();
meshes = reader.getMeshList()
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
ntheta=6
nphi=12

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

L=1.
L_side=1.
if filename=="thin" or filename=="thin_packed" or filename=="thin_2000" or filename=="thin_packed_2000":
    L_side=0.05

if filename=="test" or filename=="test100" or filename=="test40x40" or filename=="test80x80":
    right=2
    left=4
    top=6
    bot=5
elif filename=="tilted" or filename=="thin" or filename=="thin_packed" or filename=="thin_2000" or filename=="thin_packed_2000":
    top=4
    left=5
    bot=6
    right=7
else:
    print "Wrong file name"
    sys.exit(0)


######################################################################
#This section sets the boundary conditions.  Boundaries can either be#
#"temperature" or "reflecting."  For "reflecting," 1 means specular  #
#and 0 means diffuse.                                                #
######################################################################

rightbc=bcMap[right]
rightbc.bcType="temperature"
rightbc["specifiedTemperature"]=T_2
leftbc=bcMap[left]
leftbc.bcType="temperature"
leftbc["specifiedTemperature"]=T_1
topbc=bcMap[top]
topbc.bcType="reflecting"
topbc["specifiedTemperature"]=T_2
topbc["specifiedReflection"]=1.
botbc=bcMap[bot]
botbc.bcType="reflecting"
botbc["specifiedTemperature"]=T_1
botbc["specifiedReflection"]=1.

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
pmodel.advance(10000)

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

topflux=pmodel.HeatFluxIntegral(meshes[0],top)
botflux=pmodel.HeatFluxIntegral(meshes[0],bot)
leftflux=pmodel.HeatFluxIntegral(meshes[0],left)
rightflux=pmodel.HeatFluxIntegral(meshes[0],right)
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
