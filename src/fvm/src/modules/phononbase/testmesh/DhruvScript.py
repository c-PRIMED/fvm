#!/usr/bin/env python

import sys
import os

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

fileBase="./"
filename="50x50"
extension=".msh"
Kn=10
KnName='Kn_'
scale=3.1649049e-07

BZfile="Sin"
BZdisc="whocares"
levels=6
T1=300.
T2=301.
Tinit=(T1+T2)/2.
dimension=3
eVtoJoule=1.60217646e-19

def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat if it is not specified."
    sys.exit(1)

reader = FluentCase(fileBase+filename+extension)
reader.read();
meshes = reader.getMeshList()
geomFields =  fvm.models.GeomFields('geom')
metricsCalculator = fvm.models.MeshMetricsCalculatorA(geomFields,meshes)
metricsCalculator.init()

A=1.5708e-10
ntheta=1
nphi=4
tau=3.753e-12
vg=8433
omega=5.1e13
Cp=1623801.53768/eVtoJoule


#There are two way to construct the Brillioun zone.  To make a gray approximation, use constructor 1.
#   You specify the lattice constant (to get the radius of the Brillioun zone), the relaxation time, the group velocity,
#   the frequency (rad/s), the number of theta discretizations, and the number of phi discretizations.  In the grayKspace, it
#   automatically assumes a 2D physical space and therefore it only creates the upper hemisphere of the Brillioun zone.

#Constructor 2 is how you get a non gray Brillioun zone.  The first argument is the absolute path to the output file created from
#  the matlab script "MakeSiKspace.m".  The second argument should always be 3, and the last argument should always be 1.


K_space=pa.KspaceA(A,tau,vg,omega,ntheta,nphi)                                    #constructor 1
#K_space=pa.KspaceA(fileBase+BZfile+BZdisc+".txt",dimension,1)     # constructor 2



#this method calculates the diffuse thermal conductivity of the material with the given
#  wave vector space discretization.  It just evaluates the integral given by Holland.  Since we
#  are dealing with non-constant Cp, you need to specify a temperature.  Also, inside MEMOSA,
#  ALL UNITS OF ENERGY ARE IN ELECTRON VOLTS

kArray=(K_space.getHollandConductivity((T1+T2)/2)).asNumPyArray()
kArray=kArray*eVtoJoule
print "Thermal Conductivity Tensor:"
print kArray[0],' ',kArray[1],' ',kArray[2]
print kArray[3],' ',kArray[4],' ',kArray[5]
print kArray[6],' ',kArray[7],' ',kArray[8]
print 'Specific Heat:',K_space.calcSpecificHeat(300)*eVtoJoule

#this method calculates the average, max, and minimum Kn for the given kspace.
#the return value is the average Kn.
kn0=float(K_space.findKnStats(scale))
avemfp=kn0*scale
print "Average mean free path", avemfp


# ths is just a way of changing the size of the input mesh...
for mesh in meshes:
    xNodes = mesh.getNodeCoordinates().asNumPyArray()
    xNodes[:,0] *= scale
    xNodes[:,1] *= scale
    xNodes[:,2] *= scale

metricsCalculator.init()

pmacro=pext.PhononMacro(filename)

Klist=K_space.MakeList()

cmodel=pa.COMETModelA(meshes,0,geomFields,Klist,pmacro)

copts=cmodel.getOptions()
cBCs=cmodel.getBCs()
mkMap=cmodel.getMKMap();
mkMap[0]=0

FaceArea=cmodel.getWallAreaVector(meshes[0],6)
FaceAreaMag=math.sqrt(FaceArea[0]*FaceArea[0]+FaceArea[1]*FaceArea[1])
print 'Face Area Magnitude:', FaceAreaMag
print 'Scale:', scale


#Given 2 temperatures and an area vector, this method calculates the ballistic heat rate.
BallisticRate=K_space.FindBallisticHeatRate(FaceArea,T1,T2)*eVtoJoule
print "Ballistic Heat Rate: ",BallisticRate


#Options for the model
copts["initialTemperature"]=Tinit
copts.showResidual=10
copts.maxLevels=levels
copts.absTolerance=-1e-9
copts.preSweeps=1
copts.postSweeps=2
copts.withNormal=0
copts.NewtonTol=1e-5

#cBCs[4].bcType="reflecting"
cBCs[4].bcType="temperature"
cBCs[4]["specifiedTemperature"]=T1  # Left boundary
cBCs[4]["specifiedReflection"]=0.
#cBCs[5].bcType="temperature"
cBCs[5].bcType="reflecting"
cBCs[5]["specifiedReflection"]=1.
cBCs[5]["specifiedTemperature"]=T2  # Top boundary
cBCs[5]["FullyImplicit"]=0
#cBCs[6].bcType="reflecting"
cBCs[6].bcType="temperature"
cBCs[6]["specifiedReflection"]=0.
cBCs[6]["specifiedTemperature"]=T2  # Right boundary
cBCs[7].bcType="reflecting"
#cBCs[7].bcType="temperature"
cBCs[7]["specifiedTemperature"]=T1  # Bottom boundary
cBCs[7]["specifiedReflection"]=1.
cBCs[7]["FullyImplicit"]=0

print "Mesh file: "+filename
print "Scaling: "+str(scale)
print "BZ file: "+BZfile+BZdisc
print "Max Levels: ",copts.maxLevels

print "Initializing..."
cmodel.init()
print "Initialized"



cmodel.advance(0)
initialResid=cmodel.getResidual()
wall4=(cmodel.HeatFluxIntegral(meshes[0],4))*eVtoJoule             # method that calculate the heat rate on a wall
wall6=(cmodel.HeatFluxIntegral(meshes[0],6))*eVtoJoule
balance=abs(wall4+wall6)/wall4
iteration=0
print iteration,' : ',initialResid," :---------------------: ",wall4," : ",wall6," : ",balance*100.
t=0.
end=0
begin=0
residFile=open('Residual.dat','w')

resid=cmodel.getResidual()
residFile.write(str(iteration)+' '+str(resid)+' '+str(end-begin)+'\n')

total=0
step=1
balTol=.01
relTol=1.e-13
relRes=1.
div_count=0

while (balance>balTol or relRes>relTol) and iteration<total:
    begin=time.clock()
    cmodel.advance(step)
    end=time.clock()
    resid=cmodel.getResidual()    
    if resid>initialResid:
        initialResid=resid
        div_count+=1
        if div_count>8:
            print 'Divergence Detected'
            sys.exit()
    t+=(end-begin)
    wall4=(cmodel.HeatFluxIntegral(meshes[0],4))*eVtoJoule/BallisticRate
    wall6=(cmodel.HeatFluxIntegral(meshes[0],6))*eVtoJoule/BallisticRate
    wall5=cmodel.HeatFluxIntegral(meshes[0],5)*eVtoJoule/BallisticRate
    wall7=cmodel.HeatFluxIntegral(meshes[0],7)*eVtoJoule/BallisticRate
    balance=abs(wall4+wall6+wall5+wall7)/wall4*100.
    relRes=resid/initialResid 
    residFile.write(str(iteration)+' '+str(resid)+' '+str(end-begin)+'\n')
    iteration+=step
    print iteration," : ",relRes," : ",resid," : ",wall4," : ",wall6," : ",wall5," : ",wall7," : ",balance," : ",t

    
residFile.close()
resid=cmodel.getResidual()
print "Final Residual:",resid
wall4=cmodel.HeatFluxIntegral(meshes[0],4)*eVtoJoule
wall6=cmodel.HeatFluxIntegral(meshes[0],6)*eVtoJoule
wall5=cmodel.HeatFluxIntegral(meshes[0],5)*eVtoJoule
wall7=cmodel.HeatFluxIntegral(meshes[0],7)*eVtoJoule
print "Wall 4: ",wall4
print "Wall 5: ",wall5
print "Wall 6: ",wall6
print "Wall 7: ",wall7
print "Balance: ",(wall4+wall5+wall6+wall7)
print "Ballistic Ratio: ",wall4/BallisticRate
            
print "Total Solution Time:",t
print ""
print "Mesh file: "+filename
print "Scaling: "+str(scale)
print "BZ file: "+BZfile+BZdisc
print "Max Levels: ",copts.maxLevels
K_space.findKnStats(scale)


#This makes the .vtk file to be read in from Paraview

name_file=KnName+'_'+str(copts.maxLevels)+'levs'

writer = exporters.VTKWriterA(geomFields,meshes,
                              name_file+'.vtk',
                              "Phonon Transport",
                              False,0)

writer.init()
writer.writeScalarField(pmacro.temperature,"Temperature")
writer.finish()
