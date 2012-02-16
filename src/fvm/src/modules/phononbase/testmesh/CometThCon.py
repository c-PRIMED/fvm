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

fileBase=sys.argv[1]
filename=sys.argv[2]
extension=".msh"
Kn=float(sys.argv[3])
KnName='Kn_'+sys.argv[3]
initialScale=1.e-4

BZfile="SiliconIsotropicN"
BZdisc=sys.argv[4]
levels=int(sys.argv[5])
T1=float(sys.argv[6])
T2=float(sys.argv[7])
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

ntheta=10
nphi=20
tau=1e-2
vg=1000

K_space=pa.KspaceA(fileBase+BZfile+BZdisc+".txt",dimension,1)
#K_space=pa.KspaceA(1,tau,vg,4e14,ntheta,nphi)
kArray=(K_space.getHollandConductivity((T1+T2)/2)).asNumPyArray()
kArray=kArray*eVtoJoule
print "Thermal Conductivity Tensor:"
print kArray[0],' ',kArray[1],' ',kArray[2]
print kArray[3],' ',kArray[4],' ',kArray[5]
print kArray[6],' ',kArray[7],' ',kArray[8]
print 'Specific Heat:',K_space.calcSpecificHeat(300)*eVtoJoule

kn0=float(K_space.findKnStats(initialScale))
avemfp=kn0*initialScale
scale=avemfp/Kn
kn0=float(K_space.findKnStats(scale))

for mesh in meshes:
    xNodes = mesh.getNodeCoordinates().asNumPyArray()
    xNodes[:,0] *= scale
    xNodes[:,1] *= scale
    xNodes[:,2] *= scale

metricsCalculator.init()

pmacro=pext.PhononMacro(filename)

cmodel=pa.COMETModelA(meshes,0,geomFields,K_space,pmacro)

copts=cmodel.getOptions()
cBCs=cmodel.getBCs()

FaceArea=cmodel.getWallAreaVector(meshes[0],6)
FaceAreaMag=math.sqrt(FaceArea[0]*FaceArea[0]+FaceArea[1]*FaceArea[1])
print 'Face Area Magnitude:', FaceAreaMag
print 'Scale:', scale

Tinit=(T1+T2)/2

BallisticRate=K_space.FindBallisticHeatRate(FaceArea,T1,T2)*eVtoJoule
print "Ballistic Heat Rate: ",BallisticRate

copts["initialTemperature"]=Tinit
copts.showResidual=10
copts.maxLevels=levels
copts.absTolerance=-1e-9
copts.preSweeps=1
copts.postSweeps=1
copts.relFactor=1
copts.withNormal=0
copts.NewtonTol=1e-4

#cBCs[4].bcType="reflecting"
cBCs[4].bcType="temperature"
cBCs[4]["specifiedTemperature"]=T1  #left for 1d
cBCs[4]["specifiedReflection"]=0.
#cBCs[5].bcType="temperature"
cBCs[5].bcType="reflecting"
cBCs[5]["specifiedReflection"]=1.
cBCs[5]["specifiedTemperature"]=T2  # bot for .cas
cBCs[5]["FullyImplicit"]=0
#cBCs[6].bcType="reflecting"
cBCs[6].bcType="temperature"
cBCs[6]["specifiedReflection"]=0.
cBCs[6]["specifiedTemperature"]=T2  #right for 1d, top for .cas
cBCs[7].bcType="reflecting"
#cBCs[7].bcType="temperature"
cBCs[7]["specifiedTemperature"]=T1
cBCs[7]["specifiedReflection"]=1.
cBCs[7]["FullyImplicit"]=0

print "Mesh file: "+filename
print "Scaling: "+str(scale)
print "BZ file: "+BZfile+BZdisc
print "Max Levels: ",copts.maxLevels

print "Initializing..."
cmodel.init()
cmodel.setStraightLine(T1,T2)
cmodel.equilibrate()
cmodel.applyTemperatureBoundaries()
cmodel.advance(0)
initialResid=cmodel.getResidual()
wall4=(cmodel.HeatFluxIntegral(meshes[0],4))*eVtoJoule
wall6=(cmodel.HeatFluxIntegral(meshes[0],6))*eVtoJoule
balance=abs(wall4+wall6)/wall4
iteration=0
print iteration,' : ',initialResid," :---------: ",wall4," : ",wall6," : ",balance*100.
t=0.
end=0
begin=0
residFile=open('Residual.dat','w')

resid=cmodel.getResidual()
residFile.write(str(iteration)+' '+str(resid)+' '+str(end-begin)+'\n')

total=100
step=1
balTol=.01
relTol=1.e-10
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
    wall4=(cmodel.HeatFluxIntegral(meshes[0],4))*eVtoJoule
    wall6=(cmodel.HeatFluxIntegral(meshes[0],6))*eVtoJoule
    balance=abs(wall4+wall6)/wall4
    relRes=resid/initialResid 
    residFile.write(str(iteration)+' '+str(resid)+' '+str(end-begin)+'\n')
    iteration+=step
    print iteration," : ",resid," : ",relRes," : ",wall4," ",wall6," ",balance*100.

thcon=wall4*scale/FaceAreaMag/(T2-T1)
    
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
print "Thermal Conductivity: ",thcon
            
print "Total Solution Time:",t
print ""
print "Mesh file: "+filename
print "Scaling: "+str(scale)
print "BZ file: "+BZfile+BZdisc
print "Max Levels: ",copts.maxLevels
K_space.findKnStats(scale)

name_file=KnName+'_'+str(copts.maxLevels)+'levs'

writer = exporters.VTKWriterA(geomFields,meshes,
                              name_file+'.vtk',
                              "Phonon Transport",
                              False,0)

writer.init()
writer.writeScalarField(pmacro.temperature,"Temperature")
#f copts.withNormal==1:
#  writer.writeVectorField(pmacro.lam,"Lambda")
writer.finish()

resFile=open(name_file+'.res','w')
resFile.write("Final Residual: "+str(resid)+"\n")
resFile.write("Iteration Count: "+str(iteration)+"\n")
resFile.write("Total Solution Time: "+str(t)+"\n")
resFile.write("Max Levels: "+str(copts.maxLevels)+"\n")
resFile.write("Initial Guess: "+str(Tinit)+"\n")
resFile.write("Wall 4: "+str(wall4)+"\n")
resFile.write("Wall 5: "+str(wall5)+"\n")
resFile.write("Wall 6: "+str(wall6)+"\n")
resFile.write("Wall 7: "+str(wall7)+"\n")
resFile.write("Ballistic Ratio: "+str(wall4/BallisticRate)+"\n")
resFile.close()

Hfile=open('Holland.dat','w')
Hfile.write(str(kArray[0]))
Hfile.close()

preFile=open('Predicted.dat','w')
preFile.write(str(thcon))
preFile.close()

nonDimFile=open('Dimensionless.dat','w')
nonDimFile.write(str(wall4/BallisticRate))
nonDimFile.close()
