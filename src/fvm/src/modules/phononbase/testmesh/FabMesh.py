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

#print " "

fileBase="./"
filename="FabMesh_8.61e-6mX6.74e-6m"
extension=".cas"
Kn=1000.
KnName='Kn_'+str(Kn)
initialScale=1.e-4

SiBZfile="SiIsoEDIP_"
GeBZfile="GeIsoHarrison_"
BZdisc="2x8x10"
BZfine="original"
levels=4
T1=300.
T2=301.
Tinit=(T1+T2)/2.+0
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

SiKspace=pa.KspaceA(fileBase+SiBZfile+BZdisc+".txt",dimension,1)
GeKspace=pa.KspaceA(fileBase+GeBZfile+BZdisc+".txt",dimension,1)
#SiKspaceFine=pa.KspaceA(fileBase+SiBZfile+BZfine+".txt",dimension,1)
#GeKspaceFine=pa.KspaceA(fileBase+GeBZfile+BZfine+".txt",dimension,1)

SikArray=(SiKspace.getHollandConductivity((T1+T2)/2)).asNumPyArray()
SikArray=SikArray*eVtoJoule
print "Silicon Thermal Conductivity Tensor:"
print SikArray[0],' ',SikArray[1],' ',SikArray[2]
print SikArray[3],' ',SikArray[4],' ',SikArray[5]
print SikArray[6],' ',SikArray[7],' ',SikArray[8]
print 'Silicon Specific Heat:',SiKspace.calcSpecificHeat(300)*eVtoJoule

kn0=float(SiKspace.findKnStats(initialScale))
avemfp=kn0*initialScale
#print "Silicon Average mean free path", avemfp
scale=1#avemfp/Kn
kn0=float(SiKspace.findKnStats(scale))

for mesh in meshes:
    xNodes = mesh.getNodeCoordinates().asNumPyArray()
    xNodes[:,0] *= scale
    xNodes[:,1] *= scale
    xNodes[:,2] *= scale

metricsCalculator.init()

pmacro=pext.PhononMacro(filename)

Klist=SiKspace.MakeList()
GeKspace.AddToList(Klist)

cmodel=pa.COMETModelA(meshes,0,geomFields,Klist,pmacro)

copts=cmodel.getOptions()
cBCs=cmodel.getBCs()
mkMap=cmodel.getMKMap();
mkMap[0]=0
mkMap[1]=1

FaceArea=cmodel.getWallAreaVector(meshes[0],7)
FaceAreaMag=math.sqrt(FaceArea[0]*FaceArea[0]+FaceArea[1]*FaceArea[1])
#print 'Face Area Magnitude:', FaceAreaMag
print 'Scale:', scale

SiBallisticRate=SiKspace.FindBallisticHeatRate(FaceArea,T1,T2)*eVtoJoule
#print "Silicon Ballistic Heat Rate: ",SiBallisticRate


SiDOS=pa.DOS(SiKspace)
SiDOS.binMode(0,2,0,2.74e13)
SiDOS.binMode(1,2,0,7.43e13)
SiDOS.binMode(2,2,9.7e13,1.035e14)
SiDOS.binMode(3,2,7.42e13,1.035e14)
SiDOS.setDensity()

#SiDOSfine=pa.DOS(SiKspaceFine)
#SiDOSfine.binMode(0,2,0,2.74e13)
#SiDOSfine.binMode(1,2,0,7.43e13)
#SiDOSfine.binMode(2,2,9.7e13,1.035e14)
#SiDOSfine.binMode(3,2,7.42e13,1.035e14)
#SiDOSfine.setDensity()

#print " "

GekArray=(GeKspace.getHollandConductivity((T1+T2)/2)).asNumPyArray()
GekArray=GekArray*eVtoJoule
print "Germanium Thermal Conductivity Tensor:"
print GekArray[0],' ',GekArray[1],' ',GekArray[2]
print GekArray[3],' ',GekArray[4],' ',GekArray[5]
print GekArray[6],' ',GekArray[7],' ',GekArray[8]
print 'Germanium Specific Heat:',GeKspace.calcSpecificHeat(300)*eVtoJoule

kn0=float(GeKspace.findKnStats(scale))
avemfp=kn0*scale
#print "Germanium Average mean free path", avemfp

GeBallisticRate=GeKspace.FindBallisticHeatRate(FaceArea,T1,T2)*eVtoJoule
#print "Germanium Ballistic Heat Rate: ", GeBallisticRate


GeDOS=pa.DOS(GeKspace)
GeDOS.binMode(1,2,0,4.044e13)
GeDOS.binMode(0,2,0,1.5005e13)
GeDOS.binMode(2,2,5.31e13,5.66e13)
GeDOS.binMode(3,2,4.074e13,5.66e13)
GeDOS.setDensity()

#GeDOSfine=pa.DOS(GeKspaceFine)
#GeDOSfine.binMode(1,1,0,4.044e13)
#GeDOSfine.binMode(0,1,0,1.5005e13)
#GeDOSfine.binMode(2,1,5.31e13,5.66e13)
#GeDOSfine.binMode(3,1,4.074e13,5.66e13)
#GeDOSfine.setDensity()

#TtransSiGe=SiDOSfine.makeDMMtransmission(GeDOSfine,Tinit,1)
TtransSiGeDummy=SiDOS.makeDMMtransmission(GeDOS,T1,1)
transSiGe=TtransSiGeDummy.asNumPyArray()
TSiBins=SiDOS.getFreqBins()
SiKspace.setTransmission(GeKspace,TSiBins,TtransSiGeDummy)
SiMids=SiDOS.getFreqMids().asNumPyArray()

#transGeSi=GeDOSfine.makeDMMtransmission(SiDOSfine,Tinit,0)
TtransGeSiDummy=GeDOS.makeDMMtransmission(SiDOS,T1,0)
transGeSi=TtransGeSiDummy.asNumPyArray()
GeKspace.setTransmission(SiKspace,TSiBins,TtransGeSiDummy)

SiKspace.setDOS(SiDOS)
GeKspace.setDOS(GeDOS)
"""
transSiGefile=open("transSiGe","w")
for i in range(SiMids.size):
    print SiMids[i],transSiGe[i],transGeSi[i]
    transSiGefile.write(str(SiMids[i])+" "+str(transSiGe[i])+"\n")

transSiGefile.close()
"""

SiBallisticInterface=SiKspace.calcBallisticInterface(GeKspace,FaceArea,T1,T2)*eVtoJoule
#print "Si Ballistic", SiBallisticInterface
GeBallisticInterface=GeKspace.calcBallisticInterface(SiKspace,FaceArea,T1,T2)*eVtoJoule
#print "Ge Ballistic", GeBallisticInterface

#sys.exit()

copts["initialTemperature"]=Tinit
copts.showResidual=10
copts.maxLevels=levels
copts.absTolerance=-1e-9
copts.preSweeps=1
copts.postSweeps=2
copts.relFactor=1
copts.withNormal=0
copts.NewtonTol=1e-6


#cBCs[5].bcType="reflecting"
#cBCs[5].bcType="temperature"
cBCs[5]["specifiedTemperature"]=T2
cBCs[5]["specifiedReflection"]=1.
cBCs[3].bcType="temperature"
#cBCs[3].bcType="reflecting"
cBCs[3]["specifiedReflection"]=1.
cBCs[3]["specifiedTemperature"]=T1
cBCs[3]["FullyImplicit"]=0
#cBCs[7].bcType="reflecting"
cBCs[7].bcType="temperature"
cBCs[7]["specifiedReflection"]=1.
cBCs[7]["specifiedTemperature"]=T1
#cBCs[8].bcType="reflecting"
cBCs[8].bcType="temperature"
cBCs[8]["specifiedTemperature"]=T1
cBCs[8]["specifiedReflection"]=1.
cBCs[8]["FullyImplicit"]=0
#cBCs[9].bcType="reflecting"
cBCs[9].bcType="temperature"
cBCs[9]["specifiedTemperature"]=T2
cBCs[9]["specifiedReflection"]=1.
cBCs[9]["FullyImplicit"]=0
#cBCs[10].bcType="reflecting"
cBCs[10].bcType="temperature"
cBCs[10]["specifiedTemperature"]=T2
cBCs[10]["specifiedReflection"]=1.
cBCs[10]["FullyImplicit"]=0
#cBCs[11].bcType="reflecting"
cBCs[11].bcType="temperature"
cBCs[11]["specifiedTemperature"]=T2
cBCs[11]["specifiedReflection"]=1.
cBCs[11]["FullyImplicit"]=0
#cBCs[12].bcType="reflecting"
cBCs[12].bcType="temperature"
cBCs[12]["specifiedTemperature"]=T1
cBCs[12]["specifiedReflection"]=1.
cBCs[12]["FullyImplicit"]=0
#cBCs[13].bcType="reflecting"
cBCs[13].bcType="temperature"
cBCs[13]["specifiedTemperature"]=T2
cBCs[13]["specifiedReflection"]=1.
cBCs[13]["FullyImplicit"]=0
#cBCs[15].bcType="reflecting"
#cBCs[15].bcType="temperature"
cBCs[15]["specifiedTemperature"]=T1
cBCs[15]["specifiedReflection"]=1.
cBCs[15]["FullyImplicit"]=0

print "Initializing..."
cmodel.init()
print "Initialized"
cmodel.advance(0)
initialResid=cmodel.getResidual()

#Mesh 0
wall12=(cmodel.HeatFluxIntegral(meshes[0],12))*eVtoJoule
wall7=(cmodel.HeatFluxIntegral(meshes[0],7))*eVtoJoule
wall9=(cmodel.HeatFluxIntegral(meshes[0],9))*eVtoJoule
wall10=(cmodel.HeatFluxIntegral(meshes[0],10))*eVtoJoule


#Mesh 1
wall3=(cmodel.HeatFluxIntegral(meshes[1],3))*eVtoJoule
wall13=(cmodel.HeatFluxIntegral(meshes[1],13))*eVtoJoule
wall8=(cmodel.HeatFluxIntegral(meshes[1],8))*eVtoJoule
wall11=(cmodel.HeatFluxIntegral(meshes[1],11))*eVtoJoule

#Interfaces
wall5=(cmodel.HeatFluxIntegral(meshes[0],5))*eVtoJoule
wall15=(cmodel.HeatFluxIntegral(meshes[1],15))*eVtoJoule

balance=0#abs(wall5+wall9)/wall9
iteration=0
print iteration,' : ',initialResid,":",wall12+wall3," : ",wall7+wall8,":",wall9+wall13,":",wall10+wall11

t=0.
end=0
begin=0
resid=cmodel.getResidual()

total=10
step=1
balTol=.01
relTol=1.e-10
div_count=0
relRes=1.

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
            #sys.exit()
    t+=(end-begin)
    cmodel.applyTemperatureBoundaries()
    #Mesh 0
    wall12=(cmodel.HeatFluxIntegral(meshes[0],12))*eVtoJoule
    wall7=(cmodel.HeatFluxIntegral(meshes[0],7))*eVtoJoule
    wall9=(cmodel.HeatFluxIntegral(meshes[0],9))*eVtoJoule
    wall10=(cmodel.HeatFluxIntegral(meshes[0],10))*eVtoJoule
    
    #Mesh 1
    wall3=(cmodel.HeatFluxIntegral(meshes[1],3))*eVtoJoule
    wall13=(cmodel.HeatFluxIntegral(meshes[1],13))*eVtoJoule
    wall8=(cmodel.HeatFluxIntegral(meshes[1],8))*eVtoJoule
    wall11=(cmodel.HeatFluxIntegral(meshes[1],11))*eVtoJoule

    balance=wall12+wall7+wall9+wall10+wall3+wall13+wall8+wall11
    iteration+=step
    relRes=resid/initialResid
    print iteration,' : ',resid,":",relRes,":",wall12+wall3," : ",wall7+wall8,":",wall9+wall13,":",wall10+wall11,":",balance

#Mesh 0
wall12=(cmodel.HeatFluxIntegral(meshes[0],12))*eVtoJoule
wall7=(cmodel.HeatFluxIntegral(meshes[0],7))*eVtoJoule
wall9=(cmodel.HeatFluxIntegral(meshes[0],9))*eVtoJoule
wall10=(cmodel.HeatFluxIntegral(meshes[0],10))*eVtoJoule

#Mesh 1
wall3=(cmodel.HeatFluxIntegral(meshes[1],3))*eVtoJoule
wall13=(cmodel.HeatFluxIntegral(meshes[1],13))*eVtoJoule
wall8=(cmodel.HeatFluxIntegral(meshes[1],8))*eVtoJoule
wall11=(cmodel.HeatFluxIntegral(meshes[1],11))*eVtoJoule

print "Solution Time: ",t

meshlist0=fvmbaseExt.MeshList()
meshlist0.push_back(meshes[0])
meshlist1=fvmbaseExt.MeshList()
meshlist1.push_back(meshes[1])

name_file="FabMesh"

writer0 = exporters.VTKWriterA(geomFields,meshlist0,
                              name_file+'0'+'.vtk',
                              "Phonon Transport",
                              False,0)

writer0.init()
writer0.writeScalarField(pmacro.temperature,"Temperature")
for mode in range(4):
    writer0.writeScalarField(pmacro.getModeTemp(0,mode),"Mode"+str(mode)+"Mesh"+str(0))
writer0.finish()


writer1 = exporters.VTKWriterA(geomFields,meshlist1,
                              name_file+'1'+'.vtk',
                              "Phonon Transport",
                              False,0)

writer1.init()
writer1.writeScalarField(pmacro.temperature,"Temperature")
for mode in range(4):
    writer1.writeScalarField(pmacro.getModeTemp(1,mode),"Mode"+str(mode)+"Mesh"+str(1))
writer1.finish()
