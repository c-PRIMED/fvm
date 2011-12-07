#!/usr/bin/env python

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

filename="80x80Biased"
extension=".msh"

BZfile="SiliconIsotropicN"
BZdisc="60K"
dimension=3

eVtoJoule=1.60217646e-19

fileBase="/home/james/memosa/src/fvm/src/modules/phononbase/testmesh/"+filename

def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat if it is not specified."
    sys.exit(1)

reader = FluentCase(fileBase+extension)
reader.read();
meshes = reader.getMeshList()
geomFields =  fvm.models.GeomFields('geom')
metricsCalculator = fvm.models.MeshMetricsCalculatorA(geomFields,meshes)
metricsCalculator.init()

scale=1.e-6

for mesh in meshes:
    xNodes = mesh.getNodeCoordinates().asNumPyArray()
    xNodes[:,0] *= scale
    xNodes[:,1] *= scale
    xNodes[:,2] *= scale

metricsCalculator.init()

vg=.1
tau=1.
ntheta=2
nphi=3
Tinit=60.5

#K_space=pa.KspaceA(1,tau,vg,4e12,ntheta,nphi)
K_space=pa.KspaceA(BZfile+BZdisc+".txt",dimension,1)
K_space.findKnStats(scale)

pmacro=pext.PhononMacro(filename)

cmodel=pa.COMETModelA(meshes,0,geomFields,K_space,pmacro)

copts=cmodel.getOptions()
cBCs=cmodel.getBCs()

copts["initialTemperature"]=Tinit
copts.showResidual=10
copts.maxLevels=1
copts.absTolerance=-1e-9
copts.preSweeps=1
copts.postSweeps=1
copts.relFactor=1
copts.withNormal=1

cBCs[4].bcType="temperature"
cBCs[4]["specifiedTemperature"]=60  #left for 1d
cBCs[5].bcType="temperature"
#cBCs[5].bcType="reflecting"
cBCs[5]["specifiedReflection"]=1.
cBCs[5]["specifiedTemperature"]=61  # bot for .cas
#cBCs[6].bcType="reflecting"
cBCs[6].bcType="temperature"
cBCs[6]["specifiedReflection"]=1.
cBCs[6]["specifiedTemperature"]=60  #right for 1d, top for .cas
#cBCs[7].bcType="reflecting"
cBCs[7].bcType="temperature"
cBCs[7]["specifiedTemperature"]=60

print "Mesh file: "+filename
print "Scaling: "+str(scale)
print "BZ file: "+BZfile+BZdisc

if copts.withNormal==1:
    print "Including Normal Scattering"
    name_file=filename+BZdisc+str(scale)+"tempWithN.vtk"
else:
    print "Not Including Normal Scattering"
    name_file=filename+BZdisc+str(scale)+"tempWithoutN.vtk"

print "Max Levels: ",copts.maxLevels

print "Initializing..."
cmodel.init()
print "Iterating..."
t=0.
begin=time.clock()
cmodel.advance(1)
end=time.clock()
t+=(end-begin)
copts.absTolerance=1e-10

resid=cmodel.getResidual()
print "Initial Residual:",resid
wallflux1=cmodel.HeatFluxIntegral(meshes[0],4)
wallflux2=cmodel.HeatFluxIntegral(meshes[0],6)
wallflux3=cmodel.HeatFluxIntegral(meshes[0],5)
wallflux4=cmodel.HeatFluxIntegral(meshes[0],7)
print "Wall 4: ",wallflux1*eVtoJoule
print "Wall 6: ",wallflux2*eVtoJoule
print "Wall 5: ",wallflux3*eVtoJoule
print "Wall 7: ",wallflux4*eVtoJoule
print "Balance: ",(wallflux1+wallflux2+wallflux3+wallflux4)*eVtoJoule

total=30
step=5
iteration=0
balance=5
balTol=.01

while resid>copts.absTolerance and iteration<total:
    begin=time.clock()
    cmodel.advance(step)
    end=time.clock()
    t+=(end-begin)
    iteration+=step
    wallflux1=cmodel.HeatFluxIntegral(meshes[0],4)
    wallflux2=cmodel.HeatFluxIntegral(meshes[0],6)
    wallflux3=cmodel.HeatFluxIntegral(meshes[0],5)
    wallflux4=cmodel.HeatFluxIntegral(meshes[0],7)
    print ""
    print "Iteration: ",iteration
    print "Wall 4: ",wallflux1*eVtoJoule
    print "Wall 6: ",wallflux2*eVtoJoule
    print "Wall 5: ",wallflux3*eVtoJoule
    print "Wall 7: ",wallflux4*eVtoJoule
    balance=(wallflux1+wallflux2+wallflux3+wallflux4)*eVtoJoule
    print "Balance: ",balance
    resid=cmodel.getResidual()
    print "Residual: ",resid
    writer = exporters.VTKWriterA(geomFields,meshes,
                              name_file,
                              "Phonon Transport",
                              False,0)
    writer.init()
    writer.writeScalarField(pmacro.temperature,"Temperature")
    if copts.withNormal==1:
        writer.writeVectorField(pmacro.lam,"Lambda")
    writer.finish()

resid=cmodel.getResidual()
print "Final Residual:",resid
wallflux1=cmodel.HeatFluxIntegral(meshes[0],4)
wallflux2=cmodel.HeatFluxIntegral(meshes[0],6)
wallflux3=cmodel.HeatFluxIntegral(meshes[0],5)
wallflux4=cmodel.HeatFluxIntegral(meshes[0],7)
print "Wall 4: ",wallflux1*eVtoJoule
print "Wall 6: ",wallflux2*eVtoJoule
print "Wall 5: ",wallflux3*eVtoJoule
print "Wall 7: ",wallflux4*eVtoJoule
print "Balance: ",(wallflux1+wallflux2+wallflux3+wallflux4)*eVtoJoule
            
print "Total Solution Time:",t
print ""
print "Mesh file: "+filename
print "Scaling: "+str(scale)
print "BZ file: "+BZfile+BZdisc
print "Max Levels: ",copts.maxLevels
K_space.findKnStats(scale)

if copts.withNormal==1:
    print "Including Normal Scattering"
    name_file=filename+BZdisc+str(scale)+"tempWithN.vtk"
else:
    print "Not Including Normal Scattering"
    name_file=filename+BZdisc+str(scale)+"tempWithoutN.vtk"

writer.init()
writer.writeScalarField(pmacro.temperature,"Temperature")
if copts.withNormal==1:
    writer.writeVectorField(pmacro.lam,"Lambda")
writer.finish()
