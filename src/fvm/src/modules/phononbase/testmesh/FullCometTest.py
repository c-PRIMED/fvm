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

filename="1d800"
extension=".msh"

BZ="graphene_data.txt"

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

vg=.01
tau=1.
ntheta=4
nphi=8
Tinit=305

K_space=pa.KspaceA(1,tau,vg,4e12,ntheta,nphi)
#K_space=pa.KspaceA(BZ,2)

pmacro=pext.PhononMacro(filename)

cmodel=pa.COMETModelA(meshes,0,geomFields,K_space,pmacro)

copts=cmodel.getOptions()
cBCs=cmodel.getBCs()

copts["initialTemperature"]=Tinit
copts.showResidual=1
copts.maxLevels=7
copts.absTolerance=1e-9
copts.preSweeps=1
copts.postSweeps=2

cBCs[4].bcType="temperature"
cBCs[4]["specifiedTemperature"]=300  #left for 1d
#cBCs[5].bcType="temperature"
cBCs[5].bcType="reflecting"
cBCs[5]["specifiedReflection"]=1.
cBCs[5]["specifiedTemperature"]=307.5
#cBCs[6].bcType="reflecting"
cBCs[6].bcType="temperature"
cBCs[6]["specifiedReflection"]=1.
cBCs[6]["specifiedTemperature"]=310  #right for 1d
cBCs[7].bcType="reflecting"
#cBCs[7].bcType="temperature"
cBCs[7]["specifiedTemperature"]=310

print "Mesh file: "+filename
print "Knudsen #:",vg*tau
print "Max Levels: ",copts.maxLevels

print "Initializing..."
cmodel.init()
print "Iterating..."
begin=time.clock()
cmodel.advance(50)
end=time.clock()

print "Total Solution Time:",end-begin

wallflux1=cmodel.HeatFluxIntegral(meshes[0],4)
wallflux2=cmodel.HeatFluxIntegral(meshes[0],6)
wallArea=cmodel.getWallArea(meshes[0],4)

Cp=K_space.calcSpecificHeat(305)
ballistic=Cp*vg*10./4.*wallArea


print "flux1",wallflux1/ballistic
print "flux2",wallflux2/ballistic

writer = exporters.VTKWriterA(geomFields,meshes,
                              "temperatureProfile.vtk",
                              "Phonon Transport",
                              False,0)
writer.init()
writer.writeScalarField(pmacro.temperature,"Temperature")
writer.finish()
