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

filename="test40x40"
extension=".cas"

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
ntheta=5
nphi=10
Tinit=310

K_space=pa.KspaceA(1,tau,vg,4e12,ntheta,nphi)

pmacro=pext.PhononMacro(filename)

cmodel=pa.COMETModelA(meshes,0,geomFields,K_space,pmacro)

copts=cmodel.getOptions()
cBCs=cmodel.getBCs()

copts["initialTemperature"]=Tinit
copts.showResidual=1
copts.maxLevels=3
copts.absTolerance=1e-8
copts.preSweeps=0
copts.postSweeps=2

cBCs[4].bcType="temperature"
cBCs[4]["specifiedTemperature"]=300
cBCs[5].bcType="temperature"
cBCs[5]["specifiedTemperature"]=300
cBCs[6].bcType="temperature"
cBCs[6]["specifiedTemperature"]=300
cBCs[2].bcType="temperature"
cBCs[2]["specifiedTemperature"]=300

cmodel.init()
cmodel.advance(1000)

writer = exporters.VTKWriterA(geomFields,meshes,
                              "temperatureProfile.vtk",
                              "Phonon Transport",
                              False,0)
writer.init()
writer.writeScalarField(pmacro.temperature,"Temperature")
writer.finish()
