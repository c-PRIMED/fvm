#!/usr/bin/env python

#import pdb
import sys

import fvm 
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers

import time
from numpy import *

import fvm.phonon_atyped_double as pa
import fvm.phononbaseExt as pext

from FluentCase import FluentCase

fvm.set_atype('double')

fileBase="/home/james/memosa/src/fvm/src/modules/phononbase/testmesh/test"

def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat if it is not specified."
    sys.exit(1)

reader = FluentCase(fileBase+".cas")
reader.read();
meshes = reader.getMeshList()
geomFields =  fvm.models.GeomFields('geom')
metricsCalculator = fvm.models.MeshMetricsCalculatorA(geomFields,meshes)
metricsCalculator.init()

K_space=pa.KspaceA(0.1,1,1,2,4)
pmacro=pext.PhononMacro("e_dubprime")

pmodel=pa.PhononModelA(meshes,geomFields,K_space,pmacro)
popts=pmodel.getOptions()
bcMap=pmodel.getBCs()

rightbc=bcMap[2]
rightbc.bcType="temperature"
leftbc=bcMap[4]
leftbc.bcType="temperature"
topbc=bcMap[6]
topbc.bcType="temperature"
botbc=bcMap[5]
botbc.bcType="temperature"

pmodel.callBoundaryConditions()
#pmodel.printTemp()
pmodel.advance(1000)
#pmodel.printTemp()
