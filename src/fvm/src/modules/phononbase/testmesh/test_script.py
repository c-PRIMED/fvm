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

print "reading",fileBase+".cas"
reader = FluentCase(fileBase+".cas")
reader.read();
print "read mesh"
meshes = reader.getMeshList()
print "got meshlist"
geomFields =  fvm.models.GeomFields('geom')
print "got geomfields"
metricsCalculator = fvm.models.MeshMetricsCalculatorA(geomFields,meshes)
metricsCalculator.init()
print "metrics finished"

K_space=pa.KspaceA(1,1,1,1,2)
print "made k_space"
K_space.findDK3()
pmacro=pext.PhononMacro("e_dubprime")
print "made pmacro"

pmodel=pa.PhononModelA(meshes,geomFields,K_space,pmacro)
print "made pmodel"
pmodel.init()
print "initialized"
#raw_input("paused")
pmodel.callBoundaryConditions()
print "BCs finished"
