#!/usr/bin/env python
import pdb
import sys
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt
import importers
#pdb.set_trace()

atype = 'double'

if atype == 'double':
    import models_atyped_double as models
    import exporters_atyped_double as exporters
elif atype == 'tangent':
    import models_atyped_tangent_double as models
    import exporters_atyped_tangent_double as exporters


from FluentCase import FluentCase

from fvmbaseExt import vectorInt

from fvmbaseExt import VecD3

fileBase = "/home/lin/work/app-memosa/src/fvm/test/CellMark/"

          
outfile = fileBase+"-prism.dat"
    
reader = FluentCase(fileBase+"cav32.cas")

#import debug
reader.read();

meshes = reader.getMeshList()

import time
t0 = time.time()

geomFields =  models.GeomFields('geom')

metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

flowFields =  models.FlowFields('flow')

fmodel = models.FlowModelA(geomFields,flowFields,meshes)

reader.importFlowBCs(fmodel)

octree = fvmbaseExt.Octree()

mesh0 = meshes[0]

cells = mesh0.getCells()

nCells = cells.getCount()

cellCentroid = geomFields.coordinate[cells].asNumPyArray()

octree.Impl(mesh0, geomFields)


solid = fvmbaseExt.MPM()

mpmFileName = fileBase + "MPMs_Beam.dat"

#solid.setandwriteParticles(mpmFileName)

solid.Impl(mpmFileName)


option = 2

fvmbaseExt.CellMark_Impl(mesh0, geomFields, fileBase, octree, solid, option)
