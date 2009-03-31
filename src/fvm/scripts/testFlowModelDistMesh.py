#!/usr/bin/env python

import sys
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt
import importers
import time

atype = 'double'
#atype = 'tangent'

if atype == 'double':
    import models_atyped_double as models
    import exporters_atyped_double as exporters
elif atype == 'tangent':
    import models_atyped_tangent_double as models
    import exporters_atyped_tangent_double as exporters


from FluentCase import FluentCase

#fvmbaseExt.enableDebug("cdtor")

fileBase = None
numIterations = 10
#fileBase = "/home/yildirim/memosa/src/fvm/test/cav_44_tri"
#fileBase = "/home/sm/a/data/wj"

def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat if it is not specified."
    sys.exit(1)

def advance(fmodel,niter):
    for i in range(0,niter):
        try:
            fmodel.advance(1)
        except KeyboardInterrupt:
            break

# change as needed
#import debug


nmesh = 4
meshes = fvmbaseExt.MeshList()
readers = []
for n in range(0,nmesh):
   ss = "test_" + str(n) + ".cdf"
   print ss
   nc_reader = importers.NcDataReader( ss )
   thisMeshList = nc_reader.getMeshList()
   for m in  thisMeshList:
       meshes.push_back( m )

   readers.append(nc_reader)

## now create the mappers
for n in range(0,nmesh):
    readers[n].createMappers(meshes)
    

import time
t0 = time.time()

geomFields =  models.GeomFields('geom')
print "before metricsCalculator"
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)
print "after"
metricsCalculator.init()

if atype == 'tangent':
    metricsCalculator.setTangentCoords(0,7,1)
print "before flowFields " 
flowFields =  models.FlowFields('flow')

fmodel = models.FlowModelA(geomFields,flowFields,meshes)

## set bc for top to be a wall with x velocity
bc3 = fmodel.getBCMap()[3]
bc3.bcType = 'NoSlipWall'
bc3.setVar('specifiedXVelocity',1)

## set viscosity and density, this is done per mesh since each mesh has its own VC object
vcMap = fmodel.getVCMap()
for vc in vcMap.values():
    vc.setVar('density',1.0)
    vc.setVar('viscosity',1.0)


print "before BCreader " 


print "before Momentum solvers " 
momSolver = fvmbaseExt.AMG()
momSolver.relativeTolerance = 1e-1
momSolver.nMaxIterations = 20
momSolver.maxCoarseLevels=20
momSolver.verbosity=0

contSolver = fvmbaseExt.AMG()
#pc = fvmbaseExt.AMG()
#pc.verbosity=0
#contSolver = fvmbaseExt.BCGStab()
#contSolver.preconditioner = pc
contSolver.relativeTolerance = 1e-1
contSolver.nMaxIterations = 20
contSolver.verbosity=0
contSolver.maxCoarseLevels=20

foptions = fmodel.getOptions()

foptions.momentumLinearSolver = momSolver
foptions.pressureLinearSolver = contSolver

foptions.momentumTolerance=1e-3
foptions.continuityTolerance=1e-3
foptions.setVar("momentumURF",0.7)
foptions.setVar("pressureURF",0.3)
foptions.printNormalizedResiduals=False


fmodel.init()
#advance(fmodel,numIterations)

c0 = meshes[0].getCells()
c1 = meshes[1].getCells()
c2 = meshes[2].getCells()
c3 = meshes[3].getCells()

f0 = meshes[0].getFaces()
f1 = meshes[1].getFaces()
f2 = meshes[2].getFaces()
f3 = meshes[3].getFaces()

fc0 = meshes[0].getConnectivity(f0,c0)
fc1 = meshes[1].getConnectivity(f1,c1)
fc2 = meshes[2].getConnectivity(f2,c2)
fc3 = meshes[3].getConnectivity(f3,c3)


vol0 = geomFields.volume[c0].asNumPyArray()
vol1 = geomFields.volume[c1].asNumPyArray()
vol2 = geomFields.volume[c2].asNumPyArray()
vol3 = geomFields.volume[c3].asNumPyArray()

#import debug
fmodel.advance(100)
