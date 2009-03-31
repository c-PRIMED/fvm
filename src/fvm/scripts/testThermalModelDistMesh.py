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

metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

if atype == 'tangent':
    metricsCalculator.setTangentCoords(0,7,1)

thermalFields =  models.ThermalFields('therm')

tmodel = models.ThermalModelA(geomFields,thermalFields,meshes)


## set bc for top to be at 400
bc3 = tmodel.getBCMap()[3]
bc3.bcType = 'SpecifiedTemperature'
bc3.setVar('specifiedTemperature',400)

bc4 = tmodel.getBCMap()[4]
bc4.bcType = 'SpecifiedTemperature'
bc4.setVar('specifiedTemperature',400)

## set viscosity and density, this is done per mesh since each mesh has its own VC object
#vcMap = tmodel.getVCMap()
#for vc in vcMap.values():
#    vc.setVar('density',1.0)
#    vc.setVar('thermalConductivity',1.0)


tSolver = fvmbaseExt.AMG()
tSolver.relativeTolerance = 1e-6
tSolver.nMaxIterations = 20
tSolver.maxCoarseLevels=20
tSolver.verbosity=1

toptions = tmodel.getOptions()

toptions.linearSolver = tSolver

#import debug
tmodel.init()
tmodel.advance(1)
