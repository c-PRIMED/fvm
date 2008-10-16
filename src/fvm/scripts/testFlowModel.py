import sys
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt
import importers
import models_atyped_double as models
import exporters_atyped_double as exporters
#import models_atyped_tangent_double as models
#import exporters_atyped_tangent_double as exporters
from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")

fileBase = "/home/sm/prism-meshes/cav32"
reader = FluentCase(fileBase+".cas")

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
#fmodel.printBCs()

momSolver = fmodel.getMomentumSolver()
momSolver.relativeTolerance = 1e-1
momSolver.nMaxCycles = 20
momSolver.maxCoarseLevels=20
momSolver.verbosity=0

contSolver = fmodel.getContinuitySolver()
contSolver.relativeTolerance = 1e-1
contSolver.nMaxCycles = 20
contSolver.verbosity=0
contSolver.maxCoarseLevels=20

foptions = fmodel.getOptions()
foptions.momentumTolerance=1e-5
foptions.continuityTolerance=1e-9
foptions.setVar("momentumURF",0.95)
foptions.setVar("pressureURF",0.1)
foptions.printNormalizedResiduals=False

fmodel.init()
fmodel.advance(4000)

t1 = time.time()

print 'solution time = %f' % (t1-t0)

writer = exporters.FluentDataExporterA(reader,fileBase+"-prism.dat",False,0)
writer.init()
writer.writeScalarField(flowFields.pressure,1)
writer.writeVectorField(flowFields.velocity,111)
writer.writeScalarField(flowFields.massFlux,18)
writer.finish()
