import sys
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt
import importers
import models_atyped_double as models
import exporters_atyped_double as exporters
from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")

reader = FluentCase("/home/sm/a/data/case1.cas")

#import debug
reader.read();


meshes = reader.getMeshList()

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

flowFields =  models.FlowFields('flow')

fmodel = models.FlowModelA(geomFields,flowFields,meshes)

reader.importFlowBCs(fmodel)
#fmodel.printBCs()
fmodel.init()

momSolver = fmodel.getMomentumSolver()
momSolver.relativeTolerance = 1e-1
momSolver.nMaxCycles = 20
momSolver.verbosity = 1

contSolver = fmodel.getContinuitySolver()
contSolver.relativeTolerance = 1e-1
contSolver.nMaxCycles = 100
contSolver.verbosity=2
contSolver.maxCoarseLevels=1
fmodel.advance(1)

writer = exporters.FluentDataExporterA(reader,"test.dat",False,0)
writer.init()
writer.writeScalarField(flowFields.pressure,1)
writer.writeVectorField(flowFields.velocity,111)
writer.writeScalarField(flowFields.massFlux,18)
writer.finish()
