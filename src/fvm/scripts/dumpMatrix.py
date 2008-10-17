import sys
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt
import importers


import models_atyped_double as models
import exporters_atyped_double as exporters


from FluentCase import FluentCase

# change as needed
fileBase = "/home/sm/a/data/wj"

reader = FluentCase(fileBase+".cas")

reader.read();


meshes = reader.getMeshList()


geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()


flowFields =  models.FlowFields('flow')

fmodel = models.FlowModelA(geomFields,flowFields,meshes)

reader.importFlowBCs(fmodel)

momSolver = fmodel.getMomentumSolver()
momSolver.relativeTolerance = 1e-1
momSolver.nMaxCycles = 20
momSolver.maxCoarseLevels=20
momSolver.verbosity=0

contSolver = fmodel.getContinuitySolver()
contSolver.relativeTolerance = 1e-1
contSolver.nMaxCycles = 20
contSolver.verbosity=1
contSolver.maxCoarseLevels=20



fmodel.init()
fmodel.advance(1)

fmodel.dumpContinuityMatrix(fileBase)
