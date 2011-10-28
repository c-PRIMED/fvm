#!/usr/bin/env python
import sys
sys.setdlopenflags(0x100|0x2)

import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
from FluentCase import FluentCase
from   mpi4py import MPI
#fvmbaseExt.enableDebug("cdtor")

reader = FluentCase(sys.argv[1])

#import debug
reader.read();

meshes = reader.getMeshList()

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

nSpecies = 1
# species conditions
smodel = models.SpeciesModelA(geomFields,meshes,nSpecies)
bcmap = smodel.getBCMap(0)
vcmap = smodel.getVCMap(0)

vcRightZone = vcmap[0]
vcLeftZone = vcmap[1]

vcRightZone['massDiffusivity'] = 2.e-9
vcLeftZone['massDiffusivity'] = 10.e-9

bcLeft = bcmap[6]
bcTop1 = bcmap[8]
bcTop2 = bcmap[7]
bcBottom1 = bcmap[4]
bcBottom2 = bcmap[1]
bcRight = bcmap[5]

#Dirichlet on Left,Right
bcLeft.bcType = 'SpecifiedMassFraction'
bcRight.bcType = 'SpecifiedMassFraction'

bcLeft.setVar('specifiedMassFraction', 0.0)
bcRight.setVar('specifiedMassFraction', 1.0)

# Neumann on Bottom,Top
bcBottom1.bcType = 'SpecifiedMassFlux'
bcBottom1.setVar('specifiedMassFlux', 0.0)
bcTop1.bcType = 'SpecifiedMassFlux'
bcTop1.setVar('specifiedMassFlux', 0.0)
bcBottom2.bcType = 'SpecifiedMassFlux'
bcBottom2.setVar('specifiedMassFlux', 0.0)
bcTop2.bcType = 'SpecifiedMassFlux'
bcTop2.setVar('specifiedMassFlux', 0.0)

soptions = smodel.getOptions()
solver = fvmbaseExt.AMG()
solver.relativeTolerance = 1e-14 #solver tolerance
solver.absoluteTolerance = 1e-16 #solver tolerance
solver.nMaxIterations = 100
solver.maxCoarseLevels=30
solver.verbosity=0
soptions.linearSolver = solver

# model tolerances are only needed for non-linear or unstructured problems
#soptions.relativeTolerance=1e-16
#soptions.absoluteTolerance=1e-16
#redirect stdout to a disk file
saveout = sys.stdout           
outfile = open('compare.dat','w')
sys.stdout = outfile   

solver.redirectPrintToFile("bcs.dat")
smodel.printBCs()
solver.redirectPrintToScreen()
smodel.init()
#smodel.advance(1)

mesh = meshes[1]
heatFlux = smodel.getMassFluxIntegral(mesh,6,0)
print heatFlux
mesh = meshes[0]
heatFlux2 = smodel.getMassFluxIntegral(mesh,5,0)
print heatFlux2

#restore
outfile.flush()
outfile.close()
sys.stdout = saveout

speciesFields = smodel.getSpeciesFields(0)

writer = exporters.VTKWriterA(geomFields,meshes,"testSpeciesModel.vtk",
                                         "TestSpecies",False,0)
writer.init()
writer.writeScalarField(speciesFields.massFraction,"MassFraction")
writer.finish()
