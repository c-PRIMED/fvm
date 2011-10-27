import sys
sys.setdlopenflags(0x100|0x2)

import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")

reader = FluentCase("../test/SpeciesTest.cas")

#import debug
reader.read();

meshes = reader.getMeshList()

print "here"
geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

nSpecies = 2
smodel = models.SpeciesModelA(geomFields,meshes,nSpecies)

#Species 1 boundary conditions
SN = 0
bcmap1 = smodel.getBCMap(SN)
vcmap1 = smodel.getVCMap(SN)

vcZone1 = vcmap1[0]
vcZone1['massDiffusivity'] = 1e-6

bcLeft1 = bcmap1[4]
bcTop1 = bcmap1[6]
bcBottom1 = bcmap1[5]
bcRight1 = bcmap1[3]

#Dirichlet on Left,Right
bcLeft1.bcType = 'SpecifiedMassFraction'
bcRight1.bcType = 'SpecifiedMassFraction'
bcLeft1.setVar('specifiedMassFraction', 1.0)
bcRight1.setVar('specifiedMassFraction', 0.0)

# Neumann on Bottom, Top
bcBottom1.bcType = 'SpecifiedMassFlux'
bcBottom1.setVar('specifiedMassFlux', 0.0)
bcTop1.bcType = 'SpecifiedMassFlux'
bcTop1.setVar('specifiedMassFlux', 0.0)

#Species 2 BCs
SN = 1
bcmap2 = smodel.getBCMap(SN)
vcmap2 = smodel.getVCMap(SN)

vcZone2 = vcmap2[0]
vcZone2['massDiffusivity'] = 1e-6

bcLeft2 = bcmap2[4]
bcTop2 = bcmap2[6]
bcBottom2 = bcmap2[5]
bcRight2 = bcmap2[3]

#Dirichlet on Left,Right
bcLeft2.bcType = 'SpecifiedMassFraction'
bcRight2.bcType = 'SpecifiedMassFraction'
bcLeft2.setVar('specifiedMassFraction', 0.0)
bcRight2.setVar('specifiedMassFraction', 1.0)

# Neumann on Bottom, Top
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

#bmodel.printBCs()
smodel.init()
smodel.advance(1)

mesh = meshes[0]
# species 1 fluxes
print smodel.getMassFluxIntegral(mesh,3,0)
print smodel.getMassFluxIntegral(mesh,4,0)

# species 2 fluxes
print smodel.getMassFluxIntegral(mesh,3,1)
print smodel.getMassFluxIntegral(mesh,4,1)

# export species data
speciesFields1 = smodel.getSpeciesFields(0)
speciesFields2 = smodel.getSpeciesFields(1)
writer = exporters.VTKWriterA(geomFields,meshes,
                              "testSpeciesModel_MultSpecies.vtk",
                              "TestSpecies",False,0)
writer.init()
writer.writeScalarField(speciesFields1.massFraction,"MassFraction1")
writer.writeScalarField(speciesFields2.massFraction,"MassFraction2")
writer.finish()
