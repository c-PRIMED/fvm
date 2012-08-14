import sys
sys.setdlopenflags(0x100|0x2)

import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")

reader = FluentCase("../test/TwoMaterialTest.cas")

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

vcRightZone['massDiffusivity'] = 2.e-6
vcLeftZone['massDiffusivity'] = 10.e-6

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

timeStep = 1e6 # large timestep due to large mesh dimesions and small diffusivity
numTimeSteps = 50 #approximately steady state
numIterPerTimeStep = 1 # linear so only one linear solve needed

def advanceUnsteady(smodel,geomFields,meshes,numTimeSteps,numIterPerTimeStep):
   for i in range(0,numTimeSteps):

     speciesFields = smodel.getSpeciesFields(0)

     filename = 'TimeStep_Species'
     filename += `i`
     filename += '.vtk'
     print filename
     writer = exporters.VTKWriterA(geomFields,meshes,filename,
                                         "TestSpecies",False,0)
     writer.init()
     writer.writeScalarField(speciesFields.massFraction,"MassFraction")
     writer.finish()

     smodel.advance(numIterPerTimeStep)
     print 'advancing to time step %i' % i
     smodel.updateTime()

soptions = smodel.getOptions()
soptions.transient = True
soptions.setVar('timeStep',timeStep)
soptions.setVar('initialMassFraction0',0.75)
soptions.setVar('initialMassFraction1',0.25)

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

#smodel.printBCs()
smodel.init()
advanceUnsteady(smodel,geomFields,meshes,numTimeSteps,numIterPerTimeStep)

mesh = meshes[1]
heatFlux = smodel.getMassFluxIntegral(mesh,6,0)
print heatFlux
mesh = meshes[0]
heatFlux2 = smodel.getMassFluxIntegral(mesh,5,0)
print heatFlux2
