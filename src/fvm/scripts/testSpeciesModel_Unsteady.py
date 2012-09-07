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

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

nSpecies = 1
smodel = models.SpeciesModelA(geomFields,meshes,nSpecies)
bcmap = smodel.getBCMap(0)
vcmap = smodel.getVCMap(0)
print "here"

vcElectrode = vcmap[0]
vcElectrode['massDiffusivity'] = 1e-6

bcLeft = bcmap[4]
bcTop = bcmap[6]
bcBottom = bcmap[5]
bcRight = bcmap[3]

#Dirichlet on Left,Right
bcLeft.bcType = 'SpecifiedMassFraction'
bcRight.bcType = 'SpecifiedMassFraction'

bcLeft.setVar('specifiedMassFraction', 1.0)
bcRight.setVar('specifiedMassFraction', 0.0)

# Neumann on Bottom
bcBottom.bcType = 'SpecifiedMassFlux'
bcBottom.setVar('specifiedMassFlux', 0.0)

# chage to Neumann on Top
bcTop.bcType = 'SpecifiedMassFlux'
bcTop.setVar('specifiedMassFlux', 0.0)

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
soptions.setVar('initialMassFraction0',1.0)

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
advanceUnsteady(smodel,geomFields,meshes,numTimeSteps,numIterPerTimeStep)

mesh = meshes[0]
print smodel.getMassFluxIntegral(mesh,3,0)
print smodel.getMassFluxIntegral(mesh,4,0)
print smodel.getMassFluxIntegral(mesh,5,0)
print smodel.getMassFluxIntegral(mesh,6,0)
