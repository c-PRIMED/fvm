import sys
sys.setdlopenflags(0x100|0x2)

import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")

reader = FluentCase("/home/brad/ThermalTest3.cas")

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


vcElectrode = vcmap[0]
vcElectrode['massDiffusivity'] = 1

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

timeStep = 0.25
numTimeSteps = 200
numIterPerTimeStep = 50

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
soptions.setVar('initialMassFraction',0.0)

#bmodel.printBCs()
smodel.init()
advanceUnsteady(smodel,geomFields,meshes,numTimeSteps,numIterPerTimeStep)

mesh = meshes[0]
print smodel.getMassFluxIntegral(mesh,3,0)
print smodel.getMassFluxIntegral(mesh,4,0)
print smodel.getMassFluxIntegral(mesh,5,0)
print smodel.getMassFluxIntegral(mesh,6,0)
