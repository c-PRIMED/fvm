import sys
sys.setdlopenflags(0x100|0x2)

import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")

reader = FluentCase("/home/brad/TwoMaterialTest.cas")

#import debug
reader.read();

meshes_case = reader.getMeshList()

#add double shell to mesh between two materials
interfaceID = 9
shellmesh = meshes_case[1].createDoubleShell(interfaceID, meshes_case[0], interfaceID)
meshes = [meshes_case[0], meshes_case[1], shellmesh]

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

nSpecies = 1

smodel = models.SpeciesModelA(geomFields,meshes,nSpecies)
bcmap = smodel.getBCMap(0)
vcmap = smodel.getVCMap(0)

vcRightZone = vcmap[0]
vcLeftZone = vcmap[1]

vcRightZone['massDiffusivity'] = 1.e-9
vcLeftZone['massDiffusivity'] = 5.e-9

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
soptions.setVar('initialMassFraction',0.0)

# interface jump conditon (Phi_R = A*Phi_L + B)
soptions.setVar('A_coeff',2.0)
soptions.setVar('B_coeff',0.0)

solver = fvmbaseExt.BCGStab()
pc = fvmbaseExt.AMG()
pc.verbosity=0
solver.preconditioner = pc
solver.relativeTolerance = 1e-8
#solver.absoluteTolerance = 1e-16
solver.nMaxIterations = 20000
solver.maxCoarseLevels=30
solver.verbosity=0
soptions.linearSolver = solver


#smodel.printBCs()
smodel.init()
smodel.advance(1)

mesh = meshes[1]
massFlux = smodel.getMassFluxIntegral(mesh,6,0)
print massFlux
mesh = meshes[0]
massFlux2 = smodel.getMassFluxIntegral(mesh,5,0)
print massFlux2

speciesFields = smodel.getSpeciesFields(0)
writer = exporters.VTKWriterA(geomFields,meshes_case,"testSpeciesModel.vtk",
                                         "TestSpecies",False,0)
writer.init()
writer.writeScalarField(speciesFields.massFraction,"MassFraction")
writer.finish()
