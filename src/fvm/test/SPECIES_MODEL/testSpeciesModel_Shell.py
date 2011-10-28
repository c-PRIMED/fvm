import sys
sys.setdlopenflags(0x100|0x2)

import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")
from mpi4py import MPI

reader = FluentCase(sys.argv[1])
reader.read();
meshes_case = reader.getMeshList()

# add double shell to mesh between two materials
# When creating double shell mesh, the mesh it is created from
# is called the 'parent' mesh and the mesh that is passed as an argument
# is called the 'other' mesh
#
# The phi values at the interface between the two meshes are related as follows:
# Phi_other = A * Phi_parent + B 
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
bcLeft.setVar('specifiedMassFraction', 1.0)
bcRight.setVar('specifiedMassFraction', 0.0)

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

# interface jump conditon (Phi_R = A*Phi_L + B)
# or (Phi_other = A*Phi_parent + B)
soptions.setVar('A_coeff',0.5)
soptions.setVar('B_coeff',0.1)

amg = fvmbaseExt.AMG()

solver = fvmbaseExt.BCGStab()
pc = fvmbaseExt.JacobiSolver()
pc.verbosity=0
solver.preconditioner = pc
solver.relativeTolerance = 1e-16 #solver tolerance
solver.absoluteTolerance = 1e-16 #solver tolerance
solver.nMaxIterations = 100
solver.maxCoarseLevels=30
solver.verbosity=0
soptions.linearSolver = solver
soptions.relativeTolerance=1e-16 #model tolerance
soptions.absoluteTolerance=1e-16 #model tolerance

#redirect stdout to a disk file
saveout = sys.stdout           
outfile = open('compare.dat','w')
sys.stdout = outfile  

amg.redirectPrintToFile("bcs.dat")
smodel.printBCs()
amg.redirectPrintToScreen()
smodel.init()
#smodel.advance(1)

mesh = meshes[1]
massFlux = smodel.getMassFluxIntegral(mesh,6,0)
print massFlux
mesh = meshes[0]
massFlux2 = smodel.getMassFluxIntegral(mesh,5,0)
print massFlux2
#restore
outfile.flush()
outfile.close()
sys.stdout = saveout
speciesFields = smodel.getSpeciesFields(0)
writer = exporters.VTKWriterA(geomFields,meshes_case,"testSpeciesModel.vtk",
                                         "TestSpecies",False,0)
writer.init()
writer.writeScalarField(speciesFields.massFraction,"MassFraction")
writer.finish()
