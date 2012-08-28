import sys
sys.setdlopenflags(0x100|0x2)

import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")

reader = FluentCase("../test/RadialTransfer.cas")
reader.read();
meshes_case = reader.getMeshList()

# add double shell to mesh between two materials
# When creating double shell mesh, the mesh it is created from
# is called the 'parent' mesh and the mesh that is passed as an argument
# is called the 'other' mesh
#
# parent has to be electrolyte, other has to be solid electrode 
# for B-V equations to work

interfaceID = 4
shellmesh = meshes_case[1].createDoubleShell(interfaceID, meshes_case[0], interfaceID)
meshes = [meshes_case[0], meshes_case[1], shellmesh]

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

nSpecies = 1

smodel = models.SpeciesModelA(geomFields,meshes,nSpecies)
bcmap = smodel.getBCMap(0)
vcmap = smodel.getVCMap(0)

vcElectrode = vcmap[0] # solid electrode on right
vcElectrolyte = vcmap[1] # electrolyte on left

vcElectrode['massDiffusivity'] = 1.e-13
vcElectrolyte['massDiffusivity'] = 2.66e-9

bcLeft = bcmap[5]
bcRight = bcmap[6]

#Dirichlet on Left,Right
bcLeft.bcType = 'SpecifiedMassFraction'
bcRight.bcType = 'SpecifiedMassFraction'
bcLeft.setVar('specifiedMassFraction', 5000.0)
bcRight.setVar('specifiedMassFraction', 15000.0)

soptions = smodel.getOptions()

# A = Phi_S    B = Phi_E
#soptions.setVar('A_coeff',0.72)
#soptions.setVar('B_coeff',0.42)
soptions.ButlerVolmer = True
soptions.setVar('ButlerVolmerRRConstant',5.0e-7)
soptions.setVar('interfaceUnderRelax',1.0)
soptions.setVar('initialMassFraction0',15000.0)
soptions.setVar('initialMassFraction1',5000.0)

solver = fvmbaseExt.BCGStab()
pc = fvmbaseExt.JacobiSolver()
pc.verbosity=0
solver.preconditioner = pc
solver.relativeTolerance = 1e-14 #solver tolerance
solver.absoluteTolerance = 1e-14 #solver tolerance
solver.nMaxIterations = 100
solver.maxCoarseLevels=30
solver.verbosity=0
soptions.linearSolver = solver
soptions.relativeTolerance=1e-14 #model tolerance
soptions.absoluteTolerance=1e-14 #model tolerance


######################################################
## Potential Model
######################################################

elecFields = models.ElectricFields('elec')
emodel = models.ElectricModelA(geomFields,elecFields,meshes)
bcmap2 = emodel.getBCMap()
vcmap2 = emodel.getVCMap()

vcElectrode = vcmap2[0] # solid electrode on right
vcElectrolyte = vcmap2[1] # electrolyte on left

vcElectrode['dielectric_constant'] = 4.29e11
vcElectrolyte['dielectric_constant'] = 2.825e11

bcLeft = bcmap2[5]
bcRight = bcmap2[6]

#Dirichlet on Left,Right
bcLeft.bcType = 'SpecifiedPotential'
bcRight.bcType = 'SpecifiedPotential'
bcLeft.setVar('specifiedPotential', 0.0)
bcRight.setVar('specifiedPotential', 0.106037515)

eoptions = emodel.getOptions()
eoptions.ButlerVolmer = True
eoptions.setVar('ButlerVolmerRRConstant',5.0e-7)

# A = c_S    B = c_E
#eoptions.setVar('Interface_A_coeff',615.7)
#eoptions.setVar('Interface_B_coeff',384.3)

solver = fvmbaseExt.BCGStab()
pc = fvmbaseExt.JacobiSolver()
pc.verbosity=0
solver.preconditioner = pc
solver.relativeTolerance = 1e-14 #solver tolerance
solver.absoluteTolerance = 1e-14 #solver tolerance
solver.nMaxIterations = 100
solver.maxCoarseLevels=30
solver.verbosity=0
eoptions.electrostaticsLinearSolver = solver
eoptions.electrostaticsTolerance=1e-14 #model tolerance
eoptions.chargetransportTolerance=1e-14 #model tolerance

eoptions.chargetransport_enable = False

# solve coupled system

emodel.init()
smodel.init()

# set the species number that cooresponds to potential model (Lithium)
speciesFields = smodel.getSpeciesFields(0)

for i in range(0,10):

   # set the species concentrations for the species of interest (Lithium)
   elecFields.speciesConcentration = speciesFields.massFraction

   print "POTENTIAL MODEL"
   emodel.advance(50)

   #set the potential for all species   
   for j in range(0,nSpecies):
      sFields = smodel.getSpeciesFields(j)
      sFields.elecPotential = elecFields.potential

   print "SPECIES MODEL"
   smodel.advance(25)

mesh = meshes[1]
massFlux = smodel.getMassFluxIntegral(mesh,5,0)
massFlux2 = smodel.getMassFluxIntegral(mesh,4,0)
print massFlux
print massFlux2

mesh = meshes[0]
massFlux3 = smodel.getMassFluxIntegral(mesh,6,0)
massFlux4 = smodel.getMassFluxIntegral(mesh,4,0)
print massFlux3
print massFlux4

writer = exporters.VTKWriterA(geomFields,meshes_case,"testSpeciesModel.vtk",
                                         "TestSpecies",False,0)
writer.init()
writer.writeScalarField(speciesFields.massFraction,"MassFraction")
writer.finish()

writer3 = exporters.VTKWriterA(geomFields,meshes_case,"testElectricModel.vtk",
                                         "TestSpecies",False,0)
writer3.init()
writer3.writeScalarField(elecFields.potential,"sConc")
writer3.finish()

#writer4 = exporters.VTKWriterA(geomFields,meshes_case,"electhroughspecies.vtk",
#                                         "TestPotential",False,0)
#writer4.init()
#writer4.writeScalarField(speciesFields.elecPotential,"Potential")
#writer4.finish()
