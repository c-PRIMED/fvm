import sys
sys.setdlopenflags(0x100|0x2)

import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")

reader = FluentCase("../test/TwoMaterialTest.cas")
reader.read();
meshes_case = reader.getMeshList()

# add double shell to mesh between two materials
# When creating double shell mesh, the mesh it is created from
# is called the 'parent' mesh and the mesh that is passed as an argument
# is called the 'other' mesh
#
# parent has to be electrolyte, other has to be solid electrode 
# for B-V equations to work

interfaceID = 9
shellmesh = meshes_case[1].createDoubleShell(interfaceID, meshes_case[0], interfaceID)
meshes = [meshes_case[0], meshes_case[1], shellmesh]

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

nSpecies = 1

timeStep = 5e9 # large timestep due to large mesh dimesions and small diffusivity
numTimeSteps = 50 #approximately steady state
numIterPerTimeStep = 40 # nonlinear so multiple iterations per timestep needed

######################################################
## Species Model
######################################################

smodel = models.SpeciesModelA(geomFields,meshes,nSpecies)
bcmap = smodel.getBCMap(0)
vcmap = smodel.getVCMap(0)

vcRightZone = vcmap[0] # solid electrode
vcLeftZone = vcmap[1] # electrolyte

vcRightZone['massDiffusivity'] = 1.e-9
vcLeftZone['massDiffusivity'] = 1.e-9

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
bcRight.setVar('specifiedMassFraction', 1000.0)

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
soptions.transient = True
soptions.setVar('timeStep',timeStep)
soptions.setVar('initialMassFraction0',750)
soptions.setVar('initialMassFraction1',250)

# A = Phi_S    B = Phi_E
#soptions.setVar('A_coeff',0.72)
#soptions.setVar('B_coeff',0.42)
soptions.ButlerVolmer = True
soptions.setVar('interfaceUnderRelax',0.8)
#soptions.setVar('initialMassFraction',0.0)

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

vcRightZone = vcmap2[0] # solid electrode
vcLeftZone = vcmap2[1] # electrolyte

vcRightZone['dielectric_constant'] = 1.e10
vcLeftZone['dielectric_constant'] = 1.e10

bcLeft = bcmap2[6]
bcTop1 = bcmap2[8]
bcTop2 = bcmap2[7]
bcBottom1 = bcmap2[4]
bcBottom2 = bcmap2[1]
bcRight = bcmap2[5]

#Dirichlet on Left,Right
bcLeft.bcType = 'SpecifiedPotential'
bcRight.bcType = 'SpecifiedPotential'
bcLeft.setVar('specifiedPotential', 0.0)
bcRight.setVar('specifiedPotential', 1.138)

# Neumann on Bottom,Top
bcBottom1.bcType = 'SpecifiedPotentialFlux'
bcBottom1.setVar('specifiedPotentialFlux', 0.0)
bcTop1.bcType = 'SpecifiedPotentialFlux'
bcTop1.setVar('specifiedPotentialFlux', 0.0)
bcBottom2.bcType = 'SpecifiedPotentialFlux'
bcBottom2.setVar('specifiedPotentialFlux', 0.0)
bcTop2.bcType = 'SpecifiedPotentialFlux'
bcTop2.setVar('specifiedPotentialFlux', 0.0)

eoptions = emodel.getOptions()
eoptions.ButlerVolmer = True

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
eoptions.setVar('initialPotential',0.6)

################################################
######## Solve coupled system         ##########
################################################

def advanceUnsteady(smodel,emodel,elecFields,geomFields,meshes,numTimeSteps,numIterPerTimeStep):
   for t in range(0,numTimeSteps):

     speciesFields = smodel.getSpeciesFields(0)

     filename = 'TimeStep_Species'
     filename += `t`
     filename += '.vtk'
     writer = exporters.VTKWriterA(geomFields,meshes_case,filename,
                                         "TestBV",False,0)
     writer.init()
     writer.writeScalarField(speciesFields.massFraction,"MassFraction")
     writer.finish()

     filename2 = 'TimeStep_Electric'
     filename2 += `t`
     filename2 += '.vtk'
     writer2 = exporters.VTKWriterA(geomFields,meshes_case,filename2,
                                         "TestBV",False,0)
     writer2.init()
     writer2.writeScalarField(elecFields.potential,"sConc")
     writer2.finish()

     for i in range(0,numIterPerTimeStep):

        #set the potential for all species   
        for j in range(0,nSpecies):
           sFields = smodel.getSpeciesFields(j)
           sFields.elecPotential = elecFields.potential

        #print "SPECIES MODEL"
        smodel.advance(5)
   
        # set the species concentrations for the species of interest (Lithium)
        elecFields.speciesConcentration = speciesFields.massFraction

        #print "POTENTIAL MODEL"
        emodel.advance(2)

     print 'advancing to time step %i' % t
     smodel.updateTime()


# solve coupled system
emodel.init()
smodel.init()
advanceUnsteady(smodel,emodel,elecFields,geomFields,meshes,numTimeSteps,numIterPerTimeStep)

mesh = meshes[1]
massFlux = smodel.getMassFluxIntegral(mesh,6,0)
print massFlux
mesh = meshes[0]
massFlux2 = smodel.getMassFluxIntegral(mesh,5,0)
print massFlux2
