import sys
sys.setdlopenflags(0x100|0x2)

import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")

reader = FluentCase("../test/2D_Battery_Radius4.cas")
reader.read();
meshes_case = reader.getMeshList()

# add double shell to mesh between two materials
# When creating double shell mesh, the mesh it is created from
# is called the 'parent' mesh and the mesh that is passed as an argument
# is called the 'other' mesh
#
# parent has to be electrolyte, other has to be solid electrode 
# for B-V equations to work

interfaceID = 1
shellmesh = meshes_case[0].createDoubleShell(interfaceID, meshes_case[1], interfaceID)
meshes = [meshes_case[0], meshes_case[1], shellmesh]

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

nSpecies = 1

timeStep = 1e10 # large timestep due to large mesh dimesions and small diffusivity
numTimeSteps = 5 #approximately steady state
numIterPerTimeStep = 5 # nonlinear so multiple iterations per timestep needed

AppliedCurrent = 1.0 #in Amps
F = 98485.0 # Faraday const
D_electrode = 1.e-13
D_electrolyte = 1.e-11
k_electrode = 4.29e11
k_electrolyte = 2.825e11
FluxArea = 20.0
MassFlux = AppliedCurrent/F/FluxArea #in mol/s/m^2
PotentialFlux = AppliedCurrent/FluxArea # in A/m^2 or C/s/m^2
print MassFlux
print PotentialFlux

######################################################
## Species Model
######################################################
smodel = models.SpeciesModelA(geomFields,meshes,nSpecies)
bcmap = smodel.getBCMap(0)
vcmap = smodel.getVCMap(0)

vcRightZone = vcmap[1] # solid electrode
vcLeftZone = vcmap[0] # electrolyte

vcRightZone['massDiffusivity'] = D_electrode
vcLeftZone['massDiffusivity'] = D_electrolyte

bcLeft = bcmap[7]
bcTop = bcmap[6]
bcBottom = bcmap[4]
bcRight = bcmap[5]

#Dirichlet on Left,Right
#bcLeft.bcType = 'SpecifiedMassFraction'
#bcRight.bcType = 'SpecifiedMassFraction'
#bcLeft.setVar('specifiedMassFraction', 0.0)
#bcRight.setVar('specifiedMassFraction', 1000.0)
bcLeft.bcType = 'SpecifiedMassFlux'
bcLeft.setVar('specifiedMassFlux', 0.0)
bcRight.bcType = 'SpecifiedMassFlux'
bcRight.setVar('specifiedMassFlux', MassFlux)

# Neumann on Bottom,Top
bcBottom.bcType = 'SpecifiedMassFlux'
bcBottom.setVar('specifiedMassFlux', 0.0)
bcTop.bcType = 'SpecifiedMassFlux'
bcTop.setVar('specifiedMassFlux', 0.0)

soptions = smodel.getOptions()
soptions.transient = True
soptions.setVar('timeStep',timeStep)
soptions.setVar('initialMassFraction0',2000)
soptions.setVar('initialMassFraction1',14870)

# A = Phi_S    B = Phi_E
#soptions.setVar('A_coeff',0.72)
#soptions.setVar('B_coeff',0.42)
soptions.ButlerVolmer = True
soptions.setVar('interfaceUnderRelax',1.0)
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

vcRightZone = vcmap2[1] # solid electrode
vcLeftZone = vcmap2[0] # electrolyte

vcRightZone['dielectric_constant'] = k_electrode
vcLeftZone['dielectric_constant'] = k_electrolyte

bcLefte = bcmap2[7]
bcTope = bcmap2[6]
bcBottome = bcmap2[4]
bcRighte = bcmap2[5]

#Dirichlet on Left,Right
bcLefte.bcType = 'SpecifiedPotential'
bcLefte.setVar('specifiedPotential', 0.0)
#bcLefte.bcType = 'SpecifiedPotentialFlux'
#bcLefte.setVar('specifiedPotentialFlux', 0.0)
bcRighte.bcType = 'SpecifiedPotentialFlux'
bcRighte.setVar('specifiedPotentialFlux', PotentialFlux)
# Neumann on Bottom,Top
bcBottome.bcType = 'SpecifiedPotentialFlux'
bcBottome.setVar('specifiedPotentialFlux', 0.0)
bcTope.bcType = 'SpecifiedPotentialFlux'
bcTope.setVar('specifiedPotentialFlux', 0.0)

eoptions = emodel.getOptions()
eoptions.ButlerVolmer = True

# A = c_S    B = c_E
eoptions.setVar('Interface_A_coeff',14780)
eoptions.setVar('Interface_B_coeff',2000)

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
eoptions.setVar('initialPotential',0.1)

mesh = meshes[0]
massFlux1 = smodel.getMassFluxIntegral(mesh,1,0)
print massFlux1
print "############################################################"



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

        # set the species concentrations for the species of interest (Lithium)
        elecFields.speciesConcentration = speciesFields.massFraction

        #print "POTENTIAL MODEL"
        emodel.advance(200)

        #set the potential for all species   
        for j in range(0,nSpecies):
           sFields = smodel.getSpeciesFields(j)
           sFields.elecPotential = elecFields.potential

        #print "SPECIES MODEL"
        smodel.advance(50)

        mesh = meshes[1]
        massFlux1 = smodel.getMassFluxIntegral(mesh,1,0)
        print massFlux1
        print "############################################################"

     print 'Finished time step %i' % t
     smodel.updateTime()


# solve coupled system
emodel.init()
smodel.init()
advanceUnsteady(smodel,emodel,elecFields,geomFields,meshes,numTimeSteps,numIterPerTimeStep)

'''
mesh = meshes[1]
massFlux = smodel.getMassFluxIntegral(mesh,6,0)
print massFlux
mesh = meshes[0]
massFlux2 = smodel.getMassFluxIntegral(mesh,5,0)
print massFlux2
'''
