#!/usr/bin/env python

"""
the scipt is used to test the mesh dependency on tunneling model
it uses the new set of parameter from Sambit in Sep2010
"""
### import  modules ###
import pdb
import sys
import fvm
fvm.set_atype('double')
from math import *
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
from FluentCase import FluentCase
from fvm.fvmbaseExt import VecD3
import time
from optparse import OptionParser
from mpi4py import MPI


def usage():
    print __doc__
    sys.exit(1)

# map between fvm, tecplot, and xdmf types
etype = {
        'tri' : 1,
        'quad' : 2,
        'tetra' : 3,
        'hexa' : 4
        }
tectype = {
        'tri' : 'FETRIANGLE',
        'quad' : 'FEQUADRILATERAL',
        'tetra' : 'FETETRAHEDRON',
        'hexa' : 'FEBRICK'
        }
xtype = {
        'tri' : 'Triangle',
        'quad' : 'Quadrilateral',
        'tetra' : 'Tetrahedron',
        'hexa' : 'Hexahedron'
        }

parser = OptionParser()
parser.set_defaults(type='tri')
parser.add_option("--type", help="'tri'[default], 'quad', 'hexa', or 'tetra'")
parser.add_option("--xdmf", action='store_true', help="Dump data in xdmf")
parser.add_option("--time","-t",action='store_true',help="Print timing information.")
(options, args) = parser.parse_args()
if len(args) != 1:
    usage()

### mesh ID ###
topID = 4
botID = 5
sideID = 3
interiorID = 2
normal_direction  = 2      # z direction is the charging direction
nLevel  = 1000
nTrap = 2
### input and output ###
reader = FluentCase(args[0])

### physics parameters ###

dielectric_constant = 7.9

dielectric_thickness = 200e-9

applied_voltage = 100

dielectric_ionization = 3.0

dielectric_bandgap = 5.0

substrate_workfunction = 5.0

membrane_workfunction = 5.0

optical_dielectric_constant = 4.0

electron_trapdepth = 1.5

electron_trapdensity = 3e25

OP_temperature = 300

electron_effmass = 0.5

poole_frenkel_emission_frequency = 1e11

electron_capture_cross = 1e-22

electron_mobility = 50.0e-4

electron_saturation_velocity = 1.0e5

### run parameters ###
timeStep = 1e-9
timeScale = 1.1
numIterationsPerStep = 3
numTimeSteps = 100
globalTime = 0
globalCount = 0
saveFrequency = 10
totalChargeFile = open("totalCharges.dat", "w")
#targetChargeFile = open("targetCharges.dat","w")
#========================================================================================#

#----------------------------------------------------------#
def unsteadyAdvance(globalCount, globalTime, timeStep):
          
    if globalCount == 0:
        elecModel.calculateEquilibriumParameters()  
 
    for i in range(0, numTimeSteps):
        target = 1000-5
        (chargeSumT, chargeSumC) = calculateTotalCharges()

        saveTotalCharges(globalTime, chargeSumT/1e6, chargeSumC/1e6)
        #saveTargetCharges(globalTime, target)

        """
	if (globalCount % saveFrequency == 0):
            saveChargeProfile(globalCount)
            savePotentialProfile(globalCount)
        """
	elec_options['timeStep'] = timeStep

	try:
            elecModel.advance(numIterationsPerStep)
	except KeyboardInterrupt:
            break	  	
    	
        globalTime += timeStep
	globalCount += 1
        
        print "advaning to time %i\t %e" % (globalCount,globalTime)
        elecModel.updateTime()
        timeStep *= timeScale

#-----------------------------------------------------------#
def calculateTotalCharges():
    sumC = 0.0
    sumT = 0.0
    cells = meshes[0].getCells()
    nCells = cells.getSelfCount()
    charge = elecFields.charge[cells].asNumPyArray()
    volume = geomFields.volume[cells].asNumPyArray()
    for i in range(0, nCells):
        sumT = sumT + charge[i][0] + charge[i][1]
        sumC = sumC + charge[i][2]
    return sumT/nCells, sumC/nCells
#-----------------------------------------------------------#
def saveTotalCharges(time, sumC, sumT):
    totalChargeFile.write('%e\t%e\t%e\n' % (time, sumC, sumT))
    totalChargeFile.flush()

def saveTargetCharges(time, target):
    cells = meshes[0].getCells()
    nCells = cells.getSelfCount()
    charge = elecFields.charge[cells].asNumPyArray() 
    targetChargeFile.write('%e\t%e\t%e\t%e\n' % (time, charge[target][0], charge[target-1][0], charge[target+1][0]))
    targetChargeFile.flush()
#-----------------------------------------------------------#
def saveChargeProfile(nstep):
    print "saving charge profile"
    fileName = outputDir + str(nstep) + "_charge.dat"
    cells = meshes[0].getCells()
    nCells = cells.getSelfCount()
    charge = elecFields.charge[cells].asNumPyArray()
    file = open(fileName, "w")
    for i in range (nCells-1, 0-1, -1):
        file.write('%i\t%e\t%e\t%e\n' % (i, charge[i][0], charge[i][1], charge[i][2]))
    file.close()
def savePotentialProfile(nstep):
    print "saving potential"
    fileName = outputDir + str(nstep) + "_potential.dat"
    cells = meshes[0].getCells()
    nCells = cells.getSelfCount()
    potential = elecFields.potential[cells].asNumPyArray()
    file = open(fileName, "w")
    for i in range (0, nCells):
        file.write('%i\t%e\n' % (i, potential[i]))
    file.close()

#========================================================================================#

### read in meshes ###

reader.read()
meshes = reader.getMeshList()

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)
metricsCalculator.init()

elecFields =  models.ElectricFields('elec')

elecModel = models.ElectricModelA(geomFields,elecFields,meshes)

vcMap = elecModel.getVCMap()
dielectricID = interiorID
for i,vc in vcMap.iteritems():   
    vc.vcType = "dielectric"
    vc['dielectric_constant'] = dielectric_constant
   
# setup boundary conditions #
#--- potential BC: top = V; bot = 0; side = symmetry
#--- charges BC: zero Dirichlet BC applied everywhere 

membrane_voltage = applied_voltage

substrate_voltage = 0.0


bcMap = elecModel.getBCMap()
for i,bc in bcMap.iteritems():
    if i == topID:
        bc.bcType = "SpecifiedPotential"
        bc['specifiedPotential'] = membrane_voltage
    if i == botID:
        bc.bcType = "SpecifiedPotential"
        bc['specifiedPotential'] = substrate_voltage
    if i == sideID:
        bc.bcType = "Symmetry"

### setup initial condition ###
elec_options = elecModel.getOptions()
elec_options['initialPotential'] = 0
elec_options['initialTotalCharge'] = 0
elec_options['timeStep'] = timeStep

### setup enable options ###
elec_options.electrostatics_enable = True
elec_options.chargetransport_enable = True
elec_options.timeDiscretizationOrder = 1
elec_options.transient_enable = True
elec_options.injection_enable = True
elec_options.tunneling_enable = True
elec_options.emission_enable  = True
elec_options.capture_enable = True
elec_options.drift_enable = True
elec_options.trapbandtunneling_enable = True
elec_options.diffusion_enable = False

### setup physics constants ###
elec_constants = elecModel.getConstants()
elec_constants['dielectric_thickness'] = dielectric_thickness
elec_constants['voltage'] = applied_voltage
elec_constants['dielectric_ionization'] = dielectric_ionization
elec_constants['dielectric_bandgap'] = dielectric_bandgap
elec_constants['substrate_workfunction'] = substrate_workfunction
elec_constants['membrane_workfunction'] = membrane_workfunction
elec_constants['substrate_voltage'] = substrate_voltage
elec_constants['membrane_voltage'] = membrane_voltage
elec_constants['optical_dielectric_constant'] = optical_dielectric_constant
elec_constants['OP_temperature'] = OP_temperature
elec_constants['electron_effmass'] = electron_effmass
elec_constants['poole_frenkel_emission_frequency'] = poole_frenkel_emission_frequency
elec_constants['electron_capture_cross'] = electron_capture_cross
elec_constants['electron_mobility'] = electron_mobility
elec_constants['electron_saturation_velocity'] = electron_saturation_velocity
elec_constants['substrate_id'] = botID
elec_constants['membrane_id'] = topID
elec_constants['nLevel'] = nLevel
elec_constants['normal_direction'] = normal_direction
elec_constants['nTrap'] = nTrap


elec_constants.electron_trapdepth.push_back(electron_trapdepth);
elec_constants.electron_trapdensity.push_back(electron_trapdensity)
elec_constants.electron_trapdepth.push_back(1.5);
elec_constants.electron_trapdensity.push_back(electron_trapdensity)

### setup linear solve options ###
pPC = fvmbaseExt.AMG()
pPC.verbosity = 0
pSolver = fvmbaseExt.BCGStab()
pSolver.preconditioner = pPC
pSolver.relativeTolerance = 1e-20
pSolver.nMaxIterations = 100
pSolver.maxCoarseLevels=20
pSolver.absoluteTolerance = 1e-50
pSolver.verbosity=0
elec_options.electrostaticsLinearSolver = pSolver

cPC = fvmbaseExt.AMG()
cPC.verbosity = 0
cSolver = fvmbaseExt.BCGStab()
cSolver.preconditioner = cPC
cSolver.relativeTolerance = 1e-20
cSolver.nMaxIterations = 100
cSolver.maxCoarseLevels=20
cSolver.absoluteTolerance = 1e-50
cSolver.verbosity=0
elec_options.chargetransportLinearSolver = cSolver

elec_options.electrostaticsTolerance = 1e-20
elec_options.chargetransportTolerance = 1e-20
elec_options.printNormalizedResiduals = False


### advance loop ###
elecModel.init()

t1 = time.time()

unsteadyAdvance (globalCount, globalTime, timeStep)

t2 = time.time()


print  '\nsolution time = %f' % (t2-t1)
