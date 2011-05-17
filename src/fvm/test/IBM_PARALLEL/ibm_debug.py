#!/usr/bin/env python

### this script solve 3D cantileve pull-in by coupling 
### plate model and electrostatics and fluid via IBM


### import  modules ###
import pdb
import sys
import os
from math import *
import fvm
fvm.set_atype('double')
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
import fvm.fvmparallel as fvmparallel
from FluentCase import FluentCase
from fvm.fvmbaseExt import VecD3
import time
import optparse
from numpy import *
import ClientCouplingPlate
from mpi4py import MPI


etype = {
        'tri' : 1,
        'quad' : 2,
        'tetra' : 3,
        'hexa' : 4
        }

def dumpMPITimeProfile(solver_maxtime, solver_mintime):
    fname = "time_mpi_totalprocs%s_fluid.dat" % MPI.COMM_WORLD.Get_size()
    f = open(fname,'w')
    line = " solver_mintime    = " + str(solver_mintime[0])    + "\n" + \
           " solver_maxtime    = " + str(solver_maxtime[0])    + "\n"
    print line
    f.write(line)
    f.close()


def checkMarking(n):

    cells = fluidMeshes[0].getCells()
    nCells = cells.getCount()
    cellCoords = geomFields.coordinate[cells].asNumPyArray()

    fluidFile = open(fileBase_output + "fluidCells_" + str(n) + ".dat", "w")
    solidFile = open(fileBase_output + "solidCells_" + str(n) + ".dat", "w")
    IBFile = open(fileBase_output + "IBCells_" + str(n) + ".dat", "w")

    cellIBType = geomFields.ibType[cells].asNumPyArray()
    for c in range (0, nCells):    
        ibtype = cellIBType[c]
        if ibtype == -1:
            fluidFile.write("%e\t%e\t%e\n" % (cellCoords[c][0], cellCoords[c][1], cellCoords[c][2]))
        elif ibtype == -2:
            IBFile.write("%e\t%e\t%e\n" % (cellCoords[c][0], cellCoords[c][1], cellCoords[c][2]))
        elif ibtype == -3:
            solidFile.write("%e\t%e\t%e\n" % (cellCoords[c][0], cellCoords[c][1], cellCoords[c][2]))
        elif ibtype == -5:
            print ("%i\t%i\t%e\t%e\n" % (c,ibtype,  cellCoords[c][0], cellCoords[c][1]))

               
    fluidFile.close()
    solidFile.close()
    IBFile.close()


def writeForceData():
    forceFile.write('%e\t%e\t%e\n'   % (globalTime, eForce, fForce))
    forceFile.flush()

def saveVTK(n):
    writer = exporters.VTKWriterA(geomFields,fluidMeshes,
                                  fileBase_output + "elecfield-" + str(n) + ".vtk",
                                  "frogleg",
                                  False,0)
    writer.init()
    writer.writeScalarField(elecFields.potential,"potential")
    writer.writeVectorField(elecFields.electric_field,"potentialgradient")  
    writer.writeVectorField(flowFields.velocity,"velocity")
    writer.writeScalarField(flowFields.pressure, "pressure")
    writer.finish()


    writer3 = exporters.VTKWriterA(geomFields,solidBoundaryMeshes,
                                  fileBase_output + "beamBoundary-" + str(n) + ".vtk",
                                  "beam Boundary",
                                  False,0,True)
    writer3.init()
    #writer3.writeVectorField(flowFields.velocity,"velocity")
    #writer3.writeVectorField(flowFields.force,"flow_force")
    #writer3.writeVectorField(elecFields.force,"elec_force")
    writer3.finish()
### ========================== properties and parameters ===============================###

parser = optparse.OptionParser()
parser.add_option("--volt", type=float)
parser.add_option("--type", help="'tri'[default], 'quad', 'hexa', or 'tetra'")
parser.add_option("--time","-t",action='store_true',help="Print timing information.")
(options, args) = parser.parse_args()
applied_voltage = -options.volt
print "voltage =",applied_voltage

### beam
rho = 8800                    # density kg/m^3
E = 200e9                     # Young's modulus 
nu = 0.31                      # Poisson's ratio

### electric field
#applied_voltage = -90
dielectric_constant = 1.0
beam_thickness = 3e-6

### mesh id
fluidTop = 5
fluidBot = 3
fluidLeft = 6
fluidRight = 4
#electrode = 7


numTimeSteps = 1
globalTime = 0
globalCount = 0
timeStep = 2e-8
saveFrequency = 100
initialTransient = False

### ===================== mesh read ===============================================###

fileBase_output = "./" + str(int(-applied_voltage)) + "/"

### 3D fluid mesh
fluidReader = FluentCase(args[0])
fluidReader.read();
fluent_meshes_fluid = fluidReader.getMeshList()
#paritioning
nmesh = 1
npart = [MPI.COMM_WORLD.Get_size()]
print "options.typeFluid = ", options.type
etypeFluid = [etype[options.type]]
#partMesh constructor and setTypes
part_mesh_fluid = fvmparallel.MeshPartitioner( fluent_meshes_fluid, npart, etypeFluid );
part_mesh_fluid.setWeightType(0);
part_mesh_fluid.setNumFlag(0);
#actions
part_mesh_fluid.isDebug(0)
part_mesh_fluid.partition()
part_mesh_fluid.mesh()
fluidMeshes  = part_mesh_fluid.meshList()
if not MPI.COMM_WORLD.Get_rank():
   print "partition is done for Fluid Mesh"

### generate solid boundary mesh
### 2D plate mesh
beamReader = FluentCase(args[1])
beamReader.read()
solidMeshes =  beamReader.getMeshList()

if options.time:
    solver_start   = zeros(1, dtype='d')
    solver_end     = zeros(1, dtype='d')
    solver_time    = zeros(1, dtype='d')
    solver_maxtime = zeros(1, dtype='d')
    solver_mintime = zeros(1, dtype='d')
    solver_start[0] = MPI.Wtime()


geomFields =  models.GeomFields('geom')
fluidMetricsCalculator = models.MeshMetricsCalculatorA(geomFields,fluidMeshes)
fluidMetricsCalculator.init() 


solidBoundaryMeshes = [m.extractBoundaryMesh() for m in solidMeshes]
solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(geomFields,solidBoundaryMeshes)
solidBoundaryMetricsCalculator.init()


### electric model and boundary condition ###

elecFields =  models.ElectricFields('elec')

emodel = models.ElectricModelA(geomFields,elecFields,fluidMeshes)

bcMap = emodel.getBCMap()

### specify a potential field; use it as boundary condition
pf = fvmbaseExt.Field('potential')
fgs = fluidMeshes[0].getBoundaryFaceGroups()
for fg in fgs:
     if fg.id == fluidBot:
         faces = fg.site
         nFaces = faces.getCount()
         areaMag = geomFields.areaMag[faces]
         xf = geomFields.coordinate[faces].asNumPyArray()
         pot = areaMag.newSizedClone(nFaces)
         pota = pot.asNumPyArray()
         for f in range(0, nFaces):
            pota[f] = applied_voltage


if fluidBot in bcMap:
   bc = bcMap[fluidBot]
   bc.bcType = "SpecifiedPotential"
   #bc['specifiedPotential'] = applied_voltage
   bc['specifiedPotential'] = 100.0

for i in [fluidTop, fluidLeft, fluidRight]:
   if i in bcMap:
      bc = bcMap[i]
      bc.bcType = "Symmetry"

vcMap = emodel.getVCMap()
for i,vc in vcMap.iteritems():
    vc.vcType = "dielectric"
    vc['dielectric_constant'] = dielectric_constant    

### flow model and boundary condition ###

flowFields =  models.FlowFields('flow')
fmodel = models.FlowModelA(geomFields,flowFields,fluidMeshes)
fluidReader.importFlowBCs(fmodel, fluidMeshes)


### ================================= solvers ===================================###
### elec solver ###
#epc = fvmbaseExt.AMG()
elecSolver = fvmbaseExt.AMG()
elecSolver.smootherType = fvmbaseExt.AMG.JACOBI
#elecSolver.preconditioner = epc
elecSolver.relativeTolerance = 1e-3
elecSolver.nMaxIterations = 1000
elecSolver.maxCoarseLevels=0
elecSolver.verbosity=0

eoptions = emodel.getOptions()
eoptions.electrostaticsLinearSolver = elecSolver
eoptions.electrostaticsTolerance = 0.5e-5
eoptions.electrostatics_enable = 1
eoptions.chargetransport_enable = 0
eoptions.tunneling = 0
eoptions.ibm_enable = 1
eoptions.transient_enable = False
eoptions.printNormalizedResiduals = True

### flow solver ###
momPC = fvmbaseExt.AMG()
momPC.verbosity=0
momSolver = fvmbaseExt.BCGStab()
momSolver.preconditioner = momPC 
momSolver.relativeTolerance = 1e-3
momSolver.absoluteTolerance = 1e-50
momSolver.nMaxIterations = 1000
momSolver.maxCoarseLevels=20
momSolver.verbosity=0


contPC = fvmbaseExt.AMG()
contPC.verbosity=0
contSolver = fvmbaseExt.BCGStab()
contSolver.preconditioner = momPC 
contSolver.relativeTolerance = 1e-1
contSolver.absoluteTolerance = 1e-50
contSolver.nMaxIterations = 1000
contSolver.maxCoarseLevels=20
contSolver.verbosity=0





foptions = fmodel.getOptions()
foptions.momentumLinearSolver = momSolver
foptions.pressureLinearSolver = contSolver
foptions.momentumTolerance=1e-4
foptions.continuityTolerance=1e-4
foptions.setVar("momentumURF",0.7)
foptions.setVar("pressureURF",0.3)
foptions.transient=True
foptions.setVar("timeStep",timeStep)
foptions.printNormalizedResiduals=False
flowVCMap = fmodel.getVCMap()

for id, vc in flowVCMap.iteritems():
    vc['density'] = 1.225
    vc['viscosity'] = 1.8e-5


### initialize models and run ###
emodel.init()
fmodel.init()

ibManager = fvmbaseExt.IBManager(geomFields,
                                 solidBoundaryMeshes[0],
                                 fluidMeshes)
for mesh in solidBoundaryMeshes:
    faces = mesh.getFaces()
    areaMag = geomFields.areaMag[faces]
    faceCount = faces.getCount()
    pot = areaMag.newSizedClone(faceCount)
    pota = pot.asNumPyArray()
    pota[:] = 0.0
    elecFields.potential[faces] = pot

    area = geomFields.area[faces]
    vel = area.newSizedClone(faceCount)
    vela = vel.asNumPyArray()
    vela[:,:] = 0.0
    flowFields.velocity[faces] = vel


sbMeshFaces = solidBoundaryMeshes[0].getFaces()
#ibManager.fluidNeighborsPerIBFace = 4
ibManager.solidNeighborsPerIBFace = 4
ibManager.fluidNeighborsPerSolidFace = 6
ibManager.update()  
fluidMetricsCalculator.computeIBInterpolationMatrices(sbMeshFaces)
fluidMetricsCalculator.computeSolidInterpolationMatrices(sbMeshFaces)     
  
elecSolver.redirectPrintToFile("convergence.dat")
for i in range(0, 10):
     emodel.computeIBFacePotential(sbMeshFaces)
     emodel.advance(1)
elecSolver.redirectPrintToScreen()     


if options.time:
    solver_end[0]  = MPI.Wtime()
    solver_time[0] = solver_end[0] - solver_start[0]
    MPI.COMM_WORLD.Allreduce( [solver_time,MPI.DOUBLE], [solver_maxtime, MPI.DOUBLE], op=MPI.MAX) 
    MPI.COMM_WORLD.Allreduce( [solver_time,MPI.DOUBLE], [solver_mintime, MPI.DOUBLE], op=MPI.MIN) 
    if MPI.COMM_WORLD.Get_rank() == 0:
        dumpMPITimeProfile(solver_maxtime, solver_mintime)

