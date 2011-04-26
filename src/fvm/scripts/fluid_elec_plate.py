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
import ClientCouplingPlate
import time
import optparse
from mpi4py import MPI


etype = {
        'tri' : 1,
        'quad' : 2,
        'tetra' : 3,
        'hexa' : 4
        }


	


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
parser.add_option("--typeFluid", help="'tri'[default], 'quad', 'hexa', or 'tetra'")
parser.add_option("--typeSolid", help="'tri'[default], 'quad', 'hexa', or 'tetra'")
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
fluidTop = 4
fluidBot = 7
fluidSide = 3
fluidLeft = 6
fluidRight = 5
#electrode = 7


numTimeSteps = 2 
globalTime = 0
globalCount = 0
timeStep = 2e-8
saveFrequency = 1
initialTransient = False

### ===================== mesh read ===============================================###

fileBase_input  = "/home/yildirim/memosa/src/fvm/scripts/cantilever3D_coupling/"
fileBase_output =  "./" + str(int(-applied_voltage)) + "/"

print fileBase_input+"fluid_3D_new.cas" 
### 3D fluid mesh
fluidReader = FluentCase(fileBase_input+"fluid_3D_new.cas")
fluidReader.read();
fluent_meshes_fluid = fluidReader.getMeshList()
#paritioning
nmesh = 1
npart = [MPI.COMM_WORLD.Get_size()]
print "options.typeFluid = ", options.typeFluid
etypeFluid = [etype[options.typeFluid]]
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

geomFields =  models.GeomFields('geom')
fluidMetricsCalculator = models.MeshMetricsCalculatorA(geomFields,fluidMeshes)
fluidMetricsCalculator.init()

#### generate solid boundary mesh
#### 2D plate mesh
#fluent_meshes_solid =  beamReader.getMeshList()
##paritioning
#nmesh = 1
#npart = [1] 
#print "options.typeSolid = ", options.typeSolid
#etypeSolid = [etype[options.typeSolid]]
##partMesh constructor and setTypes
#part_mesh_solid = fvmparallel.MeshPartitioner( fluent_meshes_solid, npart, etypeSolid );
#part_mesh_solid.setWeightType(0);
#part_mesh_solid.setNumFlag(0);
##actions
#part_mesh_solid.isDebug(0)
#part_mesh_solid.partition()
#part_mesh_solid.mesh()
#solidMeshes  = part_mesh_solid.meshList()
#if not MPI.COMM_WORLD.Get_rank():
   #print "partition is done for Solid Mesh"

beamReader = FluentCase(fileBase_input+"beam_2D.cas")
beamReader.read();
solidMeshes =  beamReader.getMeshList()
solidBoundaryMeshes = [m.extrude(1, beam_thickness, True) for m in solidMeshes]
solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(geomFields,solidBoundaryMeshes)
solidBoundaryMetricsCalculator.init()


### output files 
probeFile = open(fileBase_output + "centerDisplacement.dat", "w")
forceFile = open(fileBase_output + "force.dat", "w")


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
             x = xf[f][0]
             y = xf[f][1]
             z = xf[f][2]
             if x >= 100e-6 and x<=150e-6 and y>=-15e-6 and y<=15e-6:
               pota[f] = applied_voltage
             else:
               pota[f] = 0.0
         pf[faces] = pot

if fluidBot in bcMap:
   bc = bcMap[fluidBot]
   bc.bcType = "SpecifiedPotential"
   bc['specifiedPotential'] = applied_voltage
   bc['specifiedPotential'] = pf

for i in [fluidSide, fluidTop, fluidLeft, fluidRight]:
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
epc = fvmbaseExt.AMG()
epc.verbosity=0
elecSolver = fvmbaseExt.BCGStab()
elecSolver.preconditioner = epc
elecSolver.relativeTolerance = 1e-3
elecSolver.nMaxIterations = 1000
elecSolver.maxCoarseLevels=20
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
    pota[:] = 0
    elecFields.potential[faces] = pot

    area = geomFields.area[faces]
    vel = area.newSizedClone(faceCount)
    vela = vel.asNumPyArray()
    vela[:,:] = 0.0
    flowFields.velocity[faces] = vel


sbMeshFaces = solidBoundaryMeshes[0].getFaces()
ibManager.fluidNeighborsPerIBFace = 4
ibManager.solidNeighborsPerIBFace = 4
ibManager.fluidNeighborsPerSolidFace = 6
ibManager.update()
#checkMarking(globalCount)


t1 = time.time()
  
#creating client object for coupling
fluid_to_solid = ClientCouplingPlate.ClientCoupling(solidBoundaryMeshes, geomFields, flowFields) 
#--------------Timestep Loop --------------------------#

for n in range(0, numTimeSteps):                
   
    # --------------- update IBM -------------------------#
    print "***       update IBM  at globalCount %i           ***" % globalCount            
    
    ibManager.update()
    fluidMetricsCalculator.computeIBInterpolationMatrices(sbMeshFaces)
    fluidMetricsCalculator.computeSolidInterpolationMatrices(sbMeshFaces)        

    #------------solve electrostatics--------#
    print "***    solving electric model  at globalCount %i  ***" % globalCount
    for i in range(0, 10):
        emodel.computeIBFacePotential(sbMeshFaces)
        emodel.advance(1)

    #---------------solve fluid --------------------------------------------------------------#
    print "***      solving flow model  at globalCount %i    ***" % globalCount
    for i in range(0,10):
        fmodel.computeIBFaceVelocity(sbMeshFaces)
        fmodel.advance(1)

        
    #------------update force on beam  ----------#
    print "***     update force at globalCount %i             ***" % globalCount
    
    solidMesh = solidMeshes[0]
    solidCells = solidMesh.getCells()
    nCells = solidCells.getCount()
    nSelfCells = solidCells.getSelfCount()
    
    nSBFaces = sbMeshFaces.getCount()

    if (nSBFaces != 2*nSelfCells+(nCells-nSelfCells)):
        print "the extruded solid boundary mesh has wrong face numbers!"

    fluid_to_solid.accept()
    solidBoundaryMetricsCalculator.recalculate_deform()
    fluid_to_solid.update([fmodel, emodel], [flowFields, elecFields])
    
    
    # -----------------update time --------------------------#
    fmodel.updateTime()
    globalTime += timeStep
    globalCount += 1

    if (n%saveFrequency == 0):
        saveVTK(n)
   #    checkMarking(globalCount)
 
t2 = time.time()



probeFile.close()
forceFile.close()
print  '\nsolution time = %f' % (t2-t1)
