#!/usr/bin/env python

### Structure and Electrostatics coupling --- 2D beam  ###

### import  modules ###
import pdb
import sys
import os
from math import *
#import tecplotExporter
import fvm
fvm.set_atype('double')
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
import fvm.fvmparallel as fvmparallel
from FluentCase import FluentCase
from fvm.fvmbaseExt import VecD3
import ClientCoupling
from mpi4py import MPI


etype = {
        'tri' : 1,
        'quad' : 2,
        'tetra' : 3,
        'hexa' : 4
        }

def dumpMPITimeProfile(part_mesh_maxtime, part_mesh_mintime, solver_maxtime, solver_mintime):
    fname = fileBaseOutput + "time_mpi_totalprocs%s.dat" % MPI.COMM_WORLD.Get_size()
    f = open(fname,'w')
    line = " part_mesh_mintime = " + str(part_mesh_mintime[0]) + "\n" + \
           " part_mesh_maxtime = " + str(part_mesh_maxtime[0]) + "\n" + \
           " solver_mintime    = " + str(solver_mintime[0])    + "\n" + \
           " solver_maxtime    = " + str(solver_maxtime[0])    + "\n"
    print line
    f.write(line)
    f.close()

	

import time


### ========================== define functions and class =============================###
class ProblemDefinition(object):
    pass

def writeForce(n, pd):
    procID = MPI.COMM_WORLD.Get_rank()
    forceFile = open( fileBaseOutput + "elecForce_" + str(pd.globalCount) + "_procID" + str(procID) +  ".dat", "w")
    sbMeshFaces = pd.solidBoundaryMeshes[0].getFaces()
    sbElecForce = pd.elecFields.force[sbMeshFaces].asNumPyArray()
    sbFlowForce = pd.flowFields.force[sbMeshFaces].asNumPyArray()
    xFaces = pd.geomFields.coordinate[sbMeshFaces].asNumPyArray()
    nfaces = sbMeshFaces.getCount()
    for i in range (0, nfaces):
        forceFile.write('%e %e %e %e \n ' % (xFaces[i][0], xFaces[i][1], sbElecForce[i][0], sbElecForce[i][1]))
    forceFile.close()

def writeDeformationProfile(n, pd):
    deformationFile = open(fileBaseOutput + "deformation_" + str(pd.globalCount) + ".dat", "w")
    sbMeshFaces = pd.solidBoundaryMeshes[0].getFaces()
    xFaces = pd.geomFields.coordinate[sbMeshFaces].asNumPyArray()
    nfaces = sbMeshFaces.getCount()
    for i in range (0, nfaces):
        deformationFile.write('%i %e %e %e\n ' % (i, xFaces[i][0], xFaces[i][1], xFaces[i][2]))        
    deformationFile.close()

def writeVelocityProfile(n, pd):
    velocityFile = open(fileBaseOutput + "velocity_" + str(pd.globalCount) + ".dat", "w")
    sbMeshFaces = pd.solidBoundaryMeshes[0].getFaces()
    velocity = pd.flowFields.velocity[sbMeshFaces].asNumPyArray()
    nfaces = sbMeshFaces.getCount()
    for i in range (0, nfaces):
        velocityFile.write('%i %e %e %e\n ' % (i, velocity[i][0], velocity[i][1], velocity[i][2]))        
    velocityFile.close()

def saveVTK(nstep, pd):
    if pd.initialTransient:
        return
    procID = MPI.COMM_WORLD.Get_rank() 
    writer2 = exporters.VTKWriterA(pd.geomFields,pd.fluidMeshes,
                                  fileBaseOutput + "elecfield-" + str(nstep) + "_proc"+str(procID) + ".vtk",
                                  "fix-fix beam",
                                  False,0)
    writer2.init()
    writer2.writeScalarField(pd.elecFields.potential,"potential")
    writer2.writeVectorField(pd.elecFields.electric_field,"potentialgradient")
    writer2.writeVectorField(pd.flowFields.velocity,"velocity")
    writer2.writeScalarField(pd.flowFields.pressure, "pressure")
    writer2.finish()

   
   
    writer3 = exporters.VTKWriterA(pd.geomFields,pd.solidBoundaryMeshes,
                                  fileBaseOutput+ "beamBoundary-" + str(nstep) + "_proc" + str(procID) +  ".vtk",
                                  "beam Boundary",
                                  False,0,True)
    writer3.init()
    writer3.writeVectorField(pd.flowFields.velocity,"velocity")
    writer3.writeVectorField(pd.flowFields.force,"flow_force")
    writer3.writeVectorField(pd.elecFields.force,"elec_force")
    writer3.finish()
    
    
           
def unsteadyAdvance(pd, numTimeSteps):
    
    for n in range(0, numTimeSteps):


        sbMeshFaces = pd.solidBoundaryMeshes[0].getFaces()

        #-------------update IBM----------------#
        print "***       update IBM  at  %i           ***" % n
        pd.ibManager.update()
        #checkMarking(n, pd)
        pd.ibManager.fluidNeighborsPerIBFace = 2 
        pd.ibManager.solidNeighborsPerIBFace = 4
        pd.ibManager.fluidNeighborsPerSolidFace =4 
        
        pd.fluidMetricsCalculator.computeIBInterpolationMatrices(sbMeshFaces)
        pd.fluidMetricsCalculator.computeSolidInterpolationMatrices(sbMeshFaces)     
        
        #------------solve electrostatics--------#
        print "***    solving electric model  at globalCount %i  ***" % pd.globalCount
        for i in range(0, 10):
            pd.emodel.computeIBFacePotential(sbMeshFaces)
            pd.emodel.advance(1)
                #break
        
        #------------solve fluid ---------------#
        print "***      solving flow model  at globalCount %i    ***" % pd.globalCount
        pd.fmodel.computeIBFaceVelocity(sbMeshFaces)
        for i in range(0,1):
            ibfaces = pd.fluidMeshes[0].getIBFaces()
            ibv = pd.flowFields.velocity[ibfaces].asNumPyArray()
            sbv = pd.flowFields.velocity[sbMeshFaces].asNumPyArray()
            if pd.fmodel.advance(0):
                break

        #---------------write out data -------------------------#
        # writeProbeData(pd)
        
        fluid_to_solid.accept(pd.flowFields) #get coord and boundary velocities
        pd.solidBoundaryMetricsCalculator.recalculate_deform()
	fluid_to_solid.update([pd.fmodel, pd.emodel],[pd.flowFields, pd.elecFields])
        #if (n%pd.saveFrequency == 0):
            #saveVTK(n, pd)
            #checkMarking(n, pd)
            #writeForce(n, pd)
            #writeDeformationProfile(n, pd)
            #writeVelocityProfile(n,pd)
       
        #---------------update time -------------------------#
        
        pd.fmodel.updateTime()
        
        pd.globalTime += pd.timeStep
        pd.globalCount += 1

def checkMarking(n, pd):

    cells = pd.fluidMeshes[0].getCells()
    nCells = cells.getCount()
    cellCoords = pd.geomFields.coordinate[cells].asNumPyArray()
    procID = MPI.COMM_WORLD.Get_rank()
    fluidFile = open(fileBaseOutput + "fluidCells_" + str(pd.globalCount) + "_" + str(procID)+ ".dat", "w")
    solidFile = open(fileBaseOutput + "solidCells_" + str(pd.globalCount) + "_" + str(procID)+ ".dat", "w")
    IBFile    = open(fileBaseOutput + "IBCells_" + str(pd.globalCount) + "_" + str(procID)+ ".dat", "w")

    cellIBType = pd.geomFields.ibType[cells].asNumPyArray()
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

        
### ========================== properties and parameters ===============================###


### electric field
applied_voltage = -100
dielectric_constant = 1.0


numTimeSteps = 10

### =========================== running script  ========================================###

pd = ProblemDefinition()


fileBase = "../"
fileBaseOutput = "./"

pd.globalTime = 0
pd.globalCount = 0
pd.timeStep = 1.e-8
pd.saveFrequency = 2
pd.fxSum = 0
pd.fySum = 0
pd.fzSum = 0
pd.initialTransient=False
#pd.probeIndex = 7106 #tip location at proc2 for 3processors run

### read in meshes ###
#beamFile = fileBase + 'frogleg_device_252k.cas'
#fluidFile = fileBase + 'frogleg_fluid_300k.cas'

beamFile  = fileBase + 'beam_500x10.cas'
fluidFile = fileBase + 'fluid_500x40.cas'


beamReader = FluentCase(beamFile)
beamReader.read()
#all fluid process will be aware of whole meshes
pd.solidMeshes = beamReader.getMeshList()

#partitioning
fluidReader = FluentCase(fluidFile)
fluidReader.read()
nmesh = 1
npart = [MPI.COMM_WORLD.Get_size()]
etype = [etype['quad']]
#partMesh constructor and setTypes
fluent_meshes = fluidReader.getMeshList()
part_mesh = fvmparallel.MeshPartitioner( fluent_meshes, npart, etype );
part_mesh.setWeightType(0);
part_mesh.setNumFlag(0);
#actions
part_mesh.isDebug(0)
part_mesh.partition()
part_mesh.mesh()
pd.fluidMeshes  = part_mesh.meshList()
if not MPI.COMM_WORLD.Get_rank():
   print "partition is done"

#boundaryMeshes
pd.solidBoundaryMeshes = [m.extractBoundaryMesh() for m in pd.solidMeshes]

### geometry field ###
pd.geomFields =  models.GeomFields('geom')
pd.fluidMetricsCalculator = models.MeshMetricsCalculatorA(pd.geomFields,
                                                          pd.fluidMeshes)
pd.fluidMetricsCalculator.init()

pd.solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(pd.geomFields,
                                                          pd.solidBoundaryMeshes)
pd.solidBoundaryMetricsCalculator.init()


## electric model and boundary condition ###
pd.elecFields =  models.ElectricFields('elec')
pd.emodel = models.ElectricModelA(pd.geomFields,
                                     pd.elecFields,
                                     pd.fluidMeshes)
### mesh id
fluidTop = 8
fluidBot = [5,4,3]
fluidLeft = 6
fluidRight = 7

bcMap = pd.emodel.getBCMap()

for id in [fluidBot[1]]:
   if id in bcMap:
      bc = bcMap[id]
      bc.bcType = "SpecifiedPotential"
      bc['specifiedPotential'] = applied_voltage

for id in [fluidTop, fluidRight, fluidLeft, fluidBot[0], fluidBot[2] ]:
    if id in bcMap:
       bc = bcMap[id]
       bc.bcType = "Symmetry"
 
  
vcMap = pd.emodel.getVCMap()
for i,vc in vcMap.iteritems():
    vc.vcType = "dielectric"
    vc['dielectric_constant'] = dielectric_constant


### flow model and boundary condition ###
pd.flowFields =  models.FlowFields('flow')
pd.fmodel = models.FlowModelA(pd.geomFields,pd.flowFields,pd.fluidMeshes)

fluidReader.importFlowBCs(pd.fmodel, pd.fluidMeshes)

#creating client object for coupling
fluid_to_solid = ClientCoupling.ClientCoupling(pd.solidBoundaryMeshes, pd.geomFields) 

# electric solver 
epc = fvmbaseExt.AMG()
epc.verbosity=0
elecSolver = fvmbaseExt.BCGStab()
elecSolver.preconditioner = epc
elecSolver.relativeTolerance = 1e-3
elecSolver.nMaxIterations = 1000
elecSolver.maxCoarseLevels=20
elecSolver.verbosity=0

eoptions = pd.emodel.getOptions()
eoptions.electrostaticsLinearSolver = elecSolver
eoptions.electrostaticsTolerance = 0.5e-5
eoptions.electrostatics_enable = 1
eoptions.chargetransport_enable = 0
eoptions.tunneling = 0
eoptions.ibm_enable = 1
eoptions.transient_enable = False
eoptions.printNormalizedResiduals = True

# flow solver
foptions = pd.fmodel.getOptions()
foptions.momentumTolerance=1e-4
foptions.continuityTolerance=1e-4
foptions.setVar("momentumURF",0.7)
foptions.setVar("pressureURF",0.3)
foptions.transient=True

foptions.setVar("timeStep",pd.timeStep)

flowVCMap = pd.fmodel.getVCMap()
for id, vc in flowVCMap.iteritems():
    vc['density'] = 1.225
    vc['viscosity'] = 1.8e-5
foptions.printNormalizedResiduals=False

### set up IBM ###

for mesh in pd.solidBoundaryMeshes:
    faces = mesh.getFaces()
    #potential field on boundary mesh
    areaMag = pd.geomFields.areaMag[faces]
    faceCount = faces.getCount()
    pot = areaMag.newSizedClone(faceCount)
    pota = pot.asNumPyArray()
    pota[:] = 0.0
    pd.elecFields.potential[faces] = pot
    #flow velocity field on boundary mesh
    area = pd.geomFields.area[faces]
    faceCount = faces.getCount()
    vel = area.newSizedClone(faceCount)
    vela = vel.asNumPyArray()
    vela[:,:] = 0.0
    pd.flowFields.velocity[faces] = vel
    
pd.ibManager = fvmbaseExt.IBManager(pd.geomFields,
                                    pd.solidBoundaryMeshes[0],
                                    pd.fluidMeshes)
pd.ibManager.fluidNeighborsPerIBFace = 2
pd.ibManager.solidNeighborsPerIBFace = 6
pd.ibManager.fluidNeighborsPerSolidFace = 6

pd.ibManager.update()

### initialize models and run ###

#checkMarking(0, pd)

pd.emodel.init()
pd.fmodel.init()

#initially, zero force is applied on beam
#createSolidForceBVFields(pd)

t1 = time.time()
                    
unsteadyAdvance(pd,numTimeSteps)
print "FLUID DONE"

#checkMarking(pd.globalCount, pd)

#t2 = time.time()

#print  '\nsolution time = %f' % (t2-t1)

