#! /usr/bin/env python

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
import ServerCoupling
from mpi4py import MPI



etype = {
        'tri' : 1,
        'quad' : 2,
        'tetra' : 3,
        'hexa' : 4
        }

def dumpMPITimeProfile(part_mesh_maxtime, part_mesh_mintime, solver_maxtime, solver_mintime):
    fname = fileBaseOutput+"time_mpi_totalprocs%s.dat" % MPI.COMM_WORLD.Get_size()
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

def writeProbeData(pd):
    cells = pd.solidMeshes[0].getCells()
    deformation = pd.structureFields.deformation[cells].asNumPyArray()[pd.probeIndex]
    #velocity = pd.flowFields.velocity[cells].asNumPyArray()[pd.probeIndex]
    pd.probeFile.write('%le %le %le %le\n' % (pd.globalTime, deformation[0], deformation[1], deformation[2]))
    pd.probeFile.flush()
    #pd.velocityFile.write('%le %le %le %le\n' % (pd.globalTime, velocity[0], velocity[1], velocity[2]))
    #pd.velocityFile.flush()
    #pd.forceFile.write('%le %le %le %le\n' % (pd.globalTime, pd.fxSum, pd.fySum, pd.fzSum))
    #pd.forceFile.flush()

def saveVTK(nstep, pd):
    if pd.initialTransient:
        return
    writer = exporters.VTKWriterA(pd.geomFields,pd.solidMeshes,
                                  fileBaseOutput + "beam-deformation-" + str(nstep) + ".vtk",
                                  "fix-fix beam",
                                  False,0)
    writer.init()
    writer.writeVectorField(pd.structureFields.deformation,"deformation")
    #writer.writeScalarField(pd.structureFields.sigmaXX,"sigmaXX")
    #writer.writeScalarField(pd.structureFields.sigmaXY,"sigmaXY")
    #writer.writeScalarField(pd.structureFields.sigmaYY,"sigmaYY")
    writer.finish()
    
        
           
def unsteadyAdvance(pd, numTimeSteps):
    
    for n in range(0, numTimeSteps):
        #------------solve structure-------------#
        print "***  solving structure model at globalCount %i   ***" % pd.globalCount
        for i in range(0, 3):
            converged = pd.smodel.advance(1)
            pd.dmodel.calculateNodeDisplacement()
            pd.dmodel.deformStructure()
            pd.solidMetricsCalculator.recalculate_deform()
            if converged:
                break

        #---------------write out data -------------------------#
        writeProbeData(pd)

        #if (n%pd.saveFrequency == 0):
            #saveVTK(n, pd)
            #writeForce(n, pd)
            #writeDeformationProfile(n, pd)
            #writeVelocityProfile(n,pd)
       
        #---------------update time -------------------------#
        
        solid_to_fluid.update() #first send coords and vels  to fluid side
	solid_to_fluid.accept() #then wait for fluid side to send you strees to continue
        pd.smodel.updateTime()
        pd.dmodel.updateTime()
        
        pd.globalTime += pd.timeStep
        pd.globalCount += 1
       
## ========================== properties and parameters ===============================###
rho = 8912                     # density kg/m^3
E = 200e9                     # Young's modulus 
nu = 0.31                      # Poisson's ratio

beamTop = 6
beamBot = 5
beamRight = 4
beamLeft = 3

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
pd.probeIndex = 5015 #tip index 
pd.probeFile    = open(fileBaseOutput+"tipDisplacement-se.dat", "w")
#pd.forceFile    = open(fileBaseOutput+"beamForce-se.dat", "w")
#pd.velocityFile = open(fileBaseOutput+"tipVelocity-se.dat", "w")

### read in mesh ###
beamFile = fileBase + 'beam_500x10.cas'

beamReader = FluentCase(beamFile)
beamReader.read()
solid_meshes = beamReader.getMeshList() 

nmesh = 1
npart = [MPI.COMM_WORLD.Get_size()]
etype = [etype['quad']]
#partMesh constructor and setTypes
part_mesh = fvmparallel.MeshPartitioner( solid_meshes, npart, etype );
part_mesh.setWeightType(0);
part_mesh.setNumFlag(0);
#actions
part_mesh.isCleanup(0)
part_mesh.isDebug(0)
part_mesh.partition()
part_mesh.mesh()
part_mesh.extractBoundaryMesh()
pd.solidMeshes  = part_mesh.meshList()
pd.solidBoundaryMeshes = [part_mesh.getBoundaryMesh()]
if not MPI.COMM_WORLD.Get_rank():
   print "partition done"


### geometry field ###

pd.geomFields =  models.GeomFields('geom')
pd.solidMetricsCalculator = models.MeshMetricsCalculatorA(pd.geomFields,
                                                          pd.solidMeshes)
pd.solidMetricsCalculator.init()

pd.solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(pd.geomFields,
                                                          pd.solidBoundaryMeshes)							  
pd.solidBoundaryMetricsCalculator.init()

### structure model and boundary condition ###
pd.structureFields =  models.StructureFields('structure')

pd.smodel = models.StructureModelA(pd.geomFields,
                                   pd.structureFields,
                                   pd.solidMeshes)

pd.dmodel = models.StructureDeformationModelA(pd.geomFields,
                                              pd.structureFields,
                                              pd.solidMeshes)

# device boundary condition
bcMap = pd.smodel.getBCMap()
for id in [beamLeft]:
    if id in bcMap:
      bc = bcMap[id]
      bc.bcType = 'SpecifiedDeformation'
      bc['specifiedXDeformation'] = 0.0
      bc['specifiedYDeformation'] = 0.0
      bc['specifiedZDeformation'] = 0.0

for id   in [beamTop,beamBot, beamRight]:    
    if id in bcMap:
       bc = bcMap[id]
       bc.bcType = 'SpecifiedForce'

vcMap = pd.smodel.getVCMap()
for i,vc in vcMap.iteritems():
    vc['density'] = rho
    vc['eta'] = E/(2.*(1+nu))
    vc['eta1'] = nu*E/((1+nu)*(1-2.0*nu))


solid_to_fluid = ServerCoupling.ServerCoupling(pd.smodel, pd.dmodel, pd.geomFields, pd.solidMeshes, pd.solidBoundaryMeshes)


### solver options ###

# structure solver
pc = fvmbaseExt.ILU0Solver()
pc.verbosity=0
#defSolver = fvmbaseExt.BCGStab()
defSolver = fvmbaseExt.CG()
defSolver.preconditioner = pc
defSolver.relativeTolerance = 1e-3
defSolver.nMaxIterations = 2000
defSolver.maxCoarseLevels=20
defSolver.verbosity=0


soptions = pd.smodel.getOptions()
soptions.deformationLinearSolver = defSolver
soptions.deformationTolerance=1.0e-6
soptions.setVar("deformationURF",1.0)
soptions.printNormalizedResiduals=False
soptions.transient=True

soptions.setVar("timeStep", pd.timeStep)


### initialize models and run ###

pd.smodel.init()
pd.dmodel.init()
                    
unsteadyAdvance(pd,numTimeSteps)

print "SOLID DONE"


