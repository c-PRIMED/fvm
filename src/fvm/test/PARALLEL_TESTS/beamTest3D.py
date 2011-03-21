#!/usr/bin/env python

import sys
import fvm
import fvm.fvmbaseExt as fvmbaseExt
fvm.set_atype('double')
import fvm.importers as importers
import fvm.fvmparallel as fvmparallel
import time
import math
from numpy import *
from mpi4py import MPI

from FluentCase import FluentCase
from optparse import OptionParser

etype = {
        'tri' : 1,
        'quad' : 2,
        'tetra' : 3,
        'hexa' : 4
        }

def dumpMPITimeProfile(part_mesh_maxtime, part_mesh_mintime, solver_maxtime, solver_mintime):
    fname = "time_mpi_totalprocs%s.dat" % MPI.COMM_WORLD.Get_size()
    f = open(fname,'w')
    line = " part_mesh_mintime = " + str(part_mesh_mintime[0]) + "\n" + \
        " part_mesh_maxtime = " + str(part_mesh_maxtime[0]) + "\n" + \
        " solver_mintime    = " + str(solver_mintime[0])    + "\n" + \
        " solver_maxtime    = " + str(solver_maxtime[0])    + "\n"
    print line
    f.write(line)
    f.close()

class ProblemDefinition(object):
    pass

def writeProbeData(pd, meshes):
    mesh = meshes[0]
    cells = mesh.getCells()
    deformation = pd.structureFields.deformation[cells].asNumPyArray()[pd.probeIndex]
    pd.probeFile.write('%le %le %le\n' % (pd.globalTime, deformation[0], deformation[1]))
    pd.probeFile.flush()

    pd = ProblemDefinition()

def writeDisplacement(pd,meshes):
    mesh = meshes[0]
    faceCells = mesh.getAllFaceCells()
    deformation =  pd.structureFields.deformation[mesh.getCells()].asNumPyArray()
    fileName = pd.fileBase + "/deformation.txt"
    file = open(fileName,"w")

    fgs = mesh.getBoundaryGroups()
    for fg in fgs:
        nFaces = fg.site.getCount()
        xf =  pd.geomFields.coordinate[fg.site].asNumPyArray()
        if fg.id==5:
            faceCells = mesh.getFaceCells(fg.site)
            for i in range(0,nFaces):
                x = xf[i][0]
                y = xf[i][1]
                def1 = deformation[faceCells(i,0)][1]
                file.write(" %e " % x)
                file.write(" %e " % def1)
                file.write("\n")
    file.close() 


    
def advanceUnsteady(pd, nTimeSteps,meshes):
    for i in range(0,nTimeSteps):
        writeProbeData(pd,meshes)
        converged = pd.smodel.advance(1)
        if converged:
            return
        pd.smodel.updateTime()

  
        pd.globalTime += pd.timeStep
        pd.globalCount += 1




parser = OptionParser()
parser.set_defaults(type='tri')
parser.add_option("--type", help="'tri'[default], 'quad', 'hexa', or 'tetra'")
parser.add_option("--xdmf", action='store_true', help="Dump data in xdmf")
parser.add_option("--time","-t",action='store_true',help="Print timing information.")
(options, args) = parser.parse_args()
if len(args) != 1:
    usage()

pd = ProblemDefinition()
pd.fileBase = "."

reader = FluentCase(args[0])
reader.read();
pd.meshes = reader.getMeshList()

nmesh = 1
npart = [MPI.COMM_WORLD.Get_size()]
etype = [etype[options.type]]
if not MPI.COMM_WORLD.Get_rank():
   print "parmesh is processing"
     
    

if options.time:
    # time profile for partmesh
    part_mesh_time     = zeros(1,dtype='d')
    part_mesh_start    = zeros(1, dtype='d')
    part_mesh_end      = zeros(1, dtype='d')
    part_mesh_maxtime  = zeros(1,dtype='d')
    part_mesh_mintime  = zeros(1, dtype='d')
    part_mesh_start[0] = MPI.Wtime()

#partMesh constructor and setTypes
part_mesh = fvmparallel.MeshPartitioner( pd.meshes, npart, etype );
part_mesh.setWeightType(0);
part_mesh.setNumFlag(0);

#actions
part_mesh.isCleanup(1)
part_mesh.isDebug(0)
part_mesh.partition()
part_mesh.mesh()
solid_meshes = part_mesh.meshList()
reader = 0
pd.meshes = 0
if options.time:
    part_mesh_end[0] = MPI.Wtime()
    part_mesh_time[0] = part_mesh_end[0] - part_mesh_start[0]
    MPI.COMM_WORLD.Allreduce( [part_mesh_time,MPI.DOUBLE], [part_mesh_maxtime, MPI.DOUBLE], op=MPI.MAX) 
    MPI.COMM_WORLD.Allreduce( [part_mesh_time,MPI.DOUBLE], [part_mesh_mintime, MPI.DOUBLE], op=MPI.MIN) 
    solver_start   = zeros(1, dtype='d')
    solver_end     = zeros(1, dtype='d')
    solver_time    = zeros(1, dtype='d')
    solver_maxtime = zeros(1, dtype='d')
    solver_mintime = zeros(1, dtype='d')
    solver_start[0] = MPI.Wtime()


mesh = solid_meshes[0]
# cell sites
cellSites = []
cellSites.append( mesh.getCells() )

pd.geomFields =  fvm.models.GeomFields('geom')
cellCentroids = []
pd.metricsCalculator = fvm.models.MeshMetricsCalculatorA(pd.geomFields,solid_meshes)


pd.metricsCalculator.init()
cellCentroids.append( pd.geomFields.coordinate[cellSites[0]].asNumPyArray() )
volume = pd.geomFields.volume[cellSites[0]].asNumPyArray() 
ncells = cellSites[0].getCount()

#if ( MPI.COMM_WORLD.Get_rank() == 0 ):
#   for n in range(0,ncells):
#      print "n = ", n, " cellCentroidX = ", cellCentroids[0][n][0], " cellCentroidY = ", cellCentroids[0][n][1], " volume = ", volume[n]


cells = mesh.getCells()

rho = 7854.0
E = 1.82e11
nu = 0.0#0.31

pd.structureFields =  fvm.models.StructureFields('structure')
pd.smodel = fvm.models.StructureModelA(pd.geomFields,pd.structureFields,solid_meshes)

bcMap = pd.smodel.getBCMap()
print bcMap.keys()
#back 
if 4 in bcMap:
   bc = bcMap[4]
   bc.bcType = 'SpecifiedDistForce'
   bc['specifiedYDistForce']=0.0

#front 
if 3 in bcMap:
   bc = bcMap[3]
   bc.bcType = 'SpecifiedDistForce'
   bc['specifiedYDistForce']=0.0



#top 
if 5 in bcMap:
   bc = bcMap[5]
   bc.bcType = 'SpecifiedDistForce'
   bc['specifiedYDistForce']=-2.77e3

#bot 
if 6 in bcMap:
   bc = bcMap[6]
   bc.bcType = 'SpecifiedDistForce'

#left
if 8 in bcMap:
  bc = bcMap[8]
  bc.bcType = 'SpecifiedDeformation'

# right
if 7 in bcMap:
   bc = bcMap[7]
   bc.bcType = 'SpecifiedDeformation'

vcMap = pd.smodel.getVCMap()
for i,vc in vcMap.iteritems():
    vc['density'] = rho
    vc['eta'] = E/(2.*(1+nu))
    vc['eta1'] = nu*E/((1+nu)*(1-nu))



defSolver = fvmbaseExt.AMG()
defSolver.smootherType = fvmbaseExt.AMG.JACOBI
defSolver.relativeTolerance = 1e-5
defSolver.nMaxIterations = 100
defSolver.maxCoarseLevels=0
defSolver.verbosity=2

soptions = pd.smodel.getOptions()
soptions.deformationLinearSolver = defSolver
soptions.deformationTolerance=1.e-6
soptions.setVar("deformationURF",1)

pd.timeStep = 1e-7
soptions.setVar("timeStep", pd.timeStep)

#soptions.setVar("dampingCoeff",1e6)

soptions.printNormalizedResiduals=False
soptions.transient=False

pd.smodel.init()
pd.globalTime = 0
pd.globalCount = 0

pd.probeIndex = 1000

#pd.smodel.dumpMatrix("beam3200")
defSolver.redirectPrintToFile("convergence.dat")
pd.smodel.advance(1)
defSolver.redirectPrintToScreen()
if options.time:
    solver_end[0] = MPI.Wtime()
    solver_time[0] = solver_end[0] - solver_start[0]
    MPI.COMM_WORLD.Allreduce( [solver_time,MPI.DOUBLE], [solver_maxtime, MPI.DOUBLE], op=MPI.MAX) 
    MPI.COMM_WORLD.Allreduce( [solver_time,MPI.DOUBLE], [solver_mintime, MPI.DOUBLE], op=MPI.MIN) 
    if MPI.COMM_WORLD.Get_rank() == 0:
        dumpMPITimeProfile(part_mesh_maxtime, part_mesh_mintime, solver_maxtime, solver_mintime)



#pd.probeFile = open(pd.fileBase + "tipDisplacement.dat", "w")

#advanceUnsteady(pd,1000,mesh)

#writeDisplacement(pd,solid_meshes)
