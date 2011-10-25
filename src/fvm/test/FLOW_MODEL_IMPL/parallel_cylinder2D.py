#!/usr/bin/env python
### import  modules ###
"""
mpirun -np   1 python ./parallel_cylinder.py  ../test/uniform-cart_5K.cas ../test/circle_66.cas --type=quad
"""
import pdb
import sys
from math import *
import fvm
fvm.set_atype('double')
#import tecplotExporter
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
import fvm.fvmparallel as fvmparallel
from FluentCase import FluentCase
from fvm.fvmbaseExt import VecD3
import time
from mpi4py  import MPI
import optparse
from Tools import *
from ComputeForce import *
from TimeStep import *
from MeshSetup import *
from tecplotParallelFlowField import *
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

    fluidFile = open("./fluidCells_" + str(n) + "_proc" + str(MPI.COMM_WORLD.Get_rank()) + ".dat", "w")
    solidFile = open("./solidCells_" + str(n) + "_proc" + str(MPI.COMM_WORLD.Get_rank()) + ".dat", "w")
    IBFile    = open("./IBCells_" + str(n)    + "_proc" + str(MPI.COMM_WORLD.Get_rank()) + ".dat", "w")

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

def dumpMPITimeProfile(solver_maxtime, solver_mintime):
    fname = "time_mpi_totalprocs%s_fluid.dat" % MPI.COMM_WORLD.Get_size()
    f = open(fname,'w')
    line = " solver_mintime    = " + str(solver_mintime[0])    + "\n" + \
           " solver_maxtime    = " + str(solver_maxtime[0])    + "\n"
    print line
    f.write(line)
    f.close()

def saveFluidVTK(n):
    writer = exporters.VTKWriterA(geomFields,fluidMeshes,
                                  fileBase + "fluid-" + str(n) + ".vtk",
                                  "fluid",
                                  False,0)
    writer.init()
    writer.writeVectorField(flowFields.velocity,"velocity")
    writer.writeScalarField(flowFields.pressure, "pressure")
    writer.finish()


timeStep = 10
numTimeSteps = 10
saveFrequency = 2

parser = optparse.OptionParser()
parser.add_option("--volt", type=float)
parser.add_option("--type", help="'tri'[default], 'quad', 'hexa', or 'tetra'")
parser.add_option("--time","-t",action='store_true',help="Print timing information.")
(options, args) = parser.parse_args()

print args[0]
print args[1]
fluidReader = FluentCase(args[0])
solidReader = FluentCase(args[1])

fluidReader.read()
solidReader.read()

fluidMeshes0 = fluidReader.getMeshList()
solidMeshes = solidReader.getMeshList()
nodeCoord = solidMeshes[0].getNodeCoordinates().asNumPyArray()


nodeCoord[:,:] *=0.5

 #paritioning
nmesh = 1
npart = [MPI.COMM_WORLD.Get_size()]
print "options.typeFluid = ", options.type
etypeFluid = [etype[options.type]]
#partMesh constructor and setTypes
part_mesh_fluid = fvmparallel.MeshPartitioner( fluidMeshes0, npart, etypeFluid );
part_mesh_fluid.setWeightType(0);
part_mesh_fluid.setNumFlag(0);
#actions
part_mesh_fluid.isDebug(0)
part_mesh_fluid.partition()
part_mesh_fluid.mesh()
fluidMeshes  = part_mesh_fluid.meshList()
if not MPI.COMM_WORLD.Get_rank():
   print "partition is done for Fluid Mesh"

if options.time:
    solver_start   = zeros(1, dtype='d')
    solver_end     = zeros(1, dtype='d')
    solver_time    = zeros(1, dtype='d')
    solver_maxtime = zeros(1, dtype='d')
    solver_mintime = zeros(1, dtype='d')
    solver_start[0] = MPI.Wtime()

geomFields =  models.GeomFields('geom')

fluidMetricsCalculator = models.MeshMetricsCalculatorA(geomFields, fluidMeshes)
fluidMetricsCalculator.init()

solidBoundaryMeshes = [m.extractBoundaryMesh() for m in solidMeshes]
solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(geomFields,solidBoundaryMeshes)
solidBoundaryMetricsCalculator.init()

flowFields =  models.FlowFields('flow')
fmodel = models.FlowModelA(geomFields,flowFields,fluidMeshes)

fluidWalls = [3,4]
fluidInlet = [5]
fluidOutlet = [6]

circleWalls = [4]

bcMap = fmodel.getBCMap()

for id in fluidWalls:
   if id in bcMap:
      bc = bcMap[id]
      bc.bcType = 'NoSlipWall'
for id in fluidInlet:
   if id in bcMap:
      bc = bcMap[id]
      bc.bcType = 'VelocityBoundary'
      bc['specifiedXVelocity']=1
      bc['specifiedYVelocity']=0
      bc['specifiedZVelocity']=0
for id in fluidOutlet:
    if id in bcMap:
       bc = bcMap[id]
       bc.bcType = 'PressureBoundary'
vcMap = fmodel.getVCMap()
for i,vc in vcMap.iteritems():
    vc['viscosity'] = 1
    vc['density'] = 1


#fluidReader.importFlowBCs(fmodel,fluidMeshes)


#flow solvers
momPC = fvmbaseExt.AMG()
momPC.verbosity=0
momSolver = fvmbaseExt.BCGStab()
momSolver.preconditioner = momPC 
momSolver.relativeTolerance = 1e-1
momSolver.absoluteTolerance = 1e-50
momSolver.nMaxIterations = 20
momSolver.verbosity=0


contPC = fvmbaseExt.AMG()
contPC.verbosity=0
contSolver = fvmbaseExt.BCGStab()
contSolver.preconditioner = contPC 
contSolver.relativeTolerance = 1e-1
contSolver.absoluteTolerance = 1e-50
contSolver.nMaxIterations = 20
contSolver.verbosity=0

foptions = fmodel.getOptions()
foptions.momentumLinearSolver = momSolver
foptions.pressureLinearSolver = contSolver
foptions.momentumTolerance=1e-5
foptions.continuityTolerance=1e-5
foptions.setVar("momentumURF",0.7)
foptions.setVar("pressureURF",0.3)
foptions.transient=True
foptions.setVar("timeStep",timeStep)
foptions.printNormalizedResiduals=True

fmodel.init()

sbMeshFaces = solidBoundaryMeshes[0].getFaces()    
ibManager = fvmbaseExt.IBManager(geomFields,solidBoundaryMeshes[0],fluidMeshes)
faceCount = sbMeshFaces.getCount()
area = geomFields.area[sbMeshFaces]
velocity = area.newSizedClone(faceCount)
velocitya = velocity.asNumPyArray()
velocitya[:,:] = 0.0
#velocitya[:,1] = 1.0
flowFields.velocity[sbMeshFaces] = velocity
velWhole = flowFields.velocity[fluidMeshes[0].getCells()].asNumPyArray()
#velWhole[:,0] = 0.0


ibManager.solidNeighborsPerIBFace = 2

globalCount = 0
globalTime = 0.0

ibManager.update()
fluidMetricsCalculator.computeIBInterpolationMatrices(sbMeshFaces)
fluidMetricsCalculator.computeSolidInterpolationMatrices(sbMeshFaces)

#redirect stdout to a disk file
saveout = sys.stdout
outfile = open('compare.dat','w')
sys.stdout = outfile

numTimeSteps=1
for n in range(0, numTimeSteps):
    ibManager.update()
    #checkMarking(n)
    fluidMetricsCalculator.computeIBInterpolationMatrices(sbMeshFaces)
    fmodel.computeIBFaceVelocity(sbMeshFaces)
    
    for i in range(0,10):
        fmodel.computeIBFaceVelocity(sbMeshFaces)
        if fmodel.advance(1):
            break

    momPC.redirectPrintToFile("compare1.dat")
    fmodel.printBCs()
	
    #calling getPressureIntegral	     
    print "getPressureIntegral ... \n"
    print fmodel.getPressureIntegral(fluidMeshes[0],3)[0]
    print fmodel.getPressureIntegral(fluidMeshes[0],3)[1]
    print fmodel.getPressureIntegral(fluidMeshes[0],3)[2]
    print "\n"
     
    print fmodel.getPressureIntegral(fluidMeshes[0],4)[0]
    print fmodel.getPressureIntegral(fluidMeshes[0],4)[1]
    print fmodel.getPressureIntegral(fluidMeshes[0],4)[2]
    print "\n"
	    
    #calling getPressureIntegralIBFaces
    print "getPressureIntegralIBFaces ... \n"
    print fmodel.getPressureIntegralonIBFaces(fluidMeshes[0])[0]
    print fmodel.getPressureIntegralonIBFaces(fluidMeshes[0])[1]
    print fmodel.getPressureIntegralonIBFaces(fluidMeshes[0])[2]  	    
    print "\n"
    
    #calling getMomentumFluxIntegra
    print "getMomentumFLuxIntegral ... \n"
    print fmodel.getMomentumFluxIntegral(fluidMeshes[0],3)[0]
    print fmodel.getMomentumFluxIntegral(fluidMeshes[0],3)[1]
    print fmodel.getMomentumFluxIntegral(fluidMeshes[0],3)[2]    
    print "\n"
    
    print fmodel.getMomentumFluxIntegral(fluidMeshes[0],4)[0]
    print fmodel.getMomentumFluxIntegral(fluidMeshes[0],4)[1]
    print fmodel.getMomentumFluxIntegral(fluidMeshes[0],4)[2]    
    print "\n"
    
    #calling getMomentumDerivativeIntegral
    print "getMomentumDerivativeIntegral ... \n"
    print fmodel.getMomentumDerivativeIntegral(fluidMeshes[0])[0]
    print fmodel.getMomentumDerivativeIntegral(fluidMeshes[0])[1]
    print fmodel.getMomentumDerivativeIntegral(fluidMeshes[0])[2]  
    print "\n"  


    #calling getStressTensor
    print "getStressTensor ... \n"
    IBCellIDsA= fvmbaseExt.newIntArray(5)
    IBCellIDs = IBCellIDsA.asNumPyArray()
    IBCellIDs[0] = 0
    IBCellIDs[1] = 10
    IBCellIDs[2] = 100
    IBCellIDs[3] = 200
    IBCellIDs[4] = 512
    stressA = fmodel.getStressTensor(fluidMeshes[0],IBCellIDsA)
    stress = stressA.asNumPyArray()
    print stress
   
    #getTranction
    print "getTraction ... \n"
    fmodel.getTraction(fluidMeshes[0])
    tractionX = flowFields.tractionX[fluidMeshes[0].getCells()]
    print tractionX.asNumPyArray()
    print "\n"
    
   
    #calling print PressureIntegrals()
    fmodel.printPressureIntegrals()
    fmodel.printMomentumFluxIntegrals()
    fmodel.printMassFluxIntegrals()
    
    momPC.redirectPrintToScreen()
    
    
    #fmodel.updateTime()
    globalCount += 1
    globalTime += timeStep


#restore
outfile.flush()
outfile.close()
sys.stdout = saveout

#if options.time:
#    solver_end[0]  = MPI.Wtime()
#    solver_time[0] = solver_end[0] - solver_start[0]
#    MPI.COMM_WORLD.Allreduce( [solver_time,MPI.DOUBLE], [solver_maxtime, MPI.DOUBLE], op=MPI.MAX) 
#    MPI.COMM_WORLD.Allreduce( [solver_time,MPI.DOUBLE], [solver_mintime, MPI.DOUBLE], op=MPI.MIN) 
#    if MPI.COMM_WORLD.Get_rank() == 0:
#        dumpMPITimeProfile(solver_maxtime, solver_mintime)    

#dumpTecplotFile(nmesh, fluidMeshes, fluidMeshes0, options.type,flowFields,"cylinder_results.dat") 
#saveFluidVTK(n)
#pdb.set_trace()
"""
def searchPoint(mesh, probe):
    cells = mesh.getCells()
    xCells = geomFields.coordinate[cells]
    xCellsA = xCells.asNumPyArray()
    cloestPoints = fvmbaseExt.newIntArray(1)
    cloestPointsA = cloestPoints.asNumPyArray()
    target = fvmbaseExt.VecD3()
    target[0] = probe[0]
    target[1] = probe[1]
    target[2] = probe[2]
    search = fvmbaseExt.KSearchTree(xCells)
    search.findNeighbors(target, 1, cloestPoints)
    point = cloestPointsA[0]
    return point
probe = [0, -5, 0]
point = searchPoint(fluidMeshes[0], probe)
line = [point]
dy = 10 / 100.
count = 1
for n in range(0, 101):
    probe[0] = 0
    probe[1] += dy
    probe[2] = 0
    point = searchPoint(fluidMeshes[0], probe)
    if point != line[count-1]:
        line.append(point)
        count +=1


file = open('line_0.dat','w')
velocity = flowFields.velocity[fluidMeshes[0].getCells()].asNumPyArray()
for n in range(0, count):
    index = line[n]
    vel = sqrt(velocity[index][0]*velocity[index][0]+velocity[index][1]*velocity[index][1])
    file.write('%i\t%i\t%e\n'   % (n, index, vel))
file.close()
"""
