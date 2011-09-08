#!/usr/bin/env python
### import  modules ###
"""
mpirun -np   1 python ./parallel_sphere.py ./cube-125k.cas ./sphere.msh --type=hexa
mpirun -np 64  python parallel_sphere.py ./fluid_1000K.cas ./sphere.msh --type=hexa

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
    
def dumpL2Error(fluidMeshes,geomFields, flowFields):
   cellSite = fluidMeshes[0].getCells()
   velocitya = flowFields.velocity[cellSite]
   velocity  = velocitya.asNumPyArray()
   
   #analytical solution
   nCells =cellSite.getSelfCount()
   xc = geomFields.coordinate[cellSite].asNumPyArray() 
   cellType = geomFields.ibType[cellSite].asNumPyArray()

   #analytical solutions
   xvelA = zeros(nCells,float)
   yvelA = zeros(nCells,float)
   zvelA = zeros(nCells,float)   
   velMagA = zeros(nCells,float)
   velMag  = zeros(nCells,float)
   
   a = 10.0
   U0 = 0.001
   for i in range(0,nCells):
     if (cellType[i] == -1): #only for fluid not for solid and ibcells
       x = xc[i][0]
       y = xc[i][1]
       z = xc[i][2]
       r = sqrt(x*x+y*y+z*z)
       alfa = acos(z/r)
       beta = atan2(y,x)

       Ur = U0 * cos(alfa) * (1-1.5*a/r+0.5*pow(a, 3)/pow(r,3))
       Ualfa = -U0 * sin(alfa) * (1-0.75*a/r-0.25*pow(a, 3)/pow(r,3))

       xvelA[i] = Ur * sin(alfa) * cos(beta) + Ualfa * cos(alfa) * cos(beta)
       yvelA[i] = Ur * sin(alfa) * sin(beta) + Ualfa * cos(alfa) * sin(beta)
       zvelA[i] = Ur * cos(alfa) - Ualfa * sin(alfa)	 
  
       #compute maginutes
       velMagA[i] = sqrt(xvelA[i]**2 + yvelA[i]**2 + zvelA[i]**2)
       velMag[i]  = sqrt(velocity[i][0]**2 + velocity[i][1]**2 + velocity[i][2]**2)
     
   diffV = velMag - velMagA
   
   term1 = dot(diffV,diffV)
   term2 = dot(velMagA,velMagA)
   L2error = sqrt( term1.sum() / term2.sum() )
   
   f = open('L2Error.dat','w')
   f.write("%20.14f \n"%(L2error))
   
   f.close()
   
        
    
def createBVFields(geomFields,meshes):
    fx = fvmbaseExt.Field('bvx')
    fy = fvmbaseExt.Field('bvy')
    fz = fvmbaseExt.Field('bvz')

    mesh = meshes[0]
    vol = geomFields.volume[mesh.getCells()]
    
    for mesh in meshes:
        fgs = mesh.getBoundaryGroups()
        for fg in fgs:
            nFaces = fg.site.getCount()

            xvel = vol.newSizedClone(nFaces)
            yvel = vol.newSizedClone(nFaces)
            zvel = vol.newSizedClone(nFaces)

            xvela = xvel.asNumPyArray()
            yvela = yvel.asNumPyArray()
            zvela = zvel.asNumPyArray()
            
            xf = geomFields.coordinate[fg.site].asNumPyArray()

            
            a = 10.0
            U0 = 0.001
          
            
            for i in range(0,nFaces):
                x = xf[i][0]
                y = xf[i][1]
                z = xf[i][2]
                
                r = sqrt(x*x+y*y+z*z)
                alfa = acos(z/r)
                beta = atan2(y,x)

                Ur = U0 * cos(alfa) * (1-1.5*a/r+0.5*pow(a, 3)/pow(r,3))
                Ualfa = -U0 * sin(alfa) * (1-0.75*a/r-0.25*pow(a, 3)/pow(r,3))

                xvela[i] = Ur * sin(alfa) * cos(beta) + Ualfa * cos(alfa) * cos(beta)
                yvela[i] = Ur * sin(alfa) * sin(beta) + Ualfa * cos(alfa) * sin(beta)
                zvela[i] = Ur * cos(alfa) - Ualfa * sin(alfa)

            #pdb.set_trace()       

              
            fx[fg.site] = xvel
            fy[fg.site] = yvel
            fz[fg.site] = zvel
            
    return fx,fy,fz

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
#nodeCoord = solidMeshes[0].getNodeCoordinates().asNumPyArray()
#nodeCoord[:,0] += 7.5

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

fx,fy,fz = createBVFields(geomFields,fluidMeshes)


bcMap = fmodel.getBCMap()

for bc in bcMap.values():
       bc['specifiedXVelocity']=fx
       bc['specifiedYVelocity']=fy
       bc['specifiedZVelocity']=fz
       bc.bcType = 'VelocityBoundary'
#       if bcId == 9:
#          bc.bcType='PressureBoundary'
	  



#flow solvers
momPC = fvmbaseExt.AMG()
momPC.verbosity=0
momSolver = fvmbaseExt.AMG()
#momSolver.preconditioner = momPC 
momSolver.relativeTolerance = 1e-1
momSolver.absoluteTolerance = 1e-50
momSolver.nMaxIterations = 20
momSolver.verbosity=0


contPC = fvmbaseExt.AMG()
contPC.verbosity=0
contSolver = fvmbaseExt.AMG()
#contSolver.preconditioner = contPC 
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
foptions.transient=False
foptions.setVar("timeStep",timeStep)
foptions.printNormalizedResiduals=False
#flowVCMap = fmodel.getVCMap()

#for id, vc in flowVCMap.iteritems():
#    vc['density'] = 1.225
#    vc['viscosity'] = 1.8e-5

fmodel.init()

sbMeshFaces = solidBoundaryMeshes[0].getFaces()    
ibManager = fvmbaseExt.IBManager(geomFields,solidBoundaryMeshes[0],fluidMeshes)
faceCount = sbMeshFaces.getCount()
area = geomFields.area[sbMeshFaces]
velocity = area.newSizedClone(faceCount)
velocitya = velocity.asNumPyArray()
velocitya[:,:] = 0.01
flowFields.velocity[sbMeshFaces] = velocity

ibManager.fluidNeighborsPerIBFace = 3
ibManager.solidNeighborsPerIBFace = 2
ibManager.fluidNeighborsPerSolidFace = 3

globalCount = 0
globalTime = 0.0

ibManager.update()
fluidMetricsCalculator.computeIBInterpolationMatrices(sbMeshFaces)
#fluidMetricsCalculator.computeSolidInterpolationMatrices(sbMeshFaces)

numTimeSteps=1
for n in range(0, numTimeSteps):
    ibManager.update()
    fluidMetricsCalculator.computeIBInterpolationMatrices(sbMeshFaces)    
    momPC.redirectPrintToFile("convergence.dat")
    for i in range(0,10):
        fmodel.computeIBFaceVelocity(sbMeshFaces)
        if fmodel.advance(1):
            break    
    momPC.redirectPrintToScreen()
    #fmodel.updateTime()
    globalCount += 1
    globalTime += timeStep

if options.time:
    solver_end[0]  = MPI.Wtime()
    solver_time[0] = solver_end[0] - solver_start[0]
    MPI.COMM_WORLD.Allreduce( [solver_time,MPI.DOUBLE], [solver_maxtime, MPI.DOUBLE], op=MPI.MAX) 
    MPI.COMM_WORLD.Allreduce( [solver_time,MPI.DOUBLE], [solver_mintime, MPI.DOUBLE], op=MPI.MIN) 
    if MPI.COMM_WORLD.Get_rank() == 0:
        dumpMPITimeProfile(solver_maxtime, solver_mintime)    

#dumpTecplotFile(nmesh, fluidMeshes, fluidMeshes0, options.type,flowFields,"sphere3D.dat") 
#dumpL2Error(fluidMeshes,geomFields, flowFields)

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
