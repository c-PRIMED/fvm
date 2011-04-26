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
import ServerCouplingPlate
from mpi4py import MPI

etype = {
        'tri'  : 1,
        'quad' : 2,
        'tetra': 3,
        'hexa' : 4
       }



def writeProbeData():
    deformation = plateFields.deformation[solidMeshes[0].getCells()].asNumPyArray()
    maxDef = deformation.min(axis = 0)
    probeFile.write('%e\t%e\t%e\n' % (globalTime, deformation[probeIndex][2], maxDef[2]))
    probeFile.flush()

def writeForceData():
    forceFile.write('%e\t%e\t%e\n'   % (globalTime, eForce, fForce))
    forceFile.flush()

def saveVTK(n):

    writer1 = exporters.VTKWriterA(geomFields,solidMeshes,
                                  fileBase_output + "deformation-" + str(n) + ".vtk",
                                  "frogleg",
                                  False,0)
    writer1.init()
    writer1.writeVectorField(plateFields.deformation,"deformation")
    writer1.finish()
### ========================== properties and parameters ===============================###

parser = optparse.OptionParser()
parser.add_option("--volt", type=float)
parser.add_option("--type", help="'tri'[default], 'quad', 'hexa', or 'tetra'")    
(options, args) = parser.parse_args()
applied_voltage = -options.volt
print "voltage =",applied_voltage



### beam
rho = 8800                    # density kg/m^3
E = 200e9                     # Young's modulus 
nu = 0.31                      # Poisson's ratio

beam_thickness = 3e-6

beamLeft = 3
beamRight = 4
beamSide = 5


numTimeSteps = 2 
globalTime = 0
globalCount = 0
timeStep = 2e-8
saveFrequency = 1
initialTransient = False
probeIndex = 0

### ===================== mesh read ===============================================###
fileBase_input  = "/home/yildirim/memosa/src/fvm/scripts/cantilever3D_coupling/"
fileBase_output = "./" + str(int(-applied_voltage)) + "/"

### 2D plate mesh
beamReader = FluentCase(fileBase_input+"beam_2D.cas")
beamReader.read()

##paritioning
#nmesh = 1
#npart = [MPI.COMM_WORLD.Get_size()]
#print "options folud.type = ", options.type
#etype = [etype[options.type]]
##partMesh constructor and setTypes
#part_mesh = fvmparallel.MeshPartitioner( fluentMeshes, npart, etype );
#part_mesh.setWeightType(0);
#part_mesh.setNumFlag(0);
##actions
#part_mesh.isDebug(0)
#part_mesh.partition()
#part_mesh.mesh()
#solidMeshes  = part_mesh.meshList()

solidMeshes = beamReader.getMeshList()
geomFields =  models.GeomFields('geom')
solidMetricsCalculator = models.MeshMetricsCalculatorA(geomFields,solidMeshes)
solidMetricsCalculator.init()

#extruding
solidBoundaryMeshes = [m.extrude(1, beam_thickness, True) for m in solidMeshes]
solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(geomFields,solidBoundaryMeshes)
solidBoundaryMetricsCalculator.init()

### find device center
cells = solidMeshes[0].getCells()
xc = geomFields.coordinate[cells].asNumPyArray()
small = 100
probeIndex = 0
for c in range(0, cells.getCount()):
    rsqr = (xc[c][0]-150e-6)*(xc[c][0]-150e-6) + xc[c][1]*xc[c][1] +xc[c][2]*xc[c][2] 
    if rsqr < small:
        small = rsqr
        probeIndex = c

print "probe Index is %i " % probeIndex

### output files 
print fileBase_output + "centerDisplacement.dat", "w"
probeFile = open(fileBase_output + "centerDisplacement.dat", "w")
forceFile = open(fileBase_output + "force.dat", "w")

### =============================== models =====================================###

### Plate Model and boundary conditions ###
plateFields =  models.PlateFields('plate')
pmodel = models.PlateModelA(geomFields, plateFields, solidMeshes)
dmodel = models.PlateDeformationModelA(geomFields, plateFields, solidMeshes)

bcMap = pmodel.getBCMap()

for id in [beamLeft]:
    bc = bcMap[id]
    bc.bcType = 'Clamped'
    bc['specifiedXRotation']=0
    bc['specifiedYRotation']=0.
    bc['specifiedZDeformation']=0. 

for id in [beamSide, beamRight]:
    bc = bcMap[id]
    bc.bcType = 'SpecifiedTraction'

vcMap = pmodel.getVCMap()
for i,vc in vcMap.iteritems():
    vc['density'] = rho
    vc['ym'] = E
    vc['nu'] = nu
### ================================= solvers ===================================###
pmodel.printBCs()
### plate solver ###
pc = fvmbaseExt.AMG()
pc.maxCoarseLevels = 0
pc.verbosity=3
defSolver = fvmbaseExt.BCGStab()
defSolver.preconditioner = pc
defSolver.relativeTolerance = 1e-9
defSolver.absoluteTolerance = 1.e-30
defSolver.nMaxIterations = 5000
defSolver.verbosity=0

poptions = pmodel.getOptions()
poptions.deformationLinearSolver = defSolver
poptions.deformationTolerance=1.0e-3
poptions.setVar("deformationURF",1.0)
poptions.printNormalizedResiduals=True
poptions.timeDiscretizationOrder = 2
poptions.transient=True
poptions.scf = 5./6.
poptions.setVar('timeStep',timeStep)



### initialize models and run ###

pmodel.init()
dmodel.init()

solid_to_fluid = ServerCouplingPlate.ServerCoupling(pmodel, dmodel, geomFields, solidMeshes, \
                    solidBoundaryMeshes, plateFields, beam_thickness)


t1 = time.time()
  
#--------------Timestep Loop --------------------------#
for n in range(0, numTimeSteps):                
   
    #------------solve structure-------------#
    print "***  solving structure model at globalCount %i   ***" % globalCount
    for i in range (0, 3):    
        pmodel.advance(1)
    dmodel.calculateNodeDisplacement()
    dmodel.deformPlate()
    solidMetricsCalculator.recalculate_deform()

    solidBoundaryMetricsCalculator.recalculate_deform() 
   
    #coupling
    print "before  solid_to_fluid.update"
    solid_to_fluid.update() #first send coords and vels  to fluid side
    solid_to_fluid.accept() #then wait for fluid side to send you strees to continue
    
    #update time    
    pmodel.updateTime()
    dmodel.updateTime()
    globalTime += timeStep
    globalCount += 1

    #------data output-----------------------#
    writeProbeData()
    if (n%saveFrequency == 0):
        saveVTK(n)
 
t2 = time.time()



probeFile.close()
forceFile.close()
print  '\nsolution time = %f' % (t2-t1)
