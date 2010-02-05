#!/usr/bin/env python

import sys

import fvm
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.fvmparallel as fvmparallel
import time
from numpy import *
import tecplot
from mpi4py import MPI

fvm.set_atype('double')

from FluentCase import FluentCase

#fvmbaseExt.enableDebug("cdtor")

fileBase = None
numIterations = 10
#fileBase = "/home/yildirim/memosa/src/fvm/test/cav_44_tri"
fileBase = "/home/yildirim/memosa/src/fvm/test/tri_894"
#fileBase = "/home/sm/a/data/wj"

def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat if it is not specified."
    sys.exit(1)

def advance(fmodel,niter):
    for i in range(0,niter):
        try:
            fmodel.advance(1)
        except KeyboardInterrupt:
            break
	    
      

# change as needed


outfile = None
if __name__ == '__main__' and fileBase is None:
    if len(sys.argv) < 2:
        usage()
    fileBase = sys.argv[1]
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

if outfile == None:
    outfile = fileBase+"-prism.dat"

#import debug

   
reader = FluentCase(fileBase+".cas")

reader.read();
import sys

fluent_meshes = reader.getMeshList()

import time
t0 = time.time()
nmesh = MPI.COMM_WORLD.Get_size()

#print "nmesh = ", nmesh
#npart = fvmparallel.IntVector(1,nmesh)  #total of distributed meshes
#etype = fvmparallel.IntVector(1,1) #triangle

npart = [nmesh]
etype = [1]

#partMesh constructor and setTypes
part_mesh = fvmparallel.PartMesh( fluent_meshes, npart, etype );
part_mesh.setWeightType(0);
part_mesh.setNumFlag(0);

#actions
part_mesh.partition()
part_mesh.mesh()
part_mesh.mesh_debug()
meshes = part_mesh.meshList()

geomFields =  fvm.models.GeomFields('geom')
metricsCalculator = fvm.models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()
#!!!!!!!!!!!!!!!!source code (metricsMeshCalcultaro)::init(), volumeField.synLocal, coordField.syncLocal() turned off

flowFields =  fvm.models.FlowFields('flow')



fmodel = fvm.models.FlowModelA(geomFields,flowFields,meshes)

## set bc for top to be a wall with x velocity
bcMap = fmodel.getBCMap()
if 3 in bcMap:
   bc3 = fmodel.getBCMap()[3]
   bc3.bcType = 'NoSlipWall'
   bc3.setVar('specifiedXVelocity',1)
				

## set viscosity and density, this is done per mesh since each mesh has its own VC object
vcMap = fmodel.getVCMap()
for vc in vcMap.values():
    vc.setVar('density',1.0)
    vc.setVar('viscosity',1.0)


momSolver = fvmbaseExt.AMG()
momSolver.relativeTolerance = 1e-1
momSolver.nMaxIterations = 20
momSolver.maxCoarseLevels=20
momSolver.verbosity=2

contSolver = fvmbaseExt.AMG()
contSolver.relativeTolerance = 1e-1
contSolver.nMaxIterations = 20
contSolver.verbosity=2
contSolver.maxCoarseLevels=20

foptions = fmodel.getOptions()

foptions.momentumLinearSolver = momSolver
foptions.pressureLinearSolver = contSolver

foptions.momentumTolerance=1e-5
foptions.continuityTolerance=1e-5
foptions.setVar("momentumURF",0.7)
foptions.setVar("pressureURF",0.3)
foptions.printNormalizedResiduals=False

fmodel.init()

#advance(fmodel,numIterations)
cellSite = meshes[0].getCells()
velField = flowFields.velocity[ cellSite ].asNumPyArray() 
#initializing values (taking velocity as their cell Number)

#print velField
#print "mesh_id = ", MPI.COMM_WORLD.Get_rank()
#if MPI.COMM_WORLD.Get_rank() == 1:
#    print  flowFields.velocity[ cellSite ].asNumPyArray()



#import debug
fmodel.advance(100)
tecplot.dumpTecplotFile(1, meshes, flowFields)







