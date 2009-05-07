#!/usr/bin/env python

import sys
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt
import importers
import fvmparallel
import time
from numpy import *
from mpi4py  import MPI
atype = 'double'
#atype = 'tangent'

if atype == 'double':
    import models_atyped_double as models
    import exporters_atyped_double as exporters
elif atype == 'tangent':
    import models_atyped_tangent_double as models
    import exporters_atyped_tangent_double as exporters


from FluentCase import FluentCase

#fvmbaseExt.enableDebug("cdtor")

fileBase = None
numIterations = 10
#fileBase = "/home/yildirim/memosa/src/fvm/test/cav_44_tri"
fileBase = "/home/yildirim/memosa/src/fvm/test/tri_894"
#fileBase = "/home/yildirim/memosa/src/fvm/test/test_tri_100by100"
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


def dumpTecplotFile(nmesh, meshes):
  #cell sites
  cellSites = []
  for n in range(0,nmesh):
     cellSites.append( meshes[n].getCells() )
     print "cellSites[", n, "].getCount = ", cellSites[n].getCount()

  #face sites
  faceSites = []
  for n in range(0,nmesh):
     faceSites.append( meshes[n].getFaces() )

  #node sites
  nodeSites = []
  for n in range(0,nmesh):
     nodeSites.append( meshes[n].getNodes() )

  #get connectivity (faceCells)
  faceCells = []
  for n in range(0,nmesh):
     faceCells.append( meshes[n].getConnectivity( faceSites[n], cellSites[n] ) )
 
  #get connectivity ( cellNodes )
  cellNodes = []
  for n in range(0,nmesh):
     cellNodes.append( meshes[n].getCellNodes() )

  #get Volume as array
  volumes = []
  for n in range(0,nmesh):
     volumes.append( geomFields.volume[cellSites[n]].asNumPyArray() )
 
  cellCentroids =[]
  for n in range(0,nmesh):
     cellCentroids.append( geomFields.coordinate[cellSites[n]].asNumPyArray() )

  velFields = []
  for n in range(0,nmesh):
     velFields.append( thermalFields.temperature[cellSites[n]].asNumPyArray() )

  coords = []
  for n in range(0,nmesh):
     coords.append( meshes[n].getNodeCoordinates().asNumPyArray() )
     print "shape( coords[", n, "] ) = ", shape( coords[n] )	
     
  file_name = "temp_proc" + str(MPI.COMM_WORLD.Get_rank()) + ".dat"
  f = open(file_name, 'w')

  f.write("Title = \" tecplot file for 2D Cavity problem \" \n")
  f.write("variables = \"x\", \"y\", \"z\", \"velX\", \"cellCentroidY\" \n")
  for n in range(0,nmesh):
     title_name = "nmesh" + str(n)
     ncell  = cellSites[n].getSelfCount()
     nnode  = nodeSites[n].getCount()
     zone_name = "Zone T = " + "\"" + title_name +  "\"" +      \
                 " N = " + str( nodeSites[n].getCount() ) +     \
		 " E = " + str( ncell ) +  \
		 " DATAPACKING = BLOCK, VARLOCATION = ([4-5]=CELLCENTERED), " + \
		 " ZONETYPE=FETRIANGLE \n"
     f.write( zone_name )
     #write x
     for i in range(0,nnode):
          f.write(str(coords[n][i][0])+"    ")
	  if ( i % 5 == 4 ):
	     f.write("\n")
     f.write("\n")	  
     
     #write y
     for i in range(0,nnode):
          f.write(str(coords[n][i][1])+"    ")
	  if ( i % 5 == 4 ):
	     f.write("\n")
     f.write("\n")	  

     #write z
     for i in range(0,nnode):
          f.write(str(coords[n][i][2])+"    ")
	  if ( i % 5 == 4 ):
	     f.write("\n")
     f.write("\n")	  
     
     
     #write velX
     for i in range(0,ncell):
        f.write( str(velFields[n][i]) + "    ")
	if ( i % 5  == 4 ):
	    f.write("\n")
     f.write("\n")

     #write velX
     for i in range(0,ncell):
        f.write( str(cellCentroids[n][i][1]) + "    ")
	if ( i % 5  == 4 ):
	    f.write("\n")
     f.write("\n")
     	  	    
   
     #connectivity
     for i in range(0,ncell):
        nnodes_per_cell = cellNodes[n].getCount(i)
        for node in range(0,nnodes_per_cell):
	    f.write( str(cellNodes[n](i,node)+1) + "     ")
        f.write("\n")

     	    
	    
     f.write("\n")	  
     
	
  
  f.close()
# change as needed
#import debug

reader = FluentCase(fileBase+".cas")
reader.read()
fluent_meshes = reader.getMeshList()

nmesh = 1

import time
t0 = time.time()



#print "nmesh = ", nmesh
#npart = fvmparallel.IntVector(1,nmesh)  #total of distributed meshes
#etype = fvmparallel.IntVector(1,1) #triangle

npart = [MPI.COMM_WORLD.Get_size()]
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

geomFields =  models.GeomFields('geom')

metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

if atype == 'tangent':
    metricsCalculator.setTangentCoords(0,7,1)

thermalFields =  models.ThermalFields('therm')

tmodel = models.ThermalModelA(geomFields,thermalFields,meshes)


## set bc for top to be at 400
bcMap = tmodel.getBCMap()

if 3 in bcMap:
   bc3 = tmodel.getBCMap()[3]
   bc3.bcType = 'SpecifiedTemperature'
   bc3.setVar('specifiedTemperature',400)
if 4 in bcMap:
   bc4 = tmodel.getBCMap()[4]
   bc4.bcType = 'SpecifiedTemperature'
   bc4.setVar('specifiedTemperature',0)

## set viscosity and density, this is done per mesh since each mesh has its own VC object
#vcMap = tmodel.getVCMap()
#for vc in vcMap.values():
#    vc.setVar('density',1.0)
#    vc.setVar('thermalConductivity',1.0)


tSolver = fvmbaseExt.AMG()
tSolver.relativeTolerance = 1e-9
tSolver.nMaxIterations = 2000
tSolver.maxCoarseLevels=20
tSolver.verbosity=2

toptions = tmodel.getOptions()

toptions.linearSolver = tSolver

#import debug
tmodel.init()
print "nmesh = ", nmesh,  "procID = ", MPI.COMM_WORLD.Get_rank() 
tmodel.advance(1)

dumpTecplotFile( nmesh, meshes)
