#!/usr/bin/env python

import fvm
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.fvmparallel as fvmparallel
import sys, time
from numpy import *
from mpi4py  import MPI
from FluentCase import FluentCase

#fvmbaseExt.enableDebug("cdtor")

fileBase = None
numIterations = 10
fileBase = "/home/yildirim/memosa/src/fvm/test/cav_26_tri"
#fileBase = "/home/yildirim/memosa/src/fvm/test/cav_44_tri"
#fileBase = "/home/yildirim/memosa/src/fvm/test/tri_894"
#fileBase = "/home/yildirim/memosa/src/fvm/test/cav_tri_915K"
#fileBase = "/home/yildirim/memosa/src/fvm/test/cav_tri_3_66M"
#fileBase = "/home/yildirim/memosa/src/fvm/test/test_tri"

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


def dumpMPITimeProfile(part_mesh_maxtime, part_mesh_mintime, solver_maxtime, solver_mintime):

      fname = "time_mpi_totalprocs"  + str(MPI.COMM_WORLD.Get_size()) + ".dat"
      f = open(fname,'w')
      line = " part_mesh_mintime = " + str(part_mesh_mintime[0]) + "\n" + \
             " part_mesh_maxtime = " + str(part_mesh_maxtime[0]) + "\n" + \
             " solver_mintime    = " + str(solver_mintime[0])    + "\n" + \
             " solver_maxtime    = " + str(solver_maxtime[0])    + "\n"
      print line
      f.write(line)
      f.close()


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

reader = FluentCase(sys.argv[1])
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

#time profile for partmesh
part_mesh_time = zeros(1,dtype='d')
part_mesh_start = zeros(1, dtype='d')
part_mesh_end   = zeros(1, dtype='d')
part_mesh_maxtime = zeros(1,dtype='d')
part_mesh_mintime = zeros(1, dtype='d')

part_mesh_start[0] = MPI.Wtime()

part_mesh = fvmparallel.PartMesh( fluent_meshes, npart, etype );
part_mesh.setWeightType(0);
part_mesh.setNumFlag(0);

#print "nmesh = ", nmesh,  "procID = ", MPI.COMM_WORLD.Get_rank() 
#actions
part_mesh.partition()
part_mesh.mesh()
#part_mesh.mesh_debug()
meshes = part_mesh.meshList()

part_mesh_end[0] = MPI.Wtime()
part_mesh_time[0] = part_mesh_end[0] - part_mesh_start[0]
MPI.COMM_WORLD.Allreduce( [part_mesh_time,MPI.DOUBLE], [part_mesh_maxtime, MPI.DOUBLE], op=MPI.MAX) 
MPI.COMM_WORLD.Allreduce( [part_mesh_time,MPI.DOUBLE], [part_mesh_mintime, MPI.DOUBLE], op=MPI.MIN) 


target_id = int(0)
group     = fvmbaseExt.IntSet()

for i in range(0,npart[0]):
   group.insert( int(i) )	

mesh0 = meshes[0]
cellSite = mesh0.getCells()
cellSiteMerger = fvmbaseExt.StorageSiteMerger( 0, group, cellSite )
cellSiteMerger.debug_print()



