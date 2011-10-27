#!/usr/bin/env python
"""
Usage: testThermalParallel.py [options] infile
options are:
--type  'tri'[default], 'quad', 'hexa', or 'tetra'
--xdmf  Dump data in xdmf
"""

import sys
import fvm
import fvm.fvmbaseExt as fvmbaseExt
fvm.set_atype('double')
import fvm.importers as importers
import fvm.fvmparallel as fvmparallel
import time
from numpy import *
from mpi4py  import MPI
from optparse import OptionParser

from FluentCase import FluentCase

def usage():
    print __doc__
    sys.exit(1)

# map between fvm, tecplot, and xdmf types
etype = {
        'tri' : 1,
        'quad' : 2,
        'tetra' : 3,
        'hexa' : 4
        }
tectype = {
        'tri' : 'FETRIANGLE',
        'quad' : 'FEQUADRILATERAL',
        'tetra' : 'FETETRAHEDRON',
        'hexa' : 'FEBRICK'
        }
xtype = {
        'tri' : 'Triangle',
        'quad' : 'Quadrilateral',
        'tetra' : 'Tetrahedron',
        'hexa' : 'Hexahedron'
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

def dumpTecplotFile(nmesh, meshes, mtype):
  #cell sites
  cellSites = []
  for n in range(0,nmesh):
     cellSites.append( meshes[n].getCells() )
#     print "cellSites[", n, "].getCount = ", cellSites[n].getCount()

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
#     print "shape( coords[", n, "] ) = ", shape( coords[n] )	
     
  f = open("temp_proc%s.dat" % MPI.COMM_WORLD.Get_rank(), 'w')
  f.write("Title = \" tecplot file for 2D Cavity problem \" \n")
  f.write("variables = \"x\", \"y\", \"z\", \"velX\", \"cellCentroidY\" \n")
  for n in range(0,nmesh):
     title_name = "nmesh%s" % n
     ncell  = cellSites[n].getSelfCount()
     nnode  = nodeSites[n].getCount()
     f.write("Zone T = \"%s\" N = %s E = %s DATAPACKING = BLOCK, VARLOCATION = ([4-5]=CELLCENTERED), ZONETYPE=%s\n" %
             (title_name,  nodeSites[n].getCount(), ncell, tectype[mtype]))
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

def writeXdmfHeader():
    nprocs = MPI.COMM_WORLD.Get_size()
    mesh_file = open("mesh.xmf", "w")
    mesh_file.write("<?xml version='1.0' ?>\n")
    mesh_file.write("<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n")
    mesh_file.write("<Xdmf xmlns:xi='http://www.w3.org/2001/XInclude' Version='2.0'>\n")
    mesh_file.write("  <Domain>\n")
    for i in range(0,nprocs):
        mesh_file.write("    <xi:include href='mesh_proc%s.xmf' />\n" % i)
    mesh_file.write("  </Domain>\n")
    mesh_file.write("</Xdmf>\n")
    mesh_file.close()
    
def dumpXdmfFile(nmesh, meshes, mtype):
    proc = MPI.COMM_WORLD.Get_rank()
    if proc == 0:
        writeXdmfHeader()

    # cell sites
    cellSites = []
    for n in range(0,nmesh):
        cellSites.append( meshes[n].getCells() )

    # face sites
    faceSites = []
    for n in range(0,nmesh):
        faceSites.append( meshes[n].getFaces() )

    # node sites
    nodeSites = []
    for n in range(0,nmesh):
        nodeSites.append( meshes[n].getNodes() )

    # get connectivity (faceCells)
    faceCells = []
    for n in range(0,nmesh):
        faceCells.append( meshes[n].getConnectivity( faceSites[n], cellSites[n] ) )
 
    # get connectivity ( cellNodes )
    cellNodes = []
    for n in range(0,nmesh):
        cellNodes.append( meshes[n].getCellNodes() )

    # get Volume as array
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
#        print "shape( coords[", n, "] ) = ", shape( coords[n] )	
     
    f = open("mesh_proc%s.xmf" % proc, 'w')
    for n in range(0,nmesh):
        f.write("<Grid Name='Mesh%s-%s' GridType='Uniform'>\n" % (proc, n))
        ncell  = cellSites[n].getSelfCount()
        nnode  = nodeSites[n].getCount()
        f.write("  <Topology TopologyType='%s' Dimensions='%s'>\n" % (xtype[mtype], ncell))
        f.write("    <DataItem Dimensions='%s %s'>\n" % (ncell, cellNodes[0].getCount(0)))

        # connectivity (topology)
        for i in range(0,ncell):
            f.write("      ")
            nnodes_per_cell = cellNodes[n].getCount(i)
            for node in range(0, nnodes_per_cell):
                f.write("%s " % cellNodes[n](i,node))
            f.write("\n")
        f.write("    </DataItem>\n")
        f.write("  </Topology>\n")

        # Geometry
        f.write("  <Geometry Type='XYZ'>\n")
        f.write("    <DataItem Dimensions='%s 3' NumberType='Float'>\n" % nnode)
        for i in range(0,nnode):
            f.write("      %s %s %s\n" % (coords[n][i][0], coords[n][i][1], coords[n][i][2]))
        f.write("    </DataItem>\n")     
        f.write("  </Geometry>\n")

        # ATTRIBUTES

        # velX
        f.write("  <Attribute Name='velX' Center='Cell'>\n")
        f.write("    <DataItem Dimensions='%s'>" % ncell)
        for i in range(0,ncell):
            if ( i % 5  == 0 ):
                f.write("\n      ")
            f.write("%s " % velFields[n][i])
        f.write("\n")
        f.write("    </DataItem>\n")     
        f.write("  </Attribute>\n")

        # cellCentroids
        f.write("  <Attribute Name='cellCentroids' Center='Cell'>\n")
        f.write("    <DataItem Dimensions='%s'>" % ncell)
        for i in range(0,ncell):
            if ( i % 5  == 0 ):
                f.write("\n      ")
            f.write("%s " % cellCentroids[n][i][1])
        f.write("\n")
        f.write("    </DataItem>\n")     
        f.write("  </Attribute>\n")

        f.write("</Grid>\n")
    f.close()


parser = OptionParser()
parser.set_defaults(type='tri')
parser.add_option("--type", help="'tri'[default], 'quad', 'hexa', or 'tetra'")
parser.add_option("--xdmf", action='store_true', help="Dump data in xdmf")
parser.add_option("--time","-t",action='store_true',help="Print timing information.")
(options, args) = parser.parse_args()
if len(args) != 1:
    usage()

numIterations = 10
reader = FluentCase(args[0])
reader.read()
fluent_meshes = reader.getMeshList()

nmesh = 1
npart = [MPI.COMM_WORLD.Get_size()]
etype = [etype[options.type]]

if not MPI.COMM_WORLD.Get_rank():
   print "parmesh is processing"
     
    

if options.time:
    # time profile for partmesh
    part_mesh_time = zeros(1,dtype='d')
    part_mesh_start = zeros(1, dtype='d')
    part_mesh_end   = zeros(1, dtype='d')
    part_mesh_maxtime = zeros(1,dtype='d')
    part_mesh_mintime = zeros(1, dtype='d')
    part_mesh_start[0] = MPI.Wtime()
#partMesh constructor and setTypes
part_mesh = fvmparallel.MeshPartitioner( fluent_meshes, npart, etype );
part_mesh.setWeightType(0);
part_mesh.setNumFlag(0);   
#actions
part_mesh.isCleanup(0)
#part_mesh.fiedler_order("ipermutation.txt")
part_mesh.isDebug(0)
part_mesh.partition()
part_mesh.mesh()
meshes = part_mesh.meshList()
reader = 0
fluent_meshes = 0

if options.time:
    part_mesh_end[0] = MPI.Wtime()
    part_mesh_time[0] = part_mesh_end[0] - part_mesh_start[0]
    MPI.COMM_WORLD.Allreduce( [part_mesh_time,MPI.DOUBLE], [part_mesh_maxtime, MPI.DOUBLE], op=MPI.MAX) 
    MPI.COMM_WORLD.Allreduce( [part_mesh_time,MPI.DOUBLE], [part_mesh_mintime, MPI.DOUBLE], op=MPI.MIN) 
    solver_start   = zeros(1, dtype='d')
    solver_end     = zeros(1, dtype='d')
    solver_time    = zeros(1, dtype='d')
    solver_maxtime = zeros(1,dtype='d')
    solver_mintime = zeros(1,dtype='d')
    solver_start[0] = MPI.Wtime()

geomFields =  fvm.models.GeomFields('geom')

metricsCalculator = fvm.models.MeshMetricsCalculatorA(geomFields,meshes)
metricsCalculator.init()

flowFields = fvm.models.FlowFields('flow')
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
    vc.setVar('viscosity',0.1)
#########################
conn =meshes[0].getCellCells() 
semi_bandwidth_mom = 3 
#spike_storage_mom = fvmbaseExt.SpikeStorage(conn, semi_bandwidth_mom) 
#pcMom = fvmbaseExt.SpikeSolver(spike_storage_mom)
#pcMom = fvmbaseExt.JacobiSolver()
#pcMom = fvmbaseExt.ILU0Solver()
pcMom = fvmbaseExt.AMG()
pcMom.verbosity = 0

momSolver = fvmbaseExt.JacobiSolver()
#momSolver.smootherType = fvmbaseExt.AMG.JACOBI
#momSolver.preconditioner=pcMom
momSolver.relativeTolerance = 1e-1
momSolver.nMaxIterations = 200
momSolver.maxCoarseLevels=0
momSolver.verbosity=0
if MPI.COMM_WORLD.Get_rank() == 0:
   momSolver.verbosity = 0

#########################
semi_bandwidth_con = 3 
#spike_storage_con = fvmbaseExt.SpikeStorage(conn, semi_bandwidth_con) 
#pcCon = fvmbaseExt.SpikeSolver(spike_storage_con)
#pcCon = fvmbaseExt.JacobiSolver()
#pcCon = fvmbaseExt.ILU0Solver()
pcCon = fvmbaseExt.AMG()
pcCon.verbosity = 0
contSolver = fvmbaseExt.JacobiSolver()
#contSolver.smootherType = fvmbaseExt.AMG.JACOBI
#contSolver.preconditioner =pcCon
contSolver.relativeTolerance = 1e-1
contSolver.nMaxIterations = 200
contSolver.verbosity=0
if MPI.COMM_WORLD.Get_rank() == 0:
   contSolver.verbosity = 0
contSolver.maxCoarseLevels=0

foptions = fmodel.getOptions()

foptions.momentumLinearSolver = momSolver
foptions.pressureLinearSolver = contSolver

foptions.momentumTolerance=1e-5
foptions.continuityTolerance=1e-5
foptions.printNormalizedResiduals=False

fmodel.init()
                                      
#fmodel.dumpContinuityMatrix("cav4096K")
#advance(fmodel,numIterations)
cellSite = meshes[0].getCells()
velField = flowFields.velocity[ cellSite ].asNumPyArray() 
#initializing values (taking velocity as their cell Number)

#print velField
#print "mesh_id = ", MPI.COMM_WORLD.Get_rank()
#if MPI.COMM_WORLD.Get_rank() == 1:
#    print  flowFields.velocity[ cellSite ].asNumPyArray()



#import debug
pcCon.redirectPrintToFile("convergence.dat")
fmodel.advance(numIterations)
pcCon.redirectPrintToScreen()
#fmodel.advanceCoupled(numIterations)


if options.time:
    solver_end[0] = MPI.Wtime()
    solver_time[0] = solver_end[0] - solver_start[0]
    MPI.COMM_WORLD.Allreduce( [solver_time,MPI.DOUBLE], [solver_maxtime, MPI.DOUBLE], op=MPI.MAX) 
    MPI.COMM_WORLD.Allreduce( [solver_time,MPI.DOUBLE], [solver_mintime, MPI.DOUBLE], op=MPI.MIN) 
    if MPI.COMM_WORLD.Get_rank() == 0:
        dumpMPITimeProfile(part_mesh_maxtime, part_mesh_mintime, solver_maxtime, solver_mintime)

#dumpTecplotFile( nmesh, meshes, options.type)
#if options.xdmf:
#    dumpXdmfFile( nmesh, meshes, options.type)
