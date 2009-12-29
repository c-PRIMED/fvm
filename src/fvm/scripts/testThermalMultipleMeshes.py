#!/usr/bin/env python
import time


"""
Usage: testThermalParallel.py [options] infile
options are:
--type  'tri'[default], 'quad', 'hexa', or 'tetra'
--xdmf  Dump data in xdmf
"""


import fvm.fvmbaseExt as fvmbaseExt
import fvm
fvm.set_atype('double')
import  importers, fvmparallel, time
from numpy import *
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

def dumpTecplotFile(nmesh, meshes, mtype):
  #cell sites
  cellSites = []
  for n in range(0,nmesh):
     cellSites.append( meshes[n].getCells() )
     #print "cellSites[", n, "].getCount = ", cellSites[n].getCount()

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

  #get connectivity ( cellCells )
  cellCells = []
  for n in range(0,nmesh):
     cellCells.append( meshes[n].getCellCells() )
  
  #print cellCells for mesh
  #meshID = int(1)
  #for n in range(0, cellSites[meshID].getCount() ):
  #    print "cellCell[", n, "]   ",
  #    for i in range(0, cellCells[meshID].getCount(n)):
  #        print cellCells[meshID](n,i), "    ",
  #    print " "             

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
     
  f = open("temp_procSer.dat" , 'w')
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

    


parser = OptionParser()
parser.set_defaults(type='tri')
parser.add_option("--type", help="'tri'[default], 'quad', 'hexa', or 'tetra'")
parser.add_option("--xdmf", action='store_true', help="Dump data in xdmf")
parser.add_option("--time","-t",action='store_true',help="Print timing information.")
(options, args) = parser.parse_args()
if len(args) != 1:
    usage()

numIterations = 10
#import ddd
reader = FluentCase(args[0])
reader.read()

meshes = reader.getMeshList()
geomFields =  fvm.models.GeomFields('geom')
nmesh =4

metricsCalculator = fvm.models.MeshMetricsCalculatorA(geomFields,meshes)
metricsCalculator.init()

thermalFields =  fvm.models.ThermalFields('therm')
tmodel = fvm.models.ThermalModelA(geomFields,thermalFields,meshes)


## set bc for top to be at 400
#bcMap = tmodel.getBCMap()
#if 6 in bcMap:
#   bc6 = tmodel.getBCMap()[6]
#   bc6.bcType = 'SpecifiedTemperature'
#   bc6.setVar('specifiedTemperature',400)
#if 2 in bcMap:
#   bc2 = tmodel.getBCMap()[2]
#   bc2.bcType = 'SpecifiedTemperature'
#   bc2.setVar('specifiedTemperature',0)
#if 4 in bcMap:
#   bc4 = tmodel.getBCMap()[4]
#   bc4.bcType = 'SpecifiedTemperature'
#   bc4.setVar('specifiedTemperature',0)
#if 5 in bcMap:
#   bc5 = tmodel.getBCMap()[5]
#   bc5.bcType = 'SpecifiedTemperature'
#   bc5.setVar('specifiedTemperature',0)


## set bc for top to be at 400
bcMap = tmodel.getBCMap()
for id, bc in bcMap.iteritems():
    if id in [10, 1]:
        bc.bcType ='SpecifiedTemperature'
        bc.setVar('specifiedTemperature',400)
    else:
        bc.bcType ='SpecifiedTemperature'
        bc.setVar('specifiedTemperature',0)      
  


tSolver = fvmbaseExt.AMG()
tSolver.relativeTolerance = 1e-4
tSolver.nMaxIterations = 100000
tSolver.maxCoarseLevels=20
tSolver.verbosity=3

toptions = tmodel.getOptions()
toptions.linearSolver = tSolver

tmodel.init()

tmodel.advance(1)
dumpTecplotFile( nmesh, meshes, options.type)

