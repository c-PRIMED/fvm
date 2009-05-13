#!/usr/bin/env python
import pdb
import sys
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt
import importers
from numpy import *
from mpi4py  import MPI
import time

atype = 'double'
#atype = 'tangent'

if atype == 'double':
    import models_atyped_double as models
    import exporters_atyped_double as exporters
elif atype == 'tangent':
    import models_atyped_tangent_double as models
    import exporters_atyped_tangent_double as exporters


from FluentCase import FluentCase

fileBase = sys.argv[1]
#fileBase = "/home/yildirim/memosa/src/fvm/test/"


def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-cellmark.dat if it is not specified."
    sys.exit(1)


#generate geometry and move to center (xc,yc,zc)
def particle_coordinates(mesh0, geometry, nradius, ntheta, xc, yc, zc):
   nparticles = nradius * ntheta
   px = geomFields.coordinate[mesh0.getCells()].newSizedClone( int(nparticles) )

   if  geometry == "cylinder2d":
       radius = 0.125
       dtheta = 2.0 * math.pi / float(ntheta)
       dr      = radius / float(nradius)
       theta = 0.0;
       indx = 0;
       ppx = px.asNumPyArray()
       for n in range(0, ntheta):
           theta = theta + dtheta
           r = 0.0;
           for i in range(0, nradius):
               r = r + dr;
               ppx[indx,0] = xc + r * math.cos(theta)
               ppx[indx,1] = yc + r * math.sin(theta)
               ppx[indx,2] = 0.0;
               indx = indx + 1

   return px


#tecplot particles
def tecplot_particles(mesh0, fvmParticles, geomFields, nsweep ):
     pxFVM         = fvmParticles.getCellIDs( 0)
     particleID    = pxFVM.asNumPyArray()       
     nparticles = size( particleID )
     cells = mesh0.getCells()
     xc = geomFields.coordinate[cells].asNumPyArray()

     file_name = "cavity_nsweep" +  str( nsweep ) + ".dat"
     f = open(file_name, 'w')
     zone_name = "Zone T ="  + "\"" + "sweep " + str(nsweep) + "\"" +  "\n"
     f.write( zone_name )
     for i in range(0,nparticles):
        cell_id = particleID[i]
        line = str(xc[cell_id][0]) + "       " + str(xc[cell_id][1]) + "        " +str(xc[cell_id][2]) + "\n"
        f.write( line )

     f.close()





outfile = None
if __name__ == '__main__' and fileBase is None:
    if len(sys.argv) < 2:
        usage()
    fileBase = sys.argv[1]
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

if outfile == None:
    outfile = fileBase+"-cellmark.dat"


reader = FluentCase( fileBase )

#import ddd
reader.read();

meshes = reader.getMeshList()

mesh0 = meshes[0]

geomFields =  models.GeomFields('geom')

metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

if atype == 'tangent':
    metricsCalculator.setTangentCoords(0,7,1)

flowFields =  models.FlowFields('flow')

fmodel = models.FlowModelA(geomFields,flowFields,meshes)

reader.importFlowBCs(fmodel)
fmodel.init()
import time
t0 = time.time()

solid = fvmbaseExt.MPM()

octree = fvmbaseExt.Octree() 

option = 1

nradius = 50;
ntheta  = 320;
nparticles = nradius * ntheta 

px = particle_coordinates(mesh0, "cylinder2d", nradius, ntheta, 0.5, 0.5, 0)
pv = geomFields.coordinate[mesh0.getCells()].newSizedClone( int(nparticles) )
pType =fvmbaseExt.newIntArray(nparticles)
pType.asNumPyArray()[:] = 1

print "px  = ", px.asNumPyArray() 
print "size = ", size( px.asNumPyArray() )
solid.setCoordinates(px)
solid.setVelocities(pv)
solid.setTypes(pType)
particles = solid.getParticles( nparticles )
geomFields.coordinate[particles]=px
flowFields.velocity[particles]  =pv

fvmbaseExt.CellMark_Impl(mesh0, geomFields, fileBase, octree, solid, option)

fvmParticles = fvmbaseExt.FVMParticles( meshes )
nsweep = sys.argv[2]
print "nsweep = ", nsweep
fvmParticles.setParticles( int(nsweep) )
tecplot_particles( mesh0, fvmParticles, geomFields, nsweep )

pxFVM  = fvmParticles.getCellIDs( 0)
print pxFVM.asNumPyArray()



cells = mesh0.getCells()

nCells = cells.getSelfCount()

print nCells





t1 = time.time()
if outfile != '/dev/stdout':
    print '\nsolution time = %f' % (t1-t0)

print 'solution time = %f' % (t1-t0)


