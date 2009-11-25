#!/usr/bin/env python
import fvm
fvm.set_atype('double')
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
from numpy import *
from mpi4py import MPI
import sys, time
from FluentCase import FluentCase

fileBase = ''
casefile = sys.argv[1]
nsweep = int(sys.argv[2])

def usage():
    print "Usage: %s casefile nsweep" % sys.argv[0]
    print "Where casefile is a Fluent case file."
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
               r += dr;
               ppx[indx,0] = xc + r * math.cos(theta)
               ppx[indx,1] = yc + r * math.sin(theta)
               ppx[indx,2] = 0.0;
               indx += 1

   return px


#tecplot particles
def tecplot_particles(mesh0, fvmParticles, geomFields, nsweep ):
     pxFVM         = fvmParticles.getCellIDs( 0)
     particleID    = pxFVM.asNumPyArray()       
     nparticles = size( particleID )
     cells = mesh0.getCells()
     xc = geomFields.coordinate[cells].asNumPyArray()
     f = open("cavity_nsweep%s.dat" % nsweep, 'w')
     f.write( "Zone T =\"sweep %d\"\n" % nsweep)
     for i in range(0,nparticles):
        cell = xc[particleID[i]]
        f.write("%f %f %s\n" % (cell[0], cell[1], cell[2]))
     f.close()


reader = FluentCase(casefile )
reader.read();

meshes = reader.getMeshList()
mesh0 = meshes[0]

geomFields =  fvm.models.GeomFields('geom')

metricsCalculator = fvm.models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

if fvm.atype == 'tangent':
    metricsCalculator.setTangentCoords(0,7,1)

flowFields =  fvm.models.FlowFields('flow')

fmodel = fvm.models.FlowModelA(geomFields,flowFields,meshes)

reader.importFlowBCs(fmodel)
fmodel.init()

t0 = time.time()
solid = fvmbaseExt.MPM()
octree = fvmbaseExt.Octree()
octree.Impl(mesh0, geomFields)

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

#print "nsweep = ", nsweep
fvmParticles.setParticles(nsweep)
tecplot_particles( mesh0, fvmParticles, geomFields, nsweep )

pxFVM  = fvmParticles.getCellIDs( 0)
#print pxFVM.asNumPyArray()

cells = mesh0.getCells()
nCells = cells.getSelfCount()
#print nCells

t1 = time.time()
print 'solution time = %f' % (t1-t0)


