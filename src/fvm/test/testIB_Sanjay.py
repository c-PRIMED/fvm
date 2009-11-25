#!/usr/bin/env python
import pdb
import sys, fvm
fvm.set_atype('double')
import fvm.fvmbaseExt as fvmbaseExt
from FluentCase import FluentCase

fileBase = "/home/lin/work/app-memosa/src/fvm/test/"

def advance(fmodel,particles,niter):
    for i in range(0,niter):
        try:
            fmodel.computeIBFaceVelocity(particles)
            if fmodel.advance(1):
                break
        except KeyboardInterrupt:
            break
        
def saveData(flowFields,reader,fileBase):
    writer = exporters.FluentDataExporterA(reader,fileBase+"testIB.dat",False,0)
    
    writer.init()
    writer.writeScalarField(flowFields.pressure,1)
    writer.writeVectorField(flowFields.velocity,111)
    writer.writeScalarField(flowFields.massFlux,18)
    writer.finish()

def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-cellmark.dat if it is not specified."
    sys.exit(1)

outfile = None
if __name__ == '__main__' and fileBase is None:
    if len(sys.argv) < 2:
        usage()
    fileBase = sys.argv[1]
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

if outfile == None:
    outfile = fileBase+"output.dat"


reader = FluentCase(fileBase+"cav32-new.cas")

reader.read();

meshes = reader.getMeshList()

mesh0 = meshes[0]

geomFields =  models.GeomFields('geom')

metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

if fvm.atype == 'tangent':
    metricsCalculator.setTangentCoords(0,7,1)

flowFields =  models.FlowFields('flow')

fmodel = models.FlowModelA(geomFields,flowFields,meshes)

reader.importFlowBCs(fmodel)
fmodel.init()


momSolver = fvmbaseExt.AMG()
momSolver.relativeTolerance = 1e-1
momSolver.nMaxIterations = 20
momSolver.maxCoarseLevels=20

momSolver.verbosity=0

contSolver = fvmbaseExt.AMG()
#pc = fvmbaseExt.AMG()
#pc.verbosity=0
#contSolver = fvmbaseExt.BCGStab()
#contSolver.preconditioner = pc
contSolver.relativeTolerance = 1e-1
contSolver.nMaxIterations = 20
contSolver.verbosity=0
contSolver.maxCoarseLevels=20

foptions = fmodel.getOptions()

foptions.momentumLinearSolver = momSolver
foptions.pressureLinearSolver = contSolver

import time
t0 = time.time()

#import ddd

solid = fvmbaseExt.MPM()

octree = fvmbaseExt.Octree() 

option = 2

fvmbaseExt.CellMark_Impl(mesh0, geomFields, fileBase, octree, solid, option)

cells = mesh0.getCells()

nCells = cells.getCount()

# add particles position and velocities to the main fields

particles = solid.getParticles()
px = solid.getCoordinates()
pv = solid.getVelocities()
geomFields.coordinate[particles]=px
flowFields.velocity[particles]=pv

metricsCalculator.computeIBInterpolationMatrices(particles)

advance(fmodel,particles,1)

saveData(flowFields,reader,fileBase)




