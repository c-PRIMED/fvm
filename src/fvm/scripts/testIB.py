#!/usr/bin/env python
import pdb
import sys
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt
import importers


atype = 'double'
#atype = 'tangent'

if atype == 'double':
    import models_atyped_double as models
    import exporters_atyped_double as exporters
elif atype == 'tangent':
    import models_atyped_tangent_double as models
    import exporters_atyped_tangent_double as exporters


from FluentCase import FluentCase


fileBase = "/home/sm/prism-meshes/"

def advance(fmodel,particles,niter):
    for i in range(0,niter):
        try:
            fmodel.computeIBFaceVelocity(particles)
            if fmodel.advance(1):
                break
        except KeyboardInterrupt:
            break
        
def saveData(flowFields,reader,fileBase):
    writer = exporters.FluentDataExporterA(reader,fileBase+"-prism.dat",False,0)
    
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
    outfile = fileBase+"-cellmark.dat"

caseBase = fileBase+"cav32"
reader = FluentCase(caseBase+".cas")

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

#import debug

solid = fvmbaseExt.MPM()

octree = fvmbaseExt.Octree() 

option = 2

fvmbaseExt.CellMark_Impl(mesh0, geomFields, fileBase, octree, solid, option)

cells = mesh0.getCells()

nCells = cells.getSelfCount()


# add particles position and velocities to the main fields

particles = solid.getParticles()
px = solid.getCoordinates()
pv = solid.getVelocities()
geomFields.coordinate[particles]=px
flowFields.velocity[particles]=pv

metricsCalculator.computeIBInterpolationMatrices(particles)

advance(fmodel,particles,200)

saveData(flowFields,reader,caseBase)




