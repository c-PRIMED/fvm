#!/usr/bin/env python
import pdb
import sys, time
import fvm
fvm.set_atype('double')
import fvm.fvmbaseExt as fvmbaseExt
from FluentCase import FluentCase

fileBase = "/home/lin/work/app-memosa/src/fvm/test/Grid/"

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
    
reader = FluentCase(fileBase+"cav32.cas")

#import ddd
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
t0 = time.time()
#pdg.set_trace()

grid = fvmbaseExt.Grid()

grid.Impl(mesh0, geomFields, flowFields, grid, fileBase)

grids = grid.getGrids();

metricsCalculator.computeGridInterpolationMatrices(grids)

grid.computeInterpolatedVelocity(grids, grid, mesh0, geomFields, fileBase)

t1 = time.time()
if outfile != '/dev/stdout':
    print '\nsolution time = %f' % (t1-t0)

print 'solution time = %f' % (t1-t0)


