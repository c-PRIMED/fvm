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
    
reader = FluentCase(fileBase+"new-beam-114k.cas")

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
#pdg.set_trace()

grid = fvmbaseExt.Grid()

grids = grid.getGrids();

grid.Impl(mesh0, geomFields, flowFields, grid, fileBase)


for mesh in meshes:
    fgs = mesh.getBoundaryGroups()
    for fg in fgs:
        if fg.id == 9:   #beam top
            faces = fg.site
            grid.setConnFaceToGrid(mesh0, geomFields, grid, faces)
            metricsCalculator.computeGridInterpolationMatrices(grids,faces)
            faceVel = grid.computeInterpolatedVelocity(grids, grid, mesh0, geomFields, fileBase, faces)
    
        if fg.id == 12:   #beam bot
            faces = fg.site
            grid.setConnFaceToGrid(mesh0, geomFields, grid, faces)
            metricsCalculator.computeGridInterpolationMatrices(grids,faces)
            faceVel = grid.computeInterpolatedVelocity(grids, grid, mesh0, geomFields, fileBase, faces)
            
        if fg.id == 11:   #beam tip
            faces = fg.site
            grid.setConnFaceToGrid(mesh0, geomFields, grid, faces)
            metricsCalculator.computeGridInterpolationMatrices(grids,faces)
            faceVel = grid.computeInterpolatedVelocity(grids, grid, mesh0, geomFields, fileBase, faces)
            
        if fg.id == 10:   #beam side
            faces = fg.site
            grid.setConnFaceToGrid(mesh0, geomFields, grid, faces)
            metricsCalculator.computeGridInterpolationMatrices(grids,faces)
            faceVel = grid.computeInterpolatedVelocity(grids, grid, mesh0, geomFields, fileBase, faces)
  

t1 = time.time()
if outfile != '/dev/stdout':
    print '\nsolution time = %f' % (t1-t0)

print 'solution time = %f' % (t1-t0)


