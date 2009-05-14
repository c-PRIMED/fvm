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
coordFile = fileBase + "Grid_Coord.dat"
velocityFile = fileBase + "Grid_Velocity.dat"

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
    
reader = FluentCase(fileBase+"3D-cantilever.cas")

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

grid = fvmbaseExt.Grid(geomFields, flowFields, coordFile, velocityFile)
gridNodes = grid.getNodes()


sideID = 10
topID = 9
botID = 12
tipID = 11

faceVelFile = open(fileBase + "InterpolatedFaceVelocity.dat", "w")
#import ddd
for mesh in meshes:
    fgs = mesh.getBoundaryGroups()
    for fg in fgs:
        if fg.id in [sideID, topID, botID, tipID]:
            faces = fg.site
            grid.setConnFaceToGrid(mesh0, faces)
            metricsCalculator.computeGridInterpolationMatrices(gridNodes,faces)
            fv = grid.computeInterpolatedVelocity(faces)
            faceVel = fv.asNumPyArray()
            faceX = geomFields.coordinate[faces].asNumPyArray()

            nFaces = faces.getCount()

            for f in range (0, nFaces):
                faceVelFile.write ("%e\t%e\t%e\t%e\t%e\t%e\n" %
                                   (faceX[f][0],faceX[f][1],faceX[f][2],
                                    faceVel[f][0],faceVel[f][1],faceVel[f][2]))


faceVelFile.close()
                                   
            
    
       
t1 = time.time()
if outfile != '/dev/stdout':
    print '\nsolution time = %f' % (t1-t0)

print 'solution time = %f' % (t1-t0)


