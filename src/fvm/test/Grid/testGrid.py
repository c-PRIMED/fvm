#!/usr/bin/env python
import sys
sys.setdlopenflags(0x100|0x2)
import os, fvmbaseExt, importers
from mpi4py import MPI

atype = 'double'
#atype = 'tangent'
if atype == 'double':
    import models_atyped_double as models
    import exporters_atyped_double as exporters
elif atype == 'tangent':
    import models_atyped_tangent_double as models
    import exporters_atyped_tangent_double as exporters

from FluentCase import FluentCase

fileBase = os.getcwd()
coordFile = os.path.join(fileBase, "Grid_Coord.dat")
velocityFile = os.path.join(fileBase, "Grid_Velocity.dat")
infile = os.path.join(fileBase, "3D-cantilever.cas")
outfile = '/dev/stdout'

reader = FluentCase(infile)
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


faceVelFile = open(outfile, "w")
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


