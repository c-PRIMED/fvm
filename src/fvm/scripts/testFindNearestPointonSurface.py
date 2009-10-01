#!/usr/bin/env python
import pdb
import sys
import math
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt
import importers
import models_atyped_double as models
import exporters_atyped_double as exporters
from FluentCase import FluentCase
from fvmbaseExt import VecD3
#fvmbaseExt.enableDebug("cdtor")



numIterations =20
fileBase = "/home/lin/work/app-memosa/src/fvm/verification/ElectroStatics/"

def advance(elec_model,niter):
    for i in range(0,niter):
        try:
            stopFlag=elec_model.advance(1)
            if stopFlag == 1:
                break
        except KeyboardInterrupt:
            break



reader = FluentCase(fileBase + "cav32_elec.cas")

reader.read()

#pdb.set_trace()

meshes = reader.getMeshList()

mesh0 = meshes[0]

import time
t0 = time.time()

geomFields =  models.GeomFields('geom')

metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

electronicFields =  models.ElectronicFields('elec')

elec_model = models.ElectronicModelA(geomFields,electronicFields,meshes)

fgs = mesh0.getBoundaryGroups()

testPoint = VecD3()

testPoint[0] = 0.0
testPoint[1] = 0.0
testPoint[2] = 0.0

#import ddd

octreeCells = fvmbaseExt.Octree()

octreeCells.Impl(mesh0, geomFields)

#########  search top surface  ####################

bcID = 6

octreeTop = fvmbaseExt.Octree()

octreeTop.Create(mesh0, geomFields, bcID)

nearestPoint = elec_model.findClosestPoint(testPoint, octreeTop)

for fg in fgs:
    if fg.id == bcID:
        faces = fg.site
        faceCentroid = geomFields.coordinate[faces].asNumPyArray()

fId = nearestPoint

endPoint = VecD3()

endPoint[0] = faceCentroid[fId][0]
endPoint[1] = faceCentroid[fId][1]
endPoint[2] = faceCentroid[fId][2]

print 'testPoint = %f %f %f' % (testPoint[0], testPoint[1], testPoint[2])
print 'nearestPoint ID = %i' % nearestPoint
print 'nearestPoint = %f %f %f' % (faceCentroid[fId][0],faceCentroid[fId][1],faceCentroid[fId][2])

N=10

radius = 0.13

#import ddd

pathPoints = elec_model.createPathAndDiscretize(testPoint,endPoint, N)

print 'pathpoints'
print pathPoints

pathPointValues = elec_model.pathPointsInterpolation(pathPoints, octreeCells, radius)

#########  search bot surface  ####################        
"""
bcID = 3

octreeBot = fvmbaseExt.Octree()

octreeBot.Create(mesh0, geomFields, bcID)

nearestPoint = elec_model.findClosestPoint(testPoint, octreeBot)

for fg in fgs:
    if fg.id == bcID:
        faces = fg.site
        faceCentroid = geomFields.coordinate[faces].asNumPyArray()

fId = nearestPoint

print 'testPoint = %f %f %f' % (testPoint[0], testPoint[1], testPoint[2])
print 'nearestPoint ID = %i' % nearestPoint
print 'nearestPoint = %f %f %f' % (faceCentroid[fId][0],faceCentroid[fId][1],faceCentroid[fId][2])

"""
                                   
t1 = time.time()

print 'solution time = %f' % (t1-t0)


