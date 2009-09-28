#!/usr/bin/env python
import pdb
import sys
import math
sys.setdlopenflags(0x100|0x2)
import pdb
import fvmbaseExt
import importers
from numpy import *

atype = 'double'

if atype == 'double':
    import models_atyped_double as models
    import exporters_atyped_double as exporters
elif atype == 'tangent':
    import models_atyped_tangent_double as models
    import exporters_atyped_tangent_double as exporters


from FluentCase import FluentCase
from fvmbaseExt import VecD3

atan=math.atan
cos=math.cos
sin=math.sin


numIterations =50
fileBase = "/home/lin/work/app-memosa/src/fvm/verification/IBM_convergence/"



def advance(tmodel,particles,niter):
    for i in range(0,niter):
        try:
            tmodel.computeIBFaceTemperature(particles)
            tmodel.advance(1)        
        except KeyboardInterrupt:
            break
       
outfile = fileBase+"-prism.dat"
    
reader = FluentCase(fileBase+"cav64_TScalar.cas")

mpmFileName = fileBase + "MPMs_test1.dat"

#import ddd
reader.read();

meshes = reader.getMeshList()

import time
t0 = time.time()

geomFields =  models.GeomFields('geom')

metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

thermalFields =  models.ThermalFields('therm')

tmodel = models.ThermalModelA(geomFields,thermalFields,meshes)

reader.importThermalBCs(tmodel)

bcMap = tmodel.getBCMap()

bcID = 3
if bcID in bcMap:
   bc = tmodel.getBCMap()[bcID]
   bc.bcType = 'SpecifiedTemperature'
   bc.setVar('specifiedTemperature',400)


tSolver = fvmbaseExt.AMG()
tSolver.relativeTolerance = 1e-5
tSolver.nMaxIterations = 1000
tSolver.maxCoarseLevels=20
tSolver.verbosity=0

toptions = tmodel.getOptions()

toptions.linearSolver = tSolver

#generate octree for mesh
mesh0 = meshes[0]

octree = fvmbaseExt.Octree()

octree.Impl(mesh0, geomFields)

#generate two Squares
#pdb.set_trace()

count=0   
        
temp = VecD3()
### geometry parameters ###
innerSide = 1.0
midSide = 3.0
outerSide = 4.2
### number of particles in x,y,z ###
nX=20
nY=20
nZ=1

gapX=outerSide/(nX)
gapY=outerSide/(nY)
gapZ=0

nMPM = (nX+1)*(nY+1)*nZ

position = zeros((nMPM,3))
velocity = zeros((nMPM,3))
type = zeros(nMPM,int)
temperature = zeros(nMPM,float)
alfa = atan(1.0)
for i in range(0, nX+1):
    for j in range(0,nY+1):
        for k in range(0,nZ):
            temp[0]=i*gapX-outerSide/2.
            temp[1]=j*gapY-outerSide/2.
            temp[2]=k*gapZ
            #inner square
            if temp[0]>=(-innerSide/2.0) and temp[0]<=(innerSide/2.0) and temp[1]>=(-innerSide/2.0) and temp[1]<=(innerSide/2.0):
                type[count] = 0
                if temp[0]> (innerSide/2.0-gapX) or temp[0]<(-innerSide/2.0+gapX):
                    type[count]=1
                if temp[1]> (innerSide/2.0-gapY) or temp[1]<(-innerSide/2.0+gapY):
                    type[count]=1;
                position[count][0]=temp[0]*cos(alfa)-temp[1]*sin(alfa);
                position[count][1]=temp[1]*cos(alfa)+temp[0]*sin(alfa);
                position[count][2]=temp[2];
	          
                count+=1;

            #outer square
            if( not(temp[0]>(-midSide/2.) and temp[0]<(midSide/2.) and temp[1]>(-midSide/2.) and temp[1]<(midSide/2.))):
                type[count] = 0;  
                if( temp[0]< (midSide/2.0+gapX) and temp[0]>(midSide/2.0-gapX) and temp[1]<(midSide/2.0+gapY) and temp[1]>(-midSide/2.0-gapY)):
                    type[count]=1;
                if( temp[0]> (-midSide/2.0-gapX) and temp[0]<(-midSide/2.0+gapX) and temp[1]<(midSide/2.0+gapY) and temp[1]>(-midSide/2.0-gapY)):
                    type[count]=1;
                if( temp[1]< (midSide/2.0+gapY) and temp[1]>(midSide/2.0-gapY) and temp[0]<(midSide/2.0+gapX) and temp[0]>(-midSide/2.0-gapX)):
                    type[count]=1;
                if( temp[1]> (-midSide/2.0-gapY) and temp[1]<(-midSide/2.0+gapY) and temp[0]<(midSide/2.0+gapX) and temp[0]>(-midSide/2.0-gapX)):
                    type[count]=1;
                position[count][0]=temp[0]*cos(alfa)-temp[1]*sin(alfa);
                position[count][1]=temp[1]*cos(alfa)+temp[0]*sin(alfa);
                position[count][2]=temp[2];
	          
                count+=1;

for p in range(0, count):
    velocity[p][0]=0.0;
    velocity[p][1]=0.0;
    velocity[p][2]=0.0;

for p in range(0, count):
    temperature[p]=300.0;

mpmFile = open(mpmFileName, "w")
print count
mpmFile.write("%i\n" % count)
for p in range(0, count):
    mpmFile.write("%e\t%e\t%e\n" % (position[p][0],position[p][1],position[p][2]))
for p in range(0, count):
    mpmFile.write("%e\t%e\t%e\n" % (velocity[p][0],velocity[p][1],velocity[p][2]))
for p in range(0, count):
    mpmFile.write("%i\n" % type[p])
for p in range(0, count):
    mpmFile.write("%e\n" % temperature[p])

mpmFile.close()


solid = fvmbaseExt.MPM(mpmFileName)

#mark IB cells

option = 1

fvmbaseExt.CellMark_Impl(mesh0, geomFields, fileBase, octree, solid, option)

cells = mesh0.getCells()

nCells = cells.getCount()

# add particles position and velocities to the main fields

particles = solid.getParticles()
px = solid.getCoordinates()
pt = solid.getTemperatures()


geomFields.coordinate[particles]=px
thermalFields.temperature[particles]=pt

metricsCalculator.computeIBInterpolationMatrices(particles)

tmodel.printBCs()
tmodel.init()

advance(tmodel,particles,numIterations)


t1 = time.time()


writer = exporters.FluentDataExporterA(reader,fileBase+"cav8_TScalar.dat",False,0)
writer.init()
writer.writeScalarField(thermalFields.temperature,3)
writer.finish()

