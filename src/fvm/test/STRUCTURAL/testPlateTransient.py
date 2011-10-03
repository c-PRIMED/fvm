#!/usr/bin/env python 
import sys
sys.setdlopenflags(0x100|0x2)

import fvm
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import numpy
from mpi4py import MPI

fvm.set_atype('double')

import math

if fvm.atype == 'double':
    import fvm.models_atyped_double as models
    import fvm.exporters_atyped_double as exporters
elif fvm.atype == 'tangent':
    import fvm.models_atyped_tangent_double as models
    import fvm.exporters_atyped_tangent_double as exporters

from FluentCase import FluentCase
from optparse import OptionParser

#fvmbaseExt.enableDebug("cdtor")

fileBase0 = None
fileBase1 = None
numIterations = 1
numEIterations = 1

fileBase0 = "plate_transient_"

def eadvance(fmodel,niter):
    for i in range(0,niter):
        try:
            stopFlag=fmodel.advance(1)
            if stopFlag == 1:
                break
        except KeyboardInterrupt:
            break

def advance(pmodel,geomFields,structureFields,meshes,niter,nStep,globalTime):
    for i in range(0,niter):
        try:                                                                                                                                
            eadvance(pmodel,numEIterations)
        except KeyboardInterrupt:
            break

def advanceUnsteady(pmodel,geomFields,plateFields,meshes,nTimeSteps,globalTime):
    fileName = fileBase0 + "middef.txt"
    file = open(fileName,"w")
    mesh0 = meshes[0]
    deformation =  plateFields.deformation[mesh0.getCells()].asNumPyArray()
    force = plateFields.force[mesh0.getCells()].asNumPyArray()
    for i in range(0,nTimeSteps):
        try:
#            if i==0:
#                force[:] = -10000.
#            else:
#                force[:] = -10000.
            eadvance(pmodel,numIterations)
            globalTime += timeStep
            print 'advancing to time %e at iteration %i' % (globalTime,i)
            print 'deflection = %e at iteration %i' % (deformation[819][2],i)
            file.write(" %e " % globalTime)
            file.write(" %e " % deformation[819][2])
            file.write("\n")
            pmodel.updateTime()
#            if (i%100)==0:
#                pmodel.getMoment(mesh0)
#                dumpTecplotFile(nmesh, meshes0, geomFields, options.type, i+1)
        except KeyboardInterrupt:
            break
            

# map between fvm, tecplot, and xdmf types
etype = {
    'tri' : 1,
    'quad' : 2,
    'tetra' : 3,
    'hexa' : 4
            }
tectype = {
    'tri' : 'FETRIANGLE',
    'quad' : 'FEQUADRILATERAL',
    'tetra' : 'FETETRAHEDRON',
    'hexa' : 'FEBRICK'
            }

def dumpTecplotFile(nmesh, meshes, geomFields, mtype, nStep):
    #cell sites
    cellSites = []
    for n in range(0,nmesh):
        cellSites.append( meshes[n].getCells() )
#     print "cellSites[", n, "].getCount = ", cellSites[n].getCount()

    #face sites
    faceSites = []
    for n in range(0,nmesh):
        faceSites.append( meshes[n].getFaces() )
        
    #node sites
    nodeSites = []
    for n in range(0,nmesh):
        nodeSites.append( meshes[n].getNodes() )
        
    #get connectivity (faceCells)
    faceCells = []
    for n in range(0,nmesh):
        faceCells.append( meshes[n].getConnectivity( faceSites[n], cellSites[n] ) )
        
    #get connectivity ( cellNodes )
    cellNodes = []
    for n in range(0,nmesh):
        cellNodes.append( meshes[n].getCellNodes() )
        
    #get Volume as array
    volumes = []
    for n in range(0,nmesh):
        volumes.append( geomFields.volume[cellSites[n]].asNumPyArray() )
        
    cellCentroids =[]
    for n in range(0,nmesh):
        cellCentroids.append( geomFields.coordinate[cellSites[n]].asNumPyArray() )
        
    coords = []
    for n in range(0,nmesh):
        coords.append( geomFields.coordinate[nodeSites[n]].asNumPyArray() )
#     print "shape( coords[", n, "] ) = ", shape( coords[n] )

    #deformation
    deformation = []
    for n in range(0,nmesh):
        deformation.append( plateFields.deformation[cellSites[n]].asNumPyArray() )

    #moment
    moment = []
    for n in range(0,nmesh):
        moment.append( plateFields.moment[cellSites[n]].asNumPyArray() )
                    
    f = open("./tecfsi-%0*d.dat" % (4,nStep),"w")
    f.write("Title = \" tecplot file for 2D Cavity problem \" \n")
    f.write("variables = \"x\", \"y\", \"betax\", \"betay\", \"defZ\", \"momentx\", \"momenty\",\"momentxy\", \"cellCentroidY\" \n")
    for n in range(0,nmesh):
        title_name = "nmesh%s" % n
        ncell  = cellSites[n].getSelfCount()
        nnode  = nodeSites[n].getCount()
        f.write("Zone T = \"%s\" N = %s E = %s DATAPACKING = BLOCK, VARLOCATION = ([3-9]=CELLCENTERED), ZONETYPE=%s\n" %
                (title_name,  nodeSites[n].getCount(), ncell, tectype[mtype]))
        #write x
        for i in range(0,nnode):
            f.write(str(coords[n][i][0])+"    ")
            if ( i % 5 == 4 ):
                f.write("\n")
        f.write("\n")

        #write y
        for i in range(0,nnode):
            f.write(str(coords[n][i][1])+"    ")
            if ( i % 5 == 4 ):
                f.write("\n")
        f.write("\n")

        #write betaX
        for i in range(0,ncell):
            f.write( str(deformation[n][i][0]) + "    ")
            if ( i % 5  == 4 ):
                f.write("\n")
        f.write("\n")
                
        #write betaY
        for i in range(0,ncell):
            f.write( str(deformation[n][i][1]) + "    ")
            if ( i % 5  == 4 ):
                f.write("\n")
        f.write("\n")                                                                                                

        #write defZ
        for i in range(0,ncell):
            f.write( str(deformation[n][i][2]) + "    ")
            if ( i % 5  == 4 ):
                f.write("\n")
        f.write("\n")

        #write momentX
        for i in range(0,ncell):
            f.write( str(moment[n][i][0]) + "    ")
            if ( i % 5  == 4 ):
                f.write("\n")
        f.write("\n")

        #write momentY
        for i in range(0,ncell):
            f.write( str(moment[n][i][1]) + "    ")
            if ( i % 5  == 4 ):
                f.write("\n")
        f.write("\n")

        #write momentXY
        for i in range(0,ncell):
            f.write( str(moment[n][i][2]) + "    ")
            if ( i % 5  == 4 ):
                f.write("\n")
        f.write("\n")
                                    
        #write velX
        for i in range(0,ncell):
            f.write( str(cellCentroids[n][i][1]) + "    ")
            if ( i % 5  == 4 ):
                f.write("\n")
        f.write("\n")
        
        #connectivity
        for i in range(0,ncell):
            nnodes_per_cell = cellNodes[n].getCount(i)
            for node in range(0,nnodes_per_cell):
                f.write( str(cellNodes[n](i,node)+1) + "     ")
            f.write("\n")
        f.write("\n")
    f.close()

parser = OptionParser()
parser.set_defaults(type='quad')
parser.add_option("--type", help="'quad'[default], 'tri', 'hexa', or 'tetra'")
parser.add_option("--xdmf", action='store_true', help="Dump data in xdmf")
parser.add_option("--time","-t",action='store_true',help="Print timing information.")
(options, args) = parser.parse_args()

nmesh = 1
fileBase0 = "plate_transient_"
reader0 = FluentCase(sys.argv[1])

#import debug
reader0.read();

meshes0 = reader0.getMeshList()

mesh0 = meshes0[0]
nmesh = 1

meshes = []
meshes.append(mesh0)

import time
t0 = time.time()

geomFields =  models.GeomFields('geom')
metricsCalculator0 = models.MeshMetricsCalculatorA(geomFields,meshes0)

metricsCalculator0.init()

nodes0 = mesh0.getNodes()

rho = 7854.
E = 2.*math.pow(10,11)
nu = 0.0

plateFields =  models.PlateFields('plate')
pmodel = models.PlateModelA(geomFields,plateFields,meshes0)

bcMap = pmodel.getBCMap()
for bcID in bcMap:
    if bcID == 6:
        bc = bcMap[bcID]
        bc.bcType = 'Clamped'
        bc['specifiedXRotation']=0
        bc['specifiedYRotation']=0.
        bc['specifiedZDeformation']=0.                       
    elif bcID == 5:
        bc = bcMap[bcID]
        bc.bcType = 'SpecifiedTraction'
    elif bcID ==4:
        bc = bcMap[bcID]
        bc.bcType = 'Clamped'
        bc['specifiedXRotation']=0
        bc['specifiedYRotation']=0.
        bc['specifiedZDeformation']=0.                        
    else:
        bc = bcMap[bcID]
        bc.bcType = 'SpecifiedTraction'
        
vcMap = pmodel.getVCMap()
for i,vc in vcMap.iteritems():
    vc['density'] = rho
    vc['ym'] = E
    vc['nu'] = nu

pc = fvmbaseExt.AMG()
pc.verbosity=0
defSolver = fvmbaseExt.BCGStab()
defSolver.preconditioner = pc
defSolver.relativeTolerance = 1e-9
defSolver.absoluteTolerance = 1.e-30
defSolver.nMaxIterations = 50000
defSolver.verbosity=1

poptions = pmodel.getOptions()
poptions.deformationLinearSolver = defSolver
poptions.deformationTolerance=1.0e-3
poptions.setVar("deformationURF",1.0)
poptions.printNormalizedResiduals=True
poptions.timeDiscretizationOrder = 2
poptions.transient=True
poptions.scf = 5./6.

numTimeSteps = 10
period = 1.e-5
timeStep = period/1000
globalTime=0.

# set the timesteps
poptions.setVar('timeStep',timeStep)

for mesh in meshes0:
    print
    print 'plate BCs'
    fgs = mesh.getBoundaryGroups()
    for fg in fgs:
        bc = pmodel.getBCMap()[fg.id]
        print '%i %s' %(fg.id,bc.bcType)

print '\n number of cells in mesh0 = %i' % (mesh0.getCells().getSelfCount())

xc = geomFields.coordinate[mesh0.getCells()].asNumPyArray()
cells0 = mesh0.getCells()
nCells0 = cells0.getCount()
fileName = fileBase0 + "coord.dat"
file = open(fileName,"w")
file.write("cellCoordinate\t\n")
for i in range(0,nCells0):
    file.write(" %i " % i)
    file.write(" %e " % xc[i][0])
    file.write(" %e " % xc[i][1])
    file.write(" %e " % xc[i][2])
    file.write("\n")
file.close()
                        
pmodel.init()
force = plateFields.force[mesh0.getCells()].asNumPyArray()
thickness = plateFields.thickness[mesh0.getCells()].asNumPyArray()

force[:] = -1.
thickness[:] = 2.e-6

#for i in range(numTimeSteps):
#    pmodel.advance(numIterations)
#    globalTime += timeStep
#    pmodel.updateTime()
#eadvance(pmodel,numTimeSteps)
pc.redirectPrintToFile("convergence.dat")
advanceUnsteady(pmodel,geomFields,plateFields,meshes,numTimeSteps,globalTime)
pc.redirectPrintToScreen()

xc = geomFields.coordinate[mesh0.getCells()].asNumPyArray()
deformation = plateFields.deformation[mesh0.getCells()].asNumPyArray()
pmodel.getMoment(mesh0)
moment = plateFields.moment[mesh0.getCells()].asNumPyArray()
cells0 = mesh0.getCells()
nCells0 = cells0.getSelfCount()
fileName = fileBase0 + "def.dat"
file = open(fileName,"w")
file.write("cellCoordinate\t\n")
for i in range(0,nCells0):
    file.write(" %i " % i)
    file.write(" %e " % xc[i][0])
    file.write(" %e " % xc[i][1])
    file.write(" %e " % xc[i][2])
    file.write(" %e " % moment[i][0])
    file.write(" %e " % moment[i][1])
    file.write(" %e " % moment[i][2])
    file.write(" %e " % deformation[i][2])
    file.write("\n")
file.close()                        

dumpTecplotFile(nmesh, meshes, geomFields, options.type, 0)

t1 = time.time()

print '\nsolution time = %f' % (t1-t0)
print '\n run complete '





