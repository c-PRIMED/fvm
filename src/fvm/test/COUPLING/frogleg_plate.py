#!/usr/bin/env python

### this script solve frogleg pull-in by coupling 
### plate model and electrostatics via IBM


### import  modules ###
import pdb
import sys
import os
from math import *
sys.setdlopenflags(0x100|0x2)
#import tecplotExporter
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
from FluentCase import FluentCase
from fvm.fvmbaseExt import VecD3
from mpi4py import MPI
import time

def checkMarking(n):

    cells = fluidMeshes[0].getCells()
    nCells = cells.getCount()
    cellCoords = geomFields.coordinate[cells].asNumPyArray()

    fluidFile = open(fileBase + "fluidCells_" + str(n) + ".dat", "w")
    solidFile = open(fileBase + "solidCells_" + str(n) + ".dat", "w")
    IBFile = open(fileBase + "IBCells_" + str(n) + ".dat", "w")

    cellIBType = geomFields.ibType[cells].asNumPyArray()
    for c in range (0, nCells):    
        ibtype = cellIBType[c]
        if ibtype == -1:
            fluidFile.write("%e\t%e\t%e\n" % (cellCoords[c][0], cellCoords[c][1], cellCoords[c][2]))
        elif ibtype == -2:
            IBFile.write("%e\t%e\t%e\n" % (cellCoords[c][0], cellCoords[c][1], cellCoords[c][2]))
        elif ibtype == -3:
            solidFile.write("%e\t%e\t%e\n" % (cellCoords[c][0], cellCoords[c][1], cellCoords[c][2]))
        elif ibtype == -5:
            print ("%i\t%i\t%e\t%e\n" % (c,ibtype,  cellCoords[c][0], cellCoords[c][1]))

               
    fluidFile.close()
    solidFile.close()
    IBFile.close()

def writeProbeData():
    deformation = plateFields.deformation[solidMeshes[0].getCells()].asNumPyArray()
    maxDef = deformation.min(axis = 0)
    probeFile.write('%e\t%e\t%e\n' % (globalTime, deformation[probeIndex][2], maxDef[2]))
    probeFile.flush() 

def saveVTK(n):
    writer = exporters.VTKWriterA(geomFields,fluidMeshes,
                                  fileBase + "elecfield-" + str(n) + ".vtk",
                                  "frogleg",
                                  False,0)
    writer.init()
    writer.writeScalarField(elecFields.potential,"potential")
    writer.writeVectorField(elecFields.electric_field,"potentialgradient")    
    writer.finish()

    writer1 = exporters.VTKWriterA(geomFields,solidMeshes,
                                  fileBase + "deformation-" + str(n) + ".vtk",
                                  "frogleg",
                                  False,0)
    writer1.init()
    writer1.writeVectorField(plateFields.deformation,"deformation")
    writer1.finish()

### ========================== properties and parameters ===============================###

### beam
rho = 8912                     # density kg/m^3
E = 200e9                     # Young's modulus 
nu = 0.31                      # Poisson's ratio

### electric field
applied_voltage = -100
dielectric_constant = 1.0
beam_thickness = 3e-6

### mesh id
fluidTop = 7
fluidBot = [9,3]
fluidSide = 8
left_electrode = 4
right_electrode = 5
central_electrode = 6
 
beam = [3]
anchors = [4, 5, 6, 7]

numTimeSteps = 1
globalTime = 0
globalCount = 0
timeStep = 5e-8
saveFrequency = 50
initialTransient = False
probeIndex = 50

### ===================== mesh read ===============================================###

fileBase = "../../"

### 2D plate mesh
beamReader = FluentCase(fileBase+"frogleg_2D_5018.cas")
beamReader.read();
solidMeshes = beamReader.getMeshList()
geomFields =  models.GeomFields('geom')
solidMetricsCalculator = models.MeshMetricsCalculatorA(geomFields,solidMeshes)
solidMetricsCalculator.init()

### 3D fluid mesh
fluidReader = FluentCase(fileBase+"fluid.cas")
fluidReader.read();
fluidMeshes = fluidReader.getMeshList()
fluidMetricsCalculator = models.MeshMetricsCalculatorA(geomFields,fluidMeshes)
fluidMetricsCalculator.init()

nodes = fluidMeshes[0].getNodes()
xn = fluidMeshes[0].getNodeCoordinates().asNumPyArray()
for n in range(0, nodes.getCount()):
    x = xn[n][0]
    y = xn[n][1]
    xn[n][0] = -y
    xn[n][1] = x
fluidMetricsCalculator.init()


### generate solid boundary mesh
solidBoundaryMeshes = [m.extrude(1, beam_thickness, True) for m in solidMeshes]
solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(geomFields,solidBoundaryMeshes)
solidBoundaryMetricsCalculator.init()


### find device center

cells = solidMeshes[0].getCells()
xc = geomFields.coordinate[cells].asNumPyArray()
small = 100
probeIndex = 0
for c in range(0, cells.getCount()):
    rsqr = xc[c][0]*xc[c][0] + xc[c][1]*xc[c][1]
    if rsqr < small:
        small = rsqr
        probeIndex = c



### output files 
probeFile = open(fileBase + "centerDisplacement.dat", "w")

### =============================== models =====================================###

### Plate Model and boundary conditions ###
plateFields =  models.PlateFields('plate')
pmodel = models.PlateModelA(geomFields,plateFields,solidMeshes)
dmodel = models.PlateDeformationModelA(geomFields,plateFields,solidMeshes)

bcMap = pmodel.getBCMap()

for id in anchors:
    bc = bcMap[id]
    bc.bcType = 'Clamped'
    bc['specifiedXRotation']=0
    bc['specifiedYRotation']=0.
    bc['specifiedZDeformation']=0. 
for id in beam:
    bc = bcMap[id]
    bc.bcType = 'SpecifiedTraction'

vcMap = pmodel.getVCMap()
for i,vc in vcMap.iteritems():
    vc['density'] = rho
    vc['ym'] = E
    vc['nu'] = nu

### electric model and boundary condition ###

elecFields =  models.ElectricFields('elec')

emodel = models.ElectricModelA(geomFields,elecFields,fluidMeshes)

bcMap = emodel.getBCMap()

bc = bcMap[central_electrode]
bc.bcType = "SpecifiedPotential"
bc['specifiedPotential'] = applied_voltage

bc = bcMap[left_electrode]
bc.bcType = "SpecifiedPotential"
bc['specifiedPotential'] = 0.0

bc = bcMap[right_electrode]
bc.bcType = "SpecifiedPotential"
bc['specifiedPotential'] = 0.0

bc = bcMap[fluidTop]
bc.bcType = "Symmetry"

for i in fluidBot:
    bc = bcMap[i]
    bc.bcType = "SpecifiedPotential"
    bc['specifiedPotential'] = 0.0

for i in [fluidSide]:
    bc = bcMap[i]
    bc.bcType = "Symmetry"

vcMap = emodel.getVCMap()
for i,vc in vcMap.iteritems():
    vc.vcType = "dielectric"
    vc['dielectric_constant'] = dielectric_constant    


### ================================= solvers ===================================###

### plate solver ###
pc = fvmbaseExt.AMG()
pc.verbosity=0
defSolver = fvmbaseExt.BCGStab()
defSolver.preconditioner = pc
defSolver.relativeTolerance = 1e-9
defSolver.absoluteTolerance = 1.e-30
defSolver.nMaxIterations = 50000
defSolver.verbosity=0

poptions = pmodel.getOptions()
poptions.deformationLinearSolver = defSolver
poptions.deformationTolerance=1.0e-3
poptions.setVar("deformationURF",1.0)
poptions.printNormalizedResiduals=True
poptions.timeDiscretizationOrder = 2
poptions.transient=True
poptions.scf = 5./6.
poptions.setVar('timeStep',timeStep)

### elec solver ###

epc = fvmbaseExt.AMG()
epc.verbosity=0
elecSolver = fvmbaseExt.BCGStab()
elecSolver.preconditioner = epc
elecSolver.relativeTolerance = 1e-3
elecSolver.nMaxIterations = 1000
elecSolver.maxCoarseLevels=20
elecSolver.verbosity=0

eoptions = emodel.getOptions()
eoptions.electrostaticsLinearSolver = elecSolver
eoptions.electrostaticsTolerance = 0.5e-5
eoptions.electrostatics_enable = 1
eoptions.chargetransport_enable = 0
eoptions.tunneling = 0
eoptions.ibm_enable = 1
eoptions.transient_enable = False
eoptions.printNormalizedResiduals = True


### initialize models and run ###

pmodel.init()
emodel.init()
dmodel.init()


ibManager = fvmbaseExt.IBManager(geomFields,
                                 solidBoundaryMeshes[0],
                                 fluidMeshes)
for mesh in solidBoundaryMeshes:
        faces = mesh.getFaces()
        areaMag = geomFields.areaMag[faces]
        faceCount = faces.getCount()
        pot = areaMag.newSizedClone(faceCount)
        pota = pot.asNumPyArray()
        pota[:] = 0
        elecFields.potential[faces] = pot

sbMeshFaces = solidBoundaryMeshes[0].getFaces()
ibManager.fluidNeighborsPerIBFace = 4
ibManager.solidNeighborsPerIBFace = 4
ibManager.fluidNeighborsPerSolidFace = 6
ibManager.update()
checkMarking(globalCount)



t1 = time.time()
pc.redirectPrintToFile("convergence.dat")
#--------------Timestep Loop --------------------------#

for n in range(0, numTimeSteps):                
   
    # --------------- update IBM -------------------------#
    print "***       update IBM  at globalCount %i           ***" % globalCount            
    
    ibManager.update()
    fluidMetricsCalculator.computeIBInterpolationMatrices(sbMeshFaces)
    fluidMetricsCalculator.computeSolidInterpolationMatrices(sbMeshFaces)        

    #------------solve electrostatics--------#
    print "***    solving electric model  at globalCount %i  ***" % globalCount
    for i in range(0, 10):
        emodel.computeIBFacePotential(sbMeshFaces)
        emodel.advance(1)
           
    emodel.computeSolidSurfaceForcePerUnitArea(sbMeshFaces)
    #saveVTK(n)
    #------------update force on beam  ----------#
    print "***     update force at globalCount %i             ***" % globalCount
    
    sbElecForce = elecFields.force[sbMeshFaces].asNumPyArray()
    
    solidMesh = solidMeshes[0]
    solidCells = solidMesh.getCells()
    nCells = solidCells.getCount()
    nSelfCells = solidCells.getSelfCount()
    
    nSBFaces = sbMeshFaces.getCount()

    if (nSBFaces != 2*nSelfCells+(nCells-nSelfCells)):
        print "the extruded solid boundary mesh has wrong face numbers!"

    force = plateFields.force[solidCells].asNumPyArray()
    thickness = plateFields.thickness[solidCells].asNumPyArray()
    force[:] = 0.
    thickness[:] = beam_thickness
    
    # force on interior cells
    for c in range(0, nSelfCells):
        botFaceIndex = c
        topFaceIndex = c+nSelfCells
        force[c] = sbElecForce[botFaceIndex][2] + sbElecForce[topFaceIndex][2]
        
    # force on boundary cells
    for c in range(nSelfCells, nCells):
        force[c] = sbElecForce[nSelfCells+c][2]    
    
    
    #pdb.set_trace()
    #------------solve structure-------------#
    print "***  solving structure model at globalCount %i   ***" % globalCount
    for i in range (0, 3):    
        pmodel.advance(1)
    dmodel.calculateNodeDisplacement()
    dmodel.deformPlate()
    solidMetricsCalculator.recalculate_deform()

    
    #------------update solid boundary mesh---------------#
    #solidBoundaryMeshes = [m.extrude(1, beam_thickness, True) for m in solidMeshes]
    sbNodes = solidBoundaryMeshes[0].getNodes()
    nSBNodes = sbNodes.getCount()
    nodes = solidMeshes[0].getNodes()
    if nSBNodes != nodes.getCount()*2:
        print "the extruded solid mesh has wrong node number!"
    nodeCoord = geomFields.coordinate[nodes].asNumPyArray()
    bNodeCoord = geomFields.coordinate[sbNodes].asNumPyArray()
    bMeshCoord = solidBoundaryMeshes[0].getNodeCoordinates().asNumPyArray()

    deformation = geomFields.nodeDisplacement[nodes].asNumPyArray()
    #pdb.set_trace()
    for sbn in range (0, nSBNodes/2):
        bNodeCoord[sbn][2] =  -beam_thickness/2 + nodeCoord[sbn][2]
        bMeshCoord[sbn][2] =  -beam_thickness/2 +  nodeCoord[sbn][2]
    for sbn in range (nSBNodes/2, nSBNodes):
        bNodeCoord[sbn][2] = beam_thickness/2 + nodeCoord[sbn - nSBNodes/2][2] 
        bMeshCoord[sbn][2] = beam_thickness/2 + nodeCoord[sbn - nSBNodes/2][2] 
    #pdb.set_trace()
    #solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(geomFields,solidBoundaryMeshes)
    #solidBoundaryMetricsCalculator.init()
    solidBoundaryMetricsCalculator.recalculate_deform() 
    

    # -----------------update time --------------------------#
    pmodel.updateTime()
    dmodel.updateTime()

    globalTime += timeStep
    globalCount += 1

    #------data output-----------------------#
    writeProbeData()
    
    if (n%saveFrequency == 0):
        saveVTK(n)
    #    checkMarking(globalCount)
 
t2 = time.time()
pc.redirectPrintToScreen()


probeFile.close()
print  '\nsolution time = %f' % (t2-t1)
