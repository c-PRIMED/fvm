#!/usr/bin/env python

### Structure and Electrostatics coupling --- 2D beam  ###

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
import optparse



### passing random variable ###


E = 200e9
nu = 0.3
rho = 8912


### ========================== define functions and class =============================###
class ProblemDefinition(object):
    pass

def writeProbeData(pd):
    cells = pd.solidMeshes[0].getCells()
    sbFaces = pd.solidBoundaryMeshes[0].getFaces()
    deformation = pd.structureFields.deformation[cells].asNumPyArray()[pd.probeIndex]
    velocity = pd.flowFields.velocity[sbFaces].asNumPyArray()[34219]
    pd.probeFile.write('%le %le %le %le\n' % (pd.globalTime, deformation[0], deformation[1], deformation[2]))
    pd.probeFile.flush()
    pd.velocityFile.write('%le %le %le %le\n' % (pd.globalTime, velocity[0], velocity[1], velocity[2]))
    pd.velocityFile.flush()
    pd.forceFile.write('%le %le %le %le\n' % (pd.globalTime, pd.fxSum, pd.fySum, pd.fzSum))
    pd.forceFile.flush()

            
def unsteadyAdvance(pd, numTimeSteps):
    
    for n in range(0, numTimeSteps):

        sbMeshFaces = pd.solidBoundaryMeshes[0].getFaces()

        #-------------update IBM----------------#
        print "***       update IBM  at  %i           ***" % n
        pd.ibManager.update()
        pd.ibManager.fluidNeighborsPerIBFace = 6
        pd.ibManager.solidNeighborsPerIBFace = 4
        pd.ibManager.fluidNeighborsPerSolidFace = 6
        
        pd.fluidMetricsCalculator.computeIBInterpolationMatrices(sbMeshFaces, 1)
        pd.fluidMetricsCalculator.computeSolidInterpolationMatrices(sbMeshFaces)     
        
        #------------solve electrostatics--------#
        print "***    solving electric model  at globalCount %i  ***" % pd.globalCount
        for i in range(0, 5):
            pd.emodel.computeIBFacePotential(sbMeshFaces)
            if pd.emodel.advance(1):
                break
        pd.emodel.computeSolidSurfaceForce(sbMeshFaces)
        
        #------------solve fluid ---------------#
        print "***      solving flow model  at globalCount %i    ***" % pd.globalCount
        pd.fmodel.computeIBFaceVelocity(sbMeshFaces)
        for i in range(0,0):
            ibfaces = pd.fluidMeshes[0].getIBFaces()
            ibv = pd.flowFields.velocity[ibfaces].asNumPyArray()
            sbv = pd.flowFields.velocity[sbMeshFaces].asNumPyArray()
            if pd.fmodel.advance(1):
                break
        pd.fmodel.computeSolidSurfaceForce(sbMeshFaces)
       
        #------------update solid force ----------#
        print "***     update force at globalCount %i             ***" % pd.globalCount
        sbForce = pd.flowFields.force[sbMeshFaces].asNumPyArray()        
        sbElecForce = pd.elecFields.force[sbMeshFaces].asNumPyArray()

        bcMap = pd.smodel.getBCMap()
        pd.fxSum = 0
        pd.fySum = 0
        pd.fzSum = 0

        for smesh in pd.solidMeshes:

            fgs = smesh.getAllFaceGroups()
            for fg in fgs:
                if fg.groupType != 'interior':
                    bc = bcMap[fg.id]
                    if bc.bcType == 'SpecifiedForce':
                        sfaces = fg.site
                        sbFaceIndices = sfaces.getCommonIndices(sbMeshFaces).asNumPyArray()
                        sFaceIndices = sbMeshFaces.getCommonIndices(sfaces).asNumPyArray()
                        faceCount = sfaces.getCount()
                        fx = pd.bForceXField[sfaces].asNumPyArray()
                        fy = pd.bForceYField[sfaces].asNumPyArray()
                        fz = pd.bForceZField[sfaces].asNumPyArray()
                        fx[:] = 0
                        fy[:] = 0
                        fz[:] = 0
                        
                        for f in range(0,faceCount):
                                                                                                                
                            fx[ sFaceIndices[f] ] = sbElecForce[ sbFaceIndices[f] ][0]
                            fy[ sFaceIndices[f] ] = sbElecForce[ sbFaceIndices[f] ][1]
                            fz[ sFaceIndices[f] ] = sbElecForce[ sbFaceIndices[f] ][2]
                            
                            pd.fxSum += fx[ sFaceIndices[f] ]
                            pd.fySum += fy[ sFaceIndices[f] ]
                            pd.fzSum += fz[ sFaceIndices[f] ]
        print "total force on device %e %e %e\n" % ( pd.fxSum, pd.fySum, pd.fzSum)

        #------------solve structure-------------#
        print "***  solving structure model at globalCount %i   ***" % pd.globalCount
        for i in range(0, 1):
            converged = pd.smodel.advance(1)
            pd.dmodel.calculateNodeDisplacement()
            pd.dmodel.deformStructure()
            pd.solidMetricsCalculator.recalculate_deform()
            if converged:
                break

        for i in range(0,len(pd.solidMeshes)):
            pd.dmodel.updateBoundaryMesh(pd.solidMeshes[i],
                                         pd.solidBoundaryMeshes[i],
                                         pd.flowFields.velocity,
                                         pd.timeStep)
        
        pd.solidBoundaryMetricsCalculator.recalculate_deform()
        
        #---------------write out data -------------------------#
        #writeProbeData(pd)              # write out backup file just in case hdf5 file failed

        cells = pd.solidMeshes[0].getCells()
        sbFaces = pd.solidBoundaryMeshes[0].getFaces()
        deformation = pd.structureFields.deformation[cells].asNumPyArray()[pd.probeIndex]
        velocity = pd.flowFields.velocity[sbFaces].asNumPyArray()[34219]
        
#        if (n%pd.saveFrequency == 0):            
#            dump_hdf5('displacement', deformation[1])
#            dump_hdf5('velocity', velocity[1])
#            dump_hdf5('totalforce', pd.fySum)
            
       
        #---------------update time -------------------------#
        
        pd.smodel.updateTime()
        pd.dmodel.updateTime()
        pd.fmodel.updateTime()
        
        pd.globalTime += pd.timeStep
        pd.globalCount += 1

        #-------------safety check --------------------------#
        maxDef = deformation[1]
        if maxDef > 2.0e-6:
            break
    
def createSolidForceBVFields(pd):
    pd.bForceXField = fvmbaseExt.Field('bForceX')
    pd.bForceYField = fvmbaseExt.Field('bForceY')
    pd.bForceZField = fvmbaseExt.Field('bForceZ')
    bcMap = pd.smodel.getBCMap()
    for smesh in pd.solidMeshes:
        fgs = smesh.getAllFaceGroups()
        for fg in fgs:
            if fg.groupType != 'interior':
                bc = bcMap[fg.id]
                if bc.bcType == 'SpecifiedForce':
                    sfaces = fg.site
                    faceCount = sfaces.getCount()
                    areaMag = pd.geomFields.areaMag[sfaces]
                    fx = areaMag.newSizedClone(faceCount)
                    fy = areaMag.newSizedClone(faceCount)
                    fz = areaMag.newSizedClone(faceCount)
                    fxa = fx.asNumPyArray()
                    fya = fy.asNumPyArray()
                    fza = fz.asNumPyArray()
                    fxa[:] = 0
                    fya[:] = 0
                    fza[:] = 0
                    pd.bForceXField[sfaces]=fx
                    pd.bForceYField[sfaces]=fy
                    pd.bForceZField[sfaces]=fz
                    bc['specifiedXForce'] = pd.bForceXField
                    bc['specifiedYForce'] = pd.bForceYField
                    bc['specifiedZForce'] = pd.bForceZField
        
### ========================== properties and parameters ===============================###

### beam
#rho = 8912                     # density kg/m^3
#E = 200e9                     # Young's modulus 
#nu = 0.31                      # Poisson's ratio

### electric field
applied_voltage = -230
dielectric_constant = 1.0

### mesh id

electrodeActive = [12, 11, 10]  # All three electrodes active
electrodeInactive = [9, 3]  
electrodeGround = [28,29]
electrodeBeamID = 7
gasWallID       = 9
 
beamTopID = 8
beamBotID = 7
beamSideID = 9
fixedBeam = [3,4,5,6]


numTimeSteps = 1

### =========================== running script  ========================================###

pd = ProblemDefinition()


fileBase = "./"

pd.globalTime = 0
pd.globalCount = 0
pd.timeStep = 1.e-8
pd.saveFrequency = 1
pd.fxSum = 0
pd.fySum = 0
pd.fzSum = 0
pd.initialTransient=False
pd.probeIndex = 671699
pd.probeFile = open(fileBase + "tipDisplacement-se.dat", "w")
pd.forceFile = open(fileBase + "beamForce-se.dat", "w")
pd.velocityFile = open(fileBase + "tipVelocity-se.dat", "w")

### read in meshes ###
beamFile = fileBase + 'Gen2_IBM_beam_12_1200.cas'
fluidFile = fileBase + 'Gen2_IBM_backgroundmesh_C2D2_wo_sub.cas'

fluidReader = FluentCase(sys.argv[1])
beamReader = FluentCase(sys.argv[2])

beamReader.read()
fluidReader.read()

pd.solidMeshes = beamReader.getMeshList()
pd.fluidMeshes = fluidReader.getMeshList()
pd.solidBoundaryMeshes = [m.extractBoundaryMesh() for m in pd.solidMeshes]

### geometry field ###

pd.geomFields =  models.GeomFields('geom')
pd.solidMetricsCalculator = models.MeshMetricsCalculatorA(pd.geomFields,
                                                          pd.solidMeshes)
pd.fluidMetricsCalculator = models.MeshMetricsCalculatorA(pd.geomFields,
                                                          pd.fluidMeshes)
pd.solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(pd.geomFields,
                                                          pd.solidBoundaryMeshes)

pd.solidMetricsCalculator.init()
pd.fluidMetricsCalculator.init()
pd.solidBoundaryMetricsCalculator.init()


### structure model and boundary condition ###

pd.structureFields =  models.StructureFields('structure')

pd.smodel = models.StructureModelA(pd.geomFields,
                                   pd.structureFields,
                                   pd.solidMeshes)

pd.dmodel = models.StructureDeformationModelA(pd.geomFields,
                                              pd.structureFields,
                                              pd.solidMeshes)


sbcMap = pd.smodel.getBCMap()
for id in fixedBeam:
	bc = sbcMap[id]
	bc.bcType = 'SpecifiedDeformation'
	bc['specifiedXDeformation'] = 0.0
	bc['specifiedYDeformation'] = 0.0
	bc['specifiedZDeformation'] = 0.0
    
for id in [beamTopID, beamBotID, beamSideID]:
    bc = sbcMap[id]
    bc.bcType = 'SpecifiedForce'


vcMap = pd.smodel.getVCMap()
for i,vc in vcMap.iteritems():
    vc['density'] = rho
    vc['eta'] = E/(2.*(1+nu))
    vc['eta1'] = nu*E/((1+nu)*(1-2.0*nu))


### electric model and boundary condition ###

pd.elecFields =  models.ElectricFields('elec')

pd.emodel = models.ElectricModelA(pd.geomFields,
                                     pd.elecFields,
                                     pd.fluidMeshes)

bcMap = pd.emodel.getBCMap()

for i in electrodeGround:
    bc = bcMap[i]
    bc.bcType = "SpecifiedPotential"
    bc['specifiedPotential'] = 0

for i in electrodeInactive:
    bc = bcMap[i]
    bc.bcType = "Symmetry"
    
for i in electrodeActive:
    bc = bcMap[i]
    bc.bcType = "SpecifiedPotential"
    bc['specifiedPotential'] = applied_voltage   
          
### setup mesh vcType ###

vcMap = pd.emodel.getVCMap()
for i,vc in vcMap.iteritems():
    vc.vcType = "dielectric"
    vc['dielectric_constant'] = dielectric_constant


### flow model and boundary condition ###
pd.flowFields =  models.FlowFields('flow')
pd.fmodel = models.FlowModelA(pd.geomFields,pd.flowFields,pd.fluidMeshes)

#fluidReader.importFlowBCs(pd.fmodel, pd.fluidMeshes)




### solver options ###

# structure solver
pc = fvmbaseExt.AMG()
pc.verbosity=0
#defSolver = fvmbaseExt.BCGStab()
defSolver = fvmbaseExt.BCGStab()
defSolver.preconditioner = pc
defSolver.relativeTolerance = 1e-1
defSolver.nMaxIterations = 20
defSolver.maxCoarseLevels=20
defSolver.verbosity=2


soptions = pd.smodel.getOptions()
soptions.deformationLinearSolver = defSolver
soptions.deformationTolerance=5.0e-1
soptions.setVar("deformationURF",1.0)
soptions.printNormalizedResiduals=False
soptions.transient=True

soptions.setVar("timeStep", pd.timeStep)

# electric solver 
epc = fvmbaseExt.AMG()
epc.verbosity=0
elecSolver = fvmbaseExt.BCGStab()
elecSolver.preconditioner = epc
elecSolver.relativeTolerance = 1e-3
elecSolver.nMaxIterations = 1000
elecSolver.maxCoarseLevels=20
elecSolver.verbosity=0

eoptions = pd.emodel.getOptions()
eoptions.electrostaticsLinearSolver = elecSolver
eoptions.electrostaticsTolerance = 0.5e-5
eoptions.electrostatics_enable = 1
eoptions.chargetransport_enable = 0
eoptions.tunneling = 0
eoptions.ibm_enable = 1
eoptions.transient_enable = False
eoptions.printNormalizedResiduals = True

# flow solver
foptions = pd.fmodel.getOptions()
foptions.momentumTolerance=1e-4
foptions.continuityTolerance=1e-4
foptions.setVar("momentumURF",0.7)
foptions.setVar("pressureURF",0.3)
foptions.transient=True

foptions.setVar("timeStep",pd.timeStep)

flowVCMap = pd.fmodel.getVCMap()
for id, vc in flowVCMap.iteritems():
    vc['density'] = 1.225
    vc['viscosity'] = 1.8e-5
foptions.printNormalizedResiduals=False

### set up IBM ###

for mesh in pd.solidBoundaryMeshes:
    faces = mesh.getFaces()
    #potential field on boundary mesh
    areaMag = pd.geomFields.areaMag[faces]
    faceCount = faces.getCount()
    pot = areaMag.newSizedClone(faceCount)
    pota = pot.asNumPyArray()
    pota[:] = 0.0
    pd.elecFields.potential[faces] = pot
    #flow velocity field on boundary mesh
    area = pd.geomFields.area[faces]
    faceCount = faces.getCount()
    vel = area.newSizedClone(faceCount)
    vela = vel.asNumPyArray()
    vela[:,:] = 0.0
    pd.flowFields.velocity[faces] = vel
    
pd.ibManager = fvmbaseExt.IBManager(pd.geomFields,
                                    pd.solidBoundaryMeshes[0],
                                    pd.fluidMeshes)

pd.ibManager.fluidNeighborsPerIBFace = 2
pd.ibManager.solidNeighborsPerIBFace = 6
pd.ibManager.fluidNeighborsPerSolidFace = 6


### initialize models and run ###

pd.smodel.init()
pd.dmodel.init()
pd.emodel.init()
pd.fmodel.init()

#initially, zero force is applied on beam
createSolidForceBVFields(pd)

t1 = time.time()
pc.redirectPrintToFile("convergence.dat")
unsteadyAdvance(pd,numTimeSteps)
pc.redirectPrintToScreen()
#checkMarking(pd.globalCount, pd)

t2 = time.time()

print  '\nsolution time = %f' % (t2-t1)

