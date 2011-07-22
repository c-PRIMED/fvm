#!/usr/bin/env python

### Generic file for setting up steady state - electrostatics, plate model ###

import pdb
import sys
from math import *
#sys.setdlopenflags(0x100|0x2)

import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
from FluentCase import FluentCase
from fvm.fvmbaseExt import VecD3
import time
from mpi4py  import MPI
import numpy as np

class problemDescription():
 
    def __init__(self, beamMesh,  backgroundMesh,  beamThickness, initialGap):
        ## Read in 2d Solid Mesh
        self.beam_thickness = beamThickness
        self.Gap = initialGap
        beamReader = FluentCase(beamMesh)
        beamReader.read();
        print "read solid mesh"
        self.solidMeshes = beamReader.getMeshList()
        self.geomFields =  models.GeomFields('geom')
        self.solidMetricsCalculator = models.MeshMetricsCalculatorA(self.geomFields,self.solidMeshes)
        self.solidMetricsCalculator.init()
        
        ## Define plate and deformation model
        self.plateFields =  models.PlateFields('plate')
        self.pmodel = models.PlateModelA(self.geomFields,self.plateFields,self.solidMeshes)
        self.dmodel = models.PlateDeformationModelA(self.geomFields,self.plateFields,self.solidMeshes)
        bcMap = self.pmodel.getBCMap()
        
        ## Apply a default Boundary Condition
        #for i, bc in bcMap.iteritems():
            #bc.bcType = 'SpecifiedTraction'
            
        ## Read in 3d Background Mesh
        fluidReader = FluentCase(backgroundMesh)
        fluidReader.read();
        self.fluidMeshes = fluidReader.getMeshList()
        self.fluidMetricsCalculator = models.MeshMetricsCalculatorA(self.geomFields,self.fluidMeshes)
        self.fluidMetricsCalculator.init()

        ## Define electric model
        self.elecFields =  models.ElectricFields('elec')
        self.emodel = models.ElectricModelA(self.geomFields,self.elecFields,self.fluidMeshes)
        bcMap = self.emodel.getBCMap()
        ## Apply Default boundary conditions
        for i, bc in bcMap.iteritems():
            bc.bcType = "Symmetry"
         
        self.solidBoundaryMeshes = [m.extrude(1, beamThickness, True) for m in self.solidMeshes]
        self.solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(self.geomFields,self.solidBoundaryMeshes)
        self.solidBoundaryMetricsCalculator.init()   
    
    def materialPropsBeam(self, rho, E, nu):
      print "--Setting Beam Material Properties--"
      vcMap = self.pmodel.getVCMap()
      for i,vc in vcMap.iteritems():
        vc['density'] = rho
        vc['ym'] = E
        vc['nu'] = nu
        
    
    def materialPropsBackground(self, dielectricConstant):
        print "--Setting Background Material Properties--"
        vcMap = self.emodel.getVCMap()
        for i,vc in vcMap.iteritems():
            vc.vcType = "dielectric"
            vc['dielectric_constant'] = dielectricConstant

    def beamClamped(self, zoneIds):
        print "--Apply Fixed Boundary Condition--"
        bcMap = self.pmodel.getBCMap()
        for bcId in zoneIds:
            bc = bcMap[bcId]
            bc.bcType = 'Clamped'
            bc['specifiedXRotation']=0
            bc['specifiedYRotation']=0.
            bc['specifiedZDeformation']=0.    
    def beamFree(self, zoneIds):
            print "--Apply Zero Traction Boundary Condition--"
            bcMap = self.pmodel.getBCMap()
            for bcId in zoneIds:
                bc = bcMap[bcId]
                bc.bcType = 'SpecifiedTraction'
                
    def appliedVoltages(self, zoneIds, potential):
        print "--Apply Voltage to Electrodes --"
        for bcId in zoneIds:
            bcMap = self.emodel.getBCMap()
            bc = bcMap[bcId]
            bc.bcType = "SpecifiedPotential"
            bc['specifiedPotential'] = potential
    def initializeModels(self):
        print "--Solving Models--"
        self.maxdeformation = []
        self.numIterations = 20
        self.globalTime = 0
        self.globalCount = 0
        self.timeStep = 3600
        self.saveFrequency = 2
        self.saveVTKstate = True
        initialTransient = False
        self.probeIndex = 3240
        
        pc = fvmbaseExt.AMG()
        pc.verbosity=0
        defSolver = fvmbaseExt.DirectSolver()
        defSolver.preconditioner = pc
        defSolver.relativeTolerance = 1e-6
        defSolver.absoluteTolerance = 1.e-15
        defSolver.nMaxIterations = 50000
        defSolver.verbosity=1

        poptions = self.pmodel.getOptions()
        poptions.deformationLinearSolver = defSolver
        poptions.deformationTolerance=1.0e-6
        poptions.setVar("deformationURF",1.0)
        poptions.printNormalizedResiduals=True
        poptions.timeDiscretizationOrder = 2
        poptions.transient=True
        poptions.scf = 5./6.
        poptions.setVar('timeStep',self.timeStep)

        ### elec solver ###

        epc = fvmbaseExt.AMG()
        epc.verbosity=0
        elecSolver = fvmbaseExt.BCGStab()
        elecSolver.preconditioner = epc
        elecSolver.relativeTolerance = 1e-6
        elecSolver.nMaxIterations = 1000
        elecSolver.maxCoarseLevels=30
        elecSolver.verbosity=0

        eoptions = self.emodel.getOptions()
        eoptions.electrostaticsLinearSolver = elecSolver
        eoptions.electrostaticsTolerance = 0.5e-5
        eoptions.electrostatics_enable = 1
        eoptions.chargetransport_enable = 0
        eoptions.tunneling = 0
        eoptions.ibm_enable = 1
        eoptions.transient_enable = False
        eoptions.printNormalizedResiduals = True


        ### initialize models and run ###

        self.pmodel.init()
        self.emodel.init()
        self.dmodel.init() 
        
    def saveVTK(self, state):
        self.saveVTKstate = state
        
    def solveModels(self, appliedVoltage):
        maxdeformation = []
        ibManager = fvmbaseExt.IBManager(self.geomFields,
                                     self.solidBoundaryMeshes[0],
                                     self.fluidMeshes)
                                     
        for mesh in self.solidBoundaryMeshes:
            faces = mesh.getFaces()
            areaMag = self.geomFields.areaMag[faces]
            faceCount = faces.getCount()
            pot = areaMag.newSizedClone(faceCount)
            pota = pot.asNumPyArray()
            pota[:] = appliedVoltage
            self.elecFields.potential[faces] = pot
            
        sbMeshFaces = self.solidBoundaryMeshes[0].getFaces()
        ibManager.fluidNeighborsPerIBFace = 4
        ibManager.solidNeighborsPerIBFace = 4
        ibManager.fluidNeighborsPerSolidFace = 6
        t1 = time.time()
  
        traceFile = open(("tracefile-%e.dat" % appliedVoltage), "w")
        #--------------Timestep Loop --------------------------#

        for n in range(0, self.numIterations):                
            #    checkMarking(globalCount)
            # --------------- update IBM -------------------------#
            print "***       update IBM  at globalCount %i           ***" % self.globalCount            
            
            ibManager.update()
            self.fluidMetricsCalculator.computeIBInterpolationMatrices(sbMeshFaces)
            self.fluidMetricsCalculator.computeSolidInterpolationMatrices(sbMeshFaces)        

            #------------solve electrostatics--------#
            print "***    solving electric model  at globalCount %i  ***" % self.globalCount
            for i in range(0, 20):
                self.emodel.computeIBFacePotential(sbMeshFaces)
                self.emodel.advance(1)
                   
                self.emodel.computeSolidSurfaceForcePerUnitArea(sbMeshFaces)
            
            #------------update force on beam  ----------#
            print "***     update force at globalCount %i             ***" % self.globalCount
            
            sbElecForce = self.elecFields.force[sbMeshFaces].asNumPyArray()
            
            solidMesh = self.solidMeshes[0]
            solidCells = solidMesh.getCells()
            nCells = solidCells.getCount()
            nSelfCells = solidCells.getSelfCount()
            
            nSBFaces = sbMeshFaces.getCount()

            if (nSBFaces != 2*nSelfCells+(nCells-nSelfCells)):
                print "the extruded solid boundary mesh has wrong face numbers!"

            force = self.plateFields.force[solidCells].asNumPyArray()
            thickness = self.plateFields.thickness[solidCells].asNumPyArray()
            force[:] = 0.
            ### Need to correct~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            thickness[:] = self.beam_thickness
            
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
            print "***  solving structure model at globalCount %i   ***" % self.globalCount
            self.pmodel.advance(1)
            self.dmodel.calculateNodeDisplacement()
            self.dmodel.deformPlate()
            self.solidMetricsCalculator.recalculate_deform()

            
            #------------update solid boundary mesh---------------#
            #solidBoundaryMeshes = [m.extrude(1, beam_thickness, True) for m in solidMeshes]
            sbNodes = self.solidBoundaryMeshes[0].getNodes()
            nSBNodes = sbNodes.getCount()
            nodes = self.solidMeshes[0].getNodes()
            if nSBNodes != nodes.getCount()*2:
                print "the extruded solid mesh has wrong node number!"
            nodeCoord = self.geomFields.coordinate[nodes].asNumPyArray()
            bNodeCoord = self.geomFields.coordinate[sbNodes].asNumPyArray()
            bMeshCoord = self.solidBoundaryMeshes[0].getNodeCoordinates().asNumPyArray()

            self.deformation = self.geomFields.nodeDisplacement[nodes].asNumPyArray()
            #pdb.set_trace()
            for sbn in range (0, nSBNodes/2):
                bNodeCoord[sbn][2] =  -self.beam_thickness/2 + nodeCoord[sbn][2]
                bMeshCoord[sbn][2] =  -self.beam_thickness/2 +  nodeCoord[sbn][2]
            for sbn in range (nSBNodes/2, nSBNodes):
                bNodeCoord[sbn][2] = self.beam_thickness/2 + nodeCoord[sbn - nSBNodes/2][2] 
                bMeshCoord[sbn][2] = self.beam_thickness/2 + nodeCoord[sbn - nSBNodes/2][2] 
            #pdb.set_trace()
            #solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(geomFields,solidBoundaryMeshes)
            #solidBoundaryMetricsCalculator.init()
            self.solidBoundaryMetricsCalculator.recalculate_deform() 

            cell = self.solidMeshes[0].getCells()
            
            deform = self.plateFields.deformation[cell].asNumPyArray()           
            temp = np.array(deform[:])
            val = 0.e0
            
            for i in range(len(temp)):
                temp2 = temp[i]
                if (temp2[2] < val):
                    val = temp2[2]
            
            maxdeformation.append(val)            
            traceFile.write("%i %e\n" %(n, maxdeformation[n]))
            traceFile.flush()
            
            # -----------------update time --------------------------#
            self.globalTime += self.timeStep
            self.globalCount += 1

            if (n%self.saveFrequency == 0 and self.saveVTKstate):
                writer = exporters.VTKWriterA(self.geomFields,self.fluidMeshes,
                                              "elecfield-"+ str(appliedVoltage) + "V-" + str(n) + ".vtk",
                                              "fix-fix beam",
                                              False,0)
                writer.init()
                writer.writeScalarField(self.elecFields.potential,"potential")
                writer.writeVectorField(self.elecFields.electric_field,"potentialgradient")
                writer.finish()

                writer1 = exporters.VTKWriterA(self.geomFields,self.solidMeshes,
                                               "structural-" + str(appliedVoltage) + "V-" +str(n) + ".vtk",
                                               "Disk",
                                               False,0)
                writer1.init()
                writer1.writeVectorField(self.plateFields.deformation,"Deformation")
                writer1.finish()
 
            if (abs(maxdeformation[n] ) > self.Gap):
                traceFile.write("Pull In")
                traceFile.flush()
                maxdeformation[n] = -self.Gap
                break
            if (n!=0):
                print abs((maxdeformation[n]-maxdeformation[n-1])/maxdeformation[n])
                if (abs((maxdeformation[n]-maxdeformation[n-1])/maxdeformation[n]) < 1e-3):
                    print "Convergence reached"
                    break

        t2 = time.time()
        print  '\nsolution time = %f' % (t2-t1)
        traceFile.close()
        return maxdeformation[n]




