from Tools import *
from ComputeForce import *
from TimeStep import *


class Simulator():

    def __init__(self, meshes, models, ibm_update, threshold, voltage,
                 transient=True, timeStep=1e-8, timevary=True):
        
        self.transient = transient
        self.timeStep = timeStep
        self.timeOption = timevary
        self.ibm_update = ibm_update
        self.threshold = threshold
        self.voltage = voltage

        self.enableElecModel = models.enableElecModel
        self.enablePlateModel = models.enablePlateModel
        self.enableFlowModel = models.enableFlowModel
        self.enableContactModel = models.enableContactModel

        self.geomFields = meshes.geomFields        
        self.solidMeshes = meshes.solidMeshes
        self.fluidMeshes = meshes.fluidMeshes
        self.solidBoundaryMeshes = meshes.solidBoundaryMeshes
        self.solidMetricsCalculator = meshes.solidMetricsCalculator
        self.solidBoundaryMetricsCalculator = meshes.solidBoundaryMetricsCalculator
        self.solidCells = self.solidMeshes[0].getCells()
        self.nSolidSelfCells = self.solidCells.getSelfCount()
        self.nSolidCells = self.solidCells.getCount()
        self.sbMeshFaces = self.solidBoundaryMeshes[0].getFaces()
        self.beam_thickness = meshes.beam_thickness
        self.gap = meshes.gap
        
        if self.enableElecModel:
            self.emodel = models.emodel
            self.elecFields = models.elecFields            

        if self.enablePlateModel:
            self.pmodel = models.pmodel
            self.poptions = self.pmodel.getOptions()
            self.plateFields = models.plateFields
            self.dmodel = models.dmodel
            self.beamForce = self.plateFields.force[self.solidCells].asNumPyArray()
            self.beamThickness = self.plateFields.thickness[self.solidCells].asNumPyArray()
            self.beamThickness[:] = self.beam_thickness
            self.deformation = self.plateFields.deformation[self.solidCells].asNumPyArray()
            self.acceleration = self.plateFields.acceleration[self.solidCells].asNumPyArray()
            self.velocity = self.plateFields.velocity[self.solidCells].asNumPyArray()

        if self.enableFlowModel:
            self.fmodel = models.fmodel
            self.flowFields = models.flowFields           

        if self.enableContactModel:
            self.cmodel = models.cmodel
            self.contactFields = models.contactFields
                    
        self.globalCount = 0
        self.globalTime = 0.0
        self.vel = 0.0
        self.acc = 0.0
        self.elecForceSum = 0.0
        self.flowForceSum = 0.0
        self.contactForceSum = 0.0

    def run(self, bias):
        self.beamForce[:] = 0.0
        self.elecForceSum = 0.0
        self.flowForceSum = 0.0
        self.contactForceSum = 0.0
        self.cloestDistance = self.deformation.min(axis=0)[2] + self.gap
        self.voltage = bias
        print 'current applied voltage %f' % self.voltage
        print 'Marching at global count %i' % self.globalCount
        print 'Marching at global time %e' % self.globalTime
        print '----------------------------------------------------------------------'
        ### use IBM to calculate electrostatic force, damping force and contact force ###
        if self.cloestDistance > self.threshold:

            if self.enableElecModel:
                self.ibm_update()
                for i in range(0, 20):
                    self.emodel.computeIBFacePotential(self.sbMeshFaces)
                    if self.emodel.advance(1):
                        break
                self.emodel.computeSolidSurfaceForcePerUnitArea(self.sbMeshFaces)
                self.elecForce = self.elecFields.force[self.sbMeshFaces].asNumPyArray()
                for c in range(0, self.nSolidSelfCells):
                    botFaceIndex = c
                    topFaceIndex = c+self.nSolidSelfCells
                    self.beamForce[c] += self.elecForce[botFaceIndex][2] + self.elecForce[topFaceIndex][2]
                    self.elecForceSum += self.elecForce[botFaceIndex][2] + self.elecForce[topFaceIndex][2]
                for c in range(self.nSolidSelfCells, self.nSolidCells):
                    self.beamForce[c] += self.elecForce[self.nSolidSelfCells+c][2]  
                    self.elecForceSum += self.elecForce[self.nSolidSelfCells+c][2] 
         
            if self.enableFlowModel:
            	for i in range(0,30):
                    self.fmodel.computeIBFaceVelocity(self.sbMeshFaces)
                    if self.fmodel.advance(1):
                        break
                self.fmodel.computeSolidSurfaceForcePerUnitArea(self.sbMeshFaces)
                self.flowForce = self.flowFields.force[self.sbMeshFaces].asNumPyArray()
                for c in range(0, self.nSolidSelfCells):
                    botFaceIndex = c
                    topFaceIndex = c+self.nSolidSelfCells
                    self.beamForce[c] += self.flowForce[botFaceIndex][2] + self.flowForce[topFaceIndex][2]
                    self.flowForceSum += self.flowForce[botFaceIndex][2] + self.flowForce[topFaceIndex][2]
                for c in range(self.nSolidSelfCells, self.nSolidCells):
                    self.beamForce[c] += self.flowForce[self.nSolidSelfCells+c][2]  
                    self.flowForceSum += self.flowForce[self.nSolidSelfCells+c][2] 

            if self.enableContactModel:
                self.cmodel.computeSolidSurfaceForcePerUnitArea(self.sbMeshFaces)
                self.contactForce = self.contactFields.force[self.sbMeshFaces].asNumPyArray()
                for c in range(0, self.nSolidSelfCells):
                    botFaceIndex = c
                    topFaceIndex = c+self.nSolidSelfCells
                    self.beamForce[c] += self.contactForce[botFaceIndex][2] + self.contactForce[topFaceIndex][2]
                    self.contactForceSum += self.contactForce[botFaceIndex][2] + self.contactForce[topFaceIndex][2]
                for c in range(self.nSolidSelfCells, self.nSolidCells):
                    self.beamForce[c] += self.contactForce[self.nSolidSelfCells+c][2]  
                    self.contactForceSum += self.contactForce[self.nSolidSelfCells+c][2] 

        ### simplied model to calculate electrostatic force, damping force and contact force ###   
        else:       
            if self.enableElecModel:
                for c in range (0, self.nSolidCells):
                    distance = self.gap + self.deformation[c][2]
                    eforce = - computeElectrostaticForce(self.voltage, distance)
                    self.beamForce[c] += eforce
                    self.elecForceSum += eforce
            
            if self.enableFlowModel:
                for c in range (0, self.nSolidCells):
                    distance = self.gap + self.deformation[c][2]
                    vel = self.velocity[c][2]
                    width = 120e-6
                    length = 500e-6
                    Cf = computeDampCoeff(distance, width, self.beam_thickness)
                    dforce = computeDampForce(vel, Cf)
                    area = width
                    dforce /= area
                    self.beamForce[c] += dforce
                    self.flowForceSum += dforce
            
            if self.enableContactModel:
                for c in range (0, self.nSolidCells):
                    distance = self.gap + self.deformation[c][2]
                    cforce = computeContactForce(distance)
                    self.beamForce[c] += cforce
                    self.contactForceSum += cforce
            
        #print '----------------------------------------------------------------------'
        
        for i in range (0, 2):    
            self.pmodel.advance(1)
        #print 'solve beam plate deformation done!'

        self.dmodel.calculateNodeDisplacement()
        self.dmodel.deformPlate()
        self.solidMetricsCalculator.recalculate_deform()
        #print 'deform beam done!'

        if self.enableFlowModel == False: 
            self.dmodel.updateBoundaryMesh(self.solidMeshes[0], 
                                           self.solidBoundaryMeshes[0], 
                                           self.beam_thickness)
        else:
            self.dmodel.updateBoundaryMesh(self.solidMeshes[0], 
                                           self.solidBoundaryMeshes[0], 
                                           self.beam_thickness, 
                                           self.timeStep, 
                                           self.flowFields.velocity)
        self.solidBoundaryMetricsCalculator.recalculate_deform()  
        print 'update boundary mesh done!'
        
        #print '----------------------------------------------------------------------'
        
        if self.transient == True:
            self.pmodel.updateTime()
            self.dmodel.updateTime()
            if self.enableFlowModel == True:
                self.fmodel.updateTime()
            self.globalTime += self.timeStep
            #print 'update time done!'
        self.globalCount += 1

        self.acc = findMax(self.acceleration)
        self.vel = findMax(self.velocity)
        
        if self.timeOption == True:
            mints = 5.0e-8
            dr = computeTravelDistance(self.cloestDistance, self.gap)
            for c in range(0, self.nSolidCells):
                acc = fabs(self.acceleration[c])
                vel = fabs(self.velocity[c][2])
                ts = computeTimeStep(dr, vel, acc)
                if ts > 1e-15 and ts < mints:
                    mints = ts
            if mints > 0:
                self.timeStep = mints 
            self.poptions.setVar('timeStep',self.timeStep)    
            print 'time step %e' % self.timeStep  


        if self.cloestDistance < 0:
            return False

    

        
            
    
