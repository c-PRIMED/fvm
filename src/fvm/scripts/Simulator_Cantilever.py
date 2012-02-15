from Tools import *
from ComputeForce import *
from TimeStep import *
from Persistence import Persistence
from OnDemandContactModel import *

class Simulator():

    def __init__(self, meshes, models, ibm_update,
                 gap1=0, gap2=0,
                 timeStep=1e-8, minR=10e-9, maxR=50e-9, dc = 7.9):
        
        self.transient = models.pmodel.getOptions().transient
        self.timeStep = timeStep
        self.timeStepN1 = timeStep
        self.timeStepN2 = timeStep
        self.timeOption = models.pmodel.getOptions().variableTimeStep
        self.ibm_update = ibm_update        
	self.minR = minR
	self.maxR = maxR
	
        self.enableElecModel = models.enableElecModel
        self.enablePlateModel = models.enablePlateModel
        self.enableFlowModel = models.enableFlowModel
        self.enableContactModel = models.enableContactModel

        self.geomFields = meshes.geomFields        
        self.solidMeshes = meshes.solidMeshes
        self.fluidMeshes = meshes.fluidMeshes
        self.fluidMeshesNew = meshes.fluidMeshesNew

        self.solidBoundaryMeshes = meshes.solidBoundaryMeshes
        self.solidMetricsCalculator = meshes.solidMetricsCalculator
        self.solidBoundaryMetricsCalculator = meshes.solidBoundaryMetricsCalculator
        self.solidCells = self.solidMeshes[0].getCells()
        self.nSolidSelfCells = self.solidCells.getSelfCount()
        self.nSolidCells = self.solidCells.getCount()
        self.xCells = self.geomFields.coordinate[self.solidCells].asNumPyArray()
        self.sbMeshFaces = self.solidBoundaryMeshes[0].getFaces()
        self.beam_thickness = meshes.beam_thickness
        self.beam_width = meshes.beam_width
        self.beam_length = meshes.beam_length	
        self.dielectric_constant = dc
        self.dielectric_thickness = meshes.dielectric_thickness
        print "dielectric thickness in sim %e" % self.dielectric_thickness
        self.gap = meshes.gap
        self.gap1 = gap1
        self.gap2 = gap2
        self.rn  = []
        
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

        #if self.enableContactModel:
        #    self.cmodel = models.cmodel
        #    self.contactFields = models.contactFields
                    
        self.globalCount = 0
        self.globalTime = 0.0
        self.vel = 0.0
        self.acc = 0.0
        self.elecForceSum = 0.0
        self.flowForceSum = 0.0
        self.contactForceSum = 0.0

    def restart(self, restartFile):
        self.globalCount = restartFile.readAttribute('globalCount')
        self.globalTime = restartFile.readAttribute('globalTime') 
        self.timeStep = restartFile.readAttribute('timeStep')
        self.timeStepN1 = restartFile.readAttribute('timeStepN1')
        self.timeStepN2 = restartFile.readAttribute('timeStepN2')
        self.poptions.setVar('timeStep',self.timeStep)
        self.poptions.timeStepN1 = self.timeStepN1
        self.poptions.timeStepN2 = self.timeStepN2
        restartFile.close()

    def saveRestartFile(self, outDir, n):
        saveFileName = outDir + 'checkpoint_' + str(n) + '.hdf5'
        f = Persistence(saveFileName,'w')
        f.saveAttribute('globalCount', self.globalCount)
        f.saveAttribute('globalTime', self.globalTime)
        ts = self.timeStep
        tsN1 = self.poptions.timeStepN1
        tsN2 = self.poptions.timeStepN2
        f.saveAttribute('timeStep', ts)
        f.saveAttribute('timeStepN1', tsN1)
        f.saveAttribute('timeStepN2', tsN2)
        f.saveFluidMeshes(self.fluidMeshes)
        f.saveSolidMeshes(self.solidMeshes)
        f.savePlateModel(self.plateFields,self.pmodel,self.solidMeshes)
        f.saveElectricModel(self.elecFields,self.emodel,self.fluidMeshes)
        f.close()     


    def run(self, bias, switch, tag):
        self.beamForce[:] = 0.0
        self.elecForceSum = 0.0
        self.flowForceSum = 0.0
        self.contactForceSum = 0.0
        self.cloestDistance = self.deformation.min(axis=0)[2] + self.gap - self.dielectric_thickness
        self.voltage = bias
        self.switch = switch
        print '======================================================================='
        #print 'current applied voltage %f' % self.voltage
        print 'Marching at global count %i' % self.globalCount
        print 'Marching at global time %e' % self.globalTime
        print 'cloest distance %e' % self.cloestDistance
        print '----------------------------------------------------------------------'
        ### use IBM to calculate electrostatic force ###
        if switch == 1:
            if self.enableElecModel:
                self.ibm_update()
                for i in range(0, 50):
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
           
        ### simplied model to calculate electrostatic force for fix fix beam ###   
        if switch == 0:       
            if self.enableElecModel:
                for c in range (0, self.nSolidCells):
                    distance = self.gap + self.deformation[c][2] - self.dielectric_thickness
                    realBias = calculateVoltage(self.voltage, 1.0, distance, self.dielectric_constant, self.dielectric_thickness)
                    eforce = - computeElectrostaticForce(realBias, distance)
           	    self.beamForce[c] += eforce
                    self.elecForceSum += eforce

        ### simplied model to calculate electrostatic force for tilted cantilever ###   
        if switch == 2:  
            startPos = 1e-6             # actuation pad starting location
            endPos   = startPos+200e-6    # actuation pad ending location, 200e-6 is dielectric length
            if self.enableElecModel:
                for c in range (0, self.nSolidCells):
                    x = self.xCells[c][0]
                    if x > startPos and x < endPos :
                        realgap = self.gap1 + (x/self.beam_length)*(self.gap2-self.gap1)
                        distance = realgap + self.deformation[c][2] - self.dielectric_thickness
                        realBias = calculateVoltage(self.voltage, 1.0, distance, 
                                                    self.dielectric_constant, self.dielectric_thickness)
                        eforce = - computeElectrostaticForce(realBias, distance)
                        self.beamForce[c] += eforce
                        self.elecForceSum += eforce            
        ### use compact damping model to compute damp force ###    
        if self.enableFlowModel:
                for c in range (0, self.nSolidCells):
                    x = self.xCells[c][0]
                    realgap = self.gap1 + (x/self.beam_length)*(self.gap2-self.gap1)
                    distance = realgap + self.deformation[c][2]- self.dielectric_thickness
                    vel = self.velocity[c][2]
                    width = self.beam_width
                    length = self.beam_length
                    if tag == 'pullin':
                    	Cf = computeDampCoeff(distance, width, self.beam_thickness)
                    if tag == 'pullout':
                    	Cf = computeDampCoeff(distance, width, self.beam_thickness)
                    dforce = computeDampForce(vel, Cf)
                    area = width
                    dforce /= area
                    self.beamForce[c] += dforce
                    self.flowForceSum += dforce
        
        ### use new contact model to compute contact force ###    
        if self.enableContactModel:
        	if self.globalCount == 0:
        	    self.rn = generateRandomNumber(self.nSolidCells)
                for c in range (0, self.nSolidCells):
                    x = self.xCells[c][0]
                    realgap = self.gap1 + (x/self.beam_length)*(self.gap2-self.gap1)
                    distance = realgap + self.deformation[c][2]- self.dielectric_thickness
                    cforce = computeContactForce(distance)
                    #cforce = computeContactForceOnDemand(0,0,0,0, distance, c, self.rn)
                    self.beamForce[c] += cforce
                    self.contactForceSum += cforce
       
         
        self.pmodel.advance(2)        
        self.dmodel.calculateNodeDisplacement()
        self.dmodel.deformPlate()
        self.solidMetricsCalculator.recalculate_deform()
        
        nodeCoordMesh = self.solidMeshes[0].getNodeCoordinates().asNumPyArray()
        nodeCoord = self.geomFields.coordinate[self.solidMeshes[0].getNodes()].asNumPyArray()
        nodeCoordMesh[:,:] = nodeCoord[:,:]

        self.dmodel.updateBoundaryMesh(self.solidMeshes[0], 
                                           self.solidBoundaryMeshes[0], 
                                           self.beam_thickness)
        self.solidBoundaryMetricsCalculator.recalculate_deform()  
        print 'update boundary mesh done!'
        
        
        if self.transient == True:
            self.pmodel.updateTime()
            self.dmodel.updateTime()
          
            self.globalTime += self.timeStep
            #print 'update time done!'
        self.globalCount += 1

        self.acc = findMax(self.acceleration)
        self.vel = findMax(self.velocity)
        
        if self.timeOption == True:
            mints = 1e-5
            dr = computeTravelDistance(self.cloestDistance, self.gap, self.minR, self.maxR)
            for c in range(0, self.nSolidCells):
                acc = fabs(self.acceleration[c])
                vel = fabs(self.velocity[c][2])
                ts = computeTimeStep(dr, vel, acc, self.minR, self.maxR)
                if ts > 0 and ts < mints:
                    mints = ts
                #if ts > 0 and ts < 5e-10:
                #    mints = 5e-10
            
            self.timeStep = mints 
            self.poptions.setVar('timeStep',self.timeStep)    
        print 'time step %e' % self.timeStep  
	print 'cloest distance %e' % self.cloestDistance

        if tag == 'pullin' and (self.cloestDistance < 5e-9 or self.globalTime >= 1.0e-3):
            return False
        if tag == 'pullout' and self.cloestDistance > self.gap*0.99:
            return False        	

            
    
