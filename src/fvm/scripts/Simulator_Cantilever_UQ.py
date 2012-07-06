from Tools import *
from ComputeForce import *
from TimeStep import *
from Persistence import Persistence
from OnDemandContactModel import *
###################################################################################
epsilon0 = 8.854187817e-12
#H = 0.215e-18
#c1 = 1.2e3
#c2 = 0.5e9
#c3 = 1.0e3
#dc = 110e-9
#tR = 8.5e-9
ds = 30e-9
d_bar = 1e-9
ta = 2e-9

def func(x):
    dx = 0
    if x > 0:
        dx = 0
    else:
        dx = -x
    return dx

def computeContactForceUQ(x, c1, c2, c3, H, dc, tR):
    attractive_force = -H/(6*pi*(pow(x,3)+pow(ds,3))) -c3*func(x-(dc+ta))/d_bar
    repulsive_force = c1*pow((func(x-dc)/d_bar),1.5) + c2*exp(-x/tR)
    contact_force = attractive_force + repulsive_force
    return contact_force

###################################################################################
class Simulator():

    def __init__(self, meshes, models, ibm_update,
                 gap1=0, gap2=0,
                 timeStep=1e-8, minR=10e-9, maxR=200e-9, dc = 7.9):
        self.meshes = meshes
        self.models = models
        self.ibm_update = ibm_update
        
        self.transient = self.models.pmodel.getOptions().transient
        self.timeStep = timeStep
        self.timeStepN1 = timeStep
        self.timeStepN2 = timeStep
        self.timeOption = self.models.pmodel.getOptions().variableTimeStep
            
	self.minR = minR
	self.maxR = maxR
	
        self.enableElecModel = self.models.enableElecModel
        self.enablePlateModel = self.models.enablePlateModel
        self.enableFlowModel = self.models.enableFlowModel
        self.enableContactModel = self.models.enableContactModel

        self.geomFields = self.meshes.geomFields        
        self.solidMeshes = self.meshes.solidMeshes
        self.fluidMeshes = self.meshes.fluidMeshes
        self.fluidMeshesNew = self.meshes.fluidMeshesNew

        self.solidBoundaryMeshes = self.meshes.solidBoundaryMeshes
        self.solidMetricsCalculator = self.meshes.solidMetricsCalculator
        self.solidBoundaryMetricsCalculator = self.meshes.solidBoundaryMetricsCalculator
        self.solidCells = self.solidMeshes[0].getCells()
        self.nSolidSelfCells = self.solidCells.getSelfCount()
        self.nSolidCells = self.solidCells.getCount()
        self.xCellsA = self.geomFields.coordinate[self.solidCells]
        self.xCells = self.xCellsA.asNumPyArray()
        self.sbMeshFaces = self.solidBoundaryMeshes[0].getFaces()
        self.beam_thickness = self.meshes.beam_thickness
        self.beam_width = self.meshes.beam_width
        self.beam_length = self.meshes.beam_length	
        self.dielectric_constant = dc
        self.dielectric_thickness = self.meshes.dielectric_thickness
        print "dielectric thickness in sim %e" % self.dielectric_thickness
        self.gap = self.meshes.gap
        self.gap1 = gap1
        self.gap2 = gap2
        self.rn  = []
        
        if self.enableElecModel:
            self.emodel = self.models.emodel
            self.elecFields = self.models.elecFields            

        if self.enablePlateModel:
            self.pmodel = self.models.pmodel
            self.poptions = self.pmodel.getOptions()
            self.plateFields = self.models.plateFields
            self.dmodel = self.models.dmodel
            self.beamForceA = self.plateFields.force[self.solidCells]
            self.beamForce = self.beamForceA.asNumPyArray()
            self.beamThicknessA = self.plateFields.thickness[self.solidCells]
            self.beamThickness = self.beamThicknessA.asNumPyArray()
            self.beamThickness[:] = self.beam_thickness
            self.deformationA = self.plateFields.deformation[self.solidCells]
            self.deformation = self.deformationA.asNumPyArray()
            self.accelerationA = self.plateFields.acceleration[self.solidCells]
            self.acceleration = self.accelerationA.asNumPyArray()
            self.velocityA = self.plateFields.velocity[self.solidCells]
            self.velocity = self.velocityA.asNumPyArray()

        if self.enableFlowModel:
            self.fmodel = self.models.fmodel
            self.flowFields = self.models.flowFields           

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


    def run(self, bias, switch, tag, source, c1, c2, c3, H, dc, tR ):
        self.beamForce[:] = 0.0
        self.elecForceSum = 0.0
        self.flowForceSum = 0.0
        self.contactForceSum = 0.0
        self.cloestDistance = 1.0
        self.farDistance = 0
        self.source = source
        for c in range (0, self.nSolidCells):
            x = self.xCells[c][0]
            realgap = self.gap1 + (x/self.beam_length)*(self.gap2-self.gap1)
            distance = realgap + self.deformation[c][2] - self.dielectric_thickness
            if distance < self.cloestDistance:
            	self.cloestDistance = distance
            if distance > self.farDistance:
                self.farDistance = distance
            	
        self.voltage = bias
        self.switch = switch
        print '======================================================================='
        #print 'current applied voltage %f' % self.voltage
        print 'Marching at global count %i' % self.globalCount
        print 'Marching at global time %e' % self.globalTime
        print 'cloest distance %e' % self.cloestDistance
        print 'furthest distance %e' % self.farDistance
        print 'time step %e' % self.timeStep  
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
                self.elecForceA = self.elecFields.force[self.sbMeshFaces]
                self.elecForce = self.elecForceA.asNumPyArray()
                
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
                    #realBias = calculateVoltage(self.voltage, 1.0, distance, self.dielectric_constant, self.dielectric_thickness)
                    #eforce = - computeElectrostaticForce(realBias, distance)
                    
                    eforce = - computeElectrostaticForceWithCharge(self.voltage, 1.0, distance, 
                                        self.dielectric_constant, self.dielectric_thickness, 0)
           	    self.beamForce[c] += eforce
                    self.elecForceSum += eforce

        ### simplied model to calculate electrostatic force for tilted cantilever ###   
        if switch == 2:  
            startPos = 10e-6             # actuation pad starting location
            endPos   = startPos+200e-6    # actuation pad ending location, 200e-6 is dielectric length 
                  
            if self.enableElecModel:
                for c in range (0, self.nSolidCells):
                    x = self.xCells[c][0]
                    if x > startPos and x < endPos :
                        realgap = self.gap1 + (x/self.beam_length)*(self.gap2-self.gap1)
                        distance = realgap + self.deformation[c][2] - self.dielectric_thickness
                        eforce = - computeElectrostaticForceWithCharge(self.voltage, 1.0, distance, 
                                        self.dielectric_constant, self.dielectric_thickness,
                                                                       self.source)
                        self.beamForce[c] += eforce
                        self.elecForceSum += eforce
                   
        ### use compact damping model to compute damp force ###  
        dfMax = 0
        dfMin = 1e10  
        if self.enableFlowModel:
            for c in range (0, self.nSolidCells):
                x = self.xCells[c][0]
                realgap = self.gap1 + (x/self.beam_length)*(self.gap2-self.gap1)
                distance = realgap + self.deformation[c][2]- self.dielectric_thickness
                vel = self.velocity[c][2]
                width = self.beam_width
                length = self.beam_length
                Cf = computeDampCoeff(distance, width, self.beam_thickness)
                dforce = computeDampForce(vel, Cf)
                area = width
                dforce /= area
                self.beamForce[c] += dforce
                self.flowForceSum += dforce
                if dforce > dfMax:
                	dfMax = dforce
                if dforce < dfMin:
                	dfMin = dforce
        print 'damping force range %e %e' % (dfMin, dfMax)
        ### use new contact model to compute contact force ###    
        if self.enableContactModel:
            for c in range (0, self.nSolidCells):
                x = self.xCells[c][0]                  
                realgap = self.gap1 + (x/self.beam_length)*(self.gap2-self.gap1)
                distance = realgap + self.deformation[c][2]- self.dielectric_thickness
                cforce = computeContactForceUQ(distance, c1, c2, c3, H, dc, tR)
                self.beamForce[c] += cforce
                self.contactForceSum += cforce
       
         
        self.pmodel.advance(2)        
        self.dmodel.calculateNodeDisplacement()
        self.dmodel.deformPlate()
        self.solidMetricsCalculator.recalculate_deform()
        self.nodeCoordMeshA = self.solidMeshes[0].getNodeCoordinates()
        nodeCoordMesh = self.nodeCoordMeshA.asNumPyArray()
        self.nodeCoordA = self.geomFields.coordinate[self.solidMeshes[0].getNodes()]
        nodeCoord = self.nodeCoordA.asNumPyArray()
        nodeCoordMesh[:,:] = nodeCoord[:,:]

        self.dmodel.updateBoundaryMesh(self.solidMeshes[0], 
                                           self.solidBoundaryMeshes[0], 
                                           self.beam_thickness)
        self.solidBoundaryMetricsCalculator.recalculate_deform()  
              
        
        if self.transient == True:
            self.pmodel.updateTime()
            self.dmodel.updateTime()
          
            self.globalTime += self.timeStep
            
        self.globalCount += 1

        self.acc = findMax(self.acceleration)
        self.vel = findMax(self.velocity)
        
        if self.timeOption == True:
            mints = 1e-5
            dr = 0
            if self.cloestDistance > 10e-9:
            	dr = computeTravelDistance(self.cloestDistance, self.gap, self.minR, self.maxR)
            else:
                dr = computeTravelDistance(10e-9, self.gap, self.minR, self.maxR)
            
            for c in range(0, self.nSolidCells):
                acc = fabs(self.acceleration[c])
                vel = fabs(self.velocity[c][2])
                ts = computeTimeStep(dr, vel, acc, self.minR, self.maxR)
                
                if ts > 0 and ts < mints:
                    mints = ts
                if ts > 0 and ts < 1e-12:
                    mints = 1e-12
            if mints/self.timeStep < 10:
                self.timeStep = mints 
                self.poptions.setVar('timeStep',self.timeStep)    
        
	

        if tag == 'pullin' and self.cloestDistance < 0.2e-6:
            return False
        if tag == 'pullout' and self.farDistance > self.gap2*0.8:
            return False        	

            
    
