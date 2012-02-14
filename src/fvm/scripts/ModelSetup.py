import fvm.fvmbaseExt as fvmbaseExt
import fvm.models_atyped_double as models
from Persistence import Persistence

class ModelSetup:

    def __init__(self, meshes, 
                 enablePlateModel=True, enableElecModel=True, 
                 enableFlowModel=True, enableContactModel=True,
                 enableDielectricCharging=False,
                 enableCreep=False ):
        self.enablePlateModel = enablePlateModel
        self.enableElecModel = enableElecModel
        self.enableFlowModel = enableFlowModel
        self.enableContactModel = enableContactModel
        self.enableDielectricCharging = enableDielectricCharging
        self.enableCreep =  enableCreep
        self.solidMeshes = meshes.solidMeshes
        self.fluidMeshes = meshes.fluidMeshes
        self.fluidMeshesNew = meshes.fluidMeshesNew
        self.geomFields = meshes.geomFields
        self.gap = meshes.gap
        self.beam_thickness = meshes.beam_thickness
        
        self.dielectric_thickness = meshes.dielectric_thickness
        
    def createModels(self):
        if self.enablePlateModel == True:
            self.plateFields =  models.PlateFields('plate')
            self.pmodel = models.PlateModelA(self.geomFields,self.plateFields,self.solidMeshes)
            self.dmodel = models.PlateDeformationModelA(self.geomFields,self.plateFields,self.solidMeshes)
        if self.enableElecModel == True:
            self.elecFields =  models.ElectricFields('elec')
            self.emodel = models.ElectricModelA(self.geomFields,self.elecFields,self.fluidMeshesNew)
        if self.enableFlowModel == True:
            self.flowFields =  models.FlowFields('flow')
            self.fmodel = models.FlowModelA(self.geomFields,self.flowFields,self.fluidMeshes)
        #currently contact model is implemented in ComputeForce.py    
        #if self.enableContactModel == True:
        #    self.contactFields =  models.ContactFields('contact')
        #    self.cmodel = models.ContactModelA(self.geomFields,self.contactFields,self.fluidMeshesNew)
        print 'models are created'
        
    def boundaryCondition(self, tag, fgs1=[], fgs2=[], value=0):
        if tag == 'beam' or tag == 'solid':
            if self.enablePlateModel == True:
                bcMap = self.pmodel.getBCMap()
                for id in fgs1:         # fix 
                    bc = bcMap[id]
                    bc.bcType = 'Clamped'
                    bc['specifiedXRotation']=0
                    bc['specifiedYRotation']=0.
                    bc['specifiedZDeformation']=0.                       
                for id in fgs2:        # force 
                    bc = bcMap[id] 
                    bc.bcType = 'SpecifiedTraction'
            print 'solid boundary condition applied'
            
        elif tag == 'fluid':
            if self.enableElecModel == True:
                bcMap = self.emodel.getBCMap()
                for id in fgs1:         
                    bc = bcMap[id]
                    bc.bcType = "Symmetry"
                for id in fgs2:
                    bc = bcMap[id]
                    if self.enableDielectricCharging == True:
                    	bc.bcType = "SpecialDielectricBoundary"
                    else:
                        bc.bcType = "SpecifiedPotential"
                    bc['specifiedPotential'] = value
            if self.enableFlowModel == True:
                bcMap = self.fmodel.getBCMap()
                for id in fgs1:
                    bc = bcMap[id]
                    bc.bcType = 'NoSlipWall'
                for id in fgs2:
                    bc = bcMap[id]
                    bc.bcType = 'NoSlipWall'
            print 'fluid background boundary condition applied'
            
        else:
            print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
            print 'Error: wrong tag name in boundary condition: [ %s ]' % tag
            print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    """
    def materialProperties(self, rhoSolid, E, nu, epsilon, rhoFluid, viscosity):
        if self.enablePlateModel == True:
            vcMap = self.pmodel.getVCMap()
            for i,vc in vcMap.iteritems():
                vc['density'] = rhoSolid
                vc['ym'] = E
                vc['nu'] = nu
                
        if self.enableElecModel == True:
            vcMap = self.emodel.getVCMap()
            for i,vc in vcMap.iteritems():
                vc['dielectric_constant'] = 1.0                        
            elec_constants = self.emodel.getConstants()
            elec_constants['dielectric_thickness'] = dielectric_thickness
	    elec_constants['dielectric_constant'] = epsilon
	    
        if self.enableFlowModel == True:
            vcMap = self.fmodel.getVCMap()
            for i,vc in vcMap.iteritems():
                vc['viscosity'] = viscosity
                vc['density'] = rhoFluid
	#if self.enableContactModel == True:
	#    constants = self.cmodel.getConstants()
	#    constants['gap'] = self.gap
	#    constants['thickness'] = self.beam_thickness
    """
    def solidProperties(self, rho, E, nu):
        if self.enablePlateModel == True:
            vcMap = self.pmodel.getVCMap()
            for i,vc in vcMap.iteritems():
                vc['density'] = rho
                vc['ym'] = E
                vc['nu'] = nu              
        print 'solid properties applied'
        
    def fluidProperties(self, rho, viscosity):
        if self.enableFlowModel == True:
            vcMap = self.fmodel.getVCMap()
            for i,vc in vcMap.iteritems():
                vc['viscosity'] = viscosity
                vc['density'] = rho            
        print 'fluid  properties applied' 
           	
    def elecProperties(self, dc, thickness): 
        if self.enableElecModel == True:
            vcMap = self.emodel.getVCMap()
            for i,vc in vcMap.iteritems():
    		vc['dielectric_constant'] = 1.0 
                                   
            elec_constants = self.emodel.getConstants()
            elec_constants['dielectric_thickness'] = thickness
	    elec_constants['dielectric_constant'] = dc
        print 'electric properties applied'
        
    def solvers(self, timeStep=1e-8):
        if self.enablePlateModel == True:
            pc = fvmbaseExt.AMG()
            pc.verbosity=0
            defSolver = fvmbaseExt.DirectSolver()
            defSolver.preconditioner = pc
            defSolver.relativeTolerance = 1e-10
            defSolver.absoluteTolerance = 1.e-20
            defSolver.nMaxIterations = 50000
            defSolver.verbosity=0
            poptions = self.pmodel.getOptions()
            poptions.deformationLinearSolver = defSolver
            poptions.deformationTolerance=1.0e-3
            poptions.setVar("deformationURF",1.0)
            poptions.printNormalizedResiduals=True
            poptions.timeDiscretizationOrder = 1
            poptions.transient=False
            poptions.scf = 5./6.
            poptions.setVar('timeStep',timeStep)
            poptions.variableTimeStep = True            
            poptions.residualStress = False
            poptions.creep  = False
          	                 
        if self.enableElecModel == True:
            epc = fvmbaseExt.AMG()
            epc.verbosity=0
            elecSolver = fvmbaseExt.BCGStab()
            elecSolver.preconditioner = epc
            elecSolver.relativeTolerance = 1e-4
            elecSolver.absoluteTolerance = 1e-20
            elecSolver.nMaxIterations = 1000
            elecSolver.maxCoarseLevels=20
            elecSolver.verbosity=0
            eoptions = self.emodel.getOptions()
            eoptions.electrostaticsLinearSolver = elecSolver
            eoptions.electrostaticsTolerance = 1e-3
            eoptions.electrostatics_enable = 1
            eoptions.chargetransport_enable = 0
            eoptions.tunneling = 0
            eoptions.ibm_enable = 1
            eoptions.transient_enable = False
            eoptions.printNormalizedResiduals = True
            
        if self.enableFlowModel == True:  
            foptions = self.fmodel.getOptions()
            foptions.momentumTolerance=1e-3
            foptions.continuityTolerance=1e-3
            foptions.setVar("momentumURF",0.7)
            foptions.setVar("pressureURF",0.3)
            foptions.transient=False
            foptions.setVar("timeStep",timeStep)
            foptions.printNormalizedResiduals=True
        print 'solvers are established'
    
    def options(self, tag, boolValue, value):
    	if self.enablePlateModel == True:
    	    if tag == 'creep':
    	        poptions = self.pmodel.getOptions()
    	    	poptions.creep = boolValue
            	poptions.A = value/3600	
 		poptions.B = 20
 		poptions.m = 1
 		poptions.n = 0.5
 		poptions.Sy0 = 0.7e9
 	    if tag == 'residual stress':
 	    	poptions = self.pmodel.getOptions()
                poptions.residualStress = boolValue
                poptions.setVar('residualStressXX',value) 
            if tag == 'transient':
            	poptions = self.pmodel.getOptions()
            	poptions.transient = boolValue
            	poptions.setVar('timeStep',value)
            if tag == 'variable timestep':
            	poptions = self.pmodel.getOptions()
            	poptions.variableTimeStep = boolValue
            	poptions.setVar('timeStep',value)
            	
        if self.enableFlowModel == True:
            if tag == 'transient':
            	foptions = self.fmodel.getOptions()
            	foptions.transient = boolValue
            	foptions.setVar('timeStep',value)	
            		
 			        
    def initModels(self):
        if self.enablePlateModel == True:
            self.pmodel.init()
            self.dmodel.init()
        if self.enableElecModel == True:    
            self.emodel.init()
        if self.enableFlowModel == True:  
            self.fmodel.init()
        #if self.enableContactModel == True:  
        #    self.cmodel.init() 
        print 'model init done'
        
    def restart(self, restartFile):
        if self.enableElecModel == True: 
            restartFile.readElectricModel(self.elecFields, self.emodel, self.fluidMeshes)
        if self.enablePlateModel == True:
            restartFile.readPlateModel(self.plateFields,self.pmodel,self.solidMeshes)
