import fvm.fvmbaseExt as fvmbaseExt
import fvm.models_atyped_double as models


class ModelSetup:

    def __init__(self, meshes, 
                 enablePlateModel=True, enableElecModel=True, enableFlowModel=True, enableContactModel=True ):
        self.enablePlateModel = enablePlateModel
        self.enableElecModel = enableElecModel
        self.enableFlowModel = enableFlowModel
        self.enableContactModel = enableContactModel
        self.solidMeshes = meshes.solidMeshes
        self.fluidMeshes = meshes.fluidMeshes
        self.geomFields = meshes.geomFields

    def createModels(self):
        if self.enablePlateModel == True:
            self.plateFields =  models.PlateFields('plate')
            self.pmodel = models.PlateModelA(self.geomFields,self.plateFields,self.solidMeshes)
            self.dmodel = models.PlateDeformationModelA(self.geomFields,self.plateFields,self.solidMeshes)
        if self.enableElecModel == True:
            self.elecFields =  models.ElectricFields('elec')
            self.emodel = models.ElectricModelA(self.geomFields,self.elecFields,self.fluidMeshes)
        if self.enableFlowModel == True:
            self.flowFields =  models.FlowFields('flow')
            self.fmodel = models.FlowModelA(self.geomFields,self.flowFields,self.fluidMeshes)
        if self.enableContactModel == True:
            self.contactFields =  models.ContactFields('contact')
            self.cmodel = models.ContactModelA(self.geomFields,self.contactFields,self.fluidMeshes)
    
    def boundaryCondition(self, tag, fgs1, fgs2, value=0):
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

        elif tag == 'fluid':
            if self.enableElecModel == True:
                bcMap = self.emodel.getBCMap()
                for id in fgs1:         
                    bc = bcMap[id]
                    bc.bcType = "Symmetry"
                for id in fgs2:
                    bc = bcMap[id]
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
        else:
            print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
            print 'Error: wrong tag name in boundary condition: [ %s ]' % tag
            print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'

    def materialProperties(self, rhoSolid, E, nu, epsilon, rhoFluid, viscosity):
        if self.enableElecModel == True:
            vcMap = self.pmodel.getVCMap()
            for i,vc in vcMap.iteritems():
                vc['density'] = rhoSolid
                vc['ym'] = E
                vc['nu'] = nu
        if self.enableElecModel == True:
            vcMap = self.emodel.getVCMap()
            for i,vc in vcMap.iteritems():
                vc.vcType = "dielectric"
                vc['dielectric_constant'] = epsilon
        if self.enableFlowModel == True:
            vcMap = self.fmodel.getVCMap()
            for i,vc in vcMap.iteritems():
                vc['viscosity'] = viscosity
                vc['density'] = rhoFluid

    def solvers(self, transient, timeStep=1e-8):
        if self.enableElecModel == True:
            pc = fvmbaseExt.AMG()
            pc.verbosity=0
            defSolver = fvmbaseExt.DirectSolver()
            defSolver.preconditioner = pc
            defSolver.relativeTolerance = 1e-9
            defSolver.absoluteTolerance = 1.e-30
            defSolver.nMaxIterations = 50000
            defSolver.verbosity=0
            poptions = self.pmodel.getOptions()
            poptions.deformationLinearSolver = defSolver
            poptions.deformationTolerance=1.0e-3
            poptions.setVar("deformationURF",1.0)
            poptions.printNormalizedResiduals=True
            poptions.timeDiscretizationOrder = 2
            poptions.transient=transient
            poptions.scf = 5./6.
            poptions.setVar('timeStep',timeStep)
            poptions.variableTimeStep = True
        if self.enableElecModel == True:
            epc = fvmbaseExt.AMG()
            epc.verbosity=0
            elecSolver = fvmbaseExt.BCGStab()
            elecSolver.preconditioner = epc
            elecSolver.relativeTolerance = 1e-4
            elecSolver.absoluteTolerance = 1e-10
            elecSolver.nMaxIterations = 1000
            elecSolver.maxCoarseLevels=20
            elecSolver.verbosity=0
            eoptions = self.emodel.getOptions()
            eoptions.electrostaticsLinearSolver = elecSolver
            eoptions.electrostaticsTolerance = 1e-8
            eoptions.electrostatics_enable = 1
            eoptions.chargetransport_enable = 0
            eoptions.tunneling = 0
            eoptions.ibm_enable = 1
            eoptions.transient_enable = False
            eoptions.printNormalizedResiduals = True
        if self.enableFlowModel == True:  
            foptions = self.fmodel.getOptions()
            foptions.momentumTolerance=1e-4
            foptions.continuityTolerance=1e-4
            foptions.setVar("momentumURF",0.7)
            foptions.setVar("pressureURF",0.3)
            foptions.transient=transient
            foptions.setVar("timeStep",timeStep)
            foptions.printNormalizedResiduals=False
            
    def initModels(self):
        if self.enableElecModel == True:
            self.pmodel.init()
            self.dmodel.init()
        if self.enableElecModel == True:    
            self.emodel.init()
        if self.enableFlowModel == True:  
            self.fmodel.init()
        if self.enableContactModel == True:  
            self.cmodel.init() 
