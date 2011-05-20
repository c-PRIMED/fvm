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
        self.fluidMeshesNew = meshes.fluidMeshesNew
        self.geomFields = meshes.geomFields
        self.gap = meshes.gap
        self.beam_thickness = meshes.beam_thickness
        self.dielectric_thickness = meshes.dielectric_thickness
        print "dielectric thickness in model %e" % self.dielectric_thickness
        print "there are %i meshes" % len(self.fluidMeshesNew)
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
        if self.enableContactModel == True:
            self.contactFields =  models.ContactFields('contact')
            self.cmodel = models.ContactModelA(self.geomFields,self.contactFields,self.fluidMeshesNew)
    
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
        elif tag == 'elec' or tag == 'dielectric':
            for mesh in self.fluidMeshes:
                for fg in mesh.getAllFaceGroups():
                    for id in fgs1:
                        if fg.id == id:
                            fg.groupType = 'dielectric interface'
        else:
            print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
            print 'Error: wrong tag name in boundary condition: [ %s ]' % tag
            print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'

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
                vc['dielectric_constant'] = epsilon                        
            elec_constants = self.emodel.getConstants()
            elec_constants['dielectric_thickness'] = self.dielectric_thickness

        if self.enableFlowModel == True:
            vcMap = self.fmodel.getVCMap()
            for i,vc in vcMap.iteritems():
                vc['viscosity'] = viscosity
                vc['density'] = rhoFluid
	if self.enableContactModel == True:
	    constants = self.cmodel.getConstants()
	    constants['gap'] = self.gap
	    constants['thickness'] = self.beam_thickness

    def solidProperties(self, rho, E, nu):
        if self.enablePlateModel == True:
            vcMap = self.pmodel.getVCMap()
            for i,vc in vcMap.iteritems():
                vc['density'] = rho
                vc['ym'] = E
                vc['nu'] = nu
        if self.enableContactModel == True:
	    constants = self.cmodel.getConstants()
	    constants['gap'] = self.gap
	    constants['thickness'] = self.beam_thickness       

    def fluidProperties(self, rho, viscosity):
        if self.enableFlowModel == True:
            vcMap = self.fmodel.getVCMap()
            for i,vc in vcMap.iteritems():
                vc['viscosity'] = viscosity
                vc['density'] = rho

    def elecProperties(self, dielectricID, airID, substrateID, dc_SiN, dc_Air, dc_Si):
       if self.enableElecModel == True:
            vcMap = self.emodel.getVCMap()
            for mesh in self.fluidMeshesNew:
                vc = vcMap[mesh.getID()]
                i = mesh.getCellZoneID()
                if i == dielectricID:
                    vc.vcType = "dielectric"
                    vc['dielectric_constant'] = dc_SiN
                elif i == airID:
                    vc.vcType = "air"
                    vc['dielectric_constant'] = dc_Air
                elif i == substrateID:
                    vc.vcType = "substrate"
                    vc['dielectric_constant'] = dc_Si
                else:
                    print "unknow cell zone id"
                print "mesh %i has type %s" % (i, vc.vcType)
                        
            elec_constants = self.emodel.getConstants()
            elec_constants['dielectric_thickness'] = self.dielectric_thickness

 
    def solvers(self, transient, timeStep=1e-8):
        if self.enablePlateModel == True:
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
            eoptions.electrostaticsTolerance = 1e-3
            eoptions.electrostatics_enable = 1
            eoptions.chargetransport_enable = 0
            eoptions.tunneling = 0
            eoptions.ibm_enable = 1
            eoptions.transient_enable = False
            eoptions.printNormalizedResiduals = True
        if self.enableFlowModel == True:  
            foptions = self.fmodel.getOptions()
            foptions.momentumTolerance=1e-1
            foptions.continuityTolerance=1e-1
            foptions.setVar("momentumURF",0.7)
            foptions.setVar("pressureURF",0.3)
            foptions.transient=transient
            foptions.setVar("timeStep",timeStep)
            foptions.printNormalizedResiduals=True
            
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
