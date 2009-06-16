
import importers

import SchemeParser

def isConsPairOrList(v):
    return isinstance(v,SchemeParser.ConsPair) or isinstance(v,list)

def scmToPy(val):
    
    if isinstance(val,SchemeParser.ConsPair):
        if not isConsPairOrList(val[1]):        
            return [scmToPy(val[0]),scmToPy(val[1])]
        if len(val[1]) ==0:
            return scmToPy(val[0])
#        elif len(val[1][1]) == 0 and str(val[1][0])[0]=='e':
#            numStr = str(val[0]) + str(val[1][0])
#            return eval(numStr)
        else:
            l = scmToPyList(val)
            if l is not None:
                return l
            else:
                return val
    elif isinstance(val,SchemeParser.Symbol):
        if (str(val)=='#f'):
            return False
        elif (str(val)=='#t'):
            return True
        else:
            return str(val)
    else:
        return val

def scmToPyList(vars):
    
    varsList = []

    while len(vars) > 0 and isinstance(vars,SchemeParser.ConsPair):
        val = vars[0]
        if isinstance(val,SchemeParser.ConsPair):
            return None
        val = scmToPy(val)
        varsList.append(val)
        vars = vars[1]
    return varsList

def scmToDict(vars):
    
    varsDict = {}

    while len(vars) > 0:
        if not isinstance(vars[0],SchemeParser.ConsPair):
            varsDict[str(vars[0])] = scmToPy(vars[1])
            return varsDict
        if isinstance(vars[0][0],SchemeParser.ConsPair):
            return None
        key = str(vars[0][0])
        val = vars[0][1]
        val = scmToPy(val)
        varsDict[key] = val
        vars = vars[1]
    return varsDict


threadTypeToZoneType = { 1:'fluid',
                         2:'interior',
                         3:'wall',
                         4:'pressure-inlet',
                         5:'pressure-outlet',
                         7:'symmetric',
                         10:'velocity-inlet',
                         17:'solid'
                         }

def getZoneType(threadType):
    if threadType in threadTypeToZoneType:
        return threadTypeToZoneType[threadType]
    else:
        raise TypeError('invalid thread type %d' % threadType)


class FluentCase(importers.FluentReader):

    
    class FluentZone():

        def __init__(self,ID,zoneName,threadType,zoneType,varString):
            self.id = ID
            if zoneType == '':
                zoneType = getZoneType(threadType)
                
            if zoneName == '':
                zoneName = "%s_%d" % (zoneType,ID)
                
            self.zoneType = zoneType
            self.zoneName = zoneName
            if varString == "":
                self.varsDict = {}
            else:
                self.varsDict = scmToDict(SchemeParser.parse(varString))
                if 'sources?' in self.varsDict and self.varsDict['sources?']:
                    sourcekey = 'source-terms'
                    if sourcekey not in self.varsDict:
                        sourcekey = 'sources'
                    sourceVar = self.varsDict[sourcekey]
                    self.varsDict['source-terms'] = scmToDict(sourceVar)
                
        def getVar(self,v):
            return self.varsDict[v]

        def getConstantVar(self,v):
            val = self.varsDict[v]
            if not isinstance(val,list):
                return val
            if isinstance(val[0],list):
                val = val[0]
            if val[0] == 'constant':
                return val[1]
            else:
                raise ValueError(v + ' value is not constant')

        def getConstantSource(self,v):
            val = self.varsDict['source-terms'][v]
            if isinstance(val[0],list):
                val = val[0]
            if val[0] == 'constant':
                return val[1]
            else:
                raise ValueError('source ' + v + ' value is not constant')
    


    
    class Material:

        def __init__(self,materialType,props):
            self.materialType = materialType
            self.props = props

        def getPropMethod(self,p):
            val = self.props[p]
            if isinstance(val[0],list):
                val=val[0]
            return val[0]
            
        def getPropData(self,p):
            val = self.props[p]
            if isinstance(val[0],list):
                val=val[0]
            return val[1]
            
        def getConstantProp(self,p):
            val = self.props[p]
            if isinstance(val[0],list):
                val=val[0]
            if val[0] == 'constant':
                return val[1]
            else:
                raise ValueError(' %s value is not constant: %s' % (p,val) )
            
        
        def getProp(self,p,tField):
            """ allows for constant or polynomial of tField"""
            val = self.props[p]
            if isinstance(val[0],list):
                val=val[0]
            if val[0] == 'constant':
                return val[1]
            elif val[0] == 'polynomial':
                return Polynomial(scmToPy(val[1]),field=tField)
                #return PyArrayUDFPolynomial(coeffs=scmToPy(val[1]),xField=tField)
            else:
                raise ValueError(' %s value is not constant or polynomial: %s' % (p,val) )

    def getVar(self,v):
        return self.varsDict[v]

    def read(self):
        print 'reading mesh'
        self.readMesh()

        vars=SchemeParser.parse(self.getVars())
        self.varsDict = scmToDict(vars)
        if self.varsDict is None:
            raise TypeError("vars is not an association list")

        self.config = scmToDict(self.getVar('case-config'))
        self.materials = {}
        self.faceZones = {}
        self.cellZones = {}
        
        mDict = scmToDict(self.getVar('materials'))
        for n,d in mDict.iteritems():
            self.materials[n] = FluentCase.Material(materialType=scmToPy(d[0]),
                                                    props=scmToDict(d[1]))

        rFaceZones = self.getFaceZones()
        for i in rFaceZones.keys():
            rfz = rFaceZones[i]
            self.faceZones[i] = FluentCase.FluentZone(ID=rfz.ID,
                                                       zoneName=rfz.zoneName,
                                                       threadType=rfz.threadType,
                                                       zoneType=rfz.zoneType,
                                                       varString=rfz.zoneVars)
        rCellZones = self.getCellZones()
        for i in rCellZones.keys():
            rfz = rCellZones[i]
            self.cellZones[i] = FluentCase.FluentZone(ID=rfz.ID,
                                                      zoneName=rfz.zoneName,
                                                      threadType=rfz.threadType,
                                                      zoneType=rfz.zoneType,
                                                      varString=rfz.zoneVars)


    def importThermalBCs(self,tmodel):
        bcMap = tmodel.getBCMap()
        for i in bcMap.keys():
            bc = bcMap[i]
            fluentZone = self.faceZones[i]
            if fluentZone.zoneType == 'wall':
                thermalBCType = fluentZone.getVar('thermal-bc')
                if thermalBCType == 0:
                    bc.bcType = 'SpecifiedTemperature'
                    bc.setVar('specifiedTemperature',fluentZone.getConstantVar('t'))
                elif thermalBCType == 1:
                    flux= fluentZone.getConstantVar('q')
                    bc.bcType = 'SpecifiedHeatFlux'
                    bc.setVar('specifiedHeatFlux',flux)
                elif thermalBCType == 3:
                    bc.bcType = 'CoupledWall'
                else:
                    raise TypeError('thermal BCType %d not handled' % thermalBCType)
            elif fluentZone.zoneType in ['velocity-inlet','pressure-inlet',
                                         'pressure-outlet',
                                         'mass-flow-inlet',
                                         'exhaust-fan', 'intake-fan',
                                         'inlet-vent', 'outlet-vent']:
                bc.bcType = 'SpecifiedTemperature'
                if fluentZone.zoneType == 'velocity-inlet':
                    bc.setVar('specifiedTemperature', fluentZone.getConstantVar('t'))
                else:
                    bc.setVar('specifiedTemperature',fluentZone.getConstantVar('t0'))
            elif fluentZone.zoneType == 'symmetry':
                pass
            else:
                raise TypeError('invalid boundary type : ' + fluentZone.zoneType)

    def importFlowBCs(self,fmodel):
        options = fmodel.getOptions()
        options.setVar('initialXVelocity', self.getVar('x-velocity/default'))
        options.setVar('initialYVelocity', self.getVar('y-velocity/default'))
        options.setVar('initialZVelocity', self.getVar('z-velocity/default'))
        options.setVar('initialPressure', self.getVar('pressure/default'))
        options.setVar('momentumURF', self.getVar('mom/relax'))
        options.setVar('pressureURF', self.getVar('pressure/relax'))
                       
        bcMap = fmodel.getBCMap()
        for i in bcMap.keys():
            bc = bcMap[i]
            fluentZone = self.faceZones[i]
            if fluentZone.zoneType == 'wall':
                motionBCType = fluentZone.getVar('motion-bc')
                bc.bcType = 'NoSlipWall'
                if motionBCType == 0:
                    pass
                elif motionBCType == 1:
                    vmag = fluentZone.getVar('vmag')
                    bc.setVar('specifiedXVelocity',vmag*fluentZone.getConstantVar('ni'))
                    bc.setVar('specifiedYVelocity',vmag*fluentZone.getConstantVar('nj'))
                    bc.setVar('specifiedZVelocity',vmag*fluentZone.getConstantVar('nk'))
                else:
                    raise TypeError('flow BCType %d not handled' % motionBCType)
            elif fluentZone.zoneType == 'velocity-inlet':
                motionBCType = fluentZone.getVar('velocity-spec')
                bc.bcType = 'VelocityBoundary'
                if motionBCType == 0:
                    vmag = fluentZone.getVar('vmag')
                    bc.setVar('specifiedXVelocity',vmag*fluentZone.getConstantVar('ni'))
                    bc.setVar('specifiedYVelocity',vmag*fluentZone.getConstantVar('nj'))
                    bc.setVar('specifiedZVelocity',vmag*fluentZone.getConstantVar('nk'))
                if motionBCType == 1:
                    vmag = fluentZone.getVar('vmag')
                    bc.setVar('specifiedXVelocity',fluentZone.getConstantVar('u'))
                    bc.setVar('specifiedYVelocity',fluentZone.getConstantVar('v'))
                    bc.setVar('specifiedZVelocity',fluentZone.getConstantVar('w'))
                else:
                    raise TypeError('flow BCType %d not handled' % motionBCType)
            elif fluentZone.zoneType == 'pressure-outlet':
                bc.bcType = 'PressureBoundary'
                bc.setVar('specifiedPressure',fluentZone.getConstantVar('p'))
            elif fluentZone.zoneType == 'pressure-inlet':
                bc.bcType = 'PressureBoundary'
                bc.setVar('specifiedPressure',fluentZone.getConstantVar('p0'))
            elif fluentZone.zoneType == 'symmetry':
                bc.bcType = 'Symmetry'
            else:
                raise TypeError('invalid boundary type : ' + fluentZone.zoneType)

        vcMap = fmodel.getVCMap()
        for i in vcMap.keys():
            fluentZone = self.cellZones[i]
            material = self.materials[fluentZone.getVar('material')]
            vc = vcMap[i]
            if material.getPropMethod('density') == 'constant':
                vc.setVar('density',material.getConstantProp('density'))
            else:
                print 'Density method is not constant. Remember to setup ideal gas density model'
            vc.setVar('viscosity',material.getConstantProp('viscosity'))
                      
                
