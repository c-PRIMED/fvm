import h5py

class Persistence(object):
    """ Provides a wrapper around a h5py.File object with methods to
    read and write persistent data for various models"""
    
    def __init__(self,fileName,mode):
        self.hdfFile = h5py.File(fileName,mode)

    def close(self):
        self.hdfFile.close()

    def saveAttribute(self, attrName, value):
        self.hdfFile.attrs[attrName] = value
        
    def readAttribute(self, attrName):
        return self.hdfFile.attrs[attrName]
        
    def writeField(self,field,site,siteLabel):
        if site not in field:
            return
        group = self.hdfFile.get(field.getName(),None)
        if group is None:
            group = self.hdfFile.create_group(field.getName())
        a = field[site].asNumPyArray()
        group.create_dataset(siteLabel,data=a)

    def readField(self,field,site,siteLabel):
        group = self.hdfFile.get(field.getName(),None)
        if group is not None:
            ds = group.get(siteLabel,None)
            if ds is not None:
                myArray = field[site].asNumPyArray()
                myArray[:] = ds.value
                return
        #raise IndexError('%s - %s not found' % (field.getName(), siteLabel))

    def saveMeshes(self,meshes):
        group = self.hdfFile.get('coordinates',None)
        if group is None:
            group = self.hdfFile.create_group('coordinates')
        for mesh in meshes:
            meshID = mesh.getID()
            nodesLabel = 'mesh-%d:nodes' % meshID

            a = mesh.getNodeCoordinates().asNumPyArray()
            group.create_dataset(nodesLabel,data=a)

    def readMeshes(self,meshes):
        group = self.hdfFile.get('coordinates',None)
        if group is not None:
            for mesh in meshes:
                meshID = mesh.getID()
                nodesLabel = 'mesh-%d:nodes' % meshID

                ds = group.get(nodesLabel,None)
                if ds is not None:
                    a = mesh.getNodeCoordinates().asNumPyArray()
                    a[:] = ds.value
                else:
                    raise IndexError(' coordinates %s not found' % nodesLabel)

    def readWriteModelFields(self, meshes, modelFieldData, op):
        for mesh in meshes:
            meshID = mesh.getID()
            cells = mesh.getCells()
            faces = mesh.getFaces()
            nodes = mesh.getNodes()
            cellsLabel = 'mesh-%d:cells' % meshID
            facesLabel = 'mesh-%d:faces' % meshID
            nodesLabel = 'mesh-%d:nodes' % meshID

            for field in modelFieldData.get('cells',[]):
                op(field,cells,cellsLabel)

            for field in modelFieldData.get('faces',[]):
                op(field,faces,facesLabel)

            for field in modelFieldData.get('nodes',[]):
                op(field,nodes,nodesLabel)

            fgs = mesh.getBoundaryGroups()
            for fg in fgs:
                fgLabel = 'mesh-%d:faceGroup-%d' % (meshID,fg.id)
                for field in modelFieldData.get('boundaryFaces',[]):
                    op(field,fg.site,fgLabel)
                    
    def saveModelFields(self, meshes, modelFieldData):
        self.readWriteModelFields(meshes,modelFieldData, self.writeField)

    def readModelFields(self, meshes, modelFieldData):
        self.readWriteModelFields(meshes,modelFieldData, self.readField)
        
    def saveModel(self, model, modelName, modelFieldData, meshes):

        self.saveModelFields(meshes, modelFieldData)

        model_group = self.hdfFile.create_group(modelName)
        x = model.getPersistenceData()
        for k in x.keys():
            model_group.create_dataset(k,data=x[k].asNumPyArray())
            
    def readModel(self,model, modelName, modelFieldData, meshes):

        self.readModelFields(meshes, modelFieldData)

        model_group = self.hdfFile.get(modelName,None)
        if model_group is not None:
            pData = model.getPersistenceData()
            for k,v in model_group.iteritems():
                a = pData[k].asNumPyArray()
                a[:] = v.value
        model.restart()

    #---------------------------------------------------------
    # Flow Model
    #---------------------------------------------------------
        
    def getFlowModelData(self,flowFields,fmodel):
        
        foptions = fmodel.getOptions()

        modelFieldData = {}

        cellFields = [flowFields.velocity, flowFields.pressure,
                      flowFields.density, flowFields.viscosity]
        faceFields = [flowFields.massFlux, flowFields.pressure]

        boundaryFields = [flowFields.momentumFlux]
        
        if foptions.transient:
            cellFields.append(flowFields.velocityN1)
            if foptions.timeDiscretizationOrder > 1:
                cellFields.append(flowFields.velocityN2)

        modelFieldData['cells'] = cellFields
        modelFieldData['faces'] = faceFields
        modelFieldData['boundaryFaces'] = boundaryFields

        return modelFieldData
    
    def saveFlowModel(self,flowFields,model,meshes):

        modelFieldData = self.getFlowModelData(flowFields,model)
        self.saveModel(model, 'fmodel', modelFieldData, meshes)

            
    def readFlowModel(self,flowFields,model,meshes):

        modelFieldData = self.getFlowModelData(flowFields,model)
        self.readModel(model, 'fmodel', modelFieldData, meshes)

    #---------------------------------------------------------
    # Structure Model
    #---------------------------------------------------------
        
    def getStructureModelData(self, structureFields, model):
        modelFieldData = {}
        options = model.getOptions()
        
        cellFields = [structureFields.deformation]
        if options.transient:
            cellFields.append(structureFields.deformationN1)
            cellFields.append(structureFields.deformationN2)

        modelFieldData['cells'] = cellFields
        
        return modelFieldData

    def saveStructureModel(self,structureFields,model,meshes):

        modelFieldData = self.getStructureModelData(structureFields,model)
        self.saveModel(model, 'smodel', modelFieldData, meshes)

            
    def readStructureModel(self,structureFields,model,meshes):

        modelFieldData = self.getStructureModelData(structureFields,model)
        self.readModel(model, 'smodel', modelFieldData, meshes)
