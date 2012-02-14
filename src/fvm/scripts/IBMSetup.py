import fvm.fvmbaseExt as fvmbaseExt

class IBM():

    def __init__(self, meshes, models):
        self.geomFields = models.geomFields
        self.enableElecModel = models.enableElecModel
        self.enableFlowModel = models.enableFlowModel
        if self.enableElecModel == True:
            self.elecFields = models.elecFields
        else:
            print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
            print 'elecModel is OFF; Canot apply IBM on potential field'
            print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        if self.enableFlowModel == True:
            self.flowFields = models.flowFields
        else:
            print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
            print 'flowModel is OFF; Canot apply IBM on flow field'
            print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        self.solidBoundaryMeshes = meshes.solidBoundaryMeshes
        self.fluidMeshes = meshes.fluidMeshesNew
        self.sbMeshFaces = self.solidBoundaryMeshes[0].getFaces()       
        self.fluidMetricsCalculator = meshes.fluidMetricsCalculator
        self.ibManager = fvmbaseExt.IBManager(self.geomFields,
                                 self.solidBoundaryMeshes[0],
                                 self.fluidMeshes)
       
    def init(self, voltage=0, vel=0):
        faceCount = self.sbMeshFaces.getCount()
        if self.enableElecModel == True:
            areaMag = self.geomFields.areaMag[self.sbMeshFaces]
            pot = areaMag.newSizedClone(faceCount)
            pota = pot.asNumPyArray()
            pota[:] = voltage
            self.elecFields.potential[self.sbMeshFaces] = pot
        if self.enableFlowModel == True:
            area = self.geomFields.area[self.sbMeshFaces]
            velocity = area.newSizedClone(faceCount)
            velocitya = velocity.asNumPyArray()
            velocitya[:,:] = 0.
            self.flowFields.velocity[self.sbMeshFaces] = velocity

    def set(self, a, b, c):
        self.ibManager.fluidNeighborsPerIBFace = a
        self.ibManager.solidNeighborsPerIBFace = b
        self.ibManager.fluidNeighborsPerSolidFace = c

    def update(self):
        self.ibManager.update()
        self.fluidMetricsCalculator.computeIBInterpolationMatrices(self.sbMeshFaces)
        
        self.fluidMetricsCalculator.computeSolidInterpolationMatrices(self.sbMeshFaces)  
       
