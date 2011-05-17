import fvm.fvmbaseExt as fvmbaseExt
import fvm.models_atyped_double as models
from FluentCase import FluentCase


	
class MeshSetup:
    
    def __init__(self, beam, fluid, beam_thickness, gap, dielectric_thickness=0):
        self.beamCaseFile = beam
        self.fluidCaseFile = fluid
        self.beam_thickness = beam_thickness
        self.dielectric_thickness = dielectric_thickness
        self.gap = gap
        self.probeIndex = 100
        
    def read(self):
        beamReader = FluentCase(self.beamCaseFile)
        beamReader.read()
        print 'beam mesh read in. Done!'
        print '------------------------------------------------------------'
        self.solidMeshes = beamReader.getMeshList()
        self.geomFields =  models.GeomFields('geom')
        self.solidMetricsCalculator = models.MeshMetricsCalculatorA(self.geomFields,self.solidMeshes)
        self.solidMetricsCalculator.init()

        fluidReader = FluentCase(self.fluidCaseFile)
        fluidReader.read()
        print 'fluid mesh read in. Done!'
        print '------------------------------------------------------------'
        self.fluidMeshes = fluidReader.getMeshList()
        self.fluidMetricsCalculator = models.MeshMetricsCalculatorA(self.geomFields,self.fluidMeshes)
        self.fluidMetricsCalculator.init()

        self.solidBoundaryMeshes = [m.extrude(1, self.beam_thickness, True) for m in self.solidMeshes]
        self.solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(self.geomFields,self.solidBoundaryMeshes)
        self.solidBoundaryMetricsCalculator.init()

    def translate(self, tag, dx=0, dy=0, dz=0):
        if tag == 'solid' or tag == 'beam':
            for mesh in self.solidMeshes:
                nodes = mesh.getNodes()
                xNodes = mesh.getNodeCoordinates().asNumPyArray()
                xNodes[:,0] += dx
                xNodes[:,1] += dy
                xNodes[:,2] += dz
            self.solidMetricsCalculator.init()
            self.solidBoundaryMeshes = [m.extrude(1, self.beam_thickness, True) for m in self.solidMeshes]
            self.solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(self.geomFields,self.solidBoundaryMeshes)
            self.solidBoundaryMetricsCalculator.init()
            print 'translate beam mesh. Done'
            print '-------------------------------------------------------------'
        elif tag == 'fluid':
            for mesh in self.fluidMeshes:
                nodes = mesh.getNodes()
                xNodes = mesh.getNodeCoordinates().asNumPyArray()
                xNodes[:,0] += dx
                xNodes[:,1] += dy
                xNodes[:,2] += dz
            self.fluidMetricsCalculator.init()
            print 'translate fluid mesh. Done'
            print '-------------------------------------------------------------'
        else:
            print 'Error: wrong tag name in translate mesh'

    def scale(self, tag, dx=1, dy=1, dz=1):
        if tag == 'solid' or tag == 'beam':
            for mesh in self.solidMeshes:
                nodes = mesh.getNodes()
                xNodes = mesh.getNodeCoordinates().asNumPyArray()
                xNodes[:,0] *= dx
                xNodes[:,1] *= dy
                xNodes[:,2] *= dz
            self.solidMetricsCalculator.init()
            self.solidBoundaryMeshes = [m.extrude(1, self.beam_thickness, True) for m in self.solidMeshes]
            self.solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(self.geomFields,self.solidBoundaryMeshes)
            self.solidBoundaryMetricsCalculator.init()
            print 'scale beam mesh. Done'
            print '-------------------------------------------------------------'
        elif tag == 'fluid':
            for mesh in self.fluidMeshes:
                nodes = mesh.getNodes()
                xNodes = mesh.getNodeCoordinates().asNumPyArray()
                xNodes[:,0] *= dx
                xNodes[:,1] *= dy
                xNodes[:,2] *= dz
            self.fluidMetricsCalculator.init()
            print 'scale fluid mesh. Done'
            print '-------------------------------------------------------------'
        else:
            print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
            print 'Error: wrong tag name in translate mesh: [ %s ]' % tag
            print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    def summary(self):
        nSelfCells = 0
        nCells = 0
        for mesh in self.solidMeshes:
            cells = mesh.getCells()
            nSelfCells += cells.getSelfCount()
            nCells += cells.getCount()
            print '-------------------------------------------------------------'
            print 'solid mesh: number of local cells %i' % nSelfCells 
            print 'solid mesh: number of total cells %i' % nCells 
            rCells = self.geomFields.coordinate[cells].asNumPyArray()
            rMin = rCells.min(axis=0)
            rMax = rCells.max(axis=0)
            print 'solid mesh x range [ %e , %e ]' % (rMin[0], rMax[0])
            print 'solid mesh y range [ %e , %e ]' % (rMin[1], rMax[1])
            print 'solid mesh z range [ %e , %e ]' % (rMin[2], rMax[2])
            if rMin[2] != 0.0 or rMax[2] != 0.0:
                print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                print 'Error: beam mesh is not centered at zero Z direction!'
                print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
            print '-------------------------------------------------------------'
        nSelfCells = 0
        nCells = 0
        for mesh in self.fluidMeshes:
            cells = mesh.getCells()
            nSelfCells += cells.getSelfCount()
            nCells += cells.getCount()
            print '-------------------------------------------------------------'
            print 'fluid mesh: number of local cells %i' % nSelfCells 
            print 'fluid mesh: number of total cells %i' % nCells 
            rCells = self.geomFields.coordinate[cells].asNumPyArray()
            rMin = rCells.min(axis=0)
            rMax = rCells.max(axis=0)
            print 'fluid mesh x range [ %e , %e ]' % (rMin[0], rMax[0])
            print 'fluid mesh y range [ %e , %e ]' % (rMin[1], rMax[1])
            print 'fluid mesh z range [ %e , %e ]' % (rMin[2], rMax[2])
            print '--------------------------------------------------------------'
        nFaces = 0
        for mesh in self.solidBoundaryMeshes:
            faces = mesh.getFaces()
            nFaces += faces.getCount()
            print '--------------------------------------------------------------'
            print 'solid boundary mesh: number of faces %i' % nFaces
            rCells = self.geomFields.coordinate[faces].asNumPyArray()
            rMin = rCells.min(axis=0)
            rMax = rCells.max(axis=0)
            print 'solid boundary mesh x range [ %e , %e ]' % (rMin[0], rMax[0])
            print 'solid boundary mesh y range [ %e , %e ]' % (rMin[1], rMax[1])
            print 'solid boudnary mesh z range [ %e , %e ]' % (rMin[2], rMax[2])
            print '--------------------------------------------------------------'

    def searchPoint(self, mesh, probe):
    	cells = mesh.getCells()
        xCells = self.geomFields.coordinate[cells]
        xCellsA = xCells.asNumPyArray()
    	cloestPoints = fvmbaseExt.newIntArray(1)
        cloestPointsA = cloestPoints.asNumPyArray()
    	target = fvmbaseExt.VecD3()
        target[0] = probe[0]
        target[1] = probe[1]
        target[2] = probe[2]
    	search = fvmbaseExt.KSearchTree(xCells)
        search.findNeighbors(target, 1, cloestPoints)
        point = cloestPointsA[0]
       	return point
	
    def findPoint(self, tag, probe):
        result = 0
        if tag == 'fluid':
            count = 0
            for mesh in self.fluidMeshes:
                count += 1
                result = self.searchPoint(mesh, probe)
            if count > 1:
                print 'multi fluid meshes search not implemented yet! '
        if tag == 'beam' or tag == 'solid':
            count = 0
            for mesh in self.solidMeshes:
                count += 1
                result = self.searchPoint(mesh, probe)
            if count > 1:
                print 'multi solid meshes search not implemented yet! '
        return result

    
    def findLine():
        print 'not implemented yet'
