import fvm.fvmbaseExt as fvmbaseExt
import fvm.models_atyped_double as models
from FluentCase import FluentCase
from Persistence import Persistence
from math import *
	
class MeshSetup:
    
    def __init__(self, beam, fluid, beam_thickness, beam_width, beam_length, gap, dielectric_thickness=0, enableDielectric=False):
        self.beamCaseFile = beam
        self.fluidCaseFile = fluid
        self.beam_thickness = beam_thickness
        self.beam_width = beam_width
        self.beam_length = beam_length
        self.dielectric_thickness = dielectric_thickness
        self.gap = gap
        self.probeIndex = 0
        self.enableDielectric = enableDielectric

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
        self.fluidMeshesNew = self.fluidMeshes
        self.fluidMetricsCalculator = models.MeshMetricsCalculatorA(self.geomFields,self.fluidMeshes)
        self.fluidMetricsCalculator.init()

        self.solidBoundaryMeshes = [m.extrude(1, self.beam_thickness, True) for m in self.solidMeshes]
        self.solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(self.geomFields,self.solidBoundaryMeshes)
        
        self.solidBoundaryMetricsCalculator.init()

    def restart(self, restartFile):
        restartFile.readFluidMeshes(self.fluidMeshes)
        restartFile.readSolidMeshes(self.solidMeshes)

        self.solidBoundaryMeshes = [m.extrude(1, self.beam_thickness, True) for m in self.solidMeshes]
        self.solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(self.geomFields,self.solidBoundaryMeshes)
        solidNodeCoordA = self.solidMeshes[0].getNodeCoordinates()
        solidNodeCoord = solidNodeCoordA.asNumPyArray()
        solidBoundaryNodeCoordA = self.solidBoundaryMeshes[0].getNodeCoordinates()
        solidBoundaryNodeCoord = solidBoundaryNodeCoordA.asNumPyArray()
        ns = self.solidMeshes[0].getNodes().getSelfCount()
        nb = self.solidBoundaryMeshes[0].getNodes().getSelfCount()
        for n in range(0, ns):
            solidBoundaryNodeCoord[n][2] += solidNodeCoord[n][2]
            solidBoundaryNodeCoord[n+ns][2] += solidNodeCoord[n][2]
        self.solidBoundaryMetricsCalculator.init()

    def translate(self, tag, dx=0, dy=0, dz=0):
        if tag == 'solid' or tag == 'beam':
            for mesh in self.solidMeshes:
                nodes = mesh.getNodes()
                xNodesA = mesh.getNodeCoordinates()
                xNodes = xNodesA.asNumPyArray()
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
                xNodesA = mesh.getNodeCoordinates()
                xNodes = xNodesA.asNumPyArray()
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
                xNodesA = mesh.getNodeCoordinates()
                xNodes = xNodesA.asNumPyArray()
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
                xNodesA = mesh.getNodeCoordinates()
                xNodes = xNodesA.asNumPyArray()
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
    
    def rotate(self, tag, axis, alfa):
        if (tag == 'solid'  or tag == 'beam') and axis == 'y':
            for  mesh in self.solidMeshes:
                nodes = mesh.getNodes()
                xNodesA = mesh.getNodeCoordinates()
                xNodes = xNodesA.asNumPyArray()
                nNodes = nodes.getCount()
                for n in range(0, nNodes):
                    x = xNodes[n,0]
                    y = xNodes[n,2]
                    xp = x*cos(alfa) - y*sin(alfa)
                    yp = x*sin(alfa) + y*cos(alfa)
                    xNodes[n,0] = xp
                    xNodes[n,2] = yp
            self.solidMetricsCalculator.init()
            for  mesh in self.solidBoundaryMeshes:
                nodes = mesh.getNodes()
                xNodesA = mesh.getNodeCoordinates()
                xNodes = xNodesA.asNumPyArray()
                nNodes = nodes.getCount()
                for n in range(0, nNodes):
                    x = xNodes[n,0]
                    y = xNodes[n,2]
                    xp = x*cos(alfa) - y*sin(alfa)
                    yp = x*sin(alfa) + y*cos(alfa)
                    xNodes[n,0] = xp
                    xNodes[n,2] = yp
            self.solidBoundaryMetricsCalculator.init()
            print 'rotate beam mesh. Done'
            print '-------------------------------------------------------------'
        else:
            print 'not implemented yet'
                    
    def createShell(self, interfaceID):
        if self.enableDielectric == True:
            numMeshes = len(self.fluidMeshes)
            if numMeshes != 2:
                print "number of fluid meshes is wrong. Can not apply dielectric"
            else:
                self.shellMesh = self.fluidMeshes[0].createShell(interfaceID, self.fluidMeshes[1],interfaceID)
                self.fluidMeshesNew=[self.fluidMeshes[0], self.fluidMeshes[1], self.shellMesh]
            self.fluidMetricsCalculator = models.MeshMetricsCalculatorA(self.geomFields,self.fluidMeshesNew)
            self.fluidMetricsCalculator.init() 
   
    def summary(self):
        for mesh in self.solidMeshes:
            cells = mesh.getCells()
            nSelfCells = cells.getSelfCount()
            nCells = cells.getCount()
            print '-------------------------------------------------------------'
            print 'solid mesh: number of local cells %i' % nSelfCells 
            print 'solid mesh: number of total cells %i' % nCells 
            rCellsA = self.geomFields.coordinate[cells]
            rCells = rCellsA.asNumPyArray()
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
        
        for mesh in self.fluidMeshesNew:
            cells = mesh.getCells()
            nSelfCells = cells.getSelfCount()
            nCells = cells.getCount()
            print '-------------------------------------------------------------'
            print 'fluid mesh: number of local cells %i' % nSelfCells 
            print 'fluid mesh: number of total cells %i' % nCells 
            rCellsA = self.geomFields.coordinate[cells]
            rCells = rCellsA.asNumPyArray()
            rMin = rCells.min(axis=0)
            rMax = rCells.max(axis=0)
            print 'fluid mesh x range [ %e , %e ]' % (rMin[0], rMax[0])
            print 'fluid mesh y range [ %e , %e ]' % (rMin[1], rMax[1])
            print 'fluid mesh z range [ %e , %e ]' % (rMin[2], rMax[2])
            print '--------------------------------------------------------------'
          
        for mesh in self.solidBoundaryMeshes:
            faces = mesh.getFaces()
            nFaces = faces.getCount()
            print '--------------------------------------------------------------'
            print 'solid boundary mesh: number of faces %i' % nFaces
            rCellsA = self.geomFields.coordinate[cells]
            rCells = rCellsA.asNumPyArray()
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
