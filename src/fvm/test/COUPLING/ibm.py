### import  modules ###
import pdb
import sys
from math import *
sys.setdlopenflags(0x100|0x2)
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
from FluentCase import FluentCase
from fvm.fvmbaseExt import VecD3
import time
from mpi4py  import MPI
from Tools import *
from ComputeForce import *
from TimeStep import *
from Persistence import Persistence

def saveFluidVTK(n):
    writer = exporters.VTKWriterA(geomFields,fluidMeshesOrig,
                                  outDir+ "fluid-" + str(n) + ".vtk",
                                  "gen5_fluid",
                                  False,0)
    writer.init()
    writer.writeScalarField(elecFields.potential,"potential")
    writer.writeVectorField(elecFields.electric_field,"potentialgradient")  
    writer.finish()

def saveBeamBoundaryVTK(n):
    
    writer3 = exporters.VTKWriterA(geomFields,solidBoundaryMeshes,
                                   outDir + "beamBoundary-" + str(n) + ".vtk",
                                   "beam Boundary",
                                   False,0,True)
    writer3.init()
    writer3.finish()

def saveBeamVTK( n):
    writer = exporters.VTKWriterA(geomFields,solidMeshes,
                                  outDir + "beam-" + str(n) + ".vtk",
                                  "gen5_beam",
                                  False,0)
    writer.init()
    writer.writeVectorField(plateFields.deformation,"deformation")
    writer.finish()

###-------------------------------------------------------------###
applied_voltage = 10
numTimeSteps = 1

electrode = [8]
wall = [5,6,7,9]
interfaceID = 4

beamFix = [5, 6]
beamFree = [3, 4]

dielectric_thickness = 200e-9
dielectricID = -1
airID = 3
substrateID = 2
dc_Si =  11.68
dc_Air =  1.0
dc_SiN =  7.5

timeStep = 1e-8
globalTime=0.
globalCount = 0
probeIndex = 1103
saveFrequency = 1

gap = 3e-6
beam_thickness = 2e-6
rho = 8912.
E = 200e9
nu = 0.0

minR = 0.1e-9
maxR = 50e-9


fileBase = './'
outDir = './'
outFile = open(outDir + 'deformation.dat','w')
forceFile = open(outDir + 'force.dat', 'w')

restartFileName = outDir + 'restart.hdf5'
restartFileName = ''
saveFileName = outDir + 'checkpoint.hdf5'

if restartFileName != '':
    print 'I am restarting......'
    restartFile = Persistence(restartFileName,'r')
else:
    restartFile = None

solidCaseFile = fileBase + 'Gen5_beam_2D.cas'
fluidCaseFile = fileBase + 'fluid_electrode_dielectric.cas'

fluidReader = FluentCase(sys.argv[1])
solidReader = FluentCase(sys.argv[2])

fluidReader.read()
solidReader.read()

fluidMeshes = fluidReader.getMeshList()
fluidMeshesOrig = fluidMeshes
solidMeshes = solidReader.getMeshList()

if restartFile is not None:
    restartFile.readFluidMeshes(fluidMeshes)
    restartFile.readSolidMeshes(solidMeshes)

fluidMeshesOrig = fluidMeshes


geomFields =  models.GeomFields('geom')

fluidMetricsCalculator = models.MeshMetricsCalculatorA(geomFields, fluidMeshes)
solidMetricsCalculator = models.MeshMetricsCalculatorA(geomFields, solidMeshes)

fluidMetricsCalculator.init()
solidMetricsCalculator.init()

solidBoundaryMeshes = [m.extrude(1, beam_thickness, True) for m in solidMeshes]
solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(geomFields,solidBoundaryMeshes)

solidNodeCoord = solidMeshes[0].getNodeCoordinates().asNumPyArray()
solidBoundaryNodeCoord = solidBoundaryMeshes[0].getNodeCoordinates().asNumPyArray()
#pdb.set_trace()
ns = solidMeshes[0].getNodes().getSelfCount()
nb = solidBoundaryMeshes[0].getNodes().getSelfCount()
for n in range(0, ns):
    solidBoundaryNodeCoord[n][2] += solidNodeCoord[n][2]
    solidBoundaryNodeCoord[n+ns][2] += solidNodeCoord[n][2]

solidBoundaryMetricsCalculator.init()


shellMesh = fluidMeshes[0].createShell(interfaceID, fluidMeshes[1],interfaceID)
fluidMeshes=[fluidMeshes[0], fluidMeshes[1], shellMesh]
fluidMetricsCalculator = models.MeshMetricsCalculatorA(geomFields, fluidMeshes)
fluidMetricsCalculator.init() 

pc = fvmbaseExt.AMG()
pc.redirectPrintToFile("mesh_boundary.dat")
f = open("mesh_boundary.dat",'w')
for mesh in solidMeshes:
    cells = mesh.getCells()
    nSelfCells = cells.getSelfCount()
    nCells = cells.getCount()
    f.write( '-------------------------------------------------------------')
    f.write( 'solid mesh: number of local cells %i \n' % nSelfCells )
    f.write( 'solid mesh: number of total cells %i \n' % nCells )
    rCells = geomFields.coordinate[cells].asNumPyArray()
    rMin = rCells.min(axis=0)
    rMax = rCells.max(axis=0)
    f.write( 'solid mesh x range [ %e , %e ] \n' % (rMin[0], rMax[0]) )
    f.write( 'solid mesh y range [ %e , %e ] \n' % (rMin[1], rMax[1]) )
    f.write( 'solid mesh z range [ %e , %e ] \n' % (rMin[2], rMax[2]) )
    if rMin[2] != 0.0 or rMax[2] != 0.0:
        f.write( 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n' )
        f.write( 'Error: beam mesh is not centered at zero Z direction! \n' )
        f.write( 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
    f.write( '-------------------------------------------------------------')
        
for mesh in fluidMeshes:
    cells = mesh.getCells()
    nSelfCells = cells.getSelfCount()
    nCells = cells.getCount()
    f.write( '-------------------------------------------------------------\n')
    f.write( 'fluid mesh: number of local cells %i \n' % nSelfCells )
    f.write( 'fluid mesh: number of total cells %i \n' % nCells )
    rCells = geomFields.coordinate[cells].asNumPyArray()
    rMin = rCells.min(axis=0)
    rMax = rCells.max(axis=0)
    f.write( 'fluid mesh x range [ %e , %e ] \n' % (rMin[0], rMax[0]) )
    f.write( 'fluid mesh y range [ %e , %e ] \n' % (rMin[1], rMax[1]) )
    f.write( 'fluid mesh z range [ %e , %e ] \n' % (rMin[2], rMax[2]) )
    f.write( '--------------------------------------------------------------\n')
          
for mesh in solidBoundaryMeshes:
    faces = mesh.getFaces()
    nFaces = faces.getCount()
    f.write( '--------------------------------------------------------------\n')
    f.write( 'solid boundary mesh: number of faces %i \n' % nFaces)
    rCells = geomFields.coordinate[faces].asNumPyArray()
    rMin = rCells.min(axis=0)
    rMax = rCells.max(axis=0)
    f.write( 'solid boundary mesh x range [ %e , %e ] \n' % (rMin[0], rMax[0]) )
    f.write( 'solid boundary mesh y range [ %e , %e ] \n' % (rMin[1], rMax[1]) )
    f.write( 'solid boudnary mesh z range [ %e , %e ] \n' % (rMin[2], rMax[2]) )
    f.write( '--------------------------------------------------------------')
f.close()
pc.redirectPrintToScreen()

