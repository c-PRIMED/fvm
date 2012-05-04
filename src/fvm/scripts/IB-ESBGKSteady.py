#!/usr/bin/env python

### 2D cantilever structure/fluid/electrostatics IBM coupling  ###

### import  modules ###
#import pdb
import sys
import os
from math import *
sys.setdlopenflags(0x100|0x2)
#import tecplotExporter
import fvm
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters
fvm.set_atype('double')
from FluentCase import FluentCase
from fvm.fvmbaseExt import VecD3
from mpi4py import MPI
import time
import tecplotESBGKEntireDomain
import tecplotESBGK
import tecplotESBGKIB
from FluentCase import FluentCase
from optparse import OptionParser

### ========================== define functions and class =============================###
class ProblemDefinition(object):
    pass

   
def saveVTK(nstep, pd):
    if pd.initialTransient:
        return

    
    writer2 = exporters.VTKWriterA(pd.geomFields,pd.macroFields,pd.fluidMeshes,
                                   "elecfield-" + str(nstep) + ".vtk",
                                   "fix-fix beam",
                                   False,0)
    writer2.init()
    writer2.writeVectorField()
    writer2.finish()
  


def checkMarking(n, pd):

    cells = pd.fluidMeshes[0].getCells()
    nCells = cells.getCount()
    cellCoords = pd.geomFields.coordinate[cells].asNumPyArray()

    fluidFile = open("fluidCells_" + str(pd.globalCount) + ".dat", "w")
    solidFile = open("solidCells_" + str(pd.globalCount) + ".dat", "w")
    IBFile = open("IBCells_" + str(pd.globalCount) + ".dat", "w")

    cellIBType = pd.geomFields.ibType[cells].asNumPyArray()
    for c in range (0, nCells):    
        ibtype = cellIBType[c]
        if ibtype == -1:
            fluidFile.write("%e\t%e\t%e\n" % (cellCoords[c][0], cellCoords[c][1], cellCoords[c][2]))
        elif ibtype == -2:
            IBFile.write("%e\t%e\t%e\n" % (cellCoords[c][0], cellCoords[c][1], cellCoords[c][2]))
        elif ibtype == -3:
            solidFile.write("%e\t%e\t%e\n" % (cellCoords[c][0], cellCoords[c][1], cellCoords[c][2]))
        elif ibtype == -5:
            print ("%i\t%i\t%e\t%e\n" % (c,ibtype,  cellCoords[c][0], cellCoords[c][1]))

               
    fluidFile.close()
    solidFile.close()
    IBFile.close()

            
def unsteadyAdvance(pd, numTimeSteps):
    
    for n in range(0, numTimeSteps):

        sbMeshFaces = pd.solidBoundaryMeshes[0].getFaces()

        #-------------update IBM----------------#
        print "***       update IBM  at  %i           ***" % n
        pd.ibManager.update()
        pd.ibManager.fluidNeighborsPerIBFace = 4
        pd.ibManager.solidNeighborsPerIBFace = 4
        pd.ibManager.fluidNeighborsPerSolidFace = 6
        
        pd.fluidMetricsCalculator.computeIBInterpolationMatrices(sbMeshFaces)
        pd.fluidMetricsCalculator.computeSolidInterpolationMatrices(sbMeshFaces)     
         
        #------------solve fluid ---------------#
        print "***      solving flow model  at globalCount %i    ***" % pd.globalCount
        bcMap = esbgk1.getBCMap()    
        esbgk1.callBoundaryConditions()
        esbgk1.computeSolidFaceDsf(sbMeshFaces,1)       
        esbgk1.ConservationofMFSolid(sbMeshFaces)
        esbgk1.computeIBFaceDsf(sbMeshFaces,1)
        esbgk1.advance(numIterationsPerStep,sbMeshFaces)

        if (n%pd.saveFrequency == 0):
            tecplotESBGKIB.esbgkTecplotFile(pd.fluidMeshes,pd.macroFields,pd.geomFields,"quad","timestep_" + str(pd.globalCount) + ".dat")
##     To blank the IB faces in tecplot file -in the software /plot/blanking/active blank
##     when collisionfrequency is greater than 0       
       
        #---------------update time -------------------------#
        if esbgk1options.transient:    
            esbgk1.updateTime()
            esbgk1.EntropyGeneration()
        pd.globalTime += pd.timeStep
        pd.globalCount += 1
 #---------------boundary update -------------------------#
        for mesh in pd.solidMeshes:
            nodes = mesh.getNodes()
            xNodes = mesh.getNodeCoordinates().asNumPyArray()
            xNodes[:,0] += 0.0
            xNodes[:,1] += 0.0
            xNodes[:,2] += 0.0
        pd.solidMetricsCalculator.init()

        for mesh in pd.solidBoundaryMeshes:
            nodes = mesh.getNodes()
            xNodes = mesh.getNodeCoordinates().asNumPyArray()
            xNodes[:,0] += 0.0
            xNodes[:,1] += 0.0
            xNodes[:,2] += 0.0
        pd.solidBoundaryMetricsCalculator.init()       
      

        for mesh in pd.solidMeshes:
            cells = mesh.getCells()
            nSelfCells = cells.getSelfCount()
            nCells = cells.getCount()
            print '---------------------------------------------------------'
            print 'solid mesh: number of local cells %i' % nSelfCells 
            print 'solid mesh: number of total cells %i' % nCells 
            rCells = pd.geomFields.coordinate[cells].asNumPyArray()
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

        for mesh in pd.solidBoundaryMeshes:
            faces = mesh.getFaces()
            nFaces = faces.getCount()
            print '--------------------------------------------------------------'
            print 'solid boundary mesh: number of faces %i' % nFaces
            rCells = pd.geomFields.coordinate[faces].asNumPyArray()
            rMin = rCells.min(axis=0)
            rMax = rCells.max(axis=0)
            print 'solid boundary mesh x range [ %e , %e ]' % (rMin[0], rMax[0])
            print 'solid boundary mesh y range [ %e , %e ]' % (rMin[1], rMax[1])
            print 'solid boudnary mesh z range [ %e , %e ]' % (rMin[2], rMax[2])
            print '--------------------------------------------------------------'
              
               
### ========================== properties and parameters ===============================###
pd = ProblemDefinition()
### beam
rho = 8800                     # density kg/m^3
E = 200e9                     # Young's modulus 
nu = 0.31                      # Poisson's ratio


### mesh id

fluidTop = 6
fluidBot = 5
fluidLeft = 3
fluidRight = 4

beamTop = 6
beamBot = 5
beamRight = 4
beamLeft = 3

#user option
pd.ibMethod=1; #1 for Interpolation 2 for Relaxation (Relaxation is not available)
numTimeSteps = 50
numIterationsPerStep=30
pd.output_interval =5
output_Coeff=5
frequency = 133862.96;
pd.timeStep = 1.0/(frequency*100.)
relTol=1e-7;
absTol=1e-22;

cellno=8   #print coeffs for BGK and ESBGK &distribution function
fgamma=2; #0 max, 1 BGK, 2ES
method =1; # 1 constant, linear otherwise
ktrial=20;

pd.nondimlength=1;
pd.rho_init=1.174; #1atmosphere pressure at 300
T_init=300;

accomCoeff=1.0;
resint=0;
restartFileName = "test"+str(resint)+"_"+str(MPI.COMM_WORLD.Get_rank())+".hdf5" 
restartFileName = ""  

gas='Air'

if (gas == 'Air'):
    Pr=3.0/4.0;
    SpHeatRatio=7.0/5.0;
    muref= 1.7116e-5;
    mu_w=0.74;
    molecularWeight=28.9;
elif (gas == 'Argon'):
    Pr=2.0/3.0;
    SpHeatRatio=5.0/3.0;
    muref= 2.117e-5;
    mu_w=0.81;
    molecularWeight=39.9;
elif (gas == 'Nitrogen'):
    Pr=3.0/4.0;
    SpHeatRatio=7.0/5.0;
    muref= 1.781e-5;
    mu_w=0.74;
    molecularWeight=28.0;
        
    
pd.u_init=(2.0*8314.0*T_init/molecularWeight)**0.5
pd.ubeam=-2.0/pd.u_init;
tbeam=1.0
print ' timeStep ',pd.timeStep



### =========================== running script  ========================================###

#pd = ProblemDefinition()

pd.timeStepND=pd.timeStep/pd.nondimlength*pd.u_init
fileBase = "/home/ba01/u140/cpekarda/Couette/CPMesh/"
outfile = None
if __name__ == '__main__' and fileBase is None:
    if len(sys.argv) < 2:
        usage()
    fileBase = sys.argv[1]
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

if outfile == None:
    outfile = fileBase+"-tecplt.dat"

pd.globalTime = 0
pd.globalCount = 0
pd.saveFrequency = pd.output_interval
pd.initialTransient=False
pd.probeIndex = 5014
pd.probeFile = open(fileBase + "tipDisplacement.dat", "w")
pd.forceFile = open(fileBase + "beamForce.dat", "w")
pd.velocityFile = open(fileBase + "tipVelocity.dat", "w")

### read in meshes ###

beamFile = fileBase + 'new0.5dxbeam.cas'
fluidFile = fileBase + 'test2fluidNS.cas'

beamReader = FluentCase(beamFile)
fluidReader = FluentCase(fluidFile)

beamReader.read()
fluidReader.read()

pd.solidMeshes = beamReader.getMeshList()
pd.fluidMeshes = fluidReader.getMeshList() 
pd.solidBoundaryMeshes = [m.extractBoundaryMesh() for m in pd.solidMeshes]
etype = 2

### geometry field ###

pd.geomFields =  models.GeomFields('geom')
pd.solidMetricsCalculator = models.MeshMetricsCalculatorA(pd.geomFields,
                                                          pd.solidMeshes)
pd.fluidMetricsCalculator = models.MeshMetricsCalculatorA(pd.geomFields,
                                                          pd.fluidMeshes)
pd.solidBoundaryMetricsCalculator = models.MeshMetricsCalculatorA(pd.geomFields,
                                                          pd.solidBoundaryMeshes)

pd.solidMetricsCalculator.init()
pd.fluidMetricsCalculator.init()
pd.solidBoundaryMetricsCalculator.init()


if restartFileName != "":
    restartFile = Persistence(restartFileName,'r')
else:
    restartFile = None

# to test quadrature
import fvm.esbgk_atyped_double as esbgk

#facegroups
fgs = pd.fluidMeshes[0].getBoundaryFaceGroups()
for fg in fgs:
    print fg.groupType
    if  fg.id == 6:
        fg.groupType = "wall"
    if  fg.id in [3,4,5]:
	fg.groupType = "wall"
    print fg.id,"  ",fg.groupType

#Quadrature- velocity mesh
#pd.quad0=esbgk.QuadratureD(20,20,20,5.5,1.0) #cartesian
pd.quad0=esbgk.QuadratureD(14,14,14,5.5,1.0) #cartesian
#quad0=esbgk.QuadratureD(8,8,0,32,0,16) #spherical
#quad2=esbgk.QuadratureD(16,16,1,8,1,4) #gauss-hermit quadrature and 3/8th rule

#import ddd
pd.macroFields=esbgk.MacroFields('flow')


print "1"

esbgk1=esbgk.KineticModelD(pd.fluidMeshes,pd.geomFields,pd.macroFields,pd.quad0)

print "2"
esbgk1options = esbgk1.getOptions()

esbgk1options.fgamma=fgamma;  #0 max, 1 for BGK,2  ES otherwise
esbgk1options.ibMethod=pd.ibMethod;
esbgk1options.printCellNumber=cellno;
esbgk1options.NewtonsMethod_ktrial=ktrial;
esbgk1options.setVar('rho_init',pd.rho_init);
esbgk1options.setVar('T_init',T_init);

esbgk1options.Prandtl=Pr;
esbgk1options.SpHeatRatio=SpHeatRatio;
esbgk1options.molecularWeight=molecularWeight;
esbgk1options.muref=muref;
esbgk1options.mu_w=mu_w;

esbgk1options.transient =False;
esbgk1options.setVar('nonDimLt',pd.nondimlength);
esbgk1options.setVar('nonDimLx',pd.nondimlength);
esbgk1options.setVar('nonDimLy',pd.nondimlength);
esbgk1options.setVar('nonDimLz',pd.nondimlength);

#Method 1:
#tSolver = fvmbaseExt.JacobiSolver()

#Method 2:
pc = fvmbaseExt.JacobiSolver()
tSolver = fvmbaseExt.BCGStab()
tSolver.preconditioner=pc

tSolver.verbosity=0

esbgk1options.KineticLinearSolver=tSolver


## initialize macroparameters and f
esbgk1.InitializeMacroparameters()

cellSites=pd.fluidMeshes[0].getCells()
cellCentroid=pd.geomFields.coordinate[cellSites].asNumPyArray()
densityField =  pd.macroFields.density[cellSites].asNumPyArray()
nuField = pd.macroFields.collisionFrequency[cellSites].asNumPyArray()
ncells=cellSites.getSelfCount()
ncellsall=cellSites.getCount()
print ncellsall
if restartFile is None:
    if method ==1:
        esbgk1.weightedMaxwellian(0.5,0,0,1.0,1.0) #initial distribution
dsfList = []
ndir = pd.quad0.getDirCount()
print 'no. of directions ', ndir

DistFunc=esbgk1.getdsf()
for i in range(0,ndir):
    dsfList.append( DistFunc.getField(i) )

if esbgk1options.transient:
    DistFunc1=esbgk1.getdsf1()
    for i in range(0,ndir):
        dsfList.append( DistFunc1.getField(i) )
    
    if esbgk1options.timeDiscretizationOrder > 1:
        DistFunc2=esbgk1.getdsf2()
        for i in range(0,ndir):
            dsfList.append( DistFunc2.getField(i) )
        
if restartFile is not None:
    restartFile.readKineticModel(pd.macroFields,esbgk1,pd.fluidMeshes,dsfList,ndir)
    restartFile.close()
    print "read f from restart file"
       
#initialize feq

print "initializing macropr and feq"
esbgk1.ComputeMacroparameters()

if (fgamma==0):
    esbgk1.initializeMaxwellianEq()
else:
    esbgk1.EquilibriumDistributionBGK()
if(fgamma==2):
    esbgk1.EquilibriumDistributionESBGK()

#collision frequency based on Prandlt
esbgk1.ComputeCollisionfrequency()


print "initializing bcs"

bcMap = esbgk1.getBCMap()
if 4 in bcMap:
    bcTop = bcMap[4]
    bcTop.bcType = 'WallBC'
    bcTop.setVar('specifiedTemperature',1.0)
if 6 in bcMap:
    bcBot = bcMap[6]
    bcBot.bcType = 'PressureInletBC'
    bcBot.setVar('specifiedTemperature',1.0)
    bcBot.setVar('specifiedPressure',1.0)
if 3 in bcMap:
    bcLeft = bcMap[3]
    bcLeft.bcType = 'PressureInletBC'
    bcLeft.setVar('specifiedTemperature',1.0)
    bcLeft.setVar('specifiedPressure',1.0)
if 5 in bcMap:
    bcRight = bcMap[5]
    bcRight.bcType = 'PressureInletBC'
    bcRight.setVar('specifiedTemperature',1.0)
    bcRight.setVar('specifiedPressure',1.0)

vcMap=esbgk1.getVCMap()
print "initializing bcs"        
esbgk1options.setVar('timeStep',pd.timeStepND)
esbgk1options.relativeTolerance=relTol;
esbgk1options.absoluteTolerance=absTol;

pd.prevIter=0
print "initializing bcs"
### solver options ###
tecplotESBGK.esbgkTecplotFile(pd.fluidMeshes,pd.macroFields,"quad","initial.dat")

### set up IBM ###

for mesh in pd.solidBoundaryMeshes:
    faces = mesh.getFaces()
    #temperature field on boundary mesh
    areaMag = pd.geomFields.areaMag[faces]
    faceCount = faces.getCount()
    tem = areaMag.newSizedClone(faceCount)
    tema = tem.asNumPyArray()
    tema[:] = tbeam
    pd.macroFields.temperature[faces] = tem
    dens = areaMag.newSizedClone(faceCount)
    densa = dens.asNumPyArray()
    densa[:] = pd.rho_init
    pd.macroFields.density[faces] = dens
    #flow velocity field on boundary mesh
    area = pd.geomFields.area[faces]
    faceCount = faces.getCount()
    vel = area.newSizedClone(faceCount)
    vela = vel.asNumPyArray()
    vela[:,1] = pd.ubeam
    vela[:,2] = 0.
    vela[:,0] = 0.
    pd.macroFields.velocity[faces] = vel

 
    
pd.ibManager = fvmbaseExt.IBManager(pd.geomFields,
                                    pd.solidBoundaryMeshes[0],
                                    pd.fluidMeshes)

pd.ibManager.fluidNeighborsPerIBFace = 4
pd.ibManager.solidNeighborsPerIBFace = 4
pd.ibManager.fluidNeighborsPerSolidFace = 6

pd.ibManager.update()
checkMarking(0, pd)

t1 = time.time()

                    
unsteadyAdvance(pd,numTimeSteps)



t2 = time.time()

print  '\nsolution time = %f' % (t2-t1)
