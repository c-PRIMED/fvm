#!/usr/bin/env python

import sys
sys.setdlopenflags(0x100|0x2)

import fvm
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import numpy

fvm.set_atype('double')

import math

if fvm.atype == 'double':
    import fvm.models_atyped_double as models
    import fvm.exporters_atyped_double as exporters
elif fvm.atype == 'tangent':
    import fvm.models_atyped_tangent_double as models
    import fvm.exporters_atyped_tangent_double as exporters

from FluentCase import FluentCase
from optparse import OptionParser

#fvmbaseExt.enableDebug("cdtor")

fileBase0 = None
fileBase1 = None
numIterations = 100
numEIterations = 100
sPot = 160.

fileBase0 = "/scratch/prism/shankha/prism/prism3/memosa/src/fvm/test/shankha/structure/trdeform14/dbeam1"
fileBase1 = "/scratch/prism/shankha/prism/prism3/memosa/src/fvm/test/shankha/structure/trdeform14/dbeam2"

def eadvance(fmodel,niter):
    for i in range(0,niter):
        try:
            stopFlag=fmodel.advance(1)
            if stopFlag == 1:
                break
        except KeyboardInterrupt:
            break

def setDirichletCommonDisplacement(dmodel,geomFields,meshes,structureFields):
    mesh0 = meshes[0]
    nodes0 = mesh0.getNodes()
    cells0 = mesh0.getCells()
    common0 = dmodel.getCommon(nodes0).asNumPyArray()
    mesh1 = meshes[1]
    nodes1 = mesh1.getNodes()
    common1 = dmodel.getCommon(nodes1).asNumPyArray()
    length = len(common0)

    def1 = geomFields.dirichletNodeDisplacement[nodes1].asNumPyArray()
    coord0N0 =  geomFields.coordinate[nodes0].asNumPyArray()
    coord0N1 =  geomFields.coordinateN1[nodes0].asNumPyArray()
    for i in range(0,length):
        id0 = common0[i]
        id1 = common1[i]
#        print '\n coord0 is %e %e %i' %(coord0N1[id0,0],coord0N1[id0,1],length)
        def1[id1] = coord0N0[id0] - coord0N1[id0]
    
def advance(smodel,dmodel,movingMeshModel,emodel,geomFields,
            structureFields,electricFields,meshes0,meshes,niter):
    for i in range(0,niter):
        try:
            bcMap = smodel.getBCMap()
            
            bcID = 3
            if bcID in bcMap:
                bc = smodel.getBCMap()[bcID]
                bc.bcType = 'SpecifiedDistForce'
                bcEID = 5
                felec=createBVFields(geomFields,meshes,bcID,bcEID,structureFields,electricFields)
                bc['specifiedXDistForce']=0
                bc['specifiedYDistForce']=felec
                bc['specifiedZDistForce']=0
                
            sk=smodel.advance(1)
            dmodel.calculateNodeDisplacement()
            dmodel.deformStructure()
            
            setDirichletCommonDisplacement(dmodel,geomFields,meshes,structureFields)
            movingMeshModel.advance()
            metricsCalculator.recalculate_deform()
            eadvance(emodel,numEIterations)
            if(sk==1):
                break
        except KeyboardInterrupt:
            break

def advanceUnsteady(smodel,dmodel,movingMeshModel,emodel,geomFields,
                    structureFields,electricFields,meshes0,meshes,nTimeSteps,globalTime):
    fileName = fileBase0 + "middef.txt"
    file = open(fileName,"w")
    mesh0 = meshes0[0]
    deformation =  structureFields.deformation[mesh0.getCells()].asNumPyArray()
    file.write(" %e " % globalTime)
    file.write(" %e " % deformation[500][0])
    file.write(" %e " % deformation[500][1])
    file.write("\n")
    for i in range(0,nTimeSteps):
        try:
            advance(smodel,dmodel,movingMeshModel,emodel,geomFields,
                    structureFields,electricFields,meshes0,meshes,numIterations)
            globalTime += timeStep
            print 'advancing to time %e at iteration %i' % (globalTime,i)
            file.write(" %e " % globalTime)
            file.write(" %e " % deformation[500][0])
            file.write(" %e " % deformation[500][1])
            file.write("\n")
            smodel.updateTime()
        except KeyboardInterrupt:
            break

def createBVFields(geomFields,meshes,id,eid,structureFields,electricFields):
    fy = fvmbaseExt.Field('bvy')
    
    mesh0 = meshes[0]
    mesh1 = meshes[1]
    vol = geomFields.volume[mesh0.getCells()]
    

    deflection =  structureFields.deformation[mesh0.getCells()].asNumPyArray()
    fgs0 = mesh0.getBoundaryGroups()
    fgs1 = mesh1.getBoundaryGroups()
    for fg1 in fgs1:
        if fg1.id == eid:
            ewall = fg1.site
    bArea = geomFields.area[ewall].asNumPyArray().copy()
    bpflux =  electricFields.potential_flux[ewall].asNumPyArray()
    for fg in fgs0:
        if fg.id==id:
            
            faceCells = mesh.getFaceCells(fg.site)
            nFaces = fg.site.getCount()
            forceY = vol.newSizedClone(nFaces)
            forceYa = forceY.asNumPyArray()
            xf =  geomFields.coordinate[fg.site].asNumPyArray()
            
            pot_top=sPot
            pot_bot=0.0
            bSurface = -3.75e-6
            perm=8.8542e-12
            
            for i in range(0,nFaces):
                c0 = faceCells(i,0)
                gap = deflection[c0,1]-bSurface
#                dpot = (pot_top-pot_bot)/gap
                magBArea = math.sqrt(bArea[i][0]*bArea[i][0]+bArea[i][1]*bArea[i][1]+bArea[i][2]*bArea[i][2])
                dpot = bpflux[i]/magBArea
                sigmat=-perm*dpot
                felec=-(sigmat*sigmat)/(2.*perm)
                forceYa[i]=felec
#                print 'force %f %f %f %e %e %e' % (xf[i,0],xf[i,1],forceYa[i],bpflux[i],dpot,magBArea)
            fy[fg.site] = forceY
            return fy
                
# change as needed

# map between fvm, tecplot, and xdmf types
etype = {
    'tri' : 1,
    'quad' : 2,
    'tetra' : 3,
    'hexa' : 4
        }
tectype = {
    'tri' : 'FETRIANGLE',
    'quad' : 'FEQUADRILATERAL',
    'tetra' : 'FETETRAHEDRON',
    'hexa' : 'FEBRICK'
        }

def dumpTecplotFile(nmesh, meshes, geomFields, mtype):
    #cell sites
    cellSites = []
    for n in range(0,nmesh):
        cellSites.append( meshes[n].getCells() )
#     print "cellSites[", n, "].getCount = ", cellSites[n].getCount()
        
    #face sites
    faceSites = []
    for n in range(0,nmesh):
        faceSites.append( meshes[n].getFaces() )
        
    #node sites
    nodeSites = []
    for n in range(0,nmesh):
        nodeSites.append( meshes[n].getNodes() )
        
    #get connectivity (faceCells)
    faceCells = []
    for n in range(0,nmesh):
        faceCells.append( meshes[n].getConnectivity( faceSites[n], cellSites[n] ) )
        
    #get connectivity ( cellNodes )
    cellNodes = []
    for n in range(0,nmesh):
        cellNodes.append( meshes[n].getCellNodes() )
        
    #get Volume as array
    volumes = []
    for n in range(0,nmesh):
        volumes.append( geomFields.volume[cellSites[n]].asNumPyArray() )
    
    
    cellCentroids =[]
    for n in range(0,nmesh):
        cellCentroids.append( geomFields.coordinate[cellSites[n]].asNumPyArray() )
        
#    defFields = []
#    for n in range(0,nmesh):
#        defFields.append( electricFields.potential[cellSites[n]].asNumPyArray() )

#    tractionXFields = []
#    for n in range(0,nmesh):
#        tractionXFields.append( structureFields.tractionX[cellSites[n]].asNumPyArray() )
        
    coords = []
    for n in range(0,nmesh):
        coords.append( geomFields.coordinate[nodeSites[n]].asNumPyArray() )
#     print "shape( coords[", n, "] ) = ", shape( coords[n] )

    f = open("tecplot_dbeam.dat","w")
    f.write("Title = \" tecplot file for 2D Cavity problem \" \n")
    f.write("variables = \"x\", \"y\", \"z\", \"cellCentroidY\" \n")
    for n in range(0,nmesh):
        title_name = "nmesh%s" % n
        ncell  = cellSites[n].getSelfCount()
        nnode  = nodeSites[n].getCount()
        f.write("Zone T = \"%s\" N = %s E = %s DATAPACKING = BLOCK, VARLOCATION = ([4]=CELLCENTERED), ZONETYPE=%s\n" %
                (title_name,  nodeSites[n].getCount(), ncell, tectype[mtype]))
        #write x
        for i in range(0,nnode):
            f.write(str(coords[n][i][0])+"    ")
            if ( i % 5 == 4 ):
                f.write("\n")
        f.write("\n")
                
        #write y
        for i in range(0,nnode):
            f.write(str(coords[n][i][1])+"    ")
            if ( i % 5 == 4 ):
                f.write("\n")
        f.write("\n")

        #write z
        for i in range(0,nnode):
            f.write(str(coords[n][i][2])+"    ")
            if ( i % 5 == 4 ):
                f.write("\n")
        f.write("\n")

#        #write defX
#        for i in range(0,ncell):
#            f.write(str(defFields[n][i]) + "    ")
#            if ( i % 5  == 4 ):
#                f.write("\n")
#        f.write("\n")

#        #write defY
#        for i in range(0,ncell):
#            f.write(str(defFields[n][i][1]) + "    ")
#            if ( i % 5  == 4 ):
#                f.write("\n")
#        f.write("\n")

#        #write sigmaXX
#        for i in range(0,ncell):
#            f.write(str(tractionXFields[n][i][0]) + "    ")
#            if ( i % 5  == 4 ):
#                f.write("\n")
#        f.write("\n")

#        #write sigmaXY
#        for i in range(0,ncell):
#            f.write(str(tractionXFields[n][i][1]) + "    ")
#            if ( i % 5  == 4 ):
#                f.write("\n")
#        f.write("\n")
                                                    
#        #write sigmaYY
#        for i in range(0,ncell):
#            f.write(str(tractionXFields[n][i][2]) + "    ")
#            if ( i % 5  == 4 ):
#                f.write("\n")
#        f.write("\n")
                                                                    
        #write velX
        for i in range(0,ncell):
            f.write( str(cellCentroids[n][i][1]) + "    ")
            if ( i % 5  == 4 ):
                f.write("\n")
        f.write("\n")

        #connectivity
        for i in range(0,ncell):
            nnodes_per_cell = cellNodes[n].getCount(i)
            for node in range(0,nnodes_per_cell):
                f.write( str(cellNodes[n](i,node)+1) + "     ")
            f.write("\n")
        f.write("\n")
    f.close()

parser = OptionParser()
parser.set_defaults(type='quad')
parser.add_option("--type", help="'quad'[default], 'tri', 'hexa', or 'tetra'")
parser.add_option("--xdmf", action='store_true', help="Dump data in xdmf")
parser.add_option("--time","-t",action='store_true',help="Print timing information.")
(options, args) = parser.parse_args()

reader0 = FluentCase(fileBase0+".cas")
reader1 = FluentCase(fileBase1+".cas")


#import debug
reader0.read();
reader1.read();

meshes0 = reader0.getMeshList()
meshes1 = reader1.getMeshList()

#for mesh in meshes:
#    mesh.getCells().clearGatherScatterMaps()

mesh0 = meshes0[0]
mesh1 = meshes1[0]
nmesh = 2

mesh0.findCommonNodes(mesh1)

meshes = []
meshes.append(mesh0)
meshes.append(mesh1)

import time
t0 = time.time()

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

nodes0 = mesh0.getNodes()
nodes1 = mesh1.getNodes()

rho = 7854.0
E = 2.0*math.pow(10,11)
nu = 0.31

if fvm.atype == 'tangent':
    metricsCalculator.setTangentCoords(0,7,1)

flowFields =  models.FlowFields('flow')
structureFields =  models.StructureFields('structure')
electricFields =  models.ElectricFields('elec')

smodel = models.StructureModelA(geomFields,structureFields,meshes0)
dmodel = models.StructureDeformationModelA(geomFields,structureFields,meshes0)
movingMeshModel = models.MovingMeshModelA(meshes1,geomFields,flowFields)
emodel = models.ElectricModelA(geomFields,electricFields,meshes1)

movingMeshModel.init()

bcMap = smodel.getBCMap()

#left (mesh0)
bcID = 6
if bcID in bcMap:
    bc = smodel.getBCMap()[bcID]
    bc.bcType = 'SpecifiedDeformation'
    bc['specifiedXDeformation']=0
    bc['specifiedYDeformation']=0
    bc['specifiedZDeformation']=0
                    
#top (mesh0)
bcID = 5
if bcID in bcMap:
    bc = smodel.getBCMap()[bcID]
    bc.bcType = 'SpecifiedTraction'
    bc['specifiedXXTraction']=0
    bc['specifiedXYTraction']=0
    bc['specifiedXZTraction']=0
    bc['specifiedYXTraction']=0
    bc['specifiedYYTraction']=0
    bc['specifiedYZTraction']=0
    bc['specifiedZXTraction']=0
    bc['specifiedZYTraction']=0
    bc['specifiedZZTraction']=0

#right (mesh0)
bcID = 4
if bcID in bcMap:
    bc = smodel.getBCMap()[bcID]
    bc.bcType = 'SpecifiedDeformation'
    bc['specifiedXDeformation']=0
    bc['specifiedYDeformation']=0
    bc['specifiedZDeformation']=0
                    
#bottom (mesh0)
#bcID = 3
#if bcID in bcMap:
#    bc = smodel.getBCMap()[bcID]
#    felec=createBVFields(geomFields,meshes0,bcID,structureFields,electricFields)
#    bc.bcType = 'SpecifiedDistForce'
#    bc['specifiedXDistForce']=0
#    bc['specifiedYDistForce']=felec
#    bc['specifiedZDistForce']=0

f2 = open("displacementOptions.dat","w")
for mesh in meshes1:
    nodes = mesh.getNodes()
    displacementOptions = geomFields.displacementOptions[nodes].asNumPyArray()
    nodeCoordinate = geomFields.coordinate[nodes].asNumPyArray()
    nodemark = numpy.zeros(nodes.getCount())
    fgs = mesh.getAllFaceGroups()
    for fg in fgs:
        if fg.id!=0:
            if fg.id == 5:
                fgsite = fg.site
                fgn = mesh.getFaceNodes(fgsite)
                nfaces =fgsite.getCount()
                for nf in range(0,nfaces):
                    nnodes = fgn.getCount(nf)
                    for nnode in range(0,nnodes):
                        nid = fgn(nf,nnode)
                        if nodemark[nid] == 0:
                            nodemark[nid] = 1
                            displacementOptions[nid] = 1
                            f2.write('%i\t' % fg.id)
                            f2.write('%i\t' % displacementOptions[nid])
                            f2.write('%f\t' % nodeCoordinate[nid,0])
                            f2.write('%f\t' % nodeCoordinate[nid,1])
                            f2.write('%f\n' % nodeCoordinate[nid,2])
    for fg in fgs:
        if fg.id!=0:
            if fg.id == 3:
                fgsite = fg.site
                fgn = mesh.getFaceNodes(fgsite)
                nfaces =fgsite.getCount()
                for nf in range(0,nfaces):
                    nnodes = fgn.getCount(nf)
                    for nnode in range(0,nnodes):
                        nid = fgn(nf,nnode)
                        if nodemark[nid] == 0:
                            nodemark[nid] = 1
                            displacementOptions[nid] = 0
                            f2.write('%i\t' % fg.id)
                            f2.write('%i\t' % displacementOptions[nid])
                            f2.write('%f\t' % nodeCoordinate[nid,0])
                            f2.write('%f\t' % nodeCoordinate[nid,1])
                            f2.write('%f\n' % nodeCoordinate[nid,2])
    for fg in fgs:
        if fg.id!=0:
            if fg.id in (6,4):
                fgsite = fg.site
                fgn = mesh.getFaceNodes(fgsite)
                nfaces =fgsite.getCount()
                for nf in range(0,nfaces):
                    nnodes = fgn.getCount(nf)
                    for nnode in range(0,nnodes):
                        nid = fgn(nf,nnode)
                        if nodemark[nid] == 0:
                            nodemark[nid] = 1
                            displacementOptions[nid] = 2
                            f2.write('%i\t' % fg.id)
                            f2.write('%i\t' % displacementOptions[nid])
                            f2.write('%f\t' % nodeCoordinate[nid,0])
                            f2.write('%f\t' % nodeCoordinate[nid,1])
                            f2.write('%f\n' % nodeCoordinate[nid,2])
f2.close()
                        
bcElecMap = emodel.getBCMap()

#top
bcID = 5
if bcID in bcElecMap:
    bc = emodel.getBCMap()[bcID]
    bc.bcType = 'SpecifiedPotential'
    bc.setVar('specifiedPotential',sPot)
    #bc.bcType = 'SpecifiedCharge'
    #bc.setVar('specifiedCharge',400)

#bot
bcID = 3
if bcID in bcElecMap:
    bc = emodel.getBCMap()[bcID]
    bc.bcType = 'SpecifiedPotential'
    bc.setVar('specifiedPotential',0)
    #bc.bcType = 'SpecifiedCharge'
    #bc.setVar('specifiedCharge',300)

#left
bcID = 6
if bcID in bcElecMap:
    bc = emodel.getBCMap()[bcID]
#    bc.bcType = 'SpecifiedPotential'
#    bc.setVar('specifiedPotential',350)
    bc.bcType = 'SpecifiedPotentialFlux'
    bc.setVar('specifiedPotentialFlux',0)
    #bc.bcType = 'SpecifiedCharge'
    #bc.setVar('specifiedCharge',350)

#right
bcID = 4
if bcID in bcElecMap:
    bc = emodel.getBCMap()[bcID]
#    bc.bcType = 'SpecifiedPotential'
#    bc.setVar('specifiedPotential',350)
    bc.bcType = 'SpecifiedPotentialFlux'
    bc.setVar('specifiedPotentialFlux',0)
    #bc.bcType = 'SpecifiedCharge'
    #bc.setVar('specifiedCharge',350)
    

vcMap = smodel.getVCMap()
for i,vc in vcMap.iteritems():
    vc['density'] = rho
    vc['eta'] = E/(2.*(1+nu))
    vc['eta1'] = nu*E/((1+nu)*(1-1.0*nu))

pc = fvmbaseExt.AMG()
pc.verbosity=0
defSolver = fvmbaseExt.BCGStab()
defSolver.preconditioner = pc
defSolver.relativeTolerance = 1e-9
defSolver.absoluteTolerance = 1.e-30
defSolver.nMaxIterations = 6000 
defSolver.verbosity=1

elecSolver = fvmbaseExt.AMG()
elecSolver.relativeTolerance = 1e-3
elecSolver.nMaxIterations = 1000
elecSolver.maxCoarseLevels=20
elecSolver.verbosity=1

soptions = smodel.getOptions()
soptions.deformationLinearSolver = defSolver
soptions.deformationTolerance=1.0e-3
soptions.setVar("deformationURF",1.0)
soptions.printNormalizedResiduals=True
soptions.transient=True

mmmoptions = movingMeshModel.getOptions()
mmmoptions.nNodeDisplacementSweeps = 500000
mmmoptions.absTolerance = 1e-14
mmmoptions.relativeTolerance = 1e-9
mmmoptions.setVar('underrelaxation',0.4)

eoptions = emodel.getOptions()
eoptions.electrostaticsLinearSolver = elecSolver
eoptions.electrostaticsTolerance = 0.5e-5
eoptions.electrostatics = 1
eoptions.chargetransport = 0
eoptions.tunneling = 0
eoptions.ibm = 0
eoptions.transient_enable = False
eoptions.printNormalizedResiduals = True


metricsCalculator.calculateBoundaryNodeNormal()

numTimeSteps = 10000
period = 8.8043e-6
timeStep = period/1000
globalTime=0.

# set the timesteps
soptions.setVar('timeStep',timeStep)

"""
if fvm.atype=='tangent':
    vcMap = fmodel.getVCMap()
    for i,vc in vcMap.iteritems():
        print vc.getVar('viscosity')
        vc.setVar('viscosity',(1.7894e-5,1))
"""

for mesh in meshes0:
    fgs = mesh.getBoundaryGroups()
    for fg in fgs:
        bc = smodel.getBCMap()[fg.id]
        print '%i %s' %(fg.id,bc.bcType)

emodel.printBCs()
#import ddd

print '\n no of cells in mesh0 = %i' % (mesh0.getCells().getSelfCount())
print '\n no of cells in mesh1 = %i' % (mesh1.getCells().getSelfCount())

smodel.init()
dmodel.init()
emodel.init()

#set up permittivity
cells1 = mesh1.getCells()
perm = electricFields.dielectric_constant[cells1].asNumPyArray()
perm[:] = 1.0
                            
#smodel.advance(numIterations)
#dmodel.calculateNodeDisplacement()
#dmodel.deformStructure()
#setDirichletCommonDisplacement(dmodel,geomFields,meshes,structureFields)
#movingMeshModel.advance()
#metricsCalculator.recalculate()
#emodel.advance(numIterations)
eadvance(emodel,numEIterations)
#advance(smodel,dmodel,movingMeshModel,emodel,geomFields,
#        structureFields,electricFields,meshes0,meshes,numIterations)
advanceUnsteady(smodel,dmodel,movingMeshModel,emodel,geomFields,
                structureFields,electricFields,meshes0,meshes,numTimeSteps,globalTime)

#cells0 = mesh0.getCells()
#nCells0 = cells0.getCount()
#xc =  geomFields.coordinate[cells0].asNumPyArray()
#print '\n The mid point is %e %e' %(xc[500,0],xc[500,1])

#smodel.getTractionX(mesh0)
#smodel.getTractionX(mesh1)

#mnumber = 0
#for mesh in meshes:
#    mnumber = mnumber + 1
#    nodes = mesh.getNodes()
#    common = dmodel.getCommon(nodes).asNumPyArray()
#    length = len(common)
#    xf =  geomFields.coordinate[nodes].asNumPyArray()
#    for i in range(0,length):
#        id = common[i]
#        print '\n Mesh %i Node %i, x = %f, y = %f\n' %(mnumber,id,xf[id][0],xf[id][1])

fileName = fileBase1 + "dirichletNodeDisplacement.txt"
file = open(fileName,"w")
xf =  geomFields.coordinate[nodes1].asNumPyArray()
nnodes1 = nodes1.getCount()
doptions = geomFields.displacementOptions[nodes1].asNumPyArray()
dvar = geomFields.dirichletNodeDisplacement[nodes1].asNumPyArray()
for i in range(0,nnodes1):
    x = xf[i][0]
    y = xf[i][1]
    file.write(" %e " % x)
    file.write(" %e " % y)
    file.write(" %i " % doptions[i])
    file.write(" %e " % dvar[i][0])
    file.write(" %e " % dvar[i][1])
    file.write("\n")
file.close()

t1 = time.time()
print '\nsolution time = %f' % (t1-t0)

deformation =  structureFields.deformation[mesh0.getCells()].asNumPyArray()
fileName = fileBase0 + "deformation.txt"
file = open(fileName,"w")
file.write("deformation\t\n")
fgs = mesh0.getAllFaceGroups()
for fg in fgs:
    nFaces = fg.site.getCount()
    xf =  geomFields.coordinate[fg.site].asNumPyArray()
    if fg.id==3:
        faceCells = mesh0.getFaceCells(fg.site)
        for i in range(0,nFaces):
            x = xf[i][0]
            y = xf[i][1]
            def0 = deformation[faceCells(i,0)][0]
            def1 = deformation[faceCells(i,0)][1]
            def2 = deformation[faceCells(i,0)][2]
            file.write(" %e " % x)
            file.write(" %e " % y)
            file.write(" %e " % def0)
            file.write(" %e " % def1)
            file.write(" %e " % def2)
            file.write("\n")
file.close() 

dumpTecplotFile( nmesh, meshes, geomFields, options.type)

#fileName = fileBase0 + "volume0.txt"
#file = open(fileName,"w")
#cells0 = mesh0.getCells()
#xc =  geomFields.coordinate[cells0].asNumPyArray()
#nCells0 = cells0.getCount()
#vol0 = structureFields.volume0[cells0].asNumPyArray()
#for i in range(0,nCells0):
#    x = xc[i][0]
#    y = xc[i][1]
#    file.write(" %e " % x)
#    file.write(" %e " % y)
#    file.write(" %e " % vol0[i])
#    file.write("\n")
#file.close()

print '\n run complete '





