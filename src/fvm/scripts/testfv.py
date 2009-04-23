#!/usr/bin/env python

import sys
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt
import importers

atype = 'double'
#atype = 'tangent'

if atype == 'double':
    import models_atyped_double as models
    import exporters_atyped_double as exporters
elif atype == 'tangent':
    import models_atyped_tangent_double as models
    import exporters_atyped_tangent_double as exporters


from FluentCase import FluentCase

#fvmbaseExt.enableDebug("cdtor")

fileBase = None
numIterations = 200
fileBase = "/home/sm/app-memosa/src/fvm/test/"



def advance(fmodel,particles,niter):
    for i in range(0,niter):
        try:
            fmodel.computeIBFaceVelocity(particles)
            fmodel.advance(1)
        except KeyboardInterrupt:
            break

def initVelocity(geomFields,velocityField,meshes):


    mesh = meshes[0]
    vol = geomFields.volume[mesh.getCells()]

    for mesh in meshes:
        cells = mesh.getCells()
        xc = geomFields.coordinate[cells].asNumPyArray()
        vc = velocityField[cells].asNumPyArray()
        rx = xc[:,0] - 0.0
        ry = xc[:,1] - 0.0
        vx = -ry
        vy = rx
            
        vc[:,0] = vx[:]
        vc[:,1] = vy[:]
            
def createBVFields(geomFields,meshes):
    fx = fvmbaseExt.Field('bvx')
    fy = fvmbaseExt.Field('bvy')

    mesh = meshes[0]
    vol = geomFields.volume[mesh.getCells()]

    for mesh in meshes:
        fgs = mesh.getBoundaryGroups()
        for fg in fgs:
            nFaces = fg.site.getCount()
            xf = geomFields.coordinate[fg.site].asNumPyArray()

            rx = xf[:,0] - 0.0
            ry = xf[:,1] - 0.0
            vx = -ry
            vy = rx
            
            xvel = vol.newSizedClone(nFaces)
            yvel = vol.newSizedClone(nFaces)

            xvela = xvel.asNumPyArray()
            yvela = yvel.asNumPyArray()

            xvela[:] = 0#vx[:]
            yvela[:] = 0#vy[:]

            fx[fg.site] = xvel
            fy[fg.site] = yvel

    return fx,fy
            
            
outfile = fileBase+"-prism.dat"
    
reader = FluentCase(fileBase+"cav32-new.cas")

#import debug
reader.read();

meshes = reader.getMeshList()

import time
t0 = time.time()

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

if atype == 'tangent':
    metricsCalculator.setTangentCoords(0,7,1)

flowFields =  models.FlowFields('flow')

fmodel = models.FlowModelA(geomFields,flowFields,meshes)

reader.importFlowBCs(fmodel)
#fmodel.printBCs()

momSolver = fvmbaseExt.AMG()
momSolver.relativeTolerance = 0.1
momSolver.nMaxIterations = 200
momSolver.maxCoarseLevels=20
momSolver.verbosity=0

contSolver = fvmbaseExt.AMG()
#pc = fvmbaseExt.AMG()
#pc.verbosity=0
#contSolver = fvmbaseExt.BCGStab()
#contSolver.preconditioner = pc
contSolver.relativeTolerance = 1e-1
contSolver.nMaxIterations = 20
contSolver.verbosity=0
contSolver.maxCoarseLevels=20

foptions = fmodel.getOptions()

foptions.momentumLinearSolver = momSolver
foptions.pressureLinearSolver = contSolver

foptions.momentumTolerance=1e-3
foptions.continuityTolerance=1e-3
foptions.setVar("momentumURF",0.7)
foptions.setVar("pressureURF",0.3)
foptions.printNormalizedResiduals=False
foptions.correctVelocity = True

"""
if atype=='tangent':
    vcMap = fmodel.getVCMap()
    for i,vc in vcMap.iteritems():
        print vc.getVar('viscosity')
        vc.setVar('viscosity',(1.7894e-5,1))
"""



solid = fvmbaseExt.MPM()

octree = fvmbaseExt.Octree() 

option = 2
mesh0 = meshes[0]
fvmbaseExt.CellMark_Impl(mesh0, geomFields, fileBase, octree, solid, option)

cells = mesh0.getCells()

nCells = cells.getCount()

# add particles position and velocities to the main fields

particles = solid.getParticles()
px = solid.getCoordinates()
pv = solid.getVelocities()
geomFields.coordinate[particles]=px
flowFields.velocity[particles]=pv

metricsCalculator.computeIBInterpolationMatrices(particles)



fx,fy = createBVFields(geomFields,meshes)


bcMap = fmodel.getBCMap()
for bc in bcMap.values():
    bc['specifiedXVelocity']=fx
    bc['specifiedYVelocity']=fy
    bc.bcType = 'VelocityBoundary'

fmodel.init()

initVelocity(geomFields,flowFields.velocity,meshes)
fmodel.computeIBFaceVelocity(particles)
advance(fmodel,particles,numIterations)

t1 = time.time()
if outfile != '/dev/stdout':
    print '\nsolution time = %f' % (t1-t0)

writer = exporters.FluentDataExporterA(reader,fileBase+"-prism.dat",False,0)

writer.init()
writer.writeScalarField(flowFields.pressure,1)
writer.writeVectorField(flowFields.velocity,111)
writer.writeScalarField(flowFields.massFlux,18)
writer.finish()

if (atype=='tangent'):
    writer = exporters.FluentDataExporterA(reader,fileBase+"-prism-tangent.dat",False,1)
    writer.init()
    writer.writeScalarField(flowFields.pressure,1)
    writer.writeVectorField(flowFields.velocity,111)
    writer.writeScalarField(flowFields.massFlux,18)
    writer.finish()


