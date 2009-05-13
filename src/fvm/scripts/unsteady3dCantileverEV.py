  #!/usr/bin/env python
 
import sys
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt
import importers
import math

atype = 'double'
#atype = 'tangent'

if atype == 'double':
    import models_atyped_double as models
    import exporters_atyped_double as exporters
elif atype == 'tangent':
    import models_atyped_tangent_double as models
    import exporters_atyped_tangent_double as exporters
    
    
from FluentCase import FluentCase

 

numIterationsPerStep = 100
fileBase = "/home/sm/prism-meshes/3dbeam/mode2EV/new-beam-114k-p=22531"
coordFile = "/home/sm/app-memosa/src/fvm/test/Grid/Grid_Coord.dat"
velocityFile = "/home/sm/app-memosa/src/fvm/test/Grid/Grid_Velocity_Mode2.dat"

## create a Field for the velocity coeffs which will be interpolated
## from the experimental measurements and one for the interpolated
## velocities

velCoeffField = fvmbaseExt.Field('velocityCoeffs')
buField = fvmbaseExt.Field('bu')
bvField = fvmbaseExt.Field('bv')
bwField = fvmbaseExt.Field('bw')


vFile = open(fileBase + "-prism-v.xy","w")
ptopFile = open(fileBase + "-prism-pIntegral-top.xy","w")
pbotFile = open(fileBase + "-prism-pIntegral-bot.xy","w")
psumFile = open(fileBase + "-prism-pIntegral.xy","w")


def advance(fmodel,niter):
    fmodel.advance(niter)
    """
    for i in range(0,niter):
        try:
            fmodel.advance(1)
        except KeyboardInterrupt:
            break
    """
    
frequency = 759734
#frequency = 114415
timeStep = 1.0/(frequency*100.)

#timeStep = 5.0e-08

numTimeSteps = 400

globalTime=0

sideID = 10
topID = 9
botID=12

def advanceUnsteady(fmodel,meshes,globalTime,nTimeSteps):
    bcMap = fmodel.getBCMap()

    for i in range(0,nTimeSteps):
        coswt = math.cos(2.0*math.pi*frequency*globalTime)

        for mesh in meshes:

            fgs = mesh.getBoundaryGroups()
            for fg in fgs:
                if fg.id in [sideID, topID, botID]:
                    faces = fg.site
                    bc = bcMap[fg.id]

                    uvel = buField[faces].asNumPyArray()
                    vvel = bvField[faces].asNumPyArray()
                    wvel = bwField[faces].asNumPyArray()

                    vCoeffs = velCoeffField[faces].asNumPyArray()

                    uvel[:] = coswt*vCoeffs[:,0]
                    vvel[:] = coswt*vCoeffs[:,1]
                    wvel[:] = coswt*vCoeffs[:,2]
                    
                    bc['specifiedXVelocity'] = buField
                    bc['specifiedYVelocity'] = bvField
                    bc['specifiedZVelocity'] = bwField


        advance(fmodel,numIterationsPerStep)

        vFile.write("%e %e\n" % (globalTime,coswt))
        pIBot = fmodel.getPVIntegral(velCoeffField,meshes[0],botID)[2]
        pITop = fmodel.getPVIntegral(velCoeffField,meshes[0],topID)[2]
        pISide = fmodel.getPVIntegral(velCoeffField,meshes[0],sideID)[2]
        pISum = pIBot+pITop+pISide
        ptopFile.write("%e %e\n" % (globalTime,pITop))
        pbotFile.write("%e %e\n" % (globalTime,pIBot))
        psumFile.write("%e %e\n" % (globalTime,pISum))

        globalTime += timeStep
        print 'advancing to time %e' % globalTime
        fmodel.updateTime()
    
# change as needed


reader = FluentCase(fileBase+".cas")

#import debug
reader.read();


meshes = reader.getMeshList()
mesh0 = meshes[0]

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


momSolver = fvmbaseExt.AMG()
momSolver.relativeTolerance = 1e-1
#momSolver.nMaxIterations = 20
#momSolver.maxCoarseLevels=20
momSolver.verbosity=0

#contSolver = fvmbaseExt.AMG()
pc = fvmbaseExt.AMG()
pc.verbosity=0
contSolver = fvmbaseExt.BCGStab()
contSolver.preconditioner = pc
contSolver.relativeTolerance = 1e-1
contSolver.nMaxIterations = 20
contSolver.verbosity=0
#contSolver.maxCoarseLevels=20

foptions = fmodel.getOptions()
foptions.transient = True
foptions.setVar('timeStep',timeStep)

foptions.momentumLinearSolver = momSolver
foptions.pressureLinearSolver = contSolver


foptions.momentumTolerance=1e-5
foptions.continuityTolerance=1e-5
foptions.setVar("momentumURF",0.9)
foptions.setVar("pressureURF",0.1)

foptions.printNormalizedResiduals=True


## the interpolator

grid = fvmbaseExt.Grid(geomFields, flowFields, coordFile, velocityFile)
gridNodes = grid.getNodes()



## interpolate velocity coeffs

for mesh in meshes:

    ## need this to create a scalar array by cloning
    vol = geomFields.volume[mesh.getCells()]

    fgs = mesh.getBoundaryGroups()

    
    for fg in fgs:
        if fg.id in [sideID, topID, botID]:
            faces = fg.site
            grid.setConnFaceToGrid(mesh, faces)
            metricsCalculator.computeGridInterpolationMatrices(gridNodes,faces)
            faceVel = grid.computeInterpolatedVelocity(faces)

            velCoeffField[faces]=faceVel

            nFaces = fg.site.getCount()
            
            uvel = vol.newSizedClone(nFaces)
            vvel = vol.newSizedClone(nFaces)
            wvel = vol.newSizedClone(nFaces)

            buField[faces] = uvel
            bvField[faces] = vvel
            bwField[faces] = wvel
            
            
fmodel.init()

advanceUnsteady(fmodel,meshes,globalTime,numTimeSteps)

t1 = time.time()

print 'solution time = %f' % (t1-t0)

print '\n\npressure integrals\n'
fmodel.printPressureIntegrals()

print '\n\nmomentum flux integrals\n'
fmodel.printMomentumFluxIntegrals()


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

vFile.close()
ptopFile.close()
pbotFile.close()
psumFile.close()
