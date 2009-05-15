#!/usr/bin/env python

import sys
import pdb
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

#fvmbaseExt.enableDebug("cdtor")

fileBase = None
numIterationsPerStep = 100
numTimeSteps = 200
fileBase = "/home/lin/work/app-memosa/src/fvm/test/2-D-Cantilever/"

vFile = open(fileBase + "velocity-fullbeam.out","w")
pFile = open(fileBase + "pIntegral-fullbeam.out","w")
#mFile = open(fileBase + "test-momIntegral.out","w")

def usage():
    print "Usage: %s filebase" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat"
    sys.exit(1)

def advance(fmodel,niter):
    fmodel.advance(niter)
    """
    for i in range(0,niter):
        try:
            fmodel.advance(1)
        except KeyboardInterrupt:
            break
    """
    
frequency = 114415
#timeStep = 1.0/(frequency*200.)
timeStep = 5.0e-08
globalTime=0

def advanceUnsteady(fmodel,meshes,globalTime,nTimeSteps, particles):    
    #pdb.set_trace()
    px = geomFields.coordinate[particles]
    pv = flowFields.velocity[particles]
    pxArray = px.asNumPyArray()
    pvArray = pv.asNumPyArray()    
    for i in range(0,nTimeSteps):
           
        v = 0.1*math.cos(2.0*math.pi*frequency*globalTime)
        pvArray[:,0] = 0
        pvArray[:,1] = v
        pvArray[:,2] = 0
             
        fmodel.computeIBFaceVelocity(particles)

        #fmodel.computeIBandSolidVelocity(particles)

        advance(fmodel,numIterationsPerStep)

        pI = fmodel.getPressureIntegralonIBFaces(mesh0)
        #mI = fmodel.getMomentumFluxIntegralonIBFaces(mesh0)
        
        #pxArray[:,1] += v * timeStep

        #solid.setVelocities(pv)

        #solid.setCoordinates(px)

        #fvmbaseExt.CellMark_Impl(mesh0, geomFields, fileBase, octree, solid, option)
       # metricsCalculator.computeIBInterpolationMatrices(particles)
       # metricsCalculator.computeIBandSolidInterpolationMatrices(particles)

        vFile.write("%e\t%e\n" % (globalTime,v))
        pFile.write("%e\t%e\t%e\t%e\n" % (globalTime, pI[0], pI[1], pI[2]))
        #mFile.write("%e\t%e\t%e\t%e\n" % (globalTime, mI[0], mI[1], mI[2]))
        
        globalTime += timeStep
        print 'advancing to time %e' % globalTime
        fmodel.updateTime()
    
# change as needed

if __name__ == '__main__' and fileBase is None:
    if len(sys.argv) != 2:
        usage()

    fileBase = sys.argv[1]


reader = FluentCase(fileBase+"fullbeam.cas")

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
#fmodel.printBCs()

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

foptions.momentumTolerance=1e-6
foptions.continuityTolerance=1e-6
foptions.setVar("momentumURF",0.9)
foptions.setVar("pressureURF",0.1)
foptions.printNormalizedResiduals=True

"""
if atype=='tangent':
    vcMap = fmodel.getVCMap()
    for i,vc in vcMap.iteritems():
        print vc.getVar('viscosity')
        vc.setVar('viscosity',(1.7894e-5,1))
"""
#import ddd
fmodel.init()

octree = fvmbaseExt.Octree()

octree.Impl(mesh0, geomFields)



mpmFileName = fileBase + "MPMs-fullbeam.dat"

solid = fvmbaseExt.MPM(mpmFileName)

option = 1

fvmbaseExt.CellMark_Impl(mesh0, geomFields, fileBase, octree, solid, option)


particles = solid.getParticles()
px = solid.getCoordinates()
pv = solid.getVelocities()
geomFields.coordinate[particles]=px
flowFields.velocity[particles]=pv
metricsCalculator.computeIBInterpolationMatrices(particles)
metricsCalculator.computeIBandSolidInterpolationMatrices(particles)
advanceUnsteady(fmodel,meshes,globalTime,numTimeSteps, particles)

t1 = time.time()

print 'solution time = %f' % (t1-t0)


writer = exporters.FluentDataExporterA(reader,fileBase+"P=83593-fullbeam.dat",False,0)
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
pFile.close()
#mFile.close()
