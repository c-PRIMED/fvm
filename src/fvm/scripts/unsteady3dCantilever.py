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

#fvmbaseExt.enableDebug("cdtor")
 
fileBase = None
numIterationsPerStep = 100
fileBase = "/home/sm/prism-meshes/3dbeam/mode2/new-beam-114k-p=83593"
#fileBase = "/home/sm/a/data/wj"

vFile = open(fileBase + "-prism-v.xy","w")
ptopFile = open(fileBase + "-prism-pIntegral-top.xy","w")
pbotFile = open(fileBase + "-prism-pIntegral-bot.xy","w")
psumFile = open(fileBase + "-prism-pIntegral.xy","w")
tautopFile = open(fileBase + "-prism-tauIntegral-top.xy","w")
taubotFile = open(fileBase + "-prism-tauIntegral-bot.xy","w")
tausumFile = open(fileBase + "-prism-tauIntegral.xy","w")

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
    
frequency = 759734
timeStep = 1.0/(frequency*100.)

#timeStep = 5.0e-08

numTimeSteps = 600

globalTime=0

sideID = 10
topID = 9
botID=12
def advanceUnsteady(fmodel,meshes,globalTime,nTimeSteps):
    bcMap = fmodel.getBCMap()
    bcSide = bcMap[sideID]
    bcBot = bcMap[botID]
    bcTop = bcMap[topID]

    for i in range(0,nTimeSteps):
        v = 0.1*math.cos(2.0*math.pi*frequency*globalTime)
        bcSide.setVar('specifiedZVelocity',v)
        bcBot.setVar('specifiedZVelocity',v)
        bcTop.setVar('specifiedZVelocity',v)

        advance(fmodel,numIterationsPerStep)

        vFile.write("%e %e\n" % (globalTime,v))
        pIBot = fmodel.getPressureIntegral(meshes[0],botID)[2]
        pITop = fmodel.getPressureIntegral(meshes[0],topID)[2]
        pISide = fmodel.getPressureIntegral(meshes[0],sideID)[2]
        tauIBot = fmodel.getMomentumFluxIntegral(meshes[0],botID)[2]
        tauITop = fmodel.getMomentumFluxIntegral(meshes[0],topID)[2]
        tauISide = fmodel.getMomentumFluxIntegral(meshes[0],sideID)[2]
        pISum = pIBot+pITop+tauIBot+tauITop+pISide+tauISide
        ptopFile.write("%e %e\n" % (globalTime,pITop))
        pbotFile.write("%e %e\n" % (globalTime,pIBot))
        tautopFile.write("%e %e\n" % (globalTime,tauITop))
        taubotFile.write("%e %e\n" % (globalTime,tauIBot))
        psumFile.write("%e %e\n" % (globalTime,pISum))

        globalTime += timeStep
        print 'advancing to time %e' % globalTime
        fmodel.updateTime()
    
# change as needed

if __name__ == '__main__' and fileBase is None:
    if len(sys.argv) != 2:
        usage()

    fileBase = sys.argv[1]


reader = FluentCase(fileBase+".cas")

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

pcc = fvmbaseExt.AMG()
pcc.verbosity=0
cSolver = fvmbaseExt.BCGStab()
cSolver.preconditioner = pcc
#cSolver=pc
pcc.maxCoarseLevels=2
cSolver.relativeTolerance = 1e-2
cSolver.nMaxIterations = 20
cSolver.verbosity=0

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

fmodel.init()
#fmodel.advance(numIterations)
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
tautopFile.close()
taubotFile.close()
tausumFile.close()
