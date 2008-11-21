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
numIterations = 1
fileBase = "/home/sm/prism-meshes/case1/case1"
#fileBase = "/home/sm/prism-meshes/cav256"
#fileBase = "/home/sm/a/data/wj"

def usage():
    print "Usage: %s filebase" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat"
    sys.exit(1)

def advance(fmodel,niter):
    for i in range(0,niter):
        try:
            fmodel.advanceCoupled(1)
        except KeyboardInterrupt:
            break

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

pc = fvmbaseExt.AMG()
pc.verbosity=0
cSolver = fvmbaseExt.BCGStab()
cSolver.preconditioner = pc
#cSolver=pc
pc.maxCoarseLevels=0
cSolver.relativeTolerance = 1e-6
cSolver.nMaxIterations = 20
cSolver.verbosity=1

foptions = fmodel.getOptions()

foptions.coupledLinearSolver = cSolver

foptions.momentumTolerance=1e-3
foptions.continuityTolerance=1e-3
foptions.setVar("momentumURF",0.999)
foptions.setVar("velocityURF",0.99)
foptions.setVar("pressureURF",0.99)
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
advance(fmodel,numIterations)

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
