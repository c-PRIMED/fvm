#!/usr/bin/env python
# encoding: utf-8
#
# Special version of testFlowModel py
# used by automated tests. Do not change unless
# you update the expected output in the test directory.
# Martin Hunt <mmh@purdue.edu>

"""
Test FVM Flow Model

Usage: FvmTestFlowModel [options] filebase [outfilename]

Where 'filebase' is a Fluent case file and
'outfilename' defaults to 'filebase'.dat.

Options:
  --iterations n   Set number of iterations to 'n' [default 10]
  --tangent        Write tangents.
"""
import sys, time
from optparse import OptionParser
sys.setdlopenflags(0x100|0x2)
import fvmbaseExt
from FluentCase import FluentCase

def usage():
    print __doc__
    sys.exit(1)

def advance(fmodel,niter):
    for i in range(0,niter):
        try:
            fmodel.advance(1)
        except KeyboardInterrupt:
            break


parser = OptionParser()
parser.set_defaults(iterations=10,tangent=False)
parser.add_option("--tangent",
                  action="store_true", dest="tangent",
                  help="Write tangents.")
parser.add_option("--iterations", type="int",
                  dest="iterations",
                  help="Iterations.")

if len(sys.argv) < 2:
    usage()

fileBase = sys.argv[1]
if len(sys.argv) == 3:
    outfile = sys.argv[2]
else:
    outfile = fileBase+".dat"

(options, args) = parser.parse_args()

if options.tangent:
    atype = 'tangent'
    print "Tangent not implemented yet."
    # these don't seem to exist
    #import models_atyped_tangent_double as models
    #import exporters_atyped_tangent_double as exporters
    sys.exit(1)
else:
    atype = 'double'
    import models_atyped_double as models
    import exporters_atyped_double as exporters

#fvmbaseExt.enableDebug("cdtor")

print "reading",fileBase+".cas"
reader = FluentCase(fileBase+".cas")
reader.read();
meshes = reader.getMeshList()

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
momSolver.nMaxIterations = 20
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


#if atype=='tangent':
#    vcMap = fmodel.getVCMap()
#    for i,vc in vcMap.iteritems():
#        print vc.getVar('viscosity')
#        vc.setVar('viscosity',(1.7894e-5,1))


fmodel.init()
#fmodel.advance(options.iterations)
advance(fmodel,options.iterations)

t1 = time.time()
if outfile != '/dev/stdout':
    print '\nsolution time = %f' % (t1-t0)

#print '\n\npressure integrals\n'
#fmodel.printPressureIntegrals()

#print '\n\nmomentum flux integrals\n'
#fmodel.printMomentumFluxIntegrals()


writer = exporters.FluentDataExporterA(reader,outfile,False,0)
writer.init()
writer.writeScalarField(flowFields.pressure,1)
writer.writeVectorField(flowFields.velocity,111)
writer.writeScalarField(flowFields.massFlux,18)
writer.finish()

if atype == 'tangent':
    writer = exporters.FluentDataExporterA(reader,outfile,False,1)
    writer.init()
    writer.writeScalarField(flowFields.pressure,1)
    writer.writeVectorField(flowFields.velocity,111)
    writer.writeScalarField(flowFields.massFlux,18)
    writer.finish()

    
