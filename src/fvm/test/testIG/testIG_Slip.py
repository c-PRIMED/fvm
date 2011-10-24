#!/usr/bin/env python

import sys

import fvm
import fvm.fvmbaseExt as fvmbaseExt
fvm.set_atype('double')
import fvm.importers as importers


from mpi4py  import MPI


from FluentCase import FluentCase

#fvmbaseExt.enableDebug("cdtor")

fileBase = None
numIterations = 50
#fileBase = "/home/shared/app-memosa/src/fvm/test/testIG_Slip"
#fileBase = "/home/shared/app-memosa/src/fvm/test/fine2"
#fileBase="/home/aerosun/a/schigull/FVM-trial/Slip_BC/test/Slip/Idealgas/Memosa/grid800/slipbc_igc_sruti"
#fileBase = "/home/sm/a/data/wj"
fileBase="./"

def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat if it is not specified."
    sys.exit(1)

def advance(models,dmodel,nstart,niter):
    for i in range(nstart,niter):
        try:
            converged = True
            for m in models:
                if not m.advance(1):
                    converged = False
            if (i%10) == 0:
                dmodel.advance(1)
            if converged:
                break
        except KeyboardInterrupt:
            break
        if ((i+1)%1000)==0:
            saveData(flowFields,reader,fileBase,i)


def saveData(flowFields,reader,fileBase,i):
    writer = fvm.exporters.FluentDataExporterA(reader,fileBase+str(i+1)+"-noslipnew.dat",False,0)
    writer.init()
    writer.writeScalarField(flowFields.pressure,1)
    writer.writeScalarField(flowFields.density,101)
    writer.writeVectorField(flowFields.velocity,111)
    writer.writeScalarField(flowFields.massFlux,18)
    writer.finish()
        
def advancenew(models,dmodel,nstart,niter):
    for i in range(nstart,niter):
        try:
            for m in models:
                m.advance(1)
                    
            #if (i%25)== 0:
            #    dmodel.advance(1)
            
        except KeyboardInterrupt:
            break
        if ((i+1)%1000)==0:
            saveData(flowFields,reader,fileBase,i)
                    

outfile = None
if __name__ == '__main__' and fileBase is None:
    if len(sys.argv) < 2:
        usage()
    fileBase = sys.argv[1]
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

if outfile == None:
    outfile = fileBase+"-prism.dat"
    
reader = FluentCase(sys.argv[1])

#import debug
#import ddd
reader.read();

meshes = reader.getMeshList()

import time
t0 = time.time()

geomFields =  fvm.models.GeomFields('geom')
metricsCalculator = fvm.models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()


flowFields =  fvm.models.FlowFields('flow')

fmodel = fvm.models.FlowModelA(geomFields,flowFields,meshes)


reader.importFlowBCs(fmodel,meshes)


## set the walls to be slip jump
for i,bc in fmodel.getBCMap().iteritems():
    if bc.bcType == 'NoSlipWall':
        bc.bcType = 'SlipJump'
        bc.setVar('accomodationCoefficient', 0.0)

#calculation of KnudsenNumber*H=meanfreepath
#meanfreepath=mu_N2/P_N2*(pi*R_N2*T_N2*0.5)**0.5

momSolver = fvmbaseExt.AMG()
momSolver.relativeTolerance = 1e-3
momSolver.nMaxIterations = 20
momSolver.maxCoarseLevels=20
momSolver.verbosity=0

contSolver = fvmbaseExt.AMG()
#pc = fvmbaseExt.AMG()
#pc.verbosity=0
#contSolver = fvmbaseExt.BCGStab()
#contSolver.preconditioner = pc
contSolver.relativeTolerance = 1e-3
contSolver.nMaxIterations = 20
contSolver.verbosity=0
contSolver.maxCoarseLevels=20

foptions = fmodel.getOptions()

foptions.momentumLinearSolver = momSolver
foptions.pressureLinearSolver = contSolver

foptions.momentumTolerance=1e-5
foptions.continuityTolerance=1e-5
foptions.printNormalizedResiduals=False

foptions.setVar("momentumURF",0.7)
foptions.setVar("pressureURF",0.3)
foptions.setVar("initialXVelocity",6)

foptions.incompressible = False
foptions.setVar('operatingPressure', 101325.0)
foptions.setVar('operatingTemperature', 300.0)
foptions.setVar('molecularWeight',28.9645)

dmodel = fvm.models.IdealGasDensityModelA(geomFields,flowFields,meshes)

## set density  model settings
for i,vc in dmodel.getVCMap().iteritems():

    ## use local pressure + op pressure 
    vc['pressure'] = flowFields.pressure
    # use constant temperature
    vc['temperature']= foptions.getVar('operatingTemperature')
    vc['operatingPressure'] = foptions.getVar('operatingPressure')
    vc['molecularWeight']= foptions.getVar('molecularWeight')
    vc['urf']= 0.5

"""
if atype=='tangent':
    vcMap = fmodel.getVCMap()
    for i,vc in vcMap.iteritems():
        print vc.getVar('viscosity')
        vc.setVar('viscosity',(1.7894e-5,1))
"""
#import debug

fmodel.init()
dmodel.init()
initIterations=0
#advancenew([fmodel],dmodel,0,initIterations)
contSolver.redirectPrintToFile("convergence.dat") 
advance([fmodel],dmodel,initIterations,numIterations)
contSolver.redirectPrintToScreen()
#advance([fmodel],dmodel,numIterations)-old

t1 = time.time()
if outfile != '/dev/stdout':
    print '\nsolution time = %f' % (t1-t0)

writer = fvm.exporters.FluentDataExporterA(reader,fileBase+"-remove.dat",False,0)

writer.init()
writer.writeScalarField(flowFields.pressure,1)
writer.writeScalarField(flowFields.density,101)
writer.writeVectorField(flowFields.velocity,111)
writer.writeScalarField(flowFields.massFlux,18)
writer.finish()


saveData(flowFields,reader,fileBase,1)
