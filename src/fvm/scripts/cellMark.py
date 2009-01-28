#!/usr/bin/env python
import pdb
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


#fileBase = None
numIterations = 1
fileBase = "/home/linsun/Work/prism/app-memosa/src/fvm/test/cav32"



def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat if it is not specified."
    sys.exit(1)

def advance(fmodel,niter):
    for i in range(0,niter):
        try:
            fmodel.advance(1)
        except KeyboardInterrupt:
            break

# change as needed

outfile = None
if __name__ == '__main__' and fileBase is None:
    if len(sys.argv) < 2:
        usage()
    fileBase = sys.argv[1]
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

if outfile == None:
    outfile = fileBase+"-ibm.dat"
    
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

foptions.momentumLinearSolver = momSolver
foptions.pressureLinearSolver = contSolver

foptions.momentumTolerance=1e-3
foptions.continuityTolerance=1e-3
foptions.setVar("momentumURF",0.95)
foptions.setVar("pressureURF",0.1)
foptions.printNormalizedResiduals=False

"""
if atype=='tangent':
    vcMap = fmodel.getVCMap()
    for i,vc in vcMap.iteritems():
        print vc.getVar('viscosity')
        vc.setVar('viscosity',(1.7894e-5,1))
"""

fmodel.init()
#fmodel.advance(numIterations)


###############
nCell=mesh0.getNumberOfCells()
print 'nCell = %i' % nCell
#nCell=1024

#pdb.set_trace()
centerCoord0=0.5
centerCoord1=0.5
centerCoord2=0.
radius=0.2
delta=0.05

for i in range(0, nCell):
    cellCoord=mesh0.getCellCoordinate(i)
    rSqr=(centerCoord0-cellCoord[0])**2+(centerCoord1-cellCoord[1])**2+(centerCoord2-cellCoord[2])**2
    if rSqr < radius**2:
        #print i, rSqr
	mesh0.setIBTypeForCell(i,2)
       
    elif rSqr > (radius+delta)**2:
        #print i, rSqr
        mesh0.setIBTypeForCell(i,0)
                   
    else:
        #print i, rSqr
        mesh0.setIBTypeForCell(i,1)
       

#file=open('./coordCell','w')
#file0=open('./fluidCell','w')
#file1=open('./IBCell','w')
#file2=open('./solidCell','w') 
for i in range(0, nCell):
    cellIbType=mesh0.getIBTypeForCell(i)
    cellCoord=mesh0.getCellCoordinate(i)
    #file.write('%f\t%f\t%f\n' % (cellCoord[0],cellCoord[1],cellCoord[2]))
    if cellIbType == 2:
       
        #file2.write('%f\t%f\t%f\n' % (cellCoord[0],cellCoord[1],cellCoord[2]))
    elif cellIbType == 1:
        print 'cellID %i ' % i 
        print 'cellType %i\n' % cellIbType
        #file1.write('%f\t%f\t%f\n' % (cellCoord[0],cellCoord[1],cellCoord[2]))
    else:
        #file0.write('%f\t%f\t%f\n' % (cellCoord[0],cellCoord[1],cellCoord[2]))
#file0.close
#file1.close
#file2.close
#file.close

#advance(fmodel,numIterations)    

t1 = time.time()
if outfile != '/dev/stdout':
    print '\nsolution time = %f' % (t1-t0)

print 'solution time = %f' % (t1-t0)

print '\n\npressure integrals\n'
fmodel.printPressureIntegrals()

print '\n\nmomentum flux integrals\n'
fmodel.printMomentumFluxIntegrals()


writer = exporters.FluentDataExporterA(reader,fileBase+"-ibm.dat",False,0)

writer.init()
writer.writeScalarField(flowFields.pressure,1)
writer.writeVectorField(flowFields.velocity,111)
writer.writeScalarField(flowFields.massFlux,18)
writer.finish()

if (atype=='tangent'):
    writer = exporters.FluentDataExporterA(reader,fileBase+"-ibm-tangent.dat",False,1)
    writer.init()
    writer.writeScalarField(flowFields.pressure,1)
    writer.writeVectorField(flowFields.velocity,111)
    writer.writeScalarField(flowFields.massFlux,18)
    writer.finish()

    
