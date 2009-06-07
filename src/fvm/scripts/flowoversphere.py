#!/usr/bin/env python
import pdb
import sys
import math
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt
import importers

#pdb.set_trace()
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
numIterations = 3000
fileBase = "/home/lin/work/app-memosa/src/fvm/verification/flowoversphere/"


sqrt = math.sqrt
cos = math.cos
sin = math.sin
acos = math.acos
atan2 = math.atan2
pow = math.pow



def advance(fmodel,particles,niter):
    for i in range(0,niter):
        try:
            fmodel.computeIBFaceVelocity(particles)
            fmodel.advance(1)
            if (i >= niter*0.9):
                pIBotX = fmodel.getPressureIntegral(meshes[0],8)[0]
                pITopX = fmodel.getPressureIntegral(meshes[0],3)[0]
                pIBotY = fmodel.getPressureIntegral(meshes[0],8)[1]
                pITopY = fmodel.getPressureIntegral(meshes[0],3)[1]
                pIBotZ = fmodel.getPressureIntegral(meshes[0],8)[2]
                pITopZ = fmodel.getPressureIntegral(meshes[0],3)[2]

                pILeftX = fmodel.getPressureIntegral(meshes[0],6)[0]
                pIRightX = fmodel.getPressureIntegral(meshes[0],5)[0]
                pILeftY = fmodel.getPressureIntegral(meshes[0],6)[1]
                pIRightY = fmodel.getPressureIntegral(meshes[0],5)[1]
                pILeftZ = fmodel.getPressureIntegral(meshes[0],6)[2]
                pIRightZ = fmodel.getPressureIntegral(meshes[0],5)[2]

                pIFrontX = fmodel.getPressureIntegral(meshes[0],7)[0]
                pIBackX  = fmodel.getPressureIntegral(meshes[0],4)[0]
                pIFrontY = fmodel.getPressureIntegral(meshes[0],7)[1]
                pIBackY  = fmodel.getPressureIntegral(meshes[0],4)[1]
                pIFrontZ = fmodel.getPressureIntegral(meshes[0],7)[2]
                pIBackz  = fmodel.getPressureIntegral(meshes[0],4)[2]
                
                momIBotX = fmodel.getMomentumFluxIntegral(meshes[0],8)[0]
                momITopX = fmodel.getMomentumFluxIntegral(meshes[0],3)[0]
                momIBotY = fmodel.getMomentumFluxIntegral(meshes[0],8)[1]
                momITopY = fmodel.getMomentumFluxIntegral(meshes[0],3)[1]
                momIBotZ = fmodel.getMomentumFluxIntegral(meshes[0],8)[2]
                momITopZ = fmodel.getMomentumFluxIntegral(meshes[0],3)[2]

                momILeftX  = fmodel.getMomentumFluxIntegral(meshes[0],6)[0]
                momIRightX = fmodel.getMomentumFluxIntegral(meshes[0],5)[0]
                momILeftY  = fmodel.getMomentumFluxIntegral(meshes[0],6)[1]
                momIRightY = fmodel.getMomentumFluxIntegral(meshes[0],5)[1]
                momILeftZ  = fmodel.getMomentumFluxIntegral(meshes[0],6)[2]
                momIRightZ = fmodel.getMomentumFluxIntegral(meshes[0],5)[2]

                momIFrontX = fmodel.getMomentumFluxIntegral(meshes[0],7)[0]
                momIBackX  = fmodel.getMomentumFluxIntegral(meshes[0],4)[0]
                momIFrontY = fmodel.getMomentumFluxIntegral(meshes[0],7)[1]
                momIBackY  = fmodel.getMomentumFluxIntegral(meshes[0],4)[1]
                momIFrontZ = fmodel.getMomentumFluxIntegral(meshes[0],7)[2]
                momIBackZ  = fmodel.getMomentumFluxIntegral(meshes[0],4)[2]
        except KeyboardInterrupt:
            break

def initVelocity(geomFields,velocityField,meshes):


    mesh = meshes[0]
    vol = geomFields.volume[mesh.getCells()]

    for mesh in meshes:
        cells = mesh.getCells()
        xc = geomFields.coordinate[cells].asNumPyArray()
        vc = velocityField[cells].asNumPyArray()
                 
        vc[:,0] = 0.0
      
            
def createBVFields(geomFields,meshes):
    fx = fvmbaseExt.Field('bvx')
    fy = fvmbaseExt.Field('bvy')
    fz = fvmbaseExt.Field('bvz')

    mesh = meshes[0]
    vol = geomFields.volume[mesh.getCells()]
    
    for mesh in meshes:
        fgs = mesh.getBoundaryGroups()
        for fg in fgs:
            nFaces = fg.site.getCount()

            xvel = vol.newSizedClone(nFaces)
            yvel = vol.newSizedClone(nFaces)
            zvel = vol.newSizedClone(nFaces)

            xvela = xvel.asNumPyArray()
            yvela = yvel.asNumPyArray()
            zvela = zvel.asNumPyArray()
            
            xf = geomFields.coordinate[fg.site].asNumPyArray()

            
            a = 10.0
            U0 = 1.0
          
            
            for i in range(0,nFaces):
                x = xf[i][0]
                y = xf[i][1]
                z = xf[i][2]
                
                r = sqrt(x*x+y*y+z*z)
                alfa = acos(z/r)
                beta = atan2(y,x)

                Ur = U0 * cos(alfa) * (1-1.5*a/r+0.5*pow(a, 3)/pow(r,3))
                Ualfa = -U0 * sin(alfa) * (1-0.75*a/r-0.25*pow(a, 3)/pow(r,3))

                xvela[i] = Ur * sin(alfa) * cos(beta) + Ualfa * cos(alfa) * cos(beta)
                yvela[i] = Ur * sin(alfa) * sin(beta) + Ualfa * cos(alfa) * sin(beta)
                zvela[i] = Ur * cos(alfa) - Ualfa * sin(alfa)

            #pdb.set_trace()       

              
            fx[fg.site] = xvel
            fy[fg.site] = yvel
            fz[fg.site] = zvel
            
    return fx,fy,fz
            
            
outfile = fileBase+"-prism.dat"
    
reader = FluentCase(fileBase+"cube-125k-adapted-1.cas")

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

mesh0 = meshes[0]

solid = fvmbaseExt.MPM()

octree = fvmbaseExt.Octree()

octree.Impl(mesh0, geomFields)



mpmFileName = fileBase + "MPMs.dat"

#solid.setandwriteParticles(mpmFileName)

solid = fvmbaseExt.MPM(mpmFileName)

#solid.Impl(mpmFileName)

option = 2

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



fx,fy,fz = createBVFields(geomFields,meshes)


bcMap = fmodel.getBCMap()
for bc in bcMap.values():
    bc['specifiedXVelocity']=fx
    bc['specifiedYVelocity']=fy
    bc['specifiedZVelocity']=fz
    bc.bcType = 'VelocityBoundary'

fmodel.init()

initVelocity(geomFields,flowFields.velocity,meshes)
fmodel.computeIBFaceVelocity(particles)
advance(fmodel,particles,numIterations)

t1 = time.time()
if outfile != '/dev/stdout':
    print '\nsolution time = %f' % (t1-t0)

writer = exporters.FluentDataExporterA(reader,fileBase+"cube-125k-adapted-withconvection.dat",False,0)

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


