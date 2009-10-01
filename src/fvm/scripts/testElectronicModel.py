#!/usr/bin/env python
import pdb
import sys
import math
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt
import importers
import models_atyped_double as models
import exporters_atyped_double as exporters
from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")



numIterationsPerStep =20
numTimeSteps = 3
timeStep = 1.0
globalTime = 0.0

fileBase = "/home/lin/work/app-memosa/src/fvm/verification/ElectroStatics/"

def advance(fmodel,niter):
    for i in range(0,niter):
        try:
            stopFlag=fmodel.advance(1)
            if stopFlag == 1:
                break
        except KeyboardInterrupt:
            break


def advanceUnsteady(fmodel, globalTime, nTimeSteps):
    for i in range(0, nTimeSteps):
        advance(fmodel, numIterationsPerStep)
        globalTime += timeStep
        fmodel.updateTime()
        


        
reader = FluentCase(fileBase + "cav32_elec.cas")

reader.read()

#pdb.set_trace()

meshes = reader.getMeshList()

import time
t0 = time.time()

geomFields =  models.GeomFields('geom')

metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

electronicFields =  models.ElectronicFields('elec')



elec_model = models.ElectronicModelA(geomFields,electronicFields,meshes)

#reader.importElectronicModelBCs(elec_model)

bcMap = elec_model.getBCMap()

#top
bcID = 6
if bcID in bcMap:
    bc = elec_model.getBCMap()[bcID]
    bc.bcType = 'SpecifiedPotential'
    bc.setVar('specifiedPotential',300)
    #bc.bcType = 'SpecifiedCharge'
    #bc.setVar('specifiedCharge',400)
#bot
bcID = 3
if bcID in bcMap:
    bc = elec_model.getBCMap()[bcID]
    bc.bcType = 'SpecifiedPotential'
    bc.setVar('specifiedPotential',400)
    #bc.bcType = 'SpecifiedCharge'
    #bc.setVar('specifiedCharge',300)
#left
bcID = 4
if bcID in bcMap:
    bc = elec_model.getBCMap()[bcID]
    bc.bcType = 'SpecifiedPotential'
    bc.setVar('specifiedPotential',350)
    #bc.bcType = 'SpecifiedPotentialFlux'
    #bc.setVar('specifiedPotentialFlux',0)
    #bc.bcType = 'SpecifiedCharge'
    #bc.setVar('specifiedCharge',350)
#right
bcID = 5
if bcID in bcMap:
    bc = elec_model.getBCMap()[bcID]
    bc.bcType = 'SpecifiedPotential'
    bc.setVar('specifiedPotential',350)
    #bc.bcType = 'SpecifiedPotentialFlux'
    #bc.setVar('specifiedPotentialFlux',0)
    #bc.bcType = 'SpecifiedCharge'
    #bc.setVar('specifiedCharge',350)
    

elec_options = elec_model.getOptions()
elec_options.electrostatics = 1
elec_options.chargetransport = 0
elec_options.tunneling = 1
elec_options.ibm = 0

elec_model.printBCs()
import ddd
elec_model.init()

#set up permittivity
cells = meshes[0].getCells()
perm = electronicFields.permittivity[cells].asNumPyArray()
perm[:] = 1.0


advanceUnsteady(elec_model, globalTime, numTimeSteps)

t1 = time.time()

print 'solution time = %f' % (t1-t0)

writer = exporters.FluentDataExporterA(reader,fileBase+"cav32_test.dat",False,0)
writer.init()
writer.writeScalarField(electronicFields.charge,3)
writer.finish()
