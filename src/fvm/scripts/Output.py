import os
import pdb
import fvm.models_atyped_double as models
import fvm.exporters_atyped_double as exporters


class Output():
    
    def __init__(self, outputDir, probeIndex, sim):
        if os.path.isdir(outputDir) == False:
            os.mkdir(outputDir)
        self.defFile = open(outputDir + 'deformation.dat', 'a')
        self.forceFile = open(outputDir + 'force.dat', 'a') 
        self.voltageFile = open(outputDir + 'voltage.dat', 'a')
        self.sim = sim
        self.probeIndex = probeIndex
        self.outputDir = outputDir
       
    def finish(self):
        self.defFile.close()
        self.forceFile.close()
        self.voltageFile.close()


    def writeData(self):
        globalTime = self.sim.globalTime
        timeStep = self.sim.timeStep
        deformation = self.sim.deformation
        maxDef = deformation.min(axis = 0)
        self.defFile.write('%e\t%e\t%e\t' % (globalTime, timeStep,maxDef[2]))
        for i in range(0, len(self.probeIndex)):
        	self.defFile.write('%e\t' % deformation[self.probeIndex[i]][2])
        self.defFile.write('\n') 
        self.defFile.flush() 

        vel = self.sim.velocity
        acc = self.sim.acceleration
        eForce = self.sim.elecForceSum
        fForce = self.sim.flowForceSum
        cForce = self.sim.contactForceSum
        self.forceFile.write('%e\t' % globalTime)
        for i in range(0, len(self.probeIndex)):
        	self.forceFile.write('%e\t' % vel[self.probeIndex[i]][2])
        self.forceFile.write('%e\t%e\t%e\n' % (eForce, fForce, cForce))
        self.forceFile.flush()

        voltage = self.sim.voltage
        self.voltageFile.write('%e\t%e\n' % (globalTime, voltage))
        self.voltageFile.flush()

    def saveFluidVTK(self, n):
        geomFields = self.sim.geomFields
        fluidMeshes = self.sim.fluidMeshes        
        elecFields = self.sim.elecFields
        if self.sim.enableFlowModel:
            flowFields = self.sim.flowFields
        writer = exporters.VTKWriterA(geomFields,fluidMeshes,
                                      self.outputDir + "fluid-" + str(n) + ".vtk",
                                      "gen5_fluid",
                                      False,0)
        writer.init()
        writer.writeScalarField(elecFields.potential,"potential")
        writer.writeVectorField(elecFields.electric_field,"potentialgradient")  
        if self.sim.enableFlowModel:
            writer.writeVectorField(flowFields.velocity,"velocity")
            writer.writeScalarField(flowFields.pressure, "pressure")
        writer.finish()

    def saveBeamVTK(self, n):
        geomFields = self.sim.geomFields
        solidMeshes = self.sim.solidMeshes
        plateFields = self.sim.plateFields
        writer = exporters.VTKWriterA(geomFields,solidMeshes,
                                  self.outputDir + "beam-" + str(n) + ".vtk",
                                  "gen5_beam",
                                  False,0)
        writer.init()
        writer.writeVectorField(plateFields.deformation,"deformation")
        writer.writeScalarField(plateFields.force, "force")
        writer.finish()
        
    def saveBeamBoundaryVTK(self, n):
        geomFields = self.sim.geomFields
        solidBoundaryMeshes = self.sim.solidBoundaryMeshes
        
        writer3 = exporters.VTKWriterA(geomFields,solidBoundaryMeshes,
                                  self.outputDir + "beamBoundary-" + str(n) + ".vtk",
                                  "beam Boundary",
                                  False,0,True)
        writer3.init()
        #writer3.writeVectorField(flowFields.velocity,"velocity")
        #writer3.writeVectorField(flowFields.force,"flow_force")
        #writer3.writeVectorField(elecFields.force,"elec_force")
        writer3.finish()



