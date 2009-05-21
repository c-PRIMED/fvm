#1 /usr/bin/env python
import fvmbaseExt
from numpy import *
from mpi4py import MPI
from array import array
from time import *



class MPMCoupling:

     def __init__(self, MeshList, FlowModel, FlowField, GeomField, solid):
         self.meshList      = MeshList
         self.mesh          = MeshList[0]  #for now!!!!!!!!!!!!!!!!!!!!!
         self.flowModel = FlowModel
         self.flowField = FlowField
         self.geomField = GeomField
         self.solid     = solid
         self.setup()

     def setup(self):
         #recv from MPM
         self.recvDT_TAG           = int(1122)
         self.recvTime_TAG         = int(2122)
         self.recvNparticle_TAG    = int(3122)
         self.recvIsContinue_TAG   = int(4122)
         self.recvParticlesPos_TAG = int(5122)
         self.recvParticlesVel_TAG = int(6122)
         #send to MPM
         self.dtSendTAG            = int(8111)
         self.timeSendTAG          = int(8121)
         self.totParticlesSendTAG  = int(8131)
         self.coordSendTAG         = int(8141)
         self.volumeSendTAG        = int(8142)
         self.stressSendTAG        = int(8151)

         self.dt = zeros(1,float);
         self.time = zeros(1,float);
         self.nparticles = zeros(1,int);
         self.isContinue = zeros(1,int);
         self.ndim = int(3)
         self.xc     = self.geomField.coordinate[self.mesh.getCells()].asNumPyArray()
         self.volume = self.geomField.volume[self.mesh.getCells()].asNumPyArray()
         

         self.procID = MPI.COMM_WORLD.Get_rank();
         #connect to server first get parent and port name then connect to port name
         self.PARENT = MPI.Comm.Get_parent()
         assert self.PARENT != MPI.COMM_NULL
         
         self.portName =  array('c', ['\0']*MPI.MAX_PORT_NAME)
         #get port name from parent (sent by MPM)
         if  self.procID == 0:
           self.PARENT.Recv([self.portName, MPI.CHAR], source=0, tag=7777)
         #Connect this port
         self.FVM_COMM_MPM = MPI.COMM_WORLD.Connect(self.portName.tostring(), MPI.INFO_NULL, root=0) 
         #construct FVMParticles class
         self.fvmParticles = fvmbaseExt.FVMParticles( self.meshList );
         self.totParticlesFVM = zeros(1,int32);
         
     def dim(self):
         return self.ndim;

     def totalParticles(self):
         self.FVM_COMM_MPM.Recv([self.nparticles, MPI.INT], source=0, tag=self.recvNparticle_TAG )
         return self.nparticles


     def updateMPM(self, dt, time, nsweep):
          #update fluid particles particles 
         mesh0 = int(0)   
         self.fvmParticles.setParticles( nsweep );
         self.totParticlesFVM[0] = int(self.fvmParticles.getNumOfFluidParticles( mesh0 ))
         self.FVM_COMM_MPM.Send( [self.totParticlesFVM, MPI.INT], dest=0, tag=self.totParticlesSendTAG )

         self.dtFVM   = zeros(1,float)
         self.timeFVM = zeros(1,float)
         self.dtFVM[0]   = dt
         self.timeFVM[0] = time
         #time and time step of FVM to MPM side
         self.FVM_COMM_MPM.Send( [self.dtFVM  , MPI.DOUBLE], dest=0, tag=self.dtSendTAG   )
         self.FVM_COMM_MPM.Send( [self.timeFVM, MPI.DOUBLE], dest=0, tag=self.timeSendTAG )



         #get coordinate, volume and stress at FVM particles
         self.cellIDs  = self.fvmParticles.getCellIDs( mesh0 )
         self.stress = self.flowModel.getStressTensor( self.meshList[ mesh0], self.cellIDs ).asNumPyArray()

         self.particlesCoord  = zeros( (self.totParticlesFVM[0],3), float )
         self.particlesVolume = zeros( self.totParticlesFVM[0], float)
         indx = int(0)
 #        print "shape(stress)  = ", shape(self.stress)
         for i in self.cellIDs.asNumPyArray():
              self.particlesCoord[indx,:] = self.xc[i,:]
              self.particlesVolume[indx]  = self.volume[i]
              indx = indx + 1


#         print "fluid particlesPos(11)  = ", self.particlesCoord[11][0], "  ", self.particlesCoord[11][1], "  ", \
#                                             self.particlesCoord[11][2], "\n"

#         print "fluid particlesStress(11) = ", self.stress[11][0], "  ",   self.stress[11][1], "  ", \
#               self.stress[11][2], "  ", self.stress[11][3] , "  ", self.stress[11][4], "   ", self.stress[11][5], "\n"

#         print "fluid particlesVolume(11) = ", self.particlesVolume[11], "\n"

         self.FVM_COMM_MPM.Send( [self.particlesCoord, MPI.DOUBLE], dest=0, tag=self.coordSendTAG  )

         self.FVM_COMM_MPM.Send( [self.stress, MPI.DOUBLE], dest = 0, tag = self.stressSendTAG );

         self.FVM_COMM_MPM.Send( [self.particlesVolume, MPI.DOUBLE], dest=0, tag=self.volumeSendTAG );
					

     def acceptMPM(self):


         self.FVM_COMM_MPM.Recv([self.isContinue, MPI.INT], source=0, tag =  self.recvIsContinue_TAG )
         if ( self.isContinue == 1 ):
             #get number of particles first
             nparticles = int( self.totalParticles() )
             self.px    = self.geomField.coordinate[self.mesh.getCells()].newSizedClone( nparticles )
             self.pv    = self.geomField.coordinate[self.mesh.getCells()].newSizedClone( nparticles )
             self.pType = fvmbaseExt.newIntArray( nparticles )
             self.pType.asNumPyArray()[:] = 1
             self.FVM_COMM_MPM.Recv([self.dt, MPI.DOUBLE], source=0, tag=self.recvDT_TAG )
             self.FVM_COMM_MPM.Recv([self.time, MPI.DOUBLE], source=0, tag=self.recvTime_TAG )

             self.FVM_COMM_MPM.Recv([self.px.asNumPyArray(), MPI.DOUBLE], source=0, tag=self.recvParticlesPos_TAG)
             self.FVM_COMM_MPM.Recv([self.pv.asNumPyArray(), MPI.DOUBLE], source=0, tag=self.recvParticlesVel_TAG)
#             print "px (FVM) = ", self.px.asNumPyArray()[791,0:3]
#             print "pv (FVM) = ", self.pv.asNumPyArray()[791,0:3]
             self.dump_coord_vel( self.nparticles, self.px.asNumPyArray(), self.pv.asNumPyArray() )
             self.particles = self.solid.getParticles( int(self.nparticles) )
             self.geomField.coordinate[self.particles] = self.px
             self.flowField.velocity[self.particles]   = self.pv 
             self.solid.setCoordinates( self.px )
             self.solid.setVelocities ( self.pv )
             self.solid.setTypes( self.pType )
              
     def  particleSite(self):
         return self.particles 

     def getCoord(self):
         return self.coordP.copy('C')

     def getVel(self):
         return self.velP.copy('C')
     
     def getNparticles( self):
         return int(self.nparticles[0])
 
     def   dump_coord_vel(self, nparticles, px, pv):
        fx = open('px.dat','w')
        fx.write(str(nparticles))
        fx.write("\n")
        for n in range(0,nparticles):
            fx.write(str(px[n,0]) )
            fx.write("     ")
            fx.write(str(px[n,1]) )
            fx.write("     ")
            fx.write(str(px[n,2]) )
            fx.write("\n")
        fx.close()

        fv = open('pv.dat','w')
        fv.write(str(nparticles))
        fv.write("\n")
        for n in range(0,nparticles):
            fv.write(str(pv[n,0]) )
            fv.write("     ")
            fv.write(str(pv[n,1]) )
            fv.write("     ")
            fv.write(str(pv[n,2]) )
            fv.write("\n")
        fv.close()


     def __del__(self):
       self.FVM_COMM_MPM.Barrier()
       self.FVM_COMM_MPM.Disconnect()
       self.PARENT.Barrier()
       self.PARENT.Disconnect()
