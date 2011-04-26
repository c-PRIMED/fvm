#1 /usr/bin/env python
import fvm.fvmbaseExt as fvmbaseExt
from numpy import *
from mpi4py import MPI
from array import array
from time import *
import math
import time 

class ClientCoupling:

     def __init__(self, solidBoundaryMeshes, geomFields, flowField):
	 self.solidBoundaryMesh = solidBoundaryMeshes[0]
	 self.geomFields = geomFields
	 self.flowField = flowField
         self.setup()
         
     def setup(self):

         self.dtMPM = zeros(1,float);
	 self.dt    = zeros(1,float)
         self.ndim  = int(3)

         self.procID = MPI.COMM_WORLD.Get_rank();
         #connect to server first get parent and port name then connect to port name
         self.PARENT = MPI.Comm.Get_parent()
         assert self.PARENT != MPI.COMM_NULL
         self.portName = str(0) 
         #get port name from parent (sent by MPM)
         if  self.procID == 0:
           self.portName = self.PARENT.recv(source=0, tag=7777)
         #Connect this port
         self.CLIENT_COMM_SERVER = MPI.COMM_WORLD.Connect(self.portName, MPI.INFO_NULL, root=0) 
         self.remoteSize   = self.CLIENT_COMM_SERVER.Get_remote_size()
	 #get coordinate from goemField 
         self.meshcoordA = self.solidBoundaryMesh.getNodeCoordinatesPtr()
         self.meshcoord  = self.meshcoordA.asNumPyArray() 
	 self.coordA  = self.geomFields.coordinate[self.solidBoundaryMesh.getNodes()]
	 self.coord   = self.coordA.asNumPyArray()
	 self.sendCoordBuf = self.coord.copy()


	 #getting velocities from server side (sendBuff is dummy)
	 self.velA = self.flowField.velocity[self.solidBoundaryMesh.getFaces()]
	 self.vel  = self.velA.asNumPyArray()
	 self.sendVelBuf = self.vel.copy()

	 
	 #force array
	 solidBoundaryFaceSite = self.solidBoundaryMesh.getFaces()
	 nfaces = solidBoundaryFaceSite.getCount()
	 self.forceA = self.geomFields.coordinate[solidBoundaryFaceSite].newSizedClone(nfaces)
	 self.force  = self.forceA.asNumPyArray()
	 self.recvForcebuf = self.force.copy()
	 
	 if not MPI.COMM_WORLD.Get_rank():
   	     print "from ClientSide, reporting remoteSize = ", self.remoteSize
     
     def update(self, ModelList, FieldList):
        for model in ModelList:
           model.computeSolidSurfaceForce(self.solidBoundaryMesh.getFaces())
        
	self.force[:,:] = 0.0        
        for field in FieldList:
           forceA = field.force[self.solidBoundaryMesh.getFaces()]
	   force  = forceA.asNumPyArray()
	   self.force += force

        #sending stress to serverside	    
        self.CLIENT_COMM_SERVER.Allreduce( [self.force, MPI.DOUBLE], [self.recvForcebuf, MPI.DOUBLE], op=MPI.SUM)


     def accept(self):
	#getting coordinate (as summed) from server side
        self.CLIENT_COMM_SERVER.Allreduce( [self.sendCoordBuf, MPI.DOUBLE], [self.coord,MPI.DOUBLE], op=MPI.SUM)
        #sharing data with mesh coordinate
	self.meshcoord[:,:] = self.coord[:,:]
        #dump coord
	#if MPI.COMM_WORLD.Get_rank() == 0:
 	#   self.dump3D(self.coord,"solidBoundaryCoord.dat")

	self.CLIENT_COMM_SERVER.Allreduce( [self.sendVelBuf, MPI.DOUBLE],[self.vel,MPI.DOUBLE], op=MPI.SUM)

	#dump vel
	#if MPI.COMM_WORLD.Get_rank() == 0:
	#   self.dump3D(self.vel,"solidBoundaryVel.dat")


     def dump3D(self,coord, fname):
        f = open(fname,'w')
        for n in range(0,len(coord)):
	   x = str( coord[n][0] )
	   y = str( coord[n][1] )
	   z = str( coord[n][2] )
	   f.write( x + "     " + y + "       " + z + "\n")

	f.close()
     def __del__(self):
       self.CLIENT_COMM_SERVER.Barrier()
       self.CLIENT_COMM_SERVER.Disconnect()
       self.PARENT.Barrier()
       self.PARENT.Disconnect()
