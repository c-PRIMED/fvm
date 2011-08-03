#! /usr/bin/env python
from numpy import *
from mpi4py import MPI
from array import array
from time import *
import math
import time
import fvm
import fvm.fvmbaseExt as fvmbaseExt

class ServerCoupling:

    def __init__(self, pModel, dModel, geomField, solidMeshes, solidBoundaryMeshes, plateField, thickness):		
        self.plateModel = pModel
        self.deformationModel = dModel
	self.geomField  = geomField
        self.solidMesh = solidMeshes[0]
	self.solidBoundaryMesh = solidBoundaryMeshes[0]
	self.plateField = plateField
	self.thickness  = thickness
        self.setup()


    def setup(self):
         self.dt       = zeros(1,float)
         self.dtClient = zeros(1,float)
       
         self.procID = MPI.COMM_WORLD.Get_rank()
         #connect to server first get parent and port name then connect to port name
         self.PARENT = MPI.Comm.Get_parent()
         assert self.PARENT != MPI.COMM_NULL
	 assert MPI.COMM_WORLD.Get_size() == 1
         self.portName = str(0) 
         #print " ServerSide: self.portName = ", self.portName 
         #print " ServerSide: self.portName = ", self.portName.tostring();
          #get port name from parent (sent by MPM)
         if  self.procID == 0:
	   self.portName = MPI.Open_port(MPI.INFO_NULL)
           print " Server Side: self.portName = ", self.portName 
	   print "type= ", type(self.portName)
           self.PARENT.send(self.portName, dest=0, tag=7776)
         #print " Server Side: self.portName = ", self.portName.tostring();
         #Connect this port
         self.SERVER_COMM_CLIENT = MPI.COMM_WORLD.Accept(self.portName, MPI.INFO_NULL, root=0) 
         self.remoteSize   = self.SERVER_COMM_CLIENT.Get_remote_size()
        
	 #get coordinate from goemField
         self.coordA = self.solidBoundaryMesh.getNodeCoordinatesPtr()
         self.coord  = self.coordA.asNumPyArray()
	 self.coordRecvBuf = self.coord.copy() 
#         self.coordA  = self.geomFields.coordinate[self.solidBoundaryMesh.getNodes()]
#         self.coord   = self.coordA.asNumPyArray()
#         self.sendBuf = self.coord.copy()


	 #velocity field, adding face velocities
	 self.vField = fvmbaseExt.Field("velocity")
	 self.faces  = self.solidBoundaryMesh.getFaces()
	 self.nfaces = self.faces.getCount()
	 area  = self.geomField.area[self.faces]

         #initialize force to zero	
	 solidCells = self.solidMesh.getCells()
         forcePlate = self.plateField.force[solidCells].asNumPyArray() 
         forcePlate[:] = 0.0
	 thickness = self.plateField.thickness[solidCells].asNumPyArray()
         thickness[:] = self.thickness

         #ading 
         self.forceA  = area.newSizedClone(self.nfaces)
	 self.force   = self.forceA.asNumPyArray()
	 self.sendForceBuf = self.force.copy()
	 if not MPI.COMM_WORLD.Get_rank():
   	     print "from ServerSide, reporting remoteSize = ", self.remoteSize

    def  update(self):
         timeStep = self.plateModel.getOptions().getVar("timeStep")
	 print "timeStep = ", timeStep
         #get new coordinates of bMesh
	 self.deformationModel.updateBoundaryMesh(self.solidMesh, self.solidBoundaryMesh, self.thickness, \
              timeStep, self.vField) 
	 velA = self.vField[self.faces]
	 vel  = velA.asNumPyArray()
	 velBuf = vel.copy()
         #self.dumpSolidBoundaryCoord(self.coord,MPI.COMM_WORLD.Get_rank())
	 #sending coordinate to client side as summing
         self.SERVER_COMM_CLIENT.Allreduce( [self.coord, MPI.DOUBLE], [self.coordRecvBuf,MPI.DOUBLE], op=MPI.SUM )

         #sending velocity to client side as summing (recvbuf is dummy, we use previous one)
         self.SERVER_COMM_CLIENT.Allreduce([vel, MPI.DOUBLE], [velBuf,MPI.DOUBLE], op=MPI.SUM)
	 #dump velcocity
	 #self.dumpSolidBoundaryVel(vel, MPI.COMM_WORLD.Get_rank())


    def  accept(self):

       self.force[:,:] = 0.0
       #receving stress to serverside	    
       self.SERVER_COMM_CLIENT.Allreduce( [self.sendForceBuf, MPI.DOUBLE], [self.force, MPI.DOUBLE], op=MPI.SUM)
       #if MPI.COMM_WORLD.Get_rank() == 0:
       #   self.dump3D(self.force,"solidBoundaryForce.dat")
      
       solidCells = self.solidMesh.getCells()
       selfCountCells = solidCells.getSelfCount()
       countCells     = solidCells.getCount()
       forcePlate = self.plateField.force[solidCells].asNumPyArray() 
       forcePlate[:] = 0.0

       #force on interior cells    # force on interior cells
       for c in range(0, selfCountCells):
          botFaceIndx = c
          topFaceIndx = c + selfCountCells
          forcePlate[c] = self.force[botFaceIndx][2] + self.force[topFaceIndx][2]
       # force on boundary cells
       for c in range(selfCountCells, countCells):
          forcePlate[c] =  self.force[selfCountCells + c][2]

    def updateTimeStep(self, timestep):
        self.timestep = zeros(1,float)
        self.timestep[0] = timestep
        mpi_proc = MPI.PROC_NULL
        if ( self.procID == 0 ):
            mpi_proc = MPI.ROOT
        self.SERVER_COMM_CLIENT.Bcast([self.timestep, MPI.DOUBLE], mpi_proc)	
 
    def dumpSolidBoundaryCoord(self,coord, rank):
        f = open('solidBoundaryCoord'+str(rank)+'.dat','w')
        for n in range(0,len(coord)):
	   x = str( coord[n][0] )
	   y = str( coord[n][1] )
	   z = str( coord[n][2] )
	   f.write( x + "     " + y + "       " + z + "\n")
    def dumpSolidBoundaryVel(self,coord, rank):
        f = open('solidBoundaryVel'+str(rank)+'.dat','w')
        for n in range(0,len(coord)):
	   x = str( coord[n][0] )
	   y = str( coord[n][1] )
	   z = str( coord[n][2] )
	   f.write( x + "     " + y + "       " + z + "\n")


	f.close()
    def dump3D(self,coord, fname):
        f = open(fname,'w')
        for n in range(0,len(coord)):
	   x = str( coord[n][0] )
	   y = str( coord[n][1] )
	   z = str( coord[n][2] )
	   f.write( x + "     " + y + "       " + z + "\n")

	f.close()

    def __del__(self):
         if ( self.procID == 0 ):
	     MPI.Close_port(self.portName)
         self.SERVER_COMM_CLIENT.Barrier()
	 self.SERVER_COMM_CLIENT.Disconnect()
	 self.PARENT.Barrier()
	 self.PARENT.Disconnect()

    
