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

    def __init__(self, sModel, dModel, geomField, solidMeshes, solidBoundaryMeshes):		
        self.solidModel = sModel
        self.deformationModel = dModel
	self.geomField  = geomField
        self.solidMesh = solidMeshes[0]
	self.solidBoundaryMesh = solidBoundaryMeshes[0]
        self.setup()


    def setup(self):
         self.dt       = zeros(1,float)
         self.dtClient = zeros(1,float)
       
         self.procID = MPI.COMM_WORLD.Get_rank()
         #connect to server first get parent and port name then connect to port name
         self.PARENT = MPI.Comm.Get_parent()
         assert self.PARENT != MPI.COMM_NULL
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
        
	 #solidBoundary mesh site
	 self.solidMesh.setNodeRepeationArrayCoupling(self.solidBoundaryMesh)
	 self.solidMesh.setCommonFacesMap(self.solidBoundaryMesh)
	 self.createSolidForceBVFields()

	 if not MPI.COMM_WORLD.Get_rank():
   	     print "from ServerSide, reporting remoteSize = ", self.remoteSize
    def  update(self):
         timeStep = self.solidModel.getOptions().getVar("timeStep")
         #get new coordinates of bMesh
	 coordA = self.solidMesh.getUpdatedNodesCoordCoupling(self.geomField, self.solidBoundaryMesh)
	 coord = coordA.asNumPyArray()
	 self.nfaces = len(coord)
         self.dumpSolidBoundaryCoord(coord,MPI.COMM_WORLD.Get_rank())
	 recvbuf = coord.copy() 
	 #sending coordinate to client side as summing
         self.SERVER_COMM_CLIENT.Allreduce( [coord, MPI.DOUBLE], [recvbuf,MPI.DOUBLE], op=MPI.SUM)

	 #now sending velocities
	 vField = fvmbaseExt.Field("velocity")
	 self.deformationModel.updateBoundaryMesh(self.solidMesh, self.solidBoundaryMesh, vField, self.solidMesh.getCommonFacesMap(),  timeStep);
	 velA = vField[self.solidBoundaryMesh.getFaces()]
	 vel  = velA.asNumPyArray()
	 #sending velocity to client side as summing (recvbuf is dummy, we use previous one)
         self.SERVER_COMM_CLIENT.Allreduce([vel, MPI.DOUBLE], [recvbuf,MPI.DOUBLE], op=MPI.SUM)


    def  accept(self):
       forceA  = self.geomField.coordinate[self.solidMesh.getCells()].newSizedClone(self.nfaces)

       force   = forceA.asNumPyArray() #deep copy, we don't want to update forceA
       force[:,:] = 0.0
       sendbuf = force.copy() 
       #receving stress to serverside	    
       self.SERVER_COMM_CLIENT.Allreduce( [sendbuf, MPI.DOUBLE], [force, MPI.DOUBLE], op=MPI.SUM)
       #if MPI.COMM_WORLD.Get_rank() == 0:
       #   self.dump3D(force,"solidBoundaryForce.dat")
       
       bcMap = self.solidModel.getBCMap()
       fgs = self.solidMesh.getBoundaryFaceGroups()
       for fg in fgs:
            bc = bcMap[fg.id]
            if bc.bcType == 'SpecifiedForce':
                sfaces = fg.site
                fxA = self.bForceXField[sfaces]
                fyA = self.bForceYField[sfaces]
                fzA = self.bForceZField[sfaces]
		#update fxA, fyB,fzA for given forceA 
		self.solidModel.updateForceOnBoundary(sfaces, forceA,  self.solidMesh.getCommonFacesMap(), \
		                                      fxA, fyA, fzA) 
 
    def dumpSolidBoundaryCoord(self,coord, rank):
        f = open('solidBoundaryCoord'+str(rank)+'.dat','w')
	print "yaziyor = ", len(coord) 
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

    def createSolidForceBVFields(self):
       self.bForceXField = fvmbaseExt.Field('bForceX')
       self.bForceYField = fvmbaseExt.Field('bForceY')
       self.bForceZField = fvmbaseExt.Field('bForceZ')
       bcMap = self.solidModel.getBCMap()
       fgs = self.solidMesh.getBoundaryFaceGroups()
       for fg in fgs:
            bc = bcMap[fg.id]
            if bc.bcType == 'SpecifiedForce':
                sfaces = fg.site
                faceCount = sfaces.getCount()
                areaMag = self.geomField.areaMag[sfaces]
                fxA = areaMag.newSizedClone(faceCount)
                fyA = areaMag.newSizedClone(faceCount)
                fzA = areaMag.newSizedClone(faceCount)
                fx = fxA.asNumPyArray()
                fy = fyA.asNumPyArray()
                fz = fzA.asNumPyArray()
                fx[:] = 0
                fy[:] = 0
                fz[:] = 0
                self.bForceXField[sfaces]=fxA
                self.bForceYField[sfaces]=fyA
                self.bForceZField[sfaces]=fzA
                bc['specifiedXForce'] = self.bForceXField
                bc['specifiedYForce'] = self.bForceYField
                bc['specifiedZForce'] = self.bForceZField

    def _del_(self):
         if ( self.procID == 0 ):
	     MPI.Close_port(self.portName)
         MPI.SERVER_COMM_CLIENT.Barrier()
	 MPI.SERVER_COMM_CLIENT.Disconnect()
	 self.PARENT.Barrier()
	 self.PARENT.Disconnect()

