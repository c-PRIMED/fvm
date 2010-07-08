#1 /usr/bin/env python
import fvm.fvmbaseExt as fvmbaseExt
from numpy import *
from mpi4py import MPI
from array import array
from time import *
import math
import time 

class MPMCoupling:

     def __init__(self, FlowModel, flowFields, MeshList, GeomField):
         self.meshList      = MeshList
         self.mesh          = MeshList[0]  #for now!!!!!!!!!!!!!!!!!!!!!
         self.geomField = GeomField
	 self.flowModel = FlowModel
	 self.flowFields = flowFields
         self.setup()
         
     def setup(self):

         self.dtMPM = zeros(1,float);
	 self.dt      = zeros(1,float)
         self.ndim = int(3)

         self.procID = MPI.COMM_WORLD.Get_rank();
         #connect to server first get parent and port name then connect to port name
         self.PARENT = MPI.Comm.Get_parent()
         assert self.PARENT != MPI.COMM_NULL
         
         self.portName =  array('c', ['\0']*MPI.MAX_PORT_NAME)
         print " self.portName = ", self.portName 
         print " self.portName = ", self.portName.tostring();
          #get port name from parent (sent by MPM)
         if  self.procID == 0:
           self.PARENT.Recv([self.portName, MPI.CHAR], source=0, tag=7777)
         print " self.portName = ", self.portName 
         print " self.portName = ", self.portName.tostring();
         #Connect this port
         self.FVM_COMM_MPM = MPI.COMM_WORLD.Connect(self.portName.tostring().rstrip('\0'), MPI.INFO_NULL, root=0) 
         
         self.remoteSize   = self.FVM_COMM_MPM.Get_remote_size()
         self.nfaces       = zeros( self.remoteSize, dtype='i' )
         self.nBndryNodes  = zeros( self.remoteSize, dtype='i' )
       
         self.couplingStep = 0
	 self.isCouplingInit = False
	 
     def couplingInit(self, dt):
	 #this should be done once, fvm and mpm time steps should be fixed  
         if (  self.isCouplingInit == False ):
	      self.dt[0] = dt
              #getting mpm time step 
              self.FVM_COMM_MPM.Bcast([self.dtMPM, MPI.DOUBLE], root=0)
	      #sending fvm time step
	      mpi_proc = MPI.PROC_NULL
	      if ( self.procID == 0 ):
	          mpi_proc = MPI.ROOT
              self.FVM_COMM_MPM.Bcast([self.dt,MPI.DOUBLE], mpi_proc)	
              self.dtMPM = self.dtMPM / 1000.0 #from "msec" to "second" conversion
	      a = int(round(self.dt/self.dtMPM))
	      b = int(round(self.dtMPM/self.dt))
	      #coupling if mod(step, couplingStep)=0
	      if ( a <= b ): 
	          self.couplingStep = b
	      else:
	          self.couplingStep = 1 
	      self.isCouplingInit = True
	      if ( MPI.COMM_WORLD.Get_rank() == 0 ):
             	 print "self.MPM = ", self.dtMPM, " self.dt = ", self.dt,  " couplingStep = ", self.couplingStep 
              
          
     def dim(self):
         return self.ndim;

     def updateMPM(self, istep):
        
         if (  istep%self.couplingStep == 0  ):
	     totCentroids = self.nfaces.sum()
             self.stress  = self.flowFields.force[self.surfaceMeshes[0].getFaces()].asNumPyArray().copy()
             recvbuf      = self.flowFields.force[self.surfaceMeshes[0].getFaces()].asNumPyArray().copy()
		     
             self.FVM_COMM_MPM.Allreduce( [self.stress, MPI.DOUBLE], [recvbuf,MPI.DOUBLE], op=MPI.SUM)

             print "face Stress(11) = ", self.stress[11][0], "  ",   self.stress[11][1], "  ", \
                   self.stress[11][2], " \n "


     def acceptMPM( self, istep):
         if (  istep%self.couplingStep == 0  ):
                #gettin nlocalfaces from MPM ( len(nlocalfaces) = remote_comm_world.size() )
                self.FVM_COMM_MPM.Allgather([None,MPI.INT],[self.nfaces,MPI.INT])
                #count
                count = self.nfaces * 4 * 3 # each face has four nodes and three coordinates
                #displ
                displ = zeros( len(count), dtype='i')
                #filling displ
                displ[0] = 0
                for i in range(1,len(count)):
                   displ[i] = displ[i-1] + count[i-1]
                #creating fvm array 
                self.faceNodesCoord =  self.geomField.coordinate[self.mesh.getCells()].newSizedClone( self.nfaces.sum()*4 )
                self.FVM_COMM_MPM.Allgatherv([None,0,0,MPI.DOUBLE],[self.faceNodesCoord.asNumPyArray(),count, displ,MPI.DOUBLE]) 
		print "self.nfaces = ", self.nfaces
                print "MPI_RANK = ", MPI.COMM_WORLD.Get_rank(), "  faceNodes coord = ", self.faceNodesCoord.asNumPyArray()                 
		self.faceNodesCoord.asNumPyArray()[:,:] = self.faceNodesCoord.asNumPyArray()[:,:] / 1000.0
                #creating meshList 
                #self.surfaceMeshes= fvmbaseExt.MeshList( fvmbaseExt.Mesh(3,99,self.faceNodesCoord) );
		self.surfaceMeshes = [ fvmbaseExt.Mesh(3,99,self.faceNodesCoord) ]
                
		#self.surfaceMesh = fvmbaseExt.Mesh(3,99,self.faceNodesCoord)
                #count
                count = self.nfaces * 3 
                #displ
                displ = zeros(len(count), dtype='i')
                #filling displ
                displ[0] = 0
                for i in range(1,len(count)):
                   displ[i] = displ[i-1] + count[i-1]
                #creating fvm array 
                self.FaceCentroidVels =  self.geomField.coordinate[self.mesh.getCells()].newSizedClone( self.nfaces.sum() )
		print "facecentroidvels = ", self.FaceCentroidVels.asNumPyArray(), "   rank = ", MPI.COMM_WORLD.Get_rank()
                self.FVM_COMM_MPM.Allgatherv([None,0,0,MPI.DOUBLE],[self.FaceCentroidVels.asNumPyArray(),count, displ,MPI.DOUBLE]) 
		if  MPI.COMM_WORLD.Get_rank() == 0:
		    self.dump_faces(self.faceNodesCoord.asNumPyArray(), self.FaceCentroidVels.asNumPyArray())

     def  getSurfaceMeshes(self):
         return self.surfaceMeshes
  
     def  getSurfaceVel(self):
         return self.FaceCentroidVels


     def dump_faces( self, faceNodesCoord, faceVel ):
        f = open('faces.dat','w')
	f.write("Title = \" MPM Faces \" \n")
        f.write("variables = \"x\", \"y\", \"z\", \"velX\", \"velY\", \"velZ\" \n")
	totFaces = self.nfaces.sum()
	f.write("#totfaces = %d\n"%(totFaces))
	indx = 0
	#loop over faces
	for n in range(0,totFaces):     
	   indx = 4 * n
	   nnode = 4 
	   title_name = str(n)
	   f.write("Zone T = \"%s\", N = %s E = %s DATAPACKING = BLOCK, VARLOCATION = ([4-6]=CELLCENTERED),  ZONETYPE=FEQUADRILATERAL\n" %  (title_name,  4, 1))   
	   #write x
	   #loop over nodes
	   for i in range(0,nnode):
               f.write(str(faceNodesCoord[indx+i][0])+"    ")
           f.write("\n")

	   #write y
	   #loop over nodes
	   for i in range(0,nnode):
               f.write(str(faceNodesCoord[indx+i][1])+"    ")
           f.write("\n")

	   #write z
	   #loop over nodes
	   for i in range(0,nnode):
               f.write(str(faceNodesCoord[indx+i][2])+"    ")
           f.write("\n")
  
           #write velX
           f.write( str(faceVel[n][0]) + "    ")
           #write velY
           f.write( str(faceVel[n][1]) + "    ")
           #write velZ
           f.write( str(faceVel[n][2]) + "    ")
           f.write("\n")
           #connectivity
           for i in range(0,nnode):
               f.write( str(i+1) + "     ")
           f.write("\n")

        f.close()   
	 
     def __del__(self):
       self.FVM_COMM_MPM.Barrier()
       self.FVM_COMM_MPM.Disconnect()
       self.PARENT.Barrier()
       self.PARENT.Disconnect()
