#! /usr/bin/env python
#important note, sometime you got solid or fluid scripts through e-mail, which
#somehow might convert filetype, then this script will throw No such file or directory
#error when it try to spawn solid or fluid scripts. The way to fix this is 
#to open scripts in vim and :set fileformat=unix and :wq
from mpi4py import MPI
from array import array
import os
import time
assert MPI.COMM_WORLD.Get_size() == 1

parent = MPI.COMM_WORLD
nprocs_solid =1  
print "ssssssssss"
solid_comm = parent.Spawn(command="./solid_plate.py",   
                      args=['--type=quad','--volt=100.0'], maxprocs=nprocs_solid, root=0)
print "solid comm"
#this is for totalview
for i in range(0,0):
   time.sleep(3000)		    
   print "time = ", i, " secs"

nprocs_fluid = 16 
fluid_comm = parent.Spawn(command="./fluid_elec_plate.py", 
                      args=['--typeFluid=hexa','--typeSolid=quad','--volt=100.0'], maxprocs=nprocs_fluid, root=0)
port=solid_comm.recv(source=0, tag=7776)
fluid_comm.send     (port, dest=0, tag=7777)

solid_comm.Barrier()
solid_comm.Disconnect()

fluid_comm.Barrier()
fluid_comm.Disconnect()


