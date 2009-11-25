#!/usr/bin/env python

import fvm
import fvm.fvmparallel as fvmparallel
import time
from numpy import *
from mpi4py  import MPI
from FluentCase import FluentCase

fileBase = None
numIterations = 10
fileBase = "/home/yildirim/memosa/src/fvm/test/tri_894"
#fileBase = "/home/yildirim/memosa/src/fvm/test/test_tri_500by500"
#fileBase = "/home/yildirim/memosa/src/fvm/test/cav32"


# change as needed

outfile = None
if __name__ == '__main__' and fileBase is None:
    if len(sys.argv) < 2:
        usage()
    fileBase = sys.argv[1]
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

if outfile == None:
    outfile = fileBase+"-prism.dat"
    
reader = FluentCase(fileBase+".cas")

#import debug
reader.read();
import sys

fluent_meshes = reader.getMeshList()

import time
t0 = time.time()
nmesh = MPI.COMM_WORLD.Get_size()

#print "nmesh = ", nmesh
#npart = fvmparallel.IntVector(1,nmesh)  #total of distributed meshes
#etype = fvmparallel.IntVector(1,1) #triangle

npart = [nmesh]
etype = [1]

part_mesh = fvmparallel.PartMesh( fluent_meshes, npart, etype );
part_mesh.setWeightType(0);
part_mesh.setNumFlag(0);

#actions
part_mesh.partition()
part_mesh.mesh()
part_mesh.debug_print()
part_mesh.mesh_debug()
meshes = part_mesh.meshList()

t1 = time.time()
#if outfile != '/dev/stdout':
#    print '\nsolution time = %f' % (t1-t0)

