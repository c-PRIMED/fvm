#!/usr/bin/env python

"""
Usage: testPartMesh.py infile
"""
import sys
sys.setdlopenflags(0x100|0x2)
import fvmbaseExt, importers, fvmparallel
from mpi4py import MPI
from FluentCase import FluentCase

numIterations = 10
def usage():
    print __doc__
    sys.exit(1)

if len(sys.argv) < 2:
    usage()
casfile = sys.argv[1]
reader = FluentCase(casfile)
reader.read()

fluent_meshes = reader.getMeshList()
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

