#!/usr/bin/env python
import time
import tecplotMesh
"""
Usage: testThermalParallel.py [options] infile
options are:
--type  'tri'[default], 'quad', 'hexa', or 'tetra'
--xdmf  Dump data in xdmf
"""

import fvm.fvmbaseExt as fvmbaseExt
import fvm
fvm.set_atype('double')
import  fvm.importers as importers
import  fvm.fvmparallel as fvmparallel
import time

from numpy import *
from mpi4py  import MPI
from optparse import OptionParser
from FluentCase import FluentCase

def usage():
    print __doc__
    sys.exit(1)

# map between fvm, tecplot, and xdmf types
etype = {
        'tri' : 1,
        'quad' : 2,
        'tetra' : 3,
        'hexa' : 4
        }
tectype = {
        'tri' : 'FETRIANGLE',
        'quad' : 'FEQUADRILATERAL',
        'tetra' : 'FETETRAHEDRON',
        'hexa' : 'FEBRICK'
        }
xtype = {
        'tri' : 'Triangle',
        'quad' : 'Quadrilateral',
        'tetra' : 'Tetrahedron',
        'hexa' : 'Hexahedron'
        }


parser = OptionParser()
parser.set_defaults(type='tri')
parser.add_option("--type", help="'tri'[default], 'quad', 'hexa', or 'tetra'")
parser.add_option("--xdmf", action='store_true', help="Dump data in xdmf")
parser.add_option("--time","-t",action='store_true',help="Print timing information.")
(options, args) = parser.parse_args()
if len(args) != 1:
    usage()

#import ddd
reader = FluentCase(args[0])
reader.read()
fluent_meshes = reader.getMeshList()

#assembling meshes
mesh_assembler = fvmbaseExt.MeshAssembler( fluent_meshes )
#mesh_assembler.debug_print()
meshes = mesh_assembler.meshList()



nmesh = 1
npart = [MPI.COMM_WORLD.Get_size()]
#etype = [etype[options.type]]
npartVec = fvmparallel.IntVector(nmesh, MPI.COMM_WORLD.Get_size())
etypeVec = fvmparallel.IntVector(nmesh, int(etype[options.type])) #triangle

#for i in range(0,nmesh):
#   npartVec.push_back( MPI.COMM_WORLD.Get_size() )
#   etypeVec.push_back( etype[options.type] )


if not MPI.COMM_WORLD.Get_rank():
   print "parmesh is processing"
     
    

if options.time:
    # time profile for partmesh
    part_mesh_time = zeros(1,dtype='d')
    part_mesh_start = zeros(1, dtype='d')
    part_mesh_end   = zeros(1, dtype='d')
    part_mesh_maxtime = zeros(1,dtype='d')
    part_mesh_mintime = zeros(1, dtype='d')
    part_mesh_start[0] = MPI.Wtime()

#partMesh constructor and setTypes
part_mesh = fvmparallel.PartMesh( meshes, npartVec, etypeVec );
part_mesh.setWeightType(2);
part_mesh.setNumFlag(0);
 
#actions
part_mesh.partition()
part_mesh.mesh()
#part_mesh.mesh_debug()
part_meshes = part_mesh.meshList()

#dismantling mesh
multi_mesher = fvmbaseExt.MeshDismantler( part_meshes )
meshes = multi_mesher.debug_face_nodes()
meshes = multi_mesher.meshList()
#tecplotMesh.dumpTecplotMesh(2, meshes, options.type)

