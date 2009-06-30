#!/usr/bin/env python

"""
Usage: testPartMesh.py [options] infile
options are:
--type  'tri'[default], 'quad', 'hexa', or 'tetra'
--xdmf  Dump data in xdmf
"""
import sys
sys.setdlopenflags(0x100|0x2)
import fvmbaseExt, importers, fvmparallel
from mpi4py import MPI
from FluentCase import FluentCase
from optparse import OptionParser

etype = {
        'tri' : 1,
        'quad' : 2,
        'tetra' : 3,
        'hexa' : 4
        }

numIterations = 10
def usage():
    print __doc__
    sys.exit(1)

parser = OptionParser()
parser.set_defaults(type='tri')
parser.add_option("--type", help="'tri'[default], 'quad', 'hexa', or 'tetra'")
parser.add_option("--xdmf", action='store_true', help="Dump data in xdmf")
(options, args) = parser.parse_args()
if len(args) != 1:
    usage()

reader = FluentCase(args[0])
reader.read()

fluent_meshes = reader.getMeshList()
nmesh = MPI.COMM_WORLD.Get_size()

npart = [nmesh]
etype = [etype[options.type]]

part_mesh = fvmparallel.PartMesh( fluent_meshes, npart, etype );
part_mesh.setWeightType(0);
part_mesh.setNumFlag(0);

#actions
part_mesh.partition()
part_mesh.mesh()
part_mesh.mesh_debug()
part_mesh.debug_print()
if options.xdmf:
    part_mesh.mesh_xdmfplot()
#meshes = part_mesh.meshList()

