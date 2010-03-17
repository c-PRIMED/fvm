#!/usr/bin/env python
import time


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
from numpy import *
from optparse import OptionParser
from mpi4py import MPI
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

numIterations = 10

reader = FluentCase(args[0])
reader.read()
#import ddd
meshes_fluent = reader.getMeshList()
mesh_assembler = fvmbaseExt.MeshAssembler( meshes_fluent )
mesh_assembler.debug_sites()


