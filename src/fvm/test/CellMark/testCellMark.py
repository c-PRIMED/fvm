#!/usr/bin/env python
"""
testCellMark.py

Usage: testCellMark.py [options] casefile mpmfile

Where 'casefile' is a Fluent case file.

Options:
  --search n   Search Option [1-4] (default 1)
"""
import fvm
fvm.set_atype('double')
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
from mpi4py import MPI
from FluentCase import FluentCase
from optparse import OptionParser

#import pdb
#pdb.set_trace()

def usage():
    print __doc__
    sys.exit(1)

# option parsing
parser = OptionParser()
parser.set_defaults(search=1)
parser.add_option("--search", type="int",
                  dest="search",
                  help="Search Option [1-4].")
(options, args) = parser.parse_args()
print "options=",options.search
if len(args) != 2 or options.search < 1 or options.search > 4:
    usage()
casefile = args[0]
mpmfile = args[1]
fileBase = ''          

reader = FluentCase(casefile)
reader.read();
meshes = reader.getMeshList()
geomFields =  fvm.models.GeomFields('geom')

metricsCalculator = fvm.models.MeshMetricsCalculatorA(geomFields,meshes)
metricsCalculator.init()

flowFields =  fvm.models.FlowFields('flow')
fmodel = fvm.models.FlowModelA(geomFields,flowFields,meshes)
reader.importFlowBCs(fmodel)

octree = fvmbaseExt.Octree()
mesh0 = meshes[0]
cells = mesh0.getCells()
nCells = cells.getCount()
cellCentroid = geomFields.coordinate[cells].asNumPyArray()
octree.Impl(mesh0, geomFields)

solid = fvmbaseExt.MPM(mpmfile)
fvmbaseExt.CellMark_Impl(mesh0, geomFields, fileBase, octree, solid, options.search)
