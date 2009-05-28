#!/bin/env python
#
# Author : Ian Curington
# 
# Example to generate a stack of layered prism elements
# Writes a "test20.xmf/h5" files with a 3D unstructured prism grid
# N vertical layers based on triangle plan grid
# containing "pressure" scalar node values,
# uses Compression for connectivity, geometry and attributes.
#

from Xdmf import *

BaseName = 'test20'

# define simple triangle mesh in horizontal plane

tri_conn   = [0,1,6, \
              1,2,6, \
              2,3,6, \
              3,4,6, \
              4,5,6, \
              5,0,6 ]

ntris = len(tri_conn)/3

tri_coords = [-2, 0, 0, \
              -1, 1, 0, \
               1, 1, 0, \
               2, 0, 0, \
               1,-1, 0, \
              -1,-1, 0, \
               0, 0, 0 ]

ncoords = len(tri_coords)/3

# prism layers, make large enough to trigger hdf5 writing
nlayers = 800 
layer_delta = 0.01

"""
Reset an existing .h5 file for truncating all or any existing contents,
by writing small junk array in 'w' mode.
"""
# Write H5 Data
arr = XdmfArray()
#arr.DebugOn() 
arr.SetNumberType(XDMF_FLOAT64_TYPE)
arr.SetShapeFromString("2 3")
arr.Generate(0.0, 1.0, 0, arr.GetNumberOfElements() - 1)
h5 = XdmfHDF()
h5.CopyType(arr)
h5.CopyShape(arr)
h5.Open( BaseName + '.h5:/_junk_', 'w')
h5.Write(arr)
h5.Close()

# XDMF init
d = XdmfDOM()

root = XdmfRoot()
root.SetDOM(d)
root.Build()

# Domain
dm = XdmfDomain()
root.Insert(dm)

# Grid
g = XdmfGrid()
g.SetName("Prism Grid")
#g.DebugOn() 

# Topology - 3D unstructured grid, prism element type
t = g.GetTopology()
t.SetTopologyType(XDMF_WEDGE)
t.SetNumberOfElements(nlayers * ntris)
#t.DebugOn() 
cc = t.GetConnectivity()

cc.SetNumberType(XDMF_INT32_TYPE)
# length of connectivity
nElem = nlayers*ntris*6
cc.SetNumberOfElements(nElem)

cc.SetHeavyDataSetName( BaseName + r'.h5:/Connectivity' )
if nElem > 1000:
    cc.SetCompression(9)

c_index = 0
for layer in range(nlayers):
   for tri in range(ntris):
      cc.SetValueFromInt64( c_index*6+0, tri_conn[3*tri+0]+layer*ncoords )
      cc.SetValueFromInt64( c_index*6+1, tri_conn[3*tri+1]+layer*ncoords )
      cc.SetValueFromInt64( c_index*6+2, tri_conn[3*tri+2]+layer*ncoords )
      cc.SetValueFromInt64( c_index*6+3, tri_conn[3*tri+0]+(layer+1)*ncoords )
      cc.SetValueFromInt64( c_index*6+4, tri_conn[3*tri+1]+(layer+1)*ncoords )
      cc.SetValueFromInt64( c_index*6+5, tri_conn[3*tri+2]+(layer+1)*ncoords )
      c_index = c_index + 1


# Geometry - origin and cell size
geo = g.GetGeometry()
geo.SetGeometryType(XDMF_GEOMETRY_XYZ)

# three explicit corners
ca = geo.GetPoints()
ca.SetNumberType(XDMF_FLOAT32_TYPE)
ca.SetHeavyDataSetName( BaseName + r'.h5:/Coordinates' )

# number of nodes * 3
nElem = (nlayers+1)*ncoords*3
ca.SetNumberOfElements(nElem)
if nElem > 1000:
    ca.SetCompression(9)

c_index = 0
for layer in range(nlayers+1):
   for node in range(ncoords):
      ca.SetValueFromFloat64( c_index*3+0, tri_coords[3*node+0] )
      ca.SetValueFromFloat64( c_index*3+1, tri_coords[3*node+1] )
      ca.SetValueFromFloat64( c_index*3+2, tri_coords[3*node+2]+layer*layer_delta )
      c_index = c_index + 1

dm.Insert(g)

# Attribute - as scalar node data
attr = XdmfAttribute()
attr.SetName("Pressure")
attr.SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_NODE);
attr.SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR);
#attr.DebugOn() 
p = attr.GetValues()
p.SetHeavyDataSetName( BaseName + r'.h5:/Pressure')

# number of nodes with values
nElem = (nlayers+1)*ncoords
p.SetNumberOfElements(nElem )
if nElem > 1000:
    p.SetCompression(9)

c_index = 0
for layer in range(nlayers+1):
   for node in range(ncoords):
      value = 10.0 + float(layer) * 2.0
      p.SetValueFromFloat64( c_index, value )
      c_index = c_index + 1


g.Insert(attr)

# Update XML and Write Values to DataItems
root.Build() # DataItems > 100 values are heavy
print d.Serialize() # prints to stdout

d.Write( BaseName + r'.xmf') # write to file


