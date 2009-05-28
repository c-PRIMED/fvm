#!/usr/bin/env python

from Xdmf import *

# Example of How to Generate Xdmf
# The Heavy Data is written separately

# Write H5 Data
arr = XdmfArray()
arr.SetNumberType(XDMF_FLOAT64_TYPE)
arr.SetShapeFromString("10 20 30")
arr.Generate(0.0, 1.0, 0, arr.GetNumberOfElements() - 1)
h5 = XdmfHDF()
h5.CopyType(arr)
h5.CopyShape(arr)
h5.Open('XdmfByHand.h5:/Mydata', 'w')
h5.Write(arr)
h5.Close()
dv = XdmfValuesHDF()
DataXml = dv.DataItemFromHDF('XdmfByHand.h5:/Mydata')
#
d = XdmfDOM() 

root = XdmfRoot()
root.SetDOM(d)
root.SetVersion(2.2) # Change the Version number because we can
root.Build()
# Information
i = XdmfInformation() # Arbitrary Name=Value Facility
i.SetName("SampleLocation")
i.SetValue("4")
root.Insert(i) # XML DOM is used as the keeper of the structure
              # Insert() creates an XML node and inserts it under
              # the parent
# Domain
dm = XdmfDomain()
root.Insert(dm)
# Grid
g = XdmfGrid()
g.SetName("Structured Grid")
# Topology
t = g.GetTopology()
t.SetTopologyType(XDMF_3DCORECTMESH)
t.GetShapeDesc().SetShapeFromString('10 20 30')
# Geometry
geo = g.GetGeometry()
geo.SetGeometryType(XDMF_GEOMETRY_ORIGIN_DXDYDZ)
geo.SetOrigin(1, 2, 3)
geo.SetDxDyDz(0.1, 0.2, 0.3)
dm.Insert(g)
# Attribute
attr = XdmfAttribute()
attr.SetName("Pressure")
attr.SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_NODE);
attr.SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR);
# Insert the raw XML
attr.SetDataXml(DataXml)
g.Insert(attr)
# Update XML and Write Values to DataItems
root.Build() # DataItems > 100 values are heavy
print d.Serialize() # prints to stdout 

d.Write('SMesh.xmf') # write to file
