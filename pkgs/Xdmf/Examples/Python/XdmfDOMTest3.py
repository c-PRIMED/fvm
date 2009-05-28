#!/bin/env python

from Xdmf import *

d = XdmfDOM()
root = d.Create()
dm = d.InsertNew(root, 'Domain')
g = d.InsertNew(dm, 'Grid')
geo = d.InsertNew(g, 'Geometry')
ds  = d.InsertNew(geo, 'DataStructure')
print d.Serialize(root)

d.Write('junk.xmf')
