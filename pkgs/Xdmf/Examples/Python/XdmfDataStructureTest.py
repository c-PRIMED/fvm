#!/bin/env python

from Xdmf import *

PointsTxt = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>


<Xdmf>
    <Domain>
        <Grid Name="Shot Points">
            <Topology Type="PolyVertex" NodesPerElement="1" NumberOfElements="5" >
            </Topology>
            <Geometry Type="XYZ">
                    <DataStructure Format="XML" DataType="Float" Precision="8"
                            Dimensions="3 2 5 3">
                            10.0        0.0         0.0
                            10.0        1.0         0.0

                            10.0        2.0         0.0
                            10.0        3.0         0.0

                            10.0        4.0         0.0
                            10.0        0.0         1.0

                            10.0        1.0         1.0
                            10.0        2.0         1.0

                            10.0        3.0         1.0
                            10.0        4.0         1.0
                            10.0        0.0         0.0
                            10.0        1.0         0.0

                            10.0        2.0         0.0
                            10.0        3.0         0.0

                            10.0        4.0         0.0
                            10.0        0.0         1.0

                            10.0        1.0         1.0
                            10.0        2.0         1.0

                            10.0        3.0         1.0
                            10.0        4.0         1.0
                            10.0        0.0         0.0
                            10.0        1.0         0.0

                            10.0        2.0         0.0
                            10.0        3.0         0.0

                            10.0        4.0         0.0
                            10.0        0.0         1.0

                            10.0        1.0         1.0
                            10.0        2.0         1.0

                            10.0        3.0         1.0
                            10.0        4.0         1.0
                    </DataStructure>
            </Geometry>
        </Grid>
    </Domain>
</Xdmf>
"""

fd = open('Points.xmf', 'w')
fd.write(PointsTxt)
fd.close()
########
dom = XdmfDOM()
dom.Parse('Points.xmf')
dm = dom.FindElement('Domain')
g = dom.FindElement('Grid', 0, dm)
geo = dom.FindElement('Geometry', 0, g)
dse = dom.FindElement('DataStructure', 0, geo)

ds = XdmfDataStructure()
ds.SetDOM(dom)
ds.SetElement(dse)
ds.UpdateInformation()
print 'Rank ', ds.GetRank()
print 'Dimensions ', ds.GetDimensions()
dd = ds.GetDataDesc()
print "Number Type ", dd.GetNumberTypeAsString()
print "Number of Elements ", dd.GetNumberOfElements()
print "Name = ", ds.GetName()
# ds.DebugOn()
ds.Update()
ds.Build()
print ds.Serialize()
print 'Values = ', ds.GetDataValues()
a = ds.GetArray()
a.SetNumberType(XDMF_FLOAT32_TYPE)
a.SetNumberOfElements(27)
# a.SetShapeFromString("5 5")
print "Shape = ", a.GetShapeAsString()
a.Generate(1, 10)
ds.DebugOn()
ds.SetFormat(XDMF_FORMAT_HDF)
ds.SetHeavyDataSetName('Jerry.h5:/NewData')
ds.Build()
print ds.Serialize()
dom.Write('NewPoints.xmf')

dsr = XdmfDataStructure()
dsr.DebugOn()
# This will clear the reference and cause I/O
dsr.SetElement(ds.GetElement())
dsr.SetDOM(dom)
dsr.Reference(ds.GetElement())
if(dsr.GetIsReference()) :
    print 'Values = ', dsr.GetDataValues()
else :
    print 'This is not a reference'


