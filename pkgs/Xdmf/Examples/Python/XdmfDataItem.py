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
                    <DataItem Format="XML" DataType="Float" Precision="8"
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
                    </DataItem>
                    <DataItem Reference="/Xdmf/Domain/Grid[@Name='Shot Points']/Geometry/DataItem[1]"/>
                    <DataItem Reference="XML">
                        /Xdmf/Domain/Grid[@Name="Shot Points"]/Geometry/DataItem[2]
                    </DataItem>
                    <DataItem Reference="/Xdmf/Domain/Grid[@Name='Shot Points']/Geometry/DataItem[1]/../DataItem[2]"/>
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
dse = dom.FindElement('DataItem', 0, geo)

ds = XdmfDataItem()
ds.DebugOn()
ds.SetDOM(dom)
ds.SetElement(dse)
ds.UpdateInformation()
ds.Update()
print 'Values = ', ds.GetDataValues()

dsre = dom.FindElement('DataItem', 1, geo)
dsr = XdmfDataItem()
dsr.SetDOM(dom)
dsr.DebugOn()
# This will clear the reference and cause I/O
dsr.SetElement(dsre)
dsr.UpdateInformation()
dsr.Update()
if(dsr.GetIsReference()) :
    print 'Getting Values'
    print 'Values = ', dsr.GetDataValues()
else :
    print 'This is not a reference'

dsre = dom.FindElement('DataItem', 3, geo)
dsr = XdmfDataItem()
dsr.SetDOM(dom)
dsr.DebugOn()
# This will clear the reference and cause I/O
dsr.SetElement(dsre)
# Cause an potential dangling reference
# dsr.SetCopyReferenceData(0)
print "UI"
dsr.UpdateInformation()
print "U"
dsr.Update()
# Cause an dangling reference by deleting refernced object
# ds = 0
if(dsr.GetIsReference()) :
    print 'Getting Values'
    print 'Values = ', dsr.GetDataValues()
else :
    print 'This is not a reference'


