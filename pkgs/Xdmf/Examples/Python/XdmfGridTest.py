#!/bin/env python

from Xdmf import *

PointsTxt = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>


<Xdmf>
    <Domain>
        <Grid Name="TestGrid">
            <Topology Type="Quadrilateral" NumberOfElements="2" >
                <DataItem Format="XML" DataType="Float"
                    Dimensions="2 4">
                    0 1 6 5
                    16 17 22 21
                </DataItem>
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
                    <DataItem Reference="/Xdmf/Domain/Grid[@Name='TestGrid']/Geometry/DataItem[1]"/>
                    <DataItem Reference="XML">
                        /Xdmf/Domain/Grid[@Name="TestGrid"]/Geometry/DataItem[2]
                    </DataItem>
                    <DataItem Reference="/Xdmf/Domain/Grid[@Name='TestGrid']/Geometry/DataItem[1]/../DataItem[2]"/>
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

te = dom.FindElementByPath('/Xdmf/Domain/Grid/Topology')
top = XdmfTopology()
top.DebugOn()
top.SetDOM(dom)
top.SetElement(te)
top.Update()
conn = top.GetConnectivity()
print 'Values = ', conn.GetValues()

ge = dom.FindElementByPath('/Xdmf/Domain/Grid/Geometry')
geo = XdmfGeometry()
geo.SetDOM(dom)
geo.SetElement(ge)
geo.Update()
points = geo.GetPoints()
print 'Geo Type = ', geo.GetGeometryTypeAsString(), ' # Points = ', geo.GetNumberOfPoints()
print 'Points = ', points.GetValues(0, 6)

