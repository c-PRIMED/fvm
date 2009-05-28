#!/bin/env python

from Xdmf import *

PointsTxt = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>


<Xdmf>
    <Domain>
        <Grid Name="Tree1" GridType="Tree">
          <Grid Name="Tree2" GridType="Tree">
            <Grid Name="FrontBack">
                <Topology Type="Quadrilateral" NumberOfElements="2" >
                    <DataItem Format="XML" DataType="Int" Dimensions="2 4">
                    0 1 2 3
                    4 5 6 7
                    </DataItem>
                </Topology>
                <Geometry Type="XYZ">
                    <DataItem Format="XML" Dimensions="2 4 3">
                    0.0    0.0    0.0
                    1.0    0.0    0.0
                    1.0    1.0    0.0
                    0.0    1.0    0.0

                    0.0    0.0    2.0
                    1.0    0.0    2.0
                    1.0    1.0    2.0
                    0.0    1.0    2.0
                    </DataItem>
                </Geometry>
            </Grid>
            <Grid Name="TopBottom">
                <Topology Type="Quadrilateral" NumberOfElements="2" >
                    <DataItem Format="XML" DataType="Int" Dimensions="2 4">
                    0 1 5 4
                    2 3 7 6
                    </DataItem>
                </Topology>
                <Geometry Type="XYZ">
                    <DataItem Format="XML" Dimensions="2 4 3">
                    0.0    0.0    0.0
                    1.0    0.0    0.0
                    1.0    1.0    0.0
                    0.0    1.0    0.0

                    0.0    0.0    2.0
                    1.0    0.0    2.0
                    1.0    1.0    2.0
                    0.0    1.0    2.0
                    </DataItem>
                </Geometry>
            </Grid>
          </Grid>
          <Grid Name="Tree3" GridType="Tree">
            <Grid Name="FrontBack">
                <Topology Type="Quadrilateral" NumberOfElements="2" >
                    <DataItem Format="XML" DataType="Int" Dimensions="2 4">
                    0 1 2 3
                    4 5 6 7
                    </DataItem>
                </Topology>
                <Geometry Type="XYZ">
                    <DataItem Format="XML" Dimensions="2 4 3">
                    0.0    0.0    0.0
                    1.0    0.0    0.0
                    1.0    1.0    0.0
                    0.0    1.0    0.0

                    0.0    0.0    2.0
                    1.0    0.0    2.0
                    1.0    1.0    2.0
                    0.0    1.0    2.0
                    </DataItem>
                </Geometry>
            </Grid>
            <Grid Name="TopBottom">
                <Topology Type="Quadrilateral" NumberOfElements="2" >
                    <DataItem Format="XML" DataType="Int" Dimensions="2 4">
                    0 1 5 4
                    2 3 7 6
                    </DataItem>
                </Topology>
                <Geometry Type="XYZ">
                    <DataItem Format="XML" Dimensions="2 4 3">
                    0.0    0.0    0.0
                    1.0    0.0    0.0
                    1.0    1.0    0.0
                    0.0    1.0    0.0

                    0.0    0.0    2.0
                    1.0    0.0    2.0
                    1.0    1.0    2.0
                    0.0    1.0    2.0
                    </DataItem>
                </Geometry>
            </Grid>
          </Grid>
        </Grid>
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

def PrintGrid(grid) :
    print 'Grid ', grid.GetName(), ' Type = ', grid.GetGridTypeAsString()
    if(grid.IsUniform()) :
        top = grid.GetTopology()
        top.DebugOn()
        conn = top.GetConnectivity()
        print 'Connectivity = ', conn.GetValues()

        geo = grid.GetGeometry()
        points = geo.GetPoints()
        print 'Geo Type = ', geo.GetGeometryTypeAsString(), ' # Points = ', geo.GetNumberOfPoints()
        print 'Points = ', points.GetValues()
    else :
        nc = grid.GetNumberOfChildren()
        for i in range(nc) :
            g = grid.GetChild(i)
            PrintGrid(g)

fd = open('Points.xmf', 'w')
fd.write(PointsTxt)
fd.close()
########
dom = XdmfDOM()
dom.Parse('Points.xmf')

ge = dom.FindElementByPath('/Xdmf/Domain/Grid[@Name="Tree1"]')
grid = XdmfGrid()
grid.SetDOM(dom)
grid.SetElement(ge)
grid.UpdateInformation()
grid.Update()
PrintGrid(grid)


