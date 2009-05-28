#!/bin/env python

from Xdmf import *

txt = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">
    <Domain>
        <Grid Name="Polyvertex" GridType="Uniform">
            <Topology TopologyType="Polyvertex"
                Dimensions="3"
                NodesPerElement="3">
                <DataItem Format="XML"
                    Dimensions="3 3"
                    NumberType="Float">
                    0 1 2
                    4 5 6
                    8 9 10
                </DataItem>
            </Topology>
            <Geometry Type="XYZ">
                <DataItem Format="XML" Dimensions="4 4 3">
                0.0    0.0    0.0
                1.0    0.0    0.0
                1.0    1.0    0.0
                0.0    1.0    0.0

                0.0    0.0    2.0
                1.0    0.0    2.0
                1.0    1.0    2.0
                0.0    1.0    2.0

                0.0    0.0    3.0
                1.0    0.0    3.0
                1.0    1.0    3.0
                0.0    1.0    3.0

                0.0    0.0    4.0
                1.0    0.0    4.0
                1.0    1.0    4.0
                0.0    1.0    4.0
                </DataItem>
            </Geometry>
            <Attribute Name="Node Centered Values" Center="Node">
                <DataItem Format="XML" Dimensions="4 4">
                100 200 300 400
                500 600 600 700
                800 900 1000 1100
                1200 1300 1400 1500
                </DataItem>
            </Attribute>
            <Attribute Name="Cell Centered Values" Center="Cell">
                <DataItem Format="XML" Dimensions="3">
                100 200 400
                </DataItem>
            </Attribute>
            <Set SetType="Node">
                <DataItem NumberType="Int" Dimensions="3" Format="XML">
                    0 2 4
                </DataItem>
                <Attribute Name="Node Centered Set Values" Center="Node">
                    <DataItem Format="XML" Dimensions="3">
                    100 200 400
                    </DataItem>
                </Attribute>
            </Set>
            <Set SetType="Face">
                <DataItem NumberType="Int" Dimensions="2" Format="XML">
                    0 1
                </DataItem>
                <DataItem NumberType="Int" Dimensions="2" Format="XML">
                    0 0
                </DataItem>
            </Set>
            <Set SetType="Edge">
                <DataItem NumberType="Int" Dimensions="4" Format="XML">
                    0 0 1 1
                </DataItem>
                <DataItem NumberType="Int" Dimensions="4" Format="XML">
                    0 0 0 0
                </DataItem>
                <DataItem NumberType="Int" Dimensions="4" Format="XML">
                    0 1 0 2
                </DataItem>
                <Attribute Name="Set Values" Center="Edge">
                    <DataItem Format="XML" Dimensions="4">
                    100 200 400 600
                    </DataItem>
                </Attribute>
                <Attribute Name="Other Set Values" Center="Edge">
                    <DataItem Format="XML" Dimensions="3">
                    1000 2000 4000 6000
                    </DataItem>
                </Attribute>
            </Set>
        </Grid>
    </Domain>
</Xdmf>
"""



d = XdmfDOM()
d.Parse(txt)

dm = d.FindElement('Domain')
gnode = d.FindElement('Grid', 0, dm)

g = XdmfGrid()
g.SetDOM(d)
g.SetElement(gnode)
g.UpdateInformation()
g.Update()


print 'Grid Sets :', g.GetNumberOfSets()
for i in range(g.GetNumberOfSets()) :
    s = g.GetSets(i)
    print 'Set #', i
    print '     Type: ', s.GetSetTypeAsString()
    print '     Size: ', s.GetSize()
    s.Update()
    ids = s.GetIds()
    print '     # Ids : ', ids.GetNumberOfElements()
    print '     Ids : ', ids.GetValues()
    cellids = s.GetCellIds()
    print '     # Cell Ids : ', cellids.GetNumberOfElements()
    print '     Cell Ids : ', cellids.GetValues()
    faceids = s.GetFaceIds()
    print '     # Face Ids : ', faceids.GetNumberOfElements()
    print '     Face Ids : ', faceids.GetValues()
    print '     #Attributes: ', s.GetNumberOfAttributes()
    for j in range(s.GetNumberOfAttributes()) :
        a = s.GetAttribute(j)
        print '         Attribute #', j
        print '         Attribute values: ', a.GetAttributeTypeAsString()
        a.Update()
        v = a.GetValues()
        print '         Attribute type: ', v.GetValues()


ns = XdmfSet()
ns.SetDOM(d)
ns.SetSetType(XDMF_SET_TYPE_NODE)
a = XdmfArray()
a.SetNumberOfElements(4)
a.SetNumberType(XDMF_INT32_TYPE)
a.SetValues(0, "100 200 300 400")
ns.SetIds(a)
g.Insert(ns)
ns.Build()
print d.Serialize()

