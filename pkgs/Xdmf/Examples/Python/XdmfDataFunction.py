#!/bin/env python

from Xdmf import *

PointsTxt = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>


<Xdmf>
    <DataItem Format="XML" DataType="Float" Precision="8"
        Dimensions="2 10">
        0 1 2 3 4 5 6 7 8 9
        100 101 102 103 104 105 106 107 108 109
    </DataItem>
    <DataItem Name="MyFunction" ItemType="Function" 
        Function="10 + $0">
        <DataItem Reference="/Xdmf/DataItem[1]" />
    </DataItem>
    <DataItem Name="MyFunction1" ItemType="Function" 
        Function="$0[5:15]">
        <DataItem Reference="/Xdmf/DataItem[1]" />
    </DataItem>
    <DataItem Name="MyFunction2" ItemType="Function" 
        Function="JOIN($0 ; $1)">
        <DataItem Reference="/Xdmf/DataItem[1]" />
        <DataItem Reference="/Xdmf/DataItem[1]" />
    </DataItem>
    <!--
        start start     stride stride   count count
    -->
    <DataItem Name="MyFunction3" ItemType="HyperSlab"  >
        <DataItem Dimensions="6" NumberType="Int" Format="XML" >
            0 0  1 2  2 5
        </DataItem>
        <DataItem Reference="/Xdmf/DataItem[1]" />
    </DataItem>
    <!--
        Parametric Coordinates
    -->
    <DataItem Name="MyFunction4" ItemType="Coordinates"  >
        <DataItem Dimensions="8" NumberType="Int" Format="XML" >
            0 0  0 9  1 0 1 9
        </DataItem>
        <DataItem Reference="/Xdmf/DataItem[1]" />
    </DataItem>
</Xdmf>
"""

fd = open('Points.xmf', 'w')
fd.write(PointsTxt)
fd.close()
########
dom = XdmfDOM()
dom.Parse('Points.xmf')
# die = dom.FindElement('DataItem', 1)
die = dom.FindElementByPath('//DataItem[@Name = "MyFunction1"]')

di = XdmfDataItem()
di.DebugOn()
di.SetDOM(dom)
di.SetElement(die)
di.UpdateInformation()
di.Update()
print 'Dims = ', di.GetArray().GetShapeAsString()
print 'Values = ', di.GetDataValues()
