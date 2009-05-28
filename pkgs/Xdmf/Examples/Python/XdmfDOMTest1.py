#!/bin/env python

from Xdmf import *
import time

NumberOfGrids = 100

# build XML
XML = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>

<Xdmf Version="2.0">
<Domain>
"""

for i in range(NumberOfGrids) :
    XML += '<Grid Name="GridNumber%04d">\n<Information>this is some text for grid #%d</Information>\n</Grid>\n' % (i + 1, i + 1)

XML += """</Domain>
</Xdmf>"""

dom = XdmfDOM()
dom.Parse(XML)
dm = dom.FindElement('Domain')
print 'dm = ', dm
print 'Domain has %d Children' % dom.GetNumberOfChildren(dm)
# print dom.Serialize(dm)
ng = dom.FindNumberOfElements('Grid', dm)
print 'DOM has %d Grids' % ng
print 'Timing Access'
start = now = time.time()
for i in range(ng) :
    g = dom.FindElement('Grid', i, dm)
    # Name = dom.GetAttribute(g, 'Name')
    Name = dom.GetAttribute(g, 'Name')
    if (Name != "GridNumber%04d" % (i + 1)) :
        print 'Error in Name'
now = time.time()
print 'Accessed %d Grid Attributes in %g Seconds' % (ng, now - start)
g = dom.GetChild(ng / 2, dm)
print 'Middle Node = ', dom.GetAttribute(g, 'Name')
g = dom.FindElement('Grid', ng / 2, dm)
print 'Middle Grid = ', dom.GetAttribute(g, 'Name')
dom.DeleteNode(g)
ng = dom.FindNumberOfElements('Grid', dm)
print 'After Deletion DOM has %d Grids' % ng
g = dom.FindElement('Grid', ng / 2, dm)
print 'Test Serialize'
print dom.Serialize(g)
print 'Test Add'
dom.InsertFromString(g, '<Information Name="Last3"> XXYYZZ </Information>')
print dom.Serialize(g)
# print 'Test Serialize'
# print dom.Serialize()
print 'Test Set Attribute an CDATA'
info = dom.FindElement('Information', 1, g)
dom.Set(info, 'CDATA', 'AABBCC')
dom.Set(info, 'NewTag', 'NewValue')
print dom.Serialize(g)
print 'Testing Information'
i = XdmfInformation()
print 'Setting DOM'
i.SetDOM(dom)
print 'Setting Element'
i.SetElement(info)
print 'Update Information'
i.UpdateInformation()
print 'Ready'
print 'info.Get("NewTag") = ', i.Get('NewTag')
i.Set('NewTag', 'BrandNewValue')
i.Set('CData', 'DDEEFF')
print dom.Serialize(g)
print 'Name = ',i.GetName()
print 'Value = ',i.GetValue()
i.Set('Value', 'New Value')
print dom.Serialize(g)
print 'Value = ',i.GetValue()
print dom.Serialize(g)
i.SetName('Jerry')
print dom.Serialize(g)
print 'element Type = ', i.GetElementType()
i.SetValue('Last Value of the Day')
print dom.Serialize(g)
i.Build()
dom.Write('Jerry.xmf')
