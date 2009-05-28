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
# Test XPath
node = dom.FindElementByPath('/Xdmf/Domain/Grid[10]')
print node
if(node) :
    print dom.Serialize(node)
node = dom.FindElementByPath('/Xdmf/Domain/Grid[11]')
print node
if(node) :
    print dom.Serialize(node)
node = dom.FindElementByPath('/Xdmf/Domain/Grid[Information]')
print node
if(node) :
    print dom.Serialize(node)
node = dom.FindElementByPath('/Xdmf/Domain/Grid[@Name = "GridNumber0020"]')
print node
if(node) :
    print dom.Serialize(node)
node = dom.FindElementByPath('/Xdmf/Domain/Grid[@Name = "GridNumber0020"]')
print node
if(node) :
    print dom.Serialize(node)
