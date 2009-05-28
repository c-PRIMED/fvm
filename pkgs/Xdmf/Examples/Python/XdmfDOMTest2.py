#!/bin/env python

from Xdmf import *

dom = XdmfDOM()
# dom.SetInputFileName('Example2.xmf')
dom.Parse('Example2.xmf')
print dom.Serialize()
