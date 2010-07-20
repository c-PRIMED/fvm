#!/usr/bin/env python

import sys
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt

#atype = 'double'
#atype = 'tangent'
atype = 'pc_3_1'

if atype == 'double':
    import models_atyped_double as models
elif atype == 'tangent':
    import models_atyped_tangent_double as models
elif atype == 'pc_3_1':
    import models_atyped_pc_3_1 as models

    
import time
t0 = time.time()


k = models.PC_3_1()
k[0] = 0.0
k[1] = 0.1
model = models.OneDConductionA(10,k)

model.solve()
x = model.getSolution()
sx = models.getStdDev(x)
sxa = sx.asNumPyArray()
xa = x.asNumPyArray()
