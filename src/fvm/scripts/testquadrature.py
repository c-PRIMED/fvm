#!/usr/bin/env python

import fvm
fvm.set_atype('double')

import fvm.fvmbaseExt as fvmbaseExt
import fvm.Quadrature as quad
import fvm.KineticModel as esbgk
import fvm.DistFunctFields as f
import fvm.MacroParameters as macropr


#cartesian
quad1=quad.Quadrature(10,12,14,5.5,1.0) 

#spherical type, constant difference
quad2=quad.Quadrature(0,10,0,16,0,32)  
quad3=quad.Quadrature(8,8,1,4,1,6)

