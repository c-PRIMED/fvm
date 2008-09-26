from enthought.traits.api import This
import numpy
from Array import ArrayBase

from PyFactory import PyFactory
from HasCImpl import HasCImpl, CMethod
from Field import Field
from Array import Array
import baseExt

from CDict import CDictTrait

from _baseExt import Vector_double_3_ as VecD3

class MultiFieldReduction(HasCImpl):

    arrays = CDictTrait('const Field*,ArrayBase*')

    CMethod("__getitem__",ArrayBase,Field)
    CMethod("__setitem__",None,Field,ArrayBase)
    CMethod("__delitem__",None,Field)
    CMethod("__iadd__",This,This)
    #CMethod("reduceSum",None,None)

    def newCopy(self):
        c = MultiFieldReduction()
        for k,a in self.arrays.iteritems():
            c[k] = a.copy()
        return c

    def __repr__(self):
        p = '['
        for  k,array in self.arrays.iteritems():
            a = array.asATypeNumPyArray()[0]
            p += k.name + " : " + str(a) + ' '
        p += ']'
        return p
    
    def __div__(self,b):
        r = MultiFieldReduction()
        for k,a in self.arrays.iteritems():
            r[k] = a.newCopy()
            a0 = a.asATypeNumPyArray()[0]
            b0 = b[k].asATypeNumPyArray()[0]
            r[k].asATypeNumPyArray()[0] = a0/(b0+1e-108)
        return r
    
    def __mul__(self,b):
        p = MultiFieldReduction()
        for k,a in self.arrays.iteritems():
            p[k] = a.newCopy()
            a0 = a.asATypeNumPyArray()[0]
            b0 = b[k].asATypeNumPyArray()[0]
            p[k].asATypeNumPyArray()[0] = a0*b0
        return p

    def isSmall(self):
        for  k,array in self.arrays.iteritems():
            a = array.asATypeNumPyArray()[0]
            if a >= 0:
                return False
        return True

    def reduceSum(self):
        #return
        asum = 0
        for k,array in self.arrays.iteritems():
            a0 = array.asATypeNumPyArray()[0]
            try:
                for i in range(0,len(a0)):
                    asum += a0[i]
            except:
                asum += a0
        
        for k,array in self.arrays.iteritems():
            av = array#.asATypeNumPyArray()
            try:
                for i in range(0,len(av[0])):
                    av[0][i] = asum
            except:
                av[0] = asum
                
PyFactory.registerCtor(MultiFieldReduction)
