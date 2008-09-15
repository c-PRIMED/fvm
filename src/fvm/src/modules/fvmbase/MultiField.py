
from enthought.traits.api import CInt, Bool, Tuple, method, Trait, This


from Array import Array, NewArray

from HasCImpl import HasCImpl, CMethod
from StorageSite import StorageSite
from HasState import HasState
from PyFactory import PyFactory
from IContainer import IContainer
from Field import Field
from Field import Field
from MultiFieldReduction import MultiFieldReduction

import baseExt

from CList import CListTrait

ArrayIndexTrait = Tuple(Field,StorageSite)

class MultiField(IContainer):

    arrays = CListTrait('ArrayBase*')
    arrayIndices = CListTrait('pair<const Field*,const StorageSite*>')

    CMethod("__setitem__",None,ArrayIndexTrait,Array)
    CMethod("__getitem__",Array,ArrayIndexTrait)
    CMethod("__delitem__",None,ArrayIndexTrait)
    CMethod("__contains__",Bool,ArrayIndexTrait)
    CMethod("__iadd__", This,This)
    CMethod("__isub__", This,This)
    CMethod("__idiv__", This,MultiFieldReduction)
    CMethod("__imul__", This,MultiFieldReduction)
    CMethod('syncLocal',None,None)
    CMethod('getOneNorm',MultiFieldReduction,None)
    CMethod('reduceSum',MultiFieldReduction,None)
    CMethod('dotWith',MultiFieldReduction,This)
    CMethod('saxpy',This,MultiFieldReduction,This)
    CMethod('msaxpy',This,MultiFieldReduction,This)
    
    method(None,Field,CListTrait('StorageSite*'),Trait(None,CInt))
    def addFieldArrays(self,field,sites,mdDim=None):
        for site in sites:
            array=field[site]
            if mdDim is not None:
                array = array[mdDim]
            self[(field,site)]=array

    def keys(self):
        return self.arrayIndices
    
    def values(self):
        return self.arrays

    def sum(self):
        return [ a.sum() for a in self.arrays]

    """def oneNorm(self):
        r = []
        for a in self.arrays:
            n = a.getOneNorm()
            r.append(n)
        return r"""

    def extract(self,keys):
        other = MultiField()
        for k in keys:
            other[k]=self[k]
            del self[k]
        return other
    
    def merge(self,other):
        nl = len(other.arrayIndices)
        for n in range(0,nl):
            k = other.arrayIndices[n]
            a = other.arrays[n]
            self[k] = a
        #for k,a in zip(other.arrayIndices,other.arrays):
        #    self.__delete__(k)
        #    self[k]=a
            
PyFactory.registerCtor(MultiField)

MultiFieldTrait = Trait(MultiField)
