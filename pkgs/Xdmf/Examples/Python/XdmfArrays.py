#!/bin/env python

from Xdmf import *

def Expression(*args) :
    e = ''
    for arg in args :
        if hasattr(arg, 'GetTagName') :
            e += arg.GetTagName() + ' '
        else :
            e += arg + ' '
    return XdmfExpr(e)


if __name__ == '__main__' :
    a1 = XdmfArray()
    a1.SetNumberType(XDMF_FLOAT32_TYPE)
    a1.SetNumberOfElements(20)
    a1.Generate(1, 20)
    a2 = XdmfArray()
    a2.SetNumberType(XDMF_INT32_TYPE)
    a2.SetNumberOfElements(5)
    a2.Generate(2, 10)
    print 'a1 Values = ', a1.GetValues()
    print 'a1[2:10] = ' + Expression(a1 , '[ 2:10 ]').GetValues()
    print 'a2 Values = ', a2.GetValues()
    print 'a1[a2] = ' + Expression(a1 , '[', a2, ']').GetValues()
    print 'a1 + a2 = ' + Expression(a1 , ' + ', a2).GetValues()
    print 'a1 * a2 = ' + Expression(a1 , ' * ', a2).GetValues()
    a2.SetNumberType(XDMF_FLOAT32_TYPE)
    a2.SetNumberOfElements(20)
    a2.Generate(21, 40)
    print 'a2 Values = ', a2.GetValues()
    print 'a1 , a2 (Interlace) = ' + Expression(a1 , ' , ', a2).GetValues()
    print 'a1 , a2, a1 (Interlace) = ' + Expression(a1 , ' , ', a2, ' , ', a1).GetValues()
    print 'a1 ; a2 (Concat) = ' + Expression(a1 , ' ; ', a2).GetValues()
    print 'where(a1 > 10) = ' + Expression('Where( ', a1 , ' > 10)').GetValues()
    print 'a2[where(a1 > 10)] = ' + Expression(a2, '[Where( ', a1 , ' > 10)]').GetValues()
