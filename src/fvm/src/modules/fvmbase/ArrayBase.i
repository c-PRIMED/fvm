
%{
#include "CException.h"
#include "ArrayBase.h"

#define PY_ARRAY_UNIQUE_SYMBOL MyNumpy
#ifndef SWIG_FILE_WITH_INIT
# define NO_IMPORT_ARRAY
#endif
#include <numpy/arrayobject.h>
  
  %}

%include "boost_shared_ptr.i"
SWIG_SHARED_PTR(ArrayPtr,ArrayBase)

class ArrayBase
{
public:
  int getDimension();
  boost::shared_ptr<ArrayBase> newSizedClone(const int size);
  
  %extend
  {
    PyObject* asNumPyArray()
    {
        int dim = self->getDimension();
  
        Py_intptr_t pyshape[10];
        int shape[10];
        self->getShape(shape);
        for(int i=0; i<10; i++) pyshape[i]=shape[i];
        PrimType primType = self->getPrimType();

        int pyArrayType;

        switch(primType)
        {
        case PRIM_TYPE_BOOL:
          pyArrayType = NPY_BOOL;
          break;
        case PRIM_TYPE_INT:
          pyArrayType = NPY_INT32;
          break;
        case PRIM_TYPE_FLOAT:
          pyArrayType = NPY_FLOAT;
          break;
        default:
          pyArrayType = NPY_DOUBLE;
          break;
        }
        return PyArray_SimpleNewFromData(dim,pyshape,pyArrayType,self->getData());
    }
  }
private:
  ArrayBase();
};



