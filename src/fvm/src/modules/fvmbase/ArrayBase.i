
%{
#include "CException.h"
#include "ArrayBase.h"

#define PY_ARRAY_UNIQUE_SYMBOL MyNumpy
#ifndef SWIG_FILE_WITH_INIT
# define NO_IMPORT_ARRAY
#endif
#include <numpy/arrayobject.h>
  
%}

%typemap(out) ArrayBase*
{
  int dim = $1->getDimension();
  
  Py_intptr_t pyshape[10];
  int shape[10];
  $1->getShape(shape);
  for(int i=0; i<10; i++) pyshape[i]=shape[i];
  PrimType primType = $1->getPrimType();

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
  cerr << "numpy array creation" << endl;
  $result = PyArray_SimpleNewFromData(dim,pyshape,pyArrayType,$1->getData());
}
