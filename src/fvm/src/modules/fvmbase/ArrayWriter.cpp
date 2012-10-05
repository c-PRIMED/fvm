// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "ArrayWriter.h"

DEFINE_HT(ArrayWriter);

void
ArrayWriter::addMethods()
{
  INHERIT_METHODS(PyCreatable);
  ADD_METHOD(ArrayWriter,write);
  ADD_METHOD(ArrayWriter,writePartial);
  ADD_METHOD(ArrayWriter,writeMasked);
}

PyReturn
ArrayWriter::write(FILE *fp)
{
  
  write(fp);
    
  return PyReturn();
}

PyReturn
ArrayWriter::writePartial(PyArgsIn args)
{
  FILE *fp = args.getFILE(0);
  const int iBeg = args.getInt(1);
  const int count = args.getInt(2);
  write(fp,iBeg,count);
    
  return PyReturn();
}

PyReturn
ArrayWriter::writeMasked(PyArgsIn args)
{
  FILE *fp = args.getFILE(0);
  const int iBeg = args.getInt(1);
  const int count = args.getInt(2);
  const Array<bool>& mask = args.getRef<Array<bool> >(3);
  write(fp,iBeg,count,&mask);
    
  return PyReturn();
}
