// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _MATRIXJML_H_
#define _MATRIXJML_H_

template<class T>
class MatrixJML
{
 public:

  MatrixJML() {}
  virtual ~MatrixJML() {}
  virtual T& getElement(const int i,const int j)=0;
  virtual void zero()=0;

};

#endif
