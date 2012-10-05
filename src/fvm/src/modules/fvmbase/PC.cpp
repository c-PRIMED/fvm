// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "PC.h"


PCSet*
createPCSet(const int order, const int dim)
{
  PCSet *pc = new PCSet(order,dim,"HG");

  cout << "The number of PC terms in an expansion is "
       << pc->GetNumberPCTerms() << endl;


  pc->SetTaylorTolerance(1.e-15);
  return pc;
}
