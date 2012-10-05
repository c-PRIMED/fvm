// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _LINEARSOLVER_H_
#define _LINEARSOLVER_H_


#include <vector>
#include "LinearSystem.h"
#include "MultiFieldReduction.h"

using namespace std;

class LinearSolver
{
public:

  LinearSolver() :
    nMaxIterations(100),
    verbosity(2),
    relativeTolerance(1e-8),
    absoluteTolerance(1e-50)
  {}

  virtual MFRPtr solve(LinearSystem & ls)=0;
  virtual void cleanup()=0;

  virtual void smooth(LinearSystem& ls) = 0;

  int nMaxIterations;
  int verbosity;
  double relativeTolerance;
  double absoluteTolerance;
};


#endif
