// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ILU0Solver_H_
#define _ILU0Solver_H_

#include <vector>
#include "LinearSystem.h"
#include "MultiFieldReduction.h"
#include "LinearSolver.h"

using namespace std;

/**
 * Solve a linear system using Jacobi iterations
 * 
 */

class ILU0Solver : public LinearSolver
{
public:
  
  ILU0Solver();
  virtual ~ILU0Solver();
  virtual MFRPtr solve(LinearSystem & ls);

  virtual void cleanup();
  virtual void smooth(LinearSystem& ls);
  
  DEFINE_TYPENAME("ILU0Solver");
private:
  void doSweeps(LinearSystem& ls, const int nSweeps);
  ILU0Solver(const ILU0Solver&);
};

#endif

