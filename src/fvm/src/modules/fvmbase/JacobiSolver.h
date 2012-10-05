// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _JACOBISOLVER_H_
#define _JACOBISOLVER_H_

#include <vector>
#include "LinearSystem.h"
#include "MultiFieldReduction.h"
#include "LinearSolver.h"

using namespace std;

/**
 * Solve a linear system using Jacobi iterations
 * 
 */

class JacobiSolver : public LinearSolver
{
public:
  
  JacobiSolver();
  virtual ~JacobiSolver();
  virtual MFRPtr solve(LinearSystem & ls);

  virtual void cleanup();
  virtual void smooth(LinearSystem& ls);
  
  DEFINE_TYPENAME("JacobiSolver");
private:
  void doSweeps(LinearSystem& ls, const int nSweeps);
  JacobiSolver(const JacobiSolver&);
};

#endif

