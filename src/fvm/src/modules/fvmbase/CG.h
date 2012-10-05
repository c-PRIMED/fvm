// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _CG_H_
#define _CG_H_

#include <vector>
#include "LinearSystem.h"
#include "MultiFieldReduction.h"
#include "LinearSolver.h"

using namespace std;

/**
 * Solve a linear system using stabilized bi conjugate-gradient method.
 * 
 */

class CG : public LinearSolver
{
public:
  
  CG();
  virtual ~CG();
  virtual MFRPtr solve(LinearSystem & ls);

  virtual void cleanup();
  virtual void smooth(LinearSystem& ls);
  
  LinearSolver *preconditioner;
private:

  CG(const CG&);
};

#endif
