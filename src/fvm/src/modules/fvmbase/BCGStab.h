// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _BCGSTab_H_
#define _BCGSTab_H_

#include <vector>
#include "LinearSystem.h"
#include "MultiFieldReduction.h"
#include "LinearSolver.h"

using namespace std;

/**
 * Solve a linear system using stabilized bi conjugate-gradient method.
 * 
 */

class BCGStab : public LinearSolver
{
public:
  
  BCGStab();
  virtual ~BCGStab();
  virtual MFRPtr solve(LinearSystem & ls);

  virtual void cleanup();
  virtual void smooth(LinearSystem& ls);
  
  int getTotalIterations() const { return _totalIterations;}
  
  DEFINE_TYPENAME("BCGStab");
  
  LinearSolver *preconditioner;
private:

  BCGStab(const BCGStab&);
  
  int _totalIterations;
};

#endif
