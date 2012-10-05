// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _DirectSolver_H_
#define _DirectSolver_H_

#include <vector>
#include <cmath>
#include "LinearSystem.h"
#include "LinearSolver.h"

#include "MultiFieldReduction.h"

using namespace std;

/**
 * Solve a linear system using UMFPack
 * 
 */

class DirectSolver : public LinearSolver
{
public:

  
  DirectSolver();
  virtual ~DirectSolver();

  virtual void cleanup();
  
  virtual MFRPtr solve(LinearSystem & ls);
  virtual void smooth(LinearSystem & ls);


private:

  DirectSolver(const DirectSolver&);

};

#endif
