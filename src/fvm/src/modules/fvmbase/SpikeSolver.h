// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _SpikeSolver_H_
#define _SpikeSolver_H_

#include <vector>
#include "LinearSystem.h"
#include "MultiFieldReduction.h"
#include "LinearSolver.h"

class SpikeStorage;
using namespace std;

/**
 * Solve a linear system using Jacobi iterations
 * 
 */

class SpikeSolver : public LinearSolver
{
public:
  
  SpikeSolver( const SpikeStorage& spike_storage);
  virtual ~SpikeSolver();
  virtual MFRPtr solve(LinearSystem & ls);

  virtual void cleanup();
  virtual void smooth(LinearSystem& ls);
  
  DEFINE_TYPENAME("SpikeSolver");
private:
  void doSweeps(LinearSystem& ls, const int nSweeps);
  SpikeSolver(const SpikeSolver&);
  const SpikeStorage& _spikeStorage;
};

#endif

