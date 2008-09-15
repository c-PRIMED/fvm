#ifndef _AMG_H_
#define _AMG_H_

#include <vector>
#include "LinearSystem.h"
#include "MultiFieldReduction.h"

using namespace std;

/**
 * Solve a linear system using algebraic multigrid.
 * 
 */

class AMG
{
public:

  enum CycleType
    {
      V_CYCLE,
      W_CYCLE,
      F_CYCLE
    };
  
  AMG(LinearSystem & ls);
  void solve();


  // these parameters can be tuned.
  int nMaxCycles;
  int maxCoarseLevels;
  int nPreSweeps;
  int nPostSweeps;
  int coarseGroupSize;
  double weightRatioThreshold;
  CycleType cycleType;
  int verbosity;
  double relativeTolerance;
  double absoluteTolerance;
  
private:

  AMG(const AMG&);
  
  LinearSystem& _finestLinearSystem;
  vector<shared_ptr<LinearSystem> > _coarseLinearSystems;

  void createCoarseLevels();
  void doSweeps(const int nSweeps, const int level);
  void cycle(const CycleType cycleType, const int level);
};

#endif
