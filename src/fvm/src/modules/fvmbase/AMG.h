#ifndef _AMG_H_
#define _AMG_H_

#include <vector>
#include "LinearSystem.h"
#include "LinearSolver.h"

#include "MultiFieldReduction.h"

using namespace std;

/**
 * Solve a linear system using algebraic multigrid.
 * 
 */

class AMG : public LinearSolver
{
public:

  enum CycleType
    {
      V_CYCLE,
      W_CYCLE,
      F_CYCLE
    };
  
  AMG();
  virtual ~AMG();

  DEFINE_TYPENAME("AMG");
  
  virtual MFRPtr solve(LinearSystem & ls);
  virtual void smooth(LinearSystem & ls);

  virtual void cleanup();
  
  // these parameters can be tuned.
  int maxCoarseLevels;
  int nPreSweeps;
  int nPostSweeps;
  int coarseGroupSize;
  double weightRatioThreshold;
  CycleType cycleType;
private:

  AMG(const AMG&);
  
  LinearSystem* _finestLinearSystem;
  vector<shared_ptr<LinearSystem> > _coarseLinearSystems;

  void createCoarseLevels();
  void doSweeps(const int nSweeps, const int level);
  void cycle(const CycleType cycleType, const int level);
};

#endif
