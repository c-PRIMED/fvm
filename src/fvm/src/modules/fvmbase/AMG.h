#ifndef _AMG_H_
#define _AMG_H_

#include <vector>
#include <cmath>
#include "LinearSystem.h"
#include "LinearSolver.h"

#include "MultiFieldReduction.h"
#include "LinearSystemMerger.h"

class LinearSystemMerger;

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

  enum SmootherType
    {
      GAUSS_SEIDEL,
      JACOBI
    };
  
  AMG();
  virtual ~AMG();

  DEFINE_TYPENAME("AMG");
  
  virtual MFRPtr solve(LinearSystem & ls);
  virtual void smooth(LinearSystem & ls);


  virtual void setMergeLevelSize(int ls_size){ 
#ifdef FVM_PARALLEL
      _mergeLevelSize = ceil(  double(ls_size) / double(MPI::COMM_WORLD.Get_size()) ); 
      _isMerge = true;
#endif

#ifndef FVM_PARALLEL
      _mergeLevelSize = ls_size;
      cout << " you can not set mergeLevelSize in serial version !!!!!!!!!!! " << endl;
      abort();
#endif
}

  virtual void cleanup();
  
  // these parameters can be tuned.
  int maxCoarseLevels;
  int nPreSweeps;
  int nPostSweeps;
  int coarseGroupSize;
  double weightRatioThreshold;
  CycleType cycleType;
  SmootherType smootherType;
  
private:

  AMG(const AMG&);
  
  LinearSystem* _finestLinearSystem;
  vector<shared_ptr<LinearSystem> > _coarseLinearSystems;
  shared_ptr<LinearSystemMerger>   _mergeLS;

  void  createCoarseLevels( );
  void  doSweeps( const int nSweeps, const int level );
  void  cycle( const CycleType cycleType, const int level );
  void  flipComm();

  static int amg_indx;
 
  int _mergeLevelSize; //where 
  int _mergeLevel; 
  bool _isMerge;

#ifdef FVM_PARALLEL
  MPI::Intracomm _commTarget;
#endif


  bool _isCOMMWORLD;
};

#endif
