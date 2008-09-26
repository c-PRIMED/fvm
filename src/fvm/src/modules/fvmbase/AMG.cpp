#include "AMG.h"

AMG::AMG(LinearSystem& ls) :
  nMaxCycles(100),
  maxCoarseLevels(20),
  nPreSweeps(1),
  nPostSweeps(0),
  coarseGroupSize(2),
  weightRatioThreshold(0.65),
  cycleType(V_CYCLE),
  verbosity(2),
  relativeTolerance(1e-8),
  absoluteTolerance(1e-16),
  _finestLinearSystem(ls)
{
  createCoarseLevels();
}
    

void
AMG::doSweeps(const int nSweeps, const int level)
{
  LinearSystem& ls = (level == 0) ?
    _finestLinearSystem : *_coarseLinearSystems[level-1];

  const MultiFieldMatrix& m = ls.getMatrix();
  MultiField& delta = ls.getDelta();
  const MultiField& b = ls.getB();
  MultiField& r = ls.getResidual();

  for(int i=0; i<nSweeps; i++)
  {
      m.forwardGS(delta,b,r);
      m.reverseGS(delta,b,r);
  }
}

void
AMG::cycle(const CycleType cycleType, const int level)
{
  doSweeps(nPreSweeps,level);

  if (level < (int)_coarseLinearSystems.size())
  {
      LinearSystem& fineLS = (level == 0) ?
        _finestLinearSystem : *_coarseLinearSystems[level-1];
      LinearSystem& coarseLS = *_coarseLinearSystems[level];

      MultiFieldMatrix& fineMatrix = fineLS.getMatrix();

      fineMatrix.computeResidual(fineLS.getDelta(),fineLS.getB(),fineLS.getResidual());
      coarseLS.getB().zero();
      fineMatrix.injectResidual(fineLS.getCoarseIndex(),
                                fineLS.getResidual(),
                                coarseLS.getB());
      coarseLS.getDelta().zero();

      cycle(cycleType,level+1);

      if (cycleType == W_CYCLE)
        cycle(W_CYCLE,level+1);
      else if (cycleType == F_CYCLE)
        cycle(V_CYCLE,level+1);

      fineMatrix.correctSolution(fineLS.getCoarseIndex(),
                                 fineLS.getDelta(),
                                 coarseLS.getDelta());
  }


  doSweeps(nPostSweeps,level);
}

void
AMG::createCoarseLevels()
{
  for(int n=0; n<maxCoarseLevels; n++)
  {
      LinearSystem& fineLS = (n == 0) ?
        _finestLinearSystem : *_coarseLinearSystems[n-1];

      shared_ptr<LinearSystem>
        coarseLS(fineLS.createCoarse(coarseGroupSize,weightRatioThreshold));

      if (verbosity > 1)
        cout << "Created coarse level " << n << " of size "
             << coarseLS->getMatrix().getSize() << endl;

      _coarseLinearSystems.push_back(coarseLS);
      if (coarseLS->getMatrix().getSize() <= 3)
        break;
      
  }           
}

void
AMG::solve()
{
  const MultiFieldMatrix& finestMatrix = _finestLinearSystem.getMatrix();
  finestMatrix.computeResidual(_finestLinearSystem.getDelta(),
                               _finestLinearSystem.getB(),
                               _finestLinearSystem.getResidual());

  MFPtr rNorm0(_finestLinearSystem.getResidual().getOneNorm());

  if (verbosity >0)
    cout << "0: " << *rNorm0 << endl;
  
  if (*rNorm0 < absoluteTolerance )
    return;
  
  for(int i=1; i<nMaxCycles; i++)
  {
      cycle(cycleType,0);
      finestMatrix.computeResidual(_finestLinearSystem.getDelta(),
                                   _finestLinearSystem.getB(),
                                   _finestLinearSystem.getResidual());
      
      MFPtr rNorm(_finestLinearSystem.getResidual().getOneNorm());
      MFPtr normRatio((*rNorm)/(*rNorm0));
      if (verbosity >0)
        cout << i << ": " << *rNorm << endl;
      if (*rNorm < absoluteTolerance || *normRatio < relativeTolerance)
        break;
  }
}
