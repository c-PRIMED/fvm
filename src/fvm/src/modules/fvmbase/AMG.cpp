#include "AMG.h"

AMG::AMG() :
  maxCoarseLevels(20),
  nPreSweeps(0),
  nPostSweeps(1),
  coarseGroupSize(2),
  weightRatioThreshold(0.65),
  cycleType(V_CYCLE),
  _finestLinearSystem(0)
{
  logCtor();
}

AMG::~AMG()
{
  logDtor();
}

void
AMG::doSweeps(const int nSweeps, const int level)
{
  LinearSystem& ls = (level == 0) ?
    *_finestLinearSystem : *_coarseLinearSystems[level-1];

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
        *_finestLinearSystem : *_coarseLinearSystems[level-1];
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
  _coarseLinearSystems.clear();
  for(int n=0; n<maxCoarseLevels; n++)
  {
      LinearSystem& fineLS = (n == 0) ?
        *_finestLinearSystem : *_coarseLinearSystems[n-1];

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
AMG::cleanup()
{
  _finestLinearSystem = 0;
  _coarseLinearSystems.clear();
}

MFRPtr
AMG::solve(LinearSystem & ls)
{
  if (_finestLinearSystem != &ls)
  {
      _finestLinearSystem = &ls;
      createCoarseLevels();
  }
  _finestLinearSystem = &ls;


  const MultiFieldMatrix& finestMatrix = _finestLinearSystem->getMatrix();
  finestMatrix.computeResidual(_finestLinearSystem->getDelta(),
                               _finestLinearSystem->getB(),
                               _finestLinearSystem->getResidual());

  MFRPtr rNorm0(_finestLinearSystem->getResidual().getOneNorm());

  if (verbosity >0)
    cout << "0: " << *rNorm0 << endl;
  
  if (*rNorm0 < absoluteTolerance )
    return rNorm0;
  
  for(int i=1; i<nMaxIterations; i++)
  {
      cycle(cycleType,0);
      finestMatrix.computeResidual(_finestLinearSystem->getDelta(),
                                   _finestLinearSystem->getB(),
                                   _finestLinearSystem->getResidual());
      
      MFRPtr rNorm(_finestLinearSystem->getResidual().getOneNorm());
      MFRPtr normRatio((*rNorm)/(*rNorm0));
      if (verbosity >0)
        cout << i << ": " << *rNorm << endl;
      if (*rNorm < absoluteTolerance || *normRatio < relativeTolerance)
        break;
  }

  return rNorm0;
}

void
AMG::smooth(LinearSystem & ls)
{
  if (_finestLinearSystem != &ls)
  {
      _finestLinearSystem = &ls;
      createCoarseLevels();
  }
  _finestLinearSystem = &ls;

  cycle(cycleType,0);
}
