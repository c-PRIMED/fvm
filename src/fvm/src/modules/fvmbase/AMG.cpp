#include "AMG.h"
#include "LinearSystemMerger.h"
#include "CRConnectivity.h"
#include <set>

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

int AMG::amg_indx = 0;

AMG::AMG() :
  maxCoarseLevels(20),
  nPreSweeps(0),
  nPostSweeps(1),
  coarseGroupSize(2),
  weightRatioThreshold(0.65),
  cycleType(V_CYCLE),
  _finestLinearSystem(0),
  _mergeLevel(-1),
  _commTarget(MPI::COMM_WORLD),
  _isMerge(false),
  _isCOMMWORLD(true)
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
      coarseLS.getDelta().zero();

      fineMatrix.injectResidual(fineLS.getCoarseIndex(),
                                fineLS.getResidual(),
                                coarseLS.getB());
 
      int nextLevel = level+1;
#ifdef  FVM_PARALLEL
      if ( level == _mergeLevel ) 
      {	
          _mergeLS->gatherB();
          nextLevel++;
      }
#endif

      cycle(cycleType,nextLevel);

      if (cycleType == W_CYCLE)
        cycle(W_CYCLE,nextLevel);
      else if (cycleType == F_CYCLE)
        cycle(V_CYCLE,nextLevel);



#ifdef  FVM_PARALLEL
      if ( level == _mergeLevel ) {
          _mergeLS->scatterDelta();
      }
#endif

      fineMatrix.correctSolution(fineLS.getCoarseIndex(),
                                 fineLS.getDelta(),
                                 coarseLS.getDelta())	;


  }
  doSweeps(nPostSweeps,level);
}

void
AMG::createCoarseLevels( )
{

  _coarseLinearSystems.clear();
  for(int n=0; n<maxCoarseLevels; n++)
  {
      LinearSystem& fineLS = (n == 0) ?
        *_finestLinearSystem : *_coarseLinearSystems[n-1];
      shared_ptr<LinearSystem>
        coarseLS(fineLS.createCoarse(coarseGroupSize,weightRatioThreshold));


#ifdef FVM_PARALLEL
     _coarseLinearSystems.push_back(coarseLS);

      if (verbosity > 1)
        cout << " proc_id = " << MPI::COMM_WORLD.Get_rank() << "  Created coarse level " << n << " of size "
             << coarseLS->getMatrix().getSize( _commTarget ) << endl;

      if ( coarseLS->getMatrix().getSize(_commTarget) <= 3  ){
        cout << "procID = " << MPI::COMM_WORLD.Get_rank() << " n = " << n << endl;
        break;
      }

      if ( coarseLS->getMatrix().getSize( _commTarget ) < _mergeLevelSize &&  _mergeLevel == -1  && _isMerge ){
         cout << " asla buraya girme " << endl;
         _mergeLevel = n;
         set<int> group;
         int size = MPI::COMM_WORLD.Get_size();
         for ( int i = 0; i < size; i++ )
            group.insert(i);

         _mergeLS = shared_ptr<LinearSystemMerger> ( new LinearSystemMerger( 0, group, *coarseLS ) );
         _mergeLS->gatherMatrix();
         _coarseLinearSystems.push_back( _mergeLS->getLS() );
          n++;
          flipComm(); //this change to COMM_WORLD to one-processedGROUP 
      } 
     if ( MPI::COMM_WORLD.Get_rank() == 0 ) cout << " n = " << n << " mergeLEvel = " << _mergeLevel << endl;
#endif 

#ifndef FVM_PARALLEL

    if ( verbosity > 1 )
        cout << "Created coarse level " << n << " of size " << coarseLS->getMatrix().getSize() << endl;
    _coarseLinearSystems.push_back(coarseLS);
    if ( coarseLS->getMatrix().getSize() <= 3 )
       break;
#endif

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

#ifdef FVM_PARALLEL
  if (verbosity >0 && MPI::COMM_WORLD.Get_rank() == 0 )
    cout << "0: " << *rNorm0 << "procID = " << MPI::COMM_WORLD.Get_rank() << endl;
#endif

#ifndef FVM_PARALLEL
   if ( verbosity > 0 )
      cout << "0: " << *rNorm0 << endl;
#endif

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
#ifdef FVM_PARALLEL
      if (verbosity >0 && MPI::COMM_WORLD.Get_rank() == 0 )
        cout << i << ": " << "procID = " << MPI::COMM_WORLD.Get_rank() <<  *rNorm << endl;
#endif

#ifndef FVM_PARALLEL
      if (verbosity >0  )
        cout << i << ": " <<  *rNorm << endl;
#endif


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
      createCoarseLevels( );
  }
  _finestLinearSystem = &ls;

  cycle(cycleType,0);
}


//this will flip COMM from COMM_WORLD to single-process goruped Communicator
void
AMG::flipComm()
{

#ifdef FVM_PARALLEL
   if ( _isCOMMWORLD ){
      int color   = ( MPI::COMM_WORLD.Get_rank() == 0 );
      int key     = MPI::COMM_WORLD.Get_rank();
      _commTarget = MPI::COMM_WORLD.Split( color, key );
      _isCOMMWORLD = !_isCOMMWORLD;
   } else {
      _commTarget = MPI::COMM_WORLD;
      _isCOMMWORLD = !_isCOMMWORLD;
   }
#endif



}
