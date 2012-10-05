// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include "AMG.h"
#include "LinearSystemMerger.h"
#include "CRConnectivity.h"
#include <set>
int AMG::amg_indx = 0;

AMG::AMG() :
  maxCoarseLevels(30),
  nPreSweeps(0),
  nPostSweeps(1),
  coarseGroupSize(2),
  weightRatioThreshold(0.65),
  cycleType(V_CYCLE),
  smootherType(GAUSS_SEIDEL),
  scaleCorrections(true),
  _finestLinearSystem(0),
  _mergeLevelSize(0),
  _mergeLevel(-1),
  _isMerge(false),
  _totalIterations(0),
#ifdef FVM_PARALLEL
  _commTarget(MPI::COMM_WORLD),
#endif
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
      if (smootherType == GAUSS_SEIDEL)
      {
          m.forwardGS(delta,b,r);
          m.reverseGS(delta,b,r);
      }
      else
      {
          m.Jacobi(delta,b,r);
          m.Jacobi(delta,b,r);
      }
  }
}

void
AMG::cycle( CycleType cycleType, const int level)
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
      if ( nextLevel == _mergeLevel ) 
      {	
          _mergeLS->gatherB();
          nextLevel++;
      }
#endif
     //if ( level == 4 ) 
     //   cycleType = W_CYCLE;

     cycle(cycleType,nextLevel);

      if (cycleType == W_CYCLE)
        cycle(W_CYCLE,nextLevel);
      else if (cycleType == F_CYCLE)
        cycle(V_CYCLE,nextLevel);

#ifdef  FVM_PARALLEL
      if ( level+1 == _mergeLevel ) {
          _mergeLS->scatterDelta();
      }
#endif

      MFRPtr scale;
      if (coarseLS.isSymmetric && scaleCorrections)
      {
          const MultiField& x = coarseLS.getDelta();
          const MultiField& b = coarseLS.getB();
          const MultiFieldMatrix& A = coarseLS.getMatrix();
          MFRPtr xb =  x.dotWith(b);
          xb->reduceSum();

          MFRPtr mxb = -(*xb);
          MFRPtr xTAx = A.quadProduct(x);

          xTAx->reduceSum();
          
          scale = (*mxb)/(*xTAx);
          scale->limit(1.0, 1.0);
          
      }

      fineMatrix.correctSolution(fineLS.getCoarseIndex(),
                                 fineLS.getDelta(),
                                 scale,
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

      coarseLS->isSymmetric = fineLS.isSymmetric;

     int isContinue =  int( fineLS.getMatrix().getLocalSize() != coarseLS->getMatrix().getLocalSize() );
#ifdef FVM_PARALLEL     
     _commTarget.Allreduce(MPI::IN_PLACE, &isContinue, 1, MPI::INT, MPI::SUM);
#endif     
     if ( isContinue == 0 )
            break;
	  
#ifdef FVM_PARALLEL
     _coarseLinearSystems.push_back(coarseLS);
       
      int min_size = coarseLS->getMatrix().getMinSize( _commTarget );


      if ( verbosity > 1 && MPI::COMM_WORLD.Get_rank() == 0 )
        cout << " proc_id = " << MPI::COMM_WORLD.Get_rank() << "  Created coarse level " << n << " of size "
             << min_size  << endl;

      if ( min_size <= 3  )
        break;

      if ( coarseLS->getMatrix().getMergeSize( _commTarget ) < _mergeLevelSize &&  _mergeLevel == -1  && _isMerge ){
         _mergeLevel = n+1;
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

#endif 

#ifndef FVM_PARALLEL
    if ( coarseLS->getMatrix().getSize() <= 3 )
       break;
    if ( verbosity > 1 )
        cout << "Created coarse level " << n << " of size " << coarseLS->getMatrix().getSize() << endl;
    _coarseLinearSystems.push_back(coarseLS);
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
    cout << "0: " << *rNorm0 <<  endl;
#endif

#ifndef FVM_PARALLEL
   if ( verbosity > 0 )
      cout << "0: " << *rNorm0 << endl;
#endif

  if (*rNorm0 < absoluteTolerance )
    return rNorm0;

  for(int i=1; i<nMaxIterations; i++)
  {
      _totalIterations++;
      cycle(cycleType,0);
      finestMatrix.computeResidual(_finestLinearSystem->getDelta(),
                                   _finestLinearSystem->getB(),
                                   _finestLinearSystem->getResidual());
      MFRPtr rNorm(_finestLinearSystem->getResidual().getOneNorm());
      MFRPtr normRatio(rNorm->normalize(*rNorm0));

#ifndef FVM_PARALLEL
    if (verbosity >0  )
        cout << i << ": " <<  *rNorm << endl;


#endif

   
#ifdef FVM_PARALLEL
     if (*rNorm < absoluteTolerance || *normRatio < relativeTolerance || i == nMaxIterations-1)
        if (verbosity >0 && MPI::COMM_WORLD.Get_rank() == 0  )
        cout <<i << ": " <<   *rNorm << endl;
#endif

      if (*rNorm < absoluteTolerance || *normRatio < relativeTolerance)
         break;

  }

  _finestLinearSystem->getDelta().sync();
  
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

  {
      cycle(cycleType,0);
  }
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


//this will dump convergence history to file rather than screen
void 
AMG::redirectPrintToFile( const string& fname ) {
    m_filestr.open ( fname.c_str() );
    m_backup = cout.rdbuf();     // back up cout's streambuf
    m_psbuf = m_filestr.rdbuf();   // get file's streambuf
    cout.rdbuf(m_psbuf); 
}


void 
AMG::redirectPrintToScreen( ) {
   cout.rdbuf(m_backup);        // restore cout's original streambuf
   m_filestr.close();
}
