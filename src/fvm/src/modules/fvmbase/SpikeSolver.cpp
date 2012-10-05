// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include "SpikeSolver.h"
#include "LinearSystemMerger.h"
#include "CRConnectivity.h"
#include "SpikeStorage.h"
#include <set>


SpikeSolver::SpikeSolver(const SpikeStorage& spike_storage):
_spikeStorage(spike_storage)
{
  logCtor();
  
}

SpikeSolver::~SpikeSolver()
{
  logDtor();
}

void
SpikeSolver::doSweeps(LinearSystem& ls, const int nSweeps)
{
  const MultiFieldMatrix& m = ls.getMatrix();
  MultiField& delta = ls.getDelta();
  const MultiField& b = ls.getB();
  MultiField& r = ls.getResidual();

  for(int i=0; i<nSweeps; i++)
  {
      m.spikeSolve(delta,b,r,_spikeStorage);
  }
}


void
SpikeSolver::cleanup()
{
}

MFRPtr
SpikeSolver::solve(LinearSystem & ls)
{
  ls.getMatrix().computeResidual(ls.getDelta(),
                                 ls.getB(),
                                 ls.getResidual());
  MFRPtr rNorm0(ls.getResidual().getOneNorm());

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
      doSweeps(ls,1);

      ls.getMatrix().computeResidual(ls.getDelta(),
                                     ls.getB(),
                                     ls.getResidual());
      MFRPtr rNorm(ls.getResidual().getOneNorm());
      MFRPtr normRatio((*rNorm)/(*rNorm0));

#ifndef FVM_PARALLEL
      if (verbosity >0  )
        cout << i << ": " <<  *rNorm << endl;
#endif


#ifdef FVM_PARALLEL
     if (*rNorm < absoluteTolerance || *normRatio < relativeTolerance || i == nMaxIterations-1)
        if (verbosity >0 && MPI::COMM_WORLD.Get_rank() == 0  )
        cout <<i << ": " << "procID = " << MPI::COMM_WORLD.Get_rank() <<  *rNorm << endl;
#endif

      if (*rNorm < absoluteTolerance || *normRatio < relativeTolerance)
         break;
      
  }
  return rNorm0;
}


void
SpikeSolver::smooth(LinearSystem & ls)
{
  doSweeps(ls,1);
}
