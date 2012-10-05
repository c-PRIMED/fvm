// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifdef  FVM_PARALLEL
#include <mpi.h>
#endif


#include "CG.h"
CG::CG() :
  preconditioner(0)
{}

CG::~CG()
{
}

void
CG::cleanup()
{
  preconditioner->cleanup();
}

MFRPtr
CG::solve(LinearSystem & ls)
{
  const MultiFieldMatrix& matrix = ls.getMatrix();

  // original system is in delta form
  shared_ptr<MultiField> x(ls.getDeltaPtr());
  shared_ptr<MultiField> bOrig(ls.getBPtr());

  matrix.computeResidual(ls.getDelta(),ls.getB(),ls.getResidual());

  MFRPtr rNorm0(ls.getResidual().getOneNorm());

  shared_ptr<MultiField> r(dynamic_pointer_cast<MultiField>(ls.getResidual().newCopy()));

  shared_ptr<MultiField> p;
  shared_ptr<MultiField> z(dynamic_pointer_cast<MultiField>(x->newClone()));
  shared_ptr<MultiField> q(dynamic_pointer_cast<MultiField>(x->newClone()));

  MFRPtr rho;
  MFRPtr rhoPrev;

  MFRPtr alpha;
  MFRPtr beta;

#ifndef  FVM_PARALLEL
  if (verbosity >0)
    cout << 0 << ": " << *rNorm0 << endl;
#endif

#ifdef  FVM_PARALLEL
  if (verbosity >0 && MPI::COMM_WORLD.Get_rank() == 0)
    cout << 0 << ": " << *rNorm0 << endl;
#endif


  for(int i = 0; i<nMaxIterations; i++)
  {
      z->zero();
      ls.replaceDelta(z);
      ls.replaceB(r);

      preconditioner->smooth(ls);


      rhoPrev = rho;
      rho = r->dotWith(*z);

      rho->reduceSum();

      if (!p)
      {
          p = dynamic_pointer_cast<MultiField>(z->newCopy());
      }
      else
      {
          MFRPtr beta = (*rho) / (*rhoPrev);
          *p *= *beta;
          *p += *z;
      }

      matrix.multiply(*q,*p);

      MFRPtr ptq = p->dotWith(*q);
      ptq->reduceSum();

      alpha = (*rho)/(*ptq);

      x->msaxpy(*alpha,*p);
      r->msaxpy(*alpha,*q);

      MFRPtr rNorm = r->getOneNorm();


      MFRPtr normRatio(rNorm->normalize(*rNorm0));

#ifndef FVM_PARALLEL
      if (verbosity >0)
        cout << i+1 << ": " << *rNorm << endl;
#endif

#ifdef  FVM_PARALLEL
  if (verbosity >0 && MPI::COMM_WORLD.Get_rank() == 0)
    cout << i+1 << ": " << *rNorm << endl;
#endif

      if (*rNorm < absoluteTolerance || *normRatio < relativeTolerance)
        break;


  }
  ls.replaceDelta(x);
  ls.replaceB(bOrig);
#ifdef FVM_PARALLEL
  x->sync();
#endif

  matrix.computeResidual(ls.getDelta(),ls.getB(),ls.getResidual());
  MFRPtr rNormn(ls.getResidual().getOneNorm());

#ifndef FVM_PARALLEL
  if (verbosity >0)
    cout << "n" << ": " << *rNormn << endl;
#endif


#ifdef  FVM_PARALLEL
  if (verbosity >0 && MPI::COMM_WORLD.Get_rank() == 0)
    cout << "n" << ": " << *rNormn << endl;
#endif



  return rNorm0;
}

void
CG::smooth(LinearSystem& ls)
{
  throw CException("cannot use CG as preconditioner");
}
