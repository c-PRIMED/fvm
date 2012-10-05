// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifdef  FVM_PARALLEL
#include <mpi.h>
#endif


#include "BCGStab.h"
BCGStab::BCGStab() :
  preconditioner(0),
  _totalIterations(0)
{}

BCGStab::~BCGStab()
{
}

void
BCGStab::cleanup()
{
  preconditioner->cleanup();
}

MFRPtr
BCGStab::solve(LinearSystem & ls)
{
  const MultiFieldMatrix& matrix = ls.getMatrix();

  // original system is in delta form
  shared_ptr<MultiField> x(ls.getDeltaPtr());
  shared_ptr<MultiField> bOrig(ls.getBPtr());

  matrix.computeResidual(ls.getDelta(),ls.getB(),ls.getResidual());

  MFRPtr rNorm0(ls.getResidual().getOneNorm());

  shared_ptr<MultiField> r(dynamic_pointer_cast<MultiField>(ls.getResidual().newCopy()));
  shared_ptr<MultiField> rTilda(dynamic_pointer_cast<MultiField>(ls.getResidual().newCopy()));

  shared_ptr<MultiField> p;
  shared_ptr<MultiField> pHat(dynamic_pointer_cast<MultiField>(x->newClone()));
  shared_ptr<MultiField> v(dynamic_pointer_cast<MultiField>(x->newClone()));
  shared_ptr<MultiField> t(dynamic_pointer_cast<MultiField>(x->newClone()));

  MFRPtr rho;
  MFRPtr rhoPrev;

  MFRPtr alpha;
  MFRPtr omega;

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
      _totalIterations++;
      rhoPrev = rho;
      rho = r->dotWith(*rTilda);

      rho->reduceSum();

      if (!p)
      {
          p = dynamic_pointer_cast<MultiField>(r->newCopy());
      }
      else
      {
          MFRPtr rhoRatio = (*rho) / (*rhoPrev);
          MFRPtr alphaByOmega = (*alpha) / (*omega);
          MFRPtr beta = (*rhoRatio) * (*alphaByOmega);
          p->msaxpy(*omega,*v);
          *p *= *beta;
          *p += *r;
      }

      pHat->zero();
      ls.replaceDelta(pHat);
      ls.replaceB(p);

      preconditioner->smooth(ls);

      matrix.multiply(*v,*pHat);

      MFRPtr rtv = rTilda->dotWith(*v);
      rtv->reduceSum();

      alpha = (*rho)/(*rtv);

      x->msaxpy(*alpha,*pHat);
      r->msaxpy(*alpha,*v);

      MFRPtr rNorm = r->getOneNorm();

      if (*rNorm < absoluteTolerance)
      {
          break;
      }

      shared_ptr<MultiField> sHat = pHat;
      sHat->zero();
      ls.replaceDelta(sHat);
      ls.replaceB(r);

      preconditioner->smooth(ls);

      matrix.multiply(*t,*sHat);

      MFRPtr tdotr = t->dotWith(*r);
      MFRPtr tdott = t->dotWith(*t);

      tdotr->reduceSum();
      tdott->reduceSum();

      omega = (*tdotr) / (*tdott);

      x->msaxpy(*omega,*sHat);
      r->msaxpy(*omega,*t);

      rNorm = r->getOneNorm();

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
BCGStab::smooth(LinearSystem& ls)
{
  throw CException("cannot use BCGStab as preconditioner");
}
