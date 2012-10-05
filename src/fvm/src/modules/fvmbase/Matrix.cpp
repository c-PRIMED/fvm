// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "Matrix.h"
#include "IContainer.h"

Matrix::Matrix() 
{}

Matrix::~Matrix()
{}


int
Matrix::createCoarsening(IContainer& coarseIndex, const int groupSize,
                         const double weighRatioThreshold)
{
  throw CException("createCoarsening not implemented");
}

shared_ptr<CRConnectivity>
Matrix::createCoarseConnectivity(const IContainer& coarseIndex,
                                 const CRConnectivity& coarseToFine,
                                 const StorageSite& coarseRowSite,
                                 const StorageSite& coarseColSite)
{
  throw CException("createCoarseConnectivity not implemented");
}

shared_ptr<Matrix>
Matrix::createCoarseMatrix(const IContainer& coarseIndex,
                           const CRConnectivity& coarseToFine,
                           const CRConnectivity& coarseConnectivity)
{
  throw CException("createCoarseMatrix not implemented");
}

void
Matrix::multiply(IContainer& yB, const IContainer& xB) const
{
  throw CException("multiply not implemented");
}

void
Matrix::multiplyAndAdd(IContainer& yB, const IContainer& xB) const
{
  throw CException("multiplyAndAdd not implemented");
}

void Matrix::forwardGS(IContainer& xB, IContainer& bB,
                       IContainer& residual) const
{
  throw CException("forwardGS not implemented");
}

void Matrix::reverseGS(IContainer& xB, IContainer& bB,
                       IContainer& residual) const
{
  throw CException("reverseGS not implemented");
}

void Matrix::Jacobi(IContainer&, const IContainer&,
                       const IContainer&) const
{
  throw CException("Jacobi not implemented");
}

void Matrix::iluSolve(IContainer&, const IContainer&,
                       const IContainer&) const
{
  throw CException("iluSolve not implemented");
}
void Matrix::spikeSolve(IContainer& xB, const IContainer& bB,
                const IContainer& residual, const SpikeStorage& spike_storage) const
{
   throw CException("spikeSolve not imlemented");
}
  

void Matrix::solveBoundary(IContainer& xB, IContainer& bB,
                           IContainer& residual) const
{
  throw CException("solveBoundary not implemented");
}

void Matrix::computeResidual(const IContainer& xB, const IContainer& bB,
                             IContainer& residual) const
{
  throw CException("computeResidual not implemented");
}
/*
void Matrix::setDirichlet(const int nr)
{
  throw CException("setDirichlet not implemented");
}
*/

