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
Matrix::multiply(IContainer& yB, const IContainer& xB)
{
  throw CException("multiply not implemented");
}

void
Matrix::multiplyAndAdd(IContainer& yB, const IContainer& xB)
{
  throw CException("multiplyAndAdd not implemented");
}

void Matrix::forwardGS(IContainer& xB, IContainer& bB,
                       IContainer& residual)
{
  throw CException("forwardGS not implemented");
}

void Matrix::reverseGS(IContainer& xB, IContainer& bB,
                       IContainer& residual)
{
  throw CException("reverseGS not implemented");
}

void Matrix::solveBoundary(IContainer& xB, IContainer& bB,
                           IContainer& residual)
{
  throw CException("solveBoundary not implemented");
}

void Matrix::computeResidual(const IContainer& xB, const IContainer& bB,
                             IContainer& residual)
{
  throw CException("computeResidual not implemented");
}


