#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "misc.h"
class IContainer;
class StorageSite;
class CRConnectivity;


class Matrix
{
public:
  Matrix();
  virtual ~Matrix();

  DEFINE_TYPENAME("Matrix");

  virtual void multiply(IContainer& yB, const IContainer& xB);
  virtual void multiplyAndAdd(IContainer& yB, const IContainer& xB);
  virtual void forwardGS(IContainer& xB, IContainer& bB,
                         IContainer& residual);
  virtual void reverseGS(IContainer& xB, IContainer& bB,
                         IContainer& residual);
  virtual void solveBoundary(IContainer& xB, IContainer& bB,
                             IContainer& residual);
  virtual void computeResidual(const IContainer& xB, const IContainer& bB,
                               IContainer& residual);

  virtual int createCoarsening(IContainer& coarseIndex, const int groupSize,
                               const double weighRatioThreshold);
  
  virtual shared_ptr<CRConnectivity>
  createCoarseConnectivity(const IContainer& coarseIndex,
                           const CRConnectivity& coarseToFine,
                           const StorageSite& coarseRowSite,
                           const StorageSite& coarseColSite);
  virtual shared_ptr<Matrix>
  createCoarseMatrix(const IContainer& coarseIndex,
                     const CRConnectivity& coarseToFine,
                     const CRConnectivity& coarseConnectivity);
};


#endif
