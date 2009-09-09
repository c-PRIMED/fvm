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

  virtual void initAssembly() = 0;
  
  virtual void multiply(IContainer& yB, const IContainer& xB) const;
  virtual void multiplyAndAdd(IContainer& yB, const IContainer& xB) const;
  virtual void forwardGS(IContainer& xB, IContainer& bB,
                         IContainer& residual) const;
  virtual void reverseGS(IContainer& xB, IContainer& bB,
                         IContainer& residual) const;
  virtual void solveBoundary(IContainer& xB, IContainer& bB,
                             IContainer& residual) const;
  virtual void computeResidual(const IContainer& xB, const IContainer& bB,
                               IContainer& residual) const;

  virtual int createCoarsening(IContainer& coarseIndex, const int groupSize,
                               const double weighRatioThreshold);
  
  virtual const CRConnectivity& getConnectivity( ) const {throw;}

  
  virtual void *getDiagData() const {throw;}
  virtual void *getOffDiagData() const {throw;}
  virtual int   getDiagDataSize() const {throw;}
  virtual int   getOffDiagDataSize() const {throw;}

  
  virtual shared_ptr<CRConnectivity>
  createCoarseConnectivity(const IContainer& coarseIndex,
                           const CRConnectivity& coarseToFine,
                           const StorageSite& coarseRowSite,
                           const StorageSite& coarseColSite);
  virtual shared_ptr<Matrix>
  createCoarseMatrix(const IContainer& coarseIndex,
                     const CRConnectivity& coarseToFine,
                     const CRConnectivity& coarseConnectivity);

  virtual bool isInvertible() {return false;}
};


#endif
