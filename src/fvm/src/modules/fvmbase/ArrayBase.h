#ifndef _ARRAYBASE_H_
#define _ARRAYBASE_H_

#include "NumType.h"
 
#include "IContainer.h"

class ArrayBase : public IContainer
{
public:

  ArrayBase() {}
  virtual ~ArrayBase() {}
  
  virtual shared_ptr<ArrayBase>
  createOffsetArray(const int offset, const int size) =0;

  virtual void
  correct(const IContainer& coarseI, const IContainer& coarseIndexI,
          const int length)=0;

  virtual void
  inject(IContainer& coarseI, const IContainer& coarseIndexI,
         const int length) const=0;

  virtual void
  setSubsetFromSubset(const ArrayBase& other, const ArrayBase& fromIndices,
                      const ArrayBase& toIndices) = 0;

  virtual void copyPartial(const IContainer& oc,
                           const int iBeg, const int iEnd)  =0;
  virtual void zeroPartial(const int iBeg, const int iEnd)  =0;

  virtual shared_ptr<ArrayBase> newSizedClone(const int size) const = 0;

  virtual void copyFrom(const IContainer& a) =0;
  virtual ArrayBase& operator+=(const ArrayBase& a) =0;
  virtual ArrayBase& operator-=(const ArrayBase& a) =0;
  virtual ArrayBase& operator/=(const ArrayBase& a) =0;
  virtual ArrayBase& operator*=(const ArrayBase& a) =0;

  virtual bool operator<(const double tolerance) const=0;

  virtual void print(ostream& os) const = 0;
  
  virtual ArrayBase& saxpy(const ArrayBase& alphabase,
                           const ArrayBase& xbase) = 0;
  virtual ArrayBase& msaxpy(const ArrayBase& alphabase,
                            const ArrayBase& xbase) = 0;
  
  virtual shared_ptr<ArrayBase>  getOneNorm() const=0;
  virtual shared_ptr<ArrayBase>  dotWith(const ArrayBase& a) const=0;
};

#endif
