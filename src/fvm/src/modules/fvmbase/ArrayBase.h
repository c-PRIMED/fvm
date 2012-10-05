// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

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
          const IContainer* scaleI,
          const int length)=0;

  virtual void
  inject(IContainer& coarseI, const IContainer& coarseIndexI,
         const int length) const=0;

  virtual void
  scatter(ArrayBase& other_, const ArrayBase& indices, const int offset=0) const=0;

  virtual void
  gather(const ArrayBase& other_, const ArrayBase& indices, const int offset=0) = 0;

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

  virtual shared_ptr<ArrayBase> operator-() const {throw;}
  virtual void limit(const double min, const double max) {throw;}
  
    virtual ArrayBase& safeDivide(const ArrayBase& a) =0;
  virtual ArrayBase& normalize(const ArrayBase& a) =0;

  virtual bool operator<(const double tolerance) const=0;
  // set self's value to be a's value if the latter is higher. used
  // for residuals
  virtual void setMax(const ArrayBase& a) =0;

  virtual void print(ostream& os) const = 0;
  
  virtual ArrayBase& saxpy(const ArrayBase& alphabase,
                           const ArrayBase& xbase) = 0;
  virtual ArrayBase& msaxpy(const ArrayBase& alphabase,
                            const ArrayBase& xbase) = 0;
  
  virtual shared_ptr<ArrayBase>  getOneNorm(const int lengthToUse) const=0;
  virtual shared_ptr<ArrayBase>  dotWith(const ArrayBase& a,
                                         const int lengthToUse) const=0;
  virtual shared_ptr<ArrayBase>  reduceSum() const=0;
  virtual void setSum(const ArrayBase& sumBase) = 0;

  virtual int getDimension() const =0;
  virtual void getShape(int *shape) const = 0;
  virtual void *getData() const = 0;
  virtual int getDataSize() const  =0;
  virtual int getLength()   const  =0;
  virtual PrimType getPrimType() const = 0;
  
};

#endif
