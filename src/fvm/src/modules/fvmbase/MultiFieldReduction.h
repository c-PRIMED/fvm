// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _MULTIFIELDREDUCTION_H_
#define _MULTIFIELDREDUCTION_H_

#include "Field.h"
#include "ArrayBase.h"

class MultiFieldReduction 
{
public:
  typedef map<const Field*,shared_ptr<ArrayBase> > ArrayMap;
  
  MultiFieldReduction();
  
  virtual ~MultiFieldReduction();
  DEFINE_TYPENAME("MultiFieldReduction");

  ArrayBase& operator[](const Field&);
  const ArrayBase& operator[](const Field&) const;
 
  MultiFieldReduction& operator+=(const MultiFieldReduction& ofield);

  void addArray(const Field& aIndex, shared_ptr<ArrayBase> a);
  bool hasArray(const Field& aIndex) const;

  shared_ptr<ArrayBase> getArrayPtr(const Field&);
  
  void reduceSum();

  bool operator<(const double tolerance) const;
  shared_ptr<MultiFieldReduction> operator/(const MultiFieldReduction& o);
  shared_ptr<MultiFieldReduction> normalize(const MultiFieldReduction& o);
  shared_ptr<MultiFieldReduction> operator*(const MultiFieldReduction& o);
  shared_ptr<MultiFieldReduction> operator-() const;

  void setMax(const MultiFieldReduction& o);
  void limit(const double min, const double max);

  void print(ostream &os) const;
  void sync();

private:
  ArrayMap _arrays;
};

inline ostream& operator<<(ostream &os,
                           const MultiFieldReduction &x)
{
  x.print(os);
  return os;
}

typedef shared_ptr<MultiFieldReduction> MFRPtr;

typedef vector<shared_ptr<MultiFieldReduction> > MFReductionVector;

#endif
