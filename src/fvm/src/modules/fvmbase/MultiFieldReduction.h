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

  void reduceSum();

  bool operator<(const double tolerance) const;
  shared_ptr<MultiFieldReduction> operator/(const MultiFieldReduction& o);
  void setMax(const MultiFieldReduction& o);

  void print(ostream &os) const;

private:
  ArrayMap _arrays;
};

inline ostream& operator<<(ostream &os,
                           const MultiFieldReduction &x)
{
  x.print(os);
  return os;
}

typedef shared_ptr<MultiFieldReduction> MFPtr;

typedef vector<shared_ptr<MultiFieldReduction> > MFReductionVector;

#endif
