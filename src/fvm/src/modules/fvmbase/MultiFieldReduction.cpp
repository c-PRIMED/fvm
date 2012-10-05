// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include "MultiFieldReduction.h"
#include "Array.h"
#include "OneToOneIndexMap.h"

#include <iostream>

MultiFieldReduction::MultiFieldReduction():
  _arrays()
{
  logCtor();
}

MultiFieldReduction::~MultiFieldReduction()
{
  logDtor();
}

ArrayBase&
MultiFieldReduction::operator[](const Field& i)
{
  ArrayMap::const_iterator pos = _arrays.find(&i);
  if (pos != _arrays.end())
  {
      return *pos->second;
  }

  ostringstream e;
  e << "MultiFieldReduction::operator[] No array found";
  throw CException(e.str());
}

shared_ptr<ArrayBase>
MultiFieldReduction::getArrayPtr(const Field& i)
{
  ArrayMap::const_iterator pos = _arrays.find(&i);
  if (pos != _arrays.end())
  {
      return pos->second;
  }

  ostringstream e;
  e << "MultiFieldReduction::operator[] No array found";
  throw CException(e.str());
}

const ArrayBase&
MultiFieldReduction::operator[](const Field& i) const
{
  ArrayMap::const_iterator pos = _arrays.find(&i);
  if (pos != _arrays.end())
  {
      return *pos->second;
  }

  ostringstream e;
  e << "MultiFieldReduction::operator[] No array found";
  throw CException(e.str());
}

MultiFieldReduction&
MultiFieldReduction::operator+=(const MultiFieldReduction& ofield)
{
  for(ArrayMap::const_iterator pos=_arrays.begin();
      pos!=_arrays.end();
      ++pos)
  {
      ArrayBase& myArray = operator[](*pos->first);
      const ArrayBase& oArray = ofield[*pos->first];
      myArray += oArray;
  }
  return *this;
}

bool
MultiFieldReduction::operator<(const double tolerance) const
{
  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      if (! (*pos.second < tolerance)) return false;
  }
  return true;
}

MFRPtr
MultiFieldReduction::operator*(const MultiFieldReduction& o)
{
  MFRPtr r(new MultiFieldReduction());
  
  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      shared_ptr<ArrayBase> aptr = dynamic_pointer_cast<ArrayBase>(pos.second->newCopy());
      *aptr *= o[*pos.first];
      r->addArray(*pos.first,aptr);
  }
  return r;
}

MFRPtr
MultiFieldReduction::operator/(const MultiFieldReduction& o)
{
  MFRPtr r(new MultiFieldReduction());
  
  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      shared_ptr<ArrayBase> aptr = dynamic_pointer_cast<ArrayBase>(pos.second->newCopy());
      aptr->safeDivide(o[*pos.first]);
      r->addArray(*pos.first,aptr);
  }
  return r;
}

MFRPtr
MultiFieldReduction::operator-() const
{
  MFRPtr r(new MultiFieldReduction());
  
  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      shared_ptr<ArrayBase> aptr = -(*pos.second);
      r->addArray(*pos.first,aptr);
  }
  return r;
}

MFRPtr
MultiFieldReduction::normalize(const MultiFieldReduction& o)
{
  MFRPtr r(new MultiFieldReduction());
  
  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      shared_ptr<ArrayBase> aptr = dynamic_pointer_cast<ArrayBase>(pos.second->newCopy());
      aptr->normalize(o[*pos.first]);
      r->addArray(*pos.first,aptr);
  }
  return r;
}

void
MultiFieldReduction::setMax(const MultiFieldReduction& o)
{
  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      pos.second->setMax(o[*pos.first]);
  }
}

void
MultiFieldReduction::limit(const double min, const double max)
{
  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      pos.second->limit(min,max);
  }
}

void
MultiFieldReduction::reduceSum()
{
  shared_ptr<ArrayBase> sum;

  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      shared_ptr<ArrayBase> thisSum = pos.second->reduceSum();
      if (sum)
        *sum += *thisSum;
      else
        sum = thisSum;
  }

  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      pos.second->setSum(*sum);
  }


}

bool
MultiFieldReduction::hasArray(const Field& i) const
{
  return  _arrays.find(&i) != _arrays.end();
}

void
MultiFieldReduction::addArray(const Field& i, shared_ptr<ArrayBase> a)
{
  _arrays[&i]=a;
}

void
MultiFieldReduction::print(ostream &os) const
{
  os << "[";

  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      os << pos.first->getName() << " : ";
      pos.second->print(os);
  }

  os << "]";
}

void
MultiFieldReduction::sync()
{
#ifdef FVM_PARALLEL
   ArrayMap::const_iterator it;
   for ( it = _arrays.begin(); it != _arrays.end(); it++ ){
         ArrayBase&  myArray = *(it->second);
         int count = myArray.getDataSize() / sizeof(double); 
         MPI::COMM_WORLD.Allreduce( MPI::IN_PLACE, myArray.getData(), count, MPI::DOUBLE, MPI::SUM);
   }
#endif

}
