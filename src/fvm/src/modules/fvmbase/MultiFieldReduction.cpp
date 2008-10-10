#include "MultiFieldReduction.h"
#include "Array.h"
#include "OneToOneIndexMap.h"

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

MFPtr
MultiFieldReduction::operator/(const MultiFieldReduction& o)
{
  MFPtr r(new MultiFieldReduction());
  
  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      shared_ptr<ArrayBase> aptr = dynamic_pointer_cast<ArrayBase>(pos.second->newCopy());
      aptr->safeDivide(o[*pos.first]);
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
MultiFieldReduction::reduceSum()
{
  ArrayBase& array0 = *(_arrays.begin()->second);

  for(ArrayMap::const_iterator posj=_arrays.begin();
      posj!=_arrays.end();
      ++posj)
  {
      ArrayBase& arrayj = *posj->second;
      if (&array0 != &arrayj)
        array0 += arrayj;
  }
  
  for(ArrayMap::const_iterator posj=_arrays.begin();
      posj!=_arrays.end();
      ++posj)
  {
      ArrayBase& arrayj = *posj->second;
      if (&array0 != &arrayj)
        arrayj.copyFrom(array0);
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
