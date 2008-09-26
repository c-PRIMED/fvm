#include "Field.h"
#include "ArrayBase.h"
#include "StorageSite.h"
#include "Array.h"
#include "OneToOneIndexMap.h"

Field::Field(const string& name):
  IContainer(),
  _name(name),
  _arrays()
{
  logCtor();
}  

Field::~Field()
{
  logDtor();
}


bool
Field::hasArray(const StorageSite& s) const
{
  return ((_arrays.find(&s) != _arrays.end()) ||
          (s.getParent() && (_arrays.find(s.getParent()) != _arrays.end())));
}

const ArrayBase&
Field::operator[](const StorageSite& s) const
{
  ArrayMap::const_iterator pos = _arrays.find(&s);
  if (pos != _arrays.end())
    return *pos->second;
  return const_cast<Field&>(*this)._create(s);
}

ArrayBase&
Field::operator[](const StorageSite& s)
{
  ArrayMap::iterator pos = _arrays.find(&s);
  if (pos != _arrays.end())
    return *pos->second;
  
  return _create(s);
}

shared_ptr<ArrayBase>
Field::getArrayPtr(const StorageSite& s)
{
  if (_arrays.find(&s) == _arrays.end())
    _create(s);
  
  return _arrays[&s];
}

void
Field::addArray(const StorageSite& s, shared_ptr<ArrayBase> a)
{
  removeArray(s);
  _arrays[&s]=a;
}

void
Field::removeArray(const StorageSite& s)
{
  if (_childSitesMap.find(&s) != _childSitesMap.end())
  {
      vector<const StorageSite*>& childSites = *_childSitesMap[&s];
      foreach (const StorageSite* site,childSites)
      {
          removeArray(*site);
      }
      delete &childSites;
      _childSitesMap.erase(&s);
      
  }
  _arrays.erase(&s);
}

void
Field::removeArrays(const StorageSiteList& sites)
{
  foreach(const StorageSite* s, sites)
  {
      removeArray(*s);
  }
}

ArrayBase& 
Field::_create(const StorageSite& s)
{
  const StorageSite* parentSite = s.getParent();
  if (parentSite)
  {
      ArrayBase& parentArray = operator[](*parentSite);
      shared_ptr<ArrayBase> anew = parentArray.createOffsetArray(s.getOffset(),
                                                                 s.getCount());
      _arrays[&s]=anew; 

      if (_childSitesMap.find(parentSite) == _childSitesMap.end())
        _childSitesMap[parentSite] = new vector<const StorageSite*>();

      vector<const StorageSite*>& childSites = *_childSitesMap[parentSite];
      childSites.push_back(&s);
      return *anew;
  }
  else
  {
      ostringstream e;
      e << "Field::operator[] No array found and no creator defined for "
        << &s << endl;
      throw CException(e.str());
  }
}

shared_ptr<IContainer>
Field::newClone() const
{
  shared_ptr<Field> c(new Field(_name));
  
  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      c->addArray(*(pos.first),
                  dynamic_pointer_cast<ArrayBase>(pos.second->newClone()));
  }
  return c;
}


shared_ptr<IContainer>
Field::newCopy() const
{
  shared_ptr<Field> c(new Field(_name));
  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      c->addArray(*(pos.first),
                  dynamic_pointer_cast<ArrayBase>(pos.second->newCopy()));
  }
  return c;
}


void
Field::zero()
{
  foreach(ArrayMap::value_type& pos, _arrays)
  {
      pos.second->zero();
  }
}

Field&
Field::operator=(const Field& oField)
{
  foreach(ArrayMap::value_type& pos, _arrays)
  {
      *pos.second = oField[*pos.first];
  }
  return *this;
}

void
Field::copyFrom(const IContainer& oc)
{
  const Field& other = dynamic_cast<const Field&>(oc);
  operator=(other);
}

#if 0
void
Field::syncGather(const StorageSite& site)
{
  ArrayBase& a = *_arrays[&site];

  const StorageSite::MappersMap& mappers = site.getMappers();

  map<const StorageSite*, bool> isUptodate;

  for(StorageSite::MappersMap::const_iterator pos = mappers.begin();
      pos != mappers.end();
      ++pos)
  {
      const StorageSite& oSite = *pos->first;
      isUptodate[&oSite]=false;
  }
  
  bool allUpdated = false;

  do
  {
      allUpdated = true;
      for(StorageSite::MappersMap::const_iterator pos = mappers.begin();
          pos != mappers.end();
          ++pos)
      {
          const StorageSite& oSite = *pos->first;

          if (!isUptodate[&oSite])
          {
              if (_arrayStatus.find(&oSite) == _arrayStatus.end())
              {
                  // ostringstream e;
                  // e << "Field::operator[] No array found and no creator defined for "
                  //  << oSite.getInt("ID");
                  // throw CException(e.str());
                  continue;
              }

              const int otherArrayStatus = _arrayStatus[&oSite];
              ArrayBase *otherArray =0;
              
              if (otherArrayStatus == ARRAY_COMPUTED)
                otherArray = _arrays[&oSite];
              else if (otherArrayStatus == ARRAY_NOT_CREATED)
                otherArray = &_create(oSite);
              else
              {
                  allUpdated = false;
                  logInfo("waiting for array in Field::syncGather ");
                  if (otherArrayStatus == getThisThreadNum())
                  {
                      throw CException("Error in syncGather ");
                  }
              }

              if (otherArray)
              {
                  const Array<int>& fromIndices = pos->second->fromIndices;
                  const Array<int>& toIndices = pos->second->toIndices;
                  a.setSubsetFromSubset(*otherArray,fromIndices,toIndices);
                  isUptodate[&oSite] = true;
              }
          }
      }
  }
  while (!allUpdated);
}

#endif
