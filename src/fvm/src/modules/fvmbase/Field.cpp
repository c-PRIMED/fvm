#include "Field.h"
#include "ArrayBase.h"
#include "ArrayCreator.h"
#include "StorageSite.h"
#include "Array.h"
#include "CDict.h"
#include "OneToOneIndexMap.h"
#include "myomp.h"

Field::Field():
  IContainer(),
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
          !s.isVarNone("parent"));
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


void
Field::addArray(const StorageSite& s, ArrayBase* a)
{
  removeArray(s);
  _arrays.add(&s,a);
  //syncScatter(s);
  //syncGather(s);
}

void
Field::removeArray(const StorageSite& s)
{
  if (_childSitesMap.find(&s) != _childSitesMap.end())
  {
      vector<const StorageSite*>& childSites = *_childSitesMap[&s];
      for (vector<const StorageSite*>::iterator pos = childSites.begin();
           pos != childSites.end();
           ++pos)
        {
            removeArray(**pos);
        }
      delete &childSites;
      _childSitesMap.erase(&s);
      
  }

  _arrays.remove(&s);
}

void
Field::removeArrays(const StorageSiteCList& sites)
{
  for (StorageSiteCList::const_iterator pos = sites.begin();
       pos != sites.end();
       ++pos)
  {
      removeArray(**pos);
  }
}

ArrayBase& 
Field::_create(const StorageSite& s)
{
  const int sID = s.getInt("ID");

  if (!s.isVarNone("parent"))
  {
      const StorageSite& parentSite = s.getChildRef<StorageSite>("parent");
      ArrayBase& parent = operator[](parentSite);
      _arrayStatus[&s] = getThisThreadNum();
      ArrayBase* anew = parent.createOffsetArray(s.getInt("offset"),
                                              s.getInt("count"));
      _arrays.add(&s,anew); 
      _arrayStatus[&s] = ARRAY_COMPUTED;

      if (_childSitesMap.find(&parentSite) == _childSitesMap.end())
        _childSitesMap[&parentSite] = new vector<const StorageSite*>();

      vector<const StorageSite*>& childSites = *_childSitesMap[&parentSite];
      childSites.push_back(&s);
      return *anew;
  }
  else
  {
      ostringstream e;
      e << "Field::operator[] No array found and no creator defined for "
        << sID;
      throw CException(e.str());
  }
}

IContainer*
Field::newClone() const
{
  Field *c = new Field();
  
  for(ArrayMap::const_iterator pos = _arrays.begin();
      pos != _arrays.end();
      ++pos)
  {
      c->addArray(*(pos->first),SafeCast<ArrayBase>(pos->second->newClone()));
  }
  return c;
}


IContainer*
Field::newCopy() const
{
  Field *c = new Field();
  
  for(ArrayMap::const_iterator pos = _arrays.begin();
      pos != _arrays.end();
      ++pos)
  {
      c->addArray(*(pos->first),SafeCast<ArrayBase>(pos->second->newCopy()));
  }
  return c;
}


void
Field::zero()
{
  for(ArrayMap::iterator pos = _arrays.begin();
      pos != _arrays.end();
      ++pos)
  {
      pos->second->zero();
  }
}

void
Field::copyFrom(const IContainer& oc)
{
  const Field& ofield = SafeCast<Field>(oc);
  for(ArrayMap::iterator pos = _arrays.begin();
      pos != _arrays.end();
      ++pos)
  {
      pos->second->copyFrom(ofield[*(pos->first)]);
  }
}


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

