#ifndef _STORAGESITE_H_
#define _STORAGESITE_H_

#include <map>
#include <vector>
#include "misc.h"
#include "RLogInterface.h"

class OneToOneIndexMap;

class StorageSite
{
public:


  typedef map<const StorageSite*, shared_ptr<OneToOneIndexMap> > MappersMap;
  
  StorageSite(const int selfCount, const int nGhost=0,
              const int offset=0, const StorageSite* parent=0);
  virtual ~StorageSite();

  DEFINE_TYPENAME("StorageSite");
  
  int getCount() const {return _count;}
  int getSelfCount() const {return _selfCount;}

  void setCount(const int selfCount, const int nGhost=0)
  {
    _count = selfCount+nGhost;
    _selfCount = selfCount;
  }

  const MappersMap& getMappers() const {return _mappers;}
  MappersMap& getMappers() {return _mappers;}

  const StorageSite* const getParent() const {return _parent;}
  int getOffset() const {return _offset;}
  
private:
  StorageSite(const StorageSite&);
  int _count;
  int _selfCount;
  const int _offset;
  const StorageSite* const _parent;
  MappersMap _mappers;
};

typedef vector<const StorageSite*> StorageSiteList;

#endif
