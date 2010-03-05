#ifndef _STORAGESITE_H_
#define _STORAGESITE_H_

#include <map>
#include <vector>
#include "misc.h"
#include "RLogInterface.h"
#include "Array.h"
#include <cassert>

class OneToOneIndexMap;

class StorageSite
{
public:


  typedef map<const StorageSite*, shared_ptr<OneToOneIndexMap> > MappersMap;
  typedef map<  const StorageSite*, shared_ptr< Array<int> >  > ScatterMap;
  typedef map<  const StorageSite*, shared_ptr< Array<int> >  > GatherMap;

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
  
  void setScatterProcID( int proc_id ) { _scatterProcID = proc_id; }
  void setGatherProcID(  int proc_id ) { _gatherProcID  = proc_id; }
  void setTag( int tag) { _tag = tag;}

  const MappersMap& getMappers() const {throw;}
  MappersMap& getMappers() {throw;}

  const ScatterMap&  getScatterMap() const  { return _scatterMap;}
  const  GatherMap&   getGatherMap() const  { return  _gatherMap;}
  ScatterMap&  getScatterMap()  { return _scatterMap;}
  GatherMap &  getGatherMap ()  { return _gatherMap; }


  int getScatterProcID() const { return _scatterProcID;}
  int getGatherProcID()  const { return _gatherProcID; }
  int getTag()           const { return _tag; }

  const StorageSite* const getParent() const {return _parent;}
  int getOffset() const {return _offset;}

  //void  scatterGatherMaps( );

private:
  StorageSite(const StorageSite&);
  int _count;
  int _selfCount;
  const int _offset;
  const StorageSite* const _parent;
  //MappersMap _mappers;
  ScatterMap         _scatterMap;
  GatherMap          _gatherMap;

  int   _scatterProcID;
  int   _gatherProcID;
  int   _tag;

};

typedef vector<const StorageSite*> StorageSiteList;

#endif
