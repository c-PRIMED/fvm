// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _STORAGESITE_H_
#define _STORAGESITE_H_

#include <map>
#include <vector>
#include "misc.h"
#include "RLogInterface.h"
#include "Array.h"
#include <cassert>

class OneToOneIndexMap;
class Mesh;

class StorageSite
{
public:


  typedef map<  const StorageSite*, shared_ptr< Array<int> >  > ScatterMap;
  typedef map<  const StorageSite*, shared_ptr< Array<int> >  > GatherMap;
  typedef map<  const StorageSite*, shared_ptr< Array<int> >  > CommonMap;
  typedef map<  const StorageSite*, const Mesh* >        MeshMap;
  //finally cellcell2 keep information in global numbering
  //field or multifield 
  typedef map<  const StorageSite*, map<int,int> >   ScatterIndex;  



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
    _countLevel1 = _count;
  }


  
  void setScatterProcID( int proc_id ) { _scatterProcID = proc_id; }
  void setGatherProcID(  int proc_id ) { _gatherProcID  = proc_id; }
  void setTag( int tag) { _tag = tag;}
   
  void setMesh( const Mesh* mesh) { _meshMap[this] = mesh; }
  const Mesh& getMesh() const { return *_meshMap[this]; }

  const ScatterMap&   getScatterMap   () const { return _scatterMap;}
  const  GatherMap&   getGatherMap    () const { return  _gatherMap;}
  const CommonMap &   getCommonMap    () const { return _commonMap; }
  const ScatterIndex&  getScatterIndex () const { return _scatterIndex; }

  ScatterMap &  getScatterMap   ()  { return _scatterMap;}
  GatherMap  &  getGatherMap    ()  { return _gatherMap; }
  CommonMap  &  getCommonMap    ()  { return _commonMap; }
  ScatterIndex&  getScatterIndex ()  { return _scatterIndex; }

  //methods for Level1 layer creation
  void setCountLevel1( const int countLevel1 ) {
    _countLevel1 = countLevel1;
  }
  int getCountLevel1() const {return _countLevel1;}
  const ScatterMap&  getScatterMapLevel1() const { return _scatterMapLevel1;}
  const  GatherMap&   getGatherMapLevel1() const { return _gatherMapLevel1 ;}

  ScatterMap&  getScatterMapLevel1()  { return _scatterMapLevel1;}
  GatherMap &  getGatherMapLevel1 ()  { return _gatherMapLevel1 ;}


  void clearGatherScatterMaps();
  
  int getScatterProcID() const { return _scatterProcID;}
  int getGatherProcID()  const { return _gatherProcID; }
  int getTag()           const { return _tag; }

  const StorageSite* const getParent() const {return _parent;}
  int getOffset() const {return _offset;}


private:
  StorageSite(const StorageSite&);
  int _count;
  int _selfCount;
  const int _offset;
  const StorageSite* const _parent;

  ScatterMap         _scatterMap;
  GatherMap          _gatherMap;
  CommonMap          _commonMap;
  ScatterIndex       _scatterIndex;
  
  mutable MeshMap            _meshMap;

  int _countLevel1;
  ScatterMap         _scatterMapLevel1;
  GatherMap          _gatherMapLevel1;
  
  int   _scatterProcID;
  int   _gatherProcID;
  int   _tag;

};

typedef vector<const StorageSite*> StorageSiteList;

#endif
