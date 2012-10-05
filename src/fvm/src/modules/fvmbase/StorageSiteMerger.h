// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _STORAGESITEMERGER_H_
#define _STORAGESITEMERGER_H_
#ifdef FVM_PARALLEL

#include <mpi.h>
#include "StorageSite.h"
#include <cassert>
#include <set>

class StorageSiteMerger
{
public:
   StorageSiteMerger( int target_proc_id, const set<int>& group, const StorageSite& cell_site );
   ~StorageSiteMerger();

   shared_ptr<StorageSite>  merge();
   void  debug_print();

   int  getSelfCount()  const { return _mergeSiteSize; }
   int  getGhostCount() const { return _mergeSiteGhostSize; }
   int  getCount()      const { return _mergeSiteSize + _mergeSiteGhostSize; }

private:
  StorageSiteMerger(const StorageSiteMerger&);
  void  init();


  int _count;
  int _selfCount;
  int _offset;

  int _groupID;              
  const set<int>&  _group;
  const StorageSite& _cellSite;
  int _mergeSiteSize;
  int _mergeSiteGhostSize;

  MPI::Intracomm _comm;

  //const StorageSite*   const _parent;
  StorageSite::ScatterMap   _scatterMap;
  StorageSite::GatherMap    _gatherMap;

  int   _scatterProcID;
  int   _gatherProcID;


};


#endif
#endif
