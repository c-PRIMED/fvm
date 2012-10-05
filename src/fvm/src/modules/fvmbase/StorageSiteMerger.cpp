// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifdef FVM_PARALLEL

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "StorageSiteMerger.h"

using namespace std;

StorageSiteMerger::StorageSiteMerger(int target_proc_id, const set<int>& group, const StorageSite& cell_site )
:_groupID(target_proc_id), _group(group), _cellSite(cell_site), _mergeSiteSize(0), _mergeSiteGhostSize(0)
{
    init();
}

StorageSiteMerger::~StorageSiteMerger()
{
  
}


void
StorageSiteMerger::init()
{
   //subcommunicator
   int color = _groupID;
   int key   = MPI::COMM_WORLD.Get_rank();
   _comm = MPI::COMM_WORLD.Split( color, key );
}


shared_ptr<StorageSite>
StorageSiteMerger::merge()
{
    StorageSite::GatherMap::const_iterator   it;
    int tot_adjacent_count = 0; //this will be substracted before sending groupID

    for ( it = _cellSite.getGatherMap().begin(); it != _cellSite.getGatherMap().end(); it++){
        const StorageSite* site = it->first;
        int  gatherID = site->getGatherProcID();
        //this check if neighbour mesh id is in group list, if it is, we want to increse tot_adjacent_count
        if ( _group.count( gatherID) > 0 ){
            int interface_count =it->second->getLength();
            tot_adjacent_count += interface_count;
        }

    }

   int local_size = _cellSite.getCount()  - tot_adjacent_count; 
   int ghost_size = local_size - _cellSite.getSelfCount();
   int self_size  = _cellSite.getSelfCount();

   _comm.Reduce( &self_size , &_mergeSiteSize     , 1, MPI::INT, MPI::SUM, _groupID );
   _comm.Reduce( &ghost_size, &_mergeSiteGhostSize, 1, MPI::INT, MPI::SUM, _groupID );
   //debug_print();
   return shared_ptr<StorageSite> ( new StorageSite(_mergeSiteSize) );

} 


void
StorageSiteMerger::debug_print()
{
  shared_ptr<StorageSite> tempSite = merge();
  stringstream ss;
  ss << "proc" << MPI::COMM_WORLD.Get_rank() << "_storage_site_merger.dat";
  ofstream  debug_file( (ss.str()).c_str() );

  debug_file << " selfCount   = " << _mergeSiteSize << endl;
  debug_file << " GhostCount  = " << _mergeSiteGhostSize << endl;
  debug_file << " count       = " << _mergeSiteSize + _mergeSiteGhostSize << endl;

  debug_file.close();

}

#endif
