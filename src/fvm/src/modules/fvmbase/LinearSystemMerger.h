#ifndef _LINEARSYSTEMMERGER_H_
#define _LINEARSYSTEMMERGER_H_

#include "LinearSystem.h"
#include <cassert>
#include <set>
#include <map>
#include <mpi.h>

// variable name convection
//XXXLocal = this variable is unique to that process and only make sense in its process
//XXXGlobl = this variable is meaningful value for all processes (all process has the same value)
//XXXX     = this variable only make sense in target process 

class StorageSiteMerger;

class LinearSystemMerger
{
public:

    typedef   shared_ptr< Array<int> >     ArrayIntPtr;
    typedef   map<int,ArrayIntPtr>         ArrayIntPtrMap;
    friend class LinearSystem;

   LinearSystemMerger( int target_proc_id, const set<int>& group, LinearSystem& ls );
   ~LinearSystemMerger();


   void  debug_print();



private:
  LinearSystemMerger(const LinearSystemMerger&);
  void  init();
  void  merge();
  void  get_neigh_mesh_counts();
  void  get_scatter_cells();
  void  get_gather_cells();
  void  get_crconnectivity();
  void  get_local_to_global_map();
  void  set_merged_crconnectivity();
  void  update_gatherCells_from_scatterCells();
  void  get_ls_vectors();

  int _targetID;
  int _groupID;
  const set<int>&     _group;
  LinearSystem& _ls;
 
  shared_ptr<LinearSystem> _mergedLS;

  int _procID;
  int _totalProcs;
  int _totalInterfaces;
  int _totalScatterCells;
  int _totalScatterCellsLocal;
  int _totalGatherCells;
  int _totalGatherCellsLocal;
  int _totalCells;    //inner cells

  //mpi buffers
  ArrayIntPtr  _neighMeshCounts;                  //size = _totalProcs, access : [procid]
  map<int, ArrayIntPtr>  _scatterInterfaceCounts; //size = sum[_neightMeshCouns[procid]], access:[procid][ interface_id ]
  ArrayIntPtr  _scatterSize;                      //size = _totalProcs 
  vector< map<int, ArrayIntPtr> > _scatterCells;  //(proc_id,interface_id) = cells
  map<int, ArrayIntPtr>  _scatterInterfaceIDs;    //size = totalInterfaces, access to ids_array [procid]

  map<int, ArrayIntPtr>  _gatherInterfaceCounts;
  ArrayIntPtr  _gatherSize;                       //size = _totalProcs 
  vector< map<int, ArrayIntPtr> >  _gatherCells;  //(proc_id,interface_id) = cells
  map<int,ArrayIntPtr>   _gatherInterfaceIDs;     //size = totalInterfaces
  vector< map<int,int> > _gatherIDsLocalToGlobalMap;  //[proc][localid] = globalID
  ArrayIntPtr  _selfCounts;                       //size = totalProcs;

  ArrayIntPtr _rowLength; //CRConnecivity row_dimension,  size = _totalProcs
  ArrayIntPtr _colLength; //CRConnectivity col_dimension, size = _totalProcs
  map< int, ArrayIntPtr >  _row;   //CRConnectivity _row, size = totalCells
  map< int, ArrayIntPtr >  _col;   //CRConnectivity _col, size = ??

  map< int, ArrayIntPtr >  _localToGlobal;  //[procID][localID] = globalID
  ArrayIntPtr   _globalToProc;  //[globalID] = procID
  ArrayIntPtr   _globalToLocal; //[globaID] = local id ( _globalToProc and _globalToLocal mostly need to be togethter)

  shared_ptr< StorageSite    > _site;
  shared_ptr< StorageSiteMerger > _siteMerger;
  shared_ptr< CRConnectivity > _conn;
   

  MPI::Intracomm _comm;


};


#endif
