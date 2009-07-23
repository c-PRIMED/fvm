#ifndef _LINEARSYSTEMMERGER_H_
#define _LINEARSYSTEMMERGER_H_

#include "LinearSystem.h"
#include <cassert>
#include <set>
#include <mpi.h>


class LinearSystemMerger
{
public:

    typedef   shared_ptr< Array<int> >     ArrayIntPtr;
    friend class LinearSystem;

   LinearSystemMerger( int target_proc_id, const set<int>& group,  LinearSystem& ls );
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

  int _targetID;
  int _groupID;
  const set<int>&     _group;
  LinearSystem& _lsCoarse;
  
  
  shared_ptr<LinearSystem> _ls;
  
  int _procID;
  int _totalProcs;
  int _totalInterfaces;
  int _totalScatterCells;
  int _totalScatterCellsLocal;
  int _totalGatherCells;
  int _totalGatherCellsLocal;
  
  ArrayIntPtr  _neighMeshCounts; //each local mesh has a certain number of adjacency meshes storing total numbers
  ArrayIntPtr  _scatterInterfaceCounts;
  ArrayIntPtr  _scatterSize;     //size = _totalProcs 
  ArrayIntPtr  _scatterCells;     //size = totalScatterCells
  ArrayIntPtr  _scatterMapIDs;   //size = totalInterfaces
  ArrayIntPtr  _gatherInterfaceCounts;
  ArrayIntPtr  _gatherSize;     //size = _totalProcs 
  ArrayIntPtr  _gatherCells;     //size = totalScatterCells
  ArrayIntPtr  _gatherMapIDs;   //size = totalInterfaces

  ArrayIntPtr _rowLength; //CRConnecivity row_dimension,  size = _totalProcs
  ArrayIntPtr _colLength; //CRConnectivity col_dimension, size = _totalProcs
  ArrayIntPtr  _row;   //CRConnectivity _row, size = totalCells
  ArrayIntPtr  _col;   //CRConnectivity _col, size = ??
  shared_ptr< CRConnectivity > _conn;
  

  MPI::Intracomm _comm;


};


#endif
