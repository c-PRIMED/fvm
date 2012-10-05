// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifdef FVM_PARALLEL

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "LinearSystemMerger.h"
#include "CRConnectivity.h"
#include "StorageSiteMerger.h"


using namespace std;

LinearSystemMerger::LinearSystemMerger(int target_proc_id, const set<int>& group, LinearSystem& ls )
:_targetID(target_proc_id), _groupID(target_proc_id), _group(group), _ls(ls)
{

    init();
    merge();
}

LinearSystemMerger::~LinearSystemMerger()
{
  
}


void
LinearSystemMerger::init()
{

   //construct merged Linearsystem
  _mergeLS = shared_ptr< LinearSystem > ( new LinearSystem() );

   //subcommunicator
   int color = _groupID;
   int key   = MPI::COMM_WORLD.Get_rank();
   _comm = MPI::COMM_WORLD.Split( color, key );
   _procID   = MPI::COMM_WORLD.Get_rank();
   _totalProcs = _group.size();
   _totalInterfaces = 0;
   _totalScatterCells = 0;
   _totalScatterCellsLocal = 0;
   _totalGatherCells = 0;
   _totalGatherCellsLocal = 0;
   _totalCells = 0;
   _neighMeshCounts = ArrayIntPtr( new Array<int>(_totalProcs) ); _neighMeshCounts->zero();
   _scatterSize     = ArrayIntPtr( new Array<int>(_totalProcs) ); _scatterSize->zero();
   _gatherSize      = ArrayIntPtr( new Array<int>(_totalProcs) ); _gatherSize->zero();
   _rowLength       = ArrayIntPtr( new Array<int>(_totalProcs) ); _rowLength->zero();//CRConnecivity row_dimension,  size = _totalProcs
   _colLength       = ArrayIntPtr( new Array<int>(_totalProcs) ); _colLength->zero();
   _selfCounts      = ArrayIntPtr( new Array<int>(_totalProcs) ); _selfCounts->zero();
   _scatterCells.reserve( _totalProcs );
   _gatherCells.reserve( _totalProcs  );
   _gatherIDsLocalToGlobalMap.resize( _totalProcs );
   //get number of all neightbourhood of a local mesh
   get_neigh_mesh_counts();
   get_scatter_cells();
   get_gather_cells();
   get_crconnectivity();
   get_local_to_global_map();
   set_merged_crconnectivity();
   set_ls_vectors();
  
}


void
LinearSystemMerger::merge()
{

   
  // _comm.Reduce( &self_size , &_mergeSiteSize     , 1, MPI::INT, MPI::SUM, _groupID );
  // _comm.Reduce( &ghost_size, &_mergeSiteGhostSize, 1, MPI::INT, MPI::SUM, _groupID );
 
} 



void
LinearSystemMerger::get_neigh_mesh_counts()
{

    vector<int> siteScatterOffsetLocal;
    vector<int> siteGatherOffsetLocal;
    vector<int> scatterMapIDsLocal;
    vector<int> gatherMapIDsLocal;

    int nmesh = 0;
    const MultiField::ArrayIndexList& arrayIndices = _ls.getB().getArrayIndices();
    foreach( MultiField::ArrayIndex k, arrayIndices ){
        const StorageSite& site = *k.second;
        const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
        foreach( const StorageSite::ScatterMap::value_type& pos, scatterMap ){
            const StorageSite& site = *pos.first;
            int  neigh_proc_id  =  site.getGatherProcID();  //neig_proc_id means other mesh, so gatherProcID() will be used for this and gatherMap
            //if neigh proc belongs to the group, then ok
            if ( _group.count( neigh_proc_id ) > 0 ){
                const Array<int>& fromIndices = *pos.second;
                _totalScatterCellsLocal += fromIndices.getLength();
                siteScatterOffsetLocal.push_back( fromIndices.getLength() );
                scatterMapIDsLocal.push_back( neigh_proc_id );
                nmesh++;
            }
       }
    }


   _comm.Gather( &nmesh, 1, MPI::INT,  _neighMeshCounts->getData(), 1, MPI::INT, _targetID );
   _comm.Gather( &_totalScatterCellsLocal, 1, MPI::INT, _scatterSize->getData(), 1, MPI::INT, _targetID );


   foreach( MultiField::ArrayIndex k, arrayIndices ){
       const StorageSite& site = *k.second;
       const StorageSite::GatherMap& gatherMap = site.getGatherMap();
       foreach( const StorageSite::GatherMap::value_type& pos, gatherMap ){
         const StorageSite& site = *pos.first;
         int  neigh_proc_id  =  site.getGatherProcID();
        //if neigh proc belongs to the group, then ok
         if ( _group.count( neigh_proc_id ) > 0 ){
            const Array<int>& toIndices = *pos.second;
            _totalGatherCellsLocal += toIndices.getLength();
            siteGatherOffsetLocal.push_back( toIndices.getLength() );
            gatherMapIDsLocal.push_back( neigh_proc_id );
         }
      }
   }


   _comm.Gather( &_totalGatherCellsLocal, 1, MPI::INT, _gatherSize->getData(), 1, MPI::INT, _targetID );
     //cout << " procid = " << _procID << " _scatterSize = " <<  (*_scatterSize)[_procID] << " _gatherSize =  " << (*_gatherSize)[_procID] << endl;

     //total interfaces
     for ( int i = 0; i < _neighMeshCounts->getLength(); i++ )
        _totalInterfaces += (*_neighMeshCounts)[i];


     ArrayIntPtr scatterInterfaceCounts = ArrayIntPtr( new Array<int>( _totalInterfaces) );
     ArrayIntPtr gatherInterfaceCounts  = ArrayIntPtr( new Array<int>( _totalInterfaces) );
     scatterInterfaceCounts->zero();
     gatherInterfaceCounts->zero();


    int displs[ _totalInterfaces ];
    displs[0] = 0;

    for ( int i = 1; i < _totalInterfaces; i++ )
         displs[i] = displs[i-1] + (*_neighMeshCounts)[i-1];

    //recvcnts = _neighMeshCounts()->getData(),
    _comm.Gatherv( &siteScatterOffsetLocal[0], nmesh, MPI::INT, scatterInterfaceCounts->getData(), &(*_neighMeshCounts)[0], displs, MPI::INT, _targetID );

    //recvcnts = _neighMeshCounts()->getData(),
    _comm.Gatherv( &siteGatherOffsetLocal[0],  nmesh, MPI::INT, gatherInterfaceCounts->getData(), &(*_neighMeshCounts)[0], displs, MPI::INT, _targetID );

   //now  filling our convenient variable
   int indx = 0;
   foreach ( const set<int>::value_type proc, _group ){
       int num_interfaces = (*_neighMeshCounts)[proc];
       ArrayIntPtr scatter_counts ( new Array<int> ( num_interfaces ) );
       ArrayIntPtr gather_counts  ( new Array<int> ( num_interfaces ) );
       for ( int n = 0; n < num_interfaces; n++ ){
           (*scatter_counts)[n] = (*scatterInterfaceCounts)[indx];
           (*gather_counts)[n]  = (*gatherInterfaceCounts )[indx];
           indx++;
       }
      _scatterInterfaceCounts[proc] = scatter_counts;
      _gatherInterfaceCounts[proc]  = gather_counts;
   }

   for ( int i = 0; i < scatterInterfaceCounts->getLength(); i++ )
      _totalScatterCells += (*scatterInterfaceCounts)[i];

   for ( int i = 0; i < gatherInterfaceCounts->getLength();  i++ )
      _totalGatherCells  += (*gatherInterfaceCounts)[i];


   //recvcnts for scatter
   int scatter_recvcnts[_totalProcs];
   for ( int proc = 0; proc < _totalProcs; proc++ )
       scatter_recvcnts[proc] = _scatterInterfaceCounts[proc]->getLength();

  //recvncnts for gather
   int gather_recvcnts[_totalProcs];
   for ( int proc = 0; proc < _totalProcs; proc++ )
       gather_recvcnts[proc] = _gatherInterfaceCounts[proc]->getLength();

   ArrayIntPtr scatterMapIDs( new Array<int>( _totalInterfaces) );

    //scatter neightbour IDS
    _comm.Gatherv( &scatterMapIDsLocal[0], nmesh,  MPI::INT, scatterMapIDs->getData(), scatter_recvcnts, displs, MPI::INT, _targetID );


    ArrayIntPtr   gatherMapIDs( new Array<int>( _totalInterfaces) );
    //scatter neighbour IDS
    _comm.Gatherv( &gatherMapIDsLocal[0], nmesh,  MPI::INT, gatherMapIDs->getData(), gather_recvcnts, displs, MPI::INT, _targetID );


    indx = 0;
    foreach ( const set<int>::value_type proc, _group ){
       int num_interfaces = (*_neighMeshCounts)[proc];
       ArrayIntPtr scatterIDs ( new Array<int> ( num_interfaces ) );
       ArrayIntPtr gatherIDs  ( new Array<int> ( num_interfaces ) );
       for ( int n = 0; n < num_interfaces; n++ ){
           (*scatterIDs)[n] = (*scatterMapIDs)[indx];
           (*gatherIDs )[n] = (*gatherMapIDs )[indx];
           indx++;
       }
      _scatterInterfaceIDs[proc] = scatterIDs;
      _gatherInterfaceIDs [proc] = gatherIDs;
    }


}

//merging scatterArray
void
LinearSystemMerger::get_scatter_cells()
{ 

 
   //fill scatter cells in local for shipping
   
    ArrayIntPtr   scatterCellsLocal( new Array<int> ( _totalScatterCellsLocal  ) );
    int indx = 0;
   const MultiField::ArrayIndexList& arrayIndices = _ls.getB().getArrayIndices();
   foreach( MultiField::ArrayIndex k, arrayIndices ){
       const StorageSite& site = *k.second;
       const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
       foreach( const StorageSite::ScatterMap::value_type& pos, scatterMap ){
          const StorageSite& site = *pos.first;
          int  neigh_proc_id  =  site.getGatherProcID();
         //if neigh proc belongs to the group, then ok
          if ( _group.count( neigh_proc_id ) > 0 ){
             const Array<int>& fromIndices = *pos.second;
             for ( int i = 0; i < fromIndices.getLength(); i++) 
                  (*scatterCellsLocal)[indx++] = fromIndices[i];
          }
       }
   }

    //fill in scatterCells
    ArrayIntPtr scatterCells( new Array<int> ( _totalScatterCells) );


    int displs[ _totalProcs ];
    displs[0] = 0;
    for ( int i = 1; i < _totalProcs; i++ )
         displs[i] = displs[i-1] + (*_scatterSize)[i-1];


     //scatterSize is receive count array, 
    _comm.Gatherv( scatterCellsLocal->getData(), _totalScatterCellsLocal, MPI::INT,
                   scatterCells->getData(), &(*_scatterSize)[0], displs, MPI::INT, _targetID );


   //convenient storage in map
   indx = 0;
   foreach ( const set<int>::value_type proc, _group ){
       const Array<int>& interfaceIDs    = *_scatterInterfaceIDs[proc];
       const Array<int>& interfaceCounts = *_scatterInterfaceCounts[proc]; 
       map<int, ArrayIntPtr> cell_map;
       for ( int i = 0; i < interfaceIDs.getLength(); i++ ){
           ArrayIntPtr  cells( new Array<int> ( interfaceCounts[i] ) );
           for ( int n = 0; n < interfaceCounts[i]; n++ )
               (*cells)[n] = (*scatterCells)[indx++];
            int interface_id = interfaceIDs[i];
            cell_map[interface_id] = cells;
       }
       _scatterCells.push_back( cell_map );
   }


}


//merging scatterArray
void
LinearSystemMerger::get_gather_cells()
{ 

   //fill gather cells in local for shipping
    ArrayIntPtr  gatherCellsLocal( new Array<int> ( _totalGatherCellsLocal ) );
    int indx = 0;
    const MultiField::ArrayIndexList& arrayIndices = _ls.getB().getArrayIndices();
    foreach( MultiField::ArrayIndex k, arrayIndices ){
       const StorageSite& site = *k.second;
       const StorageSite::GatherMap& gatherMap = site.getGatherMap();
       foreach( const StorageSite::GatherMap::value_type& pos, gatherMap ){
          const StorageSite& site = *pos.first;
          int  neigh_proc_id  =  site.getGatherProcID();
          //if neigh proc belongs to the group, then ok
          if ( _group.count( neigh_proc_id ) > 0 ){
             const Array<int>& toIndices = *pos.second;
             for ( int i = 0; i < toIndices.getLength(); i++ ) 
                 (*gatherCellsLocal)[indx++] = toIndices[i];
          }
       }
    }

    //fill in gatherCells
     ArrayIntPtr gatherCells ( new Array<int> ( _totalGatherCells ) );

    int displs[ _totalProcs ];
    displs[0] = 0;
    for ( int i = 1; i < _totalProcs; i++ )
         displs[i] = displs[i-1] + (*_gatherSize)[i-1];

     //scatterSize is receive count array, 
    _comm.Gatherv( gatherCellsLocal->getData(), _totalGatherCellsLocal, MPI::INT,
                   gatherCells->getData(), &(*_gatherSize)[0], displs, MPI::INT, _targetID );

   //conveient storage in map 
   indx = 0;
   foreach ( const set<int>::value_type proc, _group ){
       const Array<int>& interfaceIDs    = *_gatherInterfaceIDs[proc];
       const Array<int>& interfaceCounts = *_gatherInterfaceCounts[proc]; 
       map<int, ArrayIntPtr> cell_map;
       for ( int i = 0; i < interfaceIDs.getLength(); i++ ){
           ArrayIntPtr  cells( new Array<int> ( interfaceCounts[i] ) );
           for ( int n = 0; n < interfaceCounts[i]; n++ )
               (*cells)[n] = (*gatherCells)[indx++];

           int interface_id = interfaceIDs[i];
           cell_map[interface_id] = cells;
       }
       _gatherCells.push_back( cell_map );
   } 



}


void
LinearSystemMerger::get_crconnectivity()
{

    const MultiField::ArrayIndexList& arrayIndices = _ls.getB().getArrayIndices();
    foreach( MultiField::ArrayIndex k, arrayIndices ){
       const CRConnectivity& conn = _ls.getMatrix().getMatrix( k, k).getConnectivity();
       const StorageSite& site = conn.getRowSite();
       int inner_cells = site.getSelfCount();
       int rowLength = conn.getRow().getLength();
       int colLength = conn.getCol().getLength();
       _comm.Gather( &rowLength  , 1, MPI::INT, _rowLength->getData() , 1, MPI::INT, _targetID );
       _comm.Gather( &colLength  , 1, MPI::INT, _colLength->getData() , 1, MPI::INT, _targetID );
       _comm.Gather( &inner_cells, 1, MPI::INT, _selfCounts->getData(), 1, MPI::INT, _targetID );
   
    }

  //compute aggloremation size
   int row_size = 0;
   int col_size = 0;
   for ( int i = 0; i < _totalProcs; i++ ){
       row_size  += (*_rowLength)[i];
       col_size  += (*_colLength)[i];
   }


   //forming arrays
    ArrayIntPtr  row( new Array<int> ( row_size ) );
    ArrayIntPtr  col( new Array<int> ( col_size ) );

    int displs_row[ _totalProcs ];
    int displs_col[ _totalProcs ];
    displs_row[0] = 0;
    displs_col[0] = 0;
    for ( int i = 1; i < _totalProcs; i++ ){
         displs_row[i] = displs_row[i-1] + (*_rowLength)[i-1];
         displs_col[i] = displs_col[i-1] + (*_colLength)[i-1];
    }


  //getting _row and _col
    foreach( MultiField::ArrayIndex k, arrayIndices ){
       const CRConnectivity& conn = _ls.getMatrix().getMatrix( k, k).getConnectivity();
       const Array<int>&  row_local = conn.getRow();
       const Array<int>&  col_local = conn.getCol();
       int rowLength = row_local.getLength();
       int colLength = col_local.getLength();
      _comm.Gatherv( row_local.getData(), rowLength, MPI::INT, row->getData(), &(*_rowLength)[0], displs_row, MPI::INT, _targetID );
      _comm.Gatherv( col_local.getData(), colLength, MPI::INT, col->getData(), &(*_colLength)[0], displs_col, MPI::INT, _targetID );
   }

  //for ( int i = 0; i < col->getLength(); i++ )
    // cout << " col(" << i << ") = " << (*col)[i] << endl;

   //copying to map data structures, procid->row, or procid->col
   int local = 0;
   foreach ( const set<int>::value_type proc_id, _group ){
     _row[proc_id] = ArrayIntPtr ( new Array<int> ( (*_rowLength)[local] ) );
     _col[proc_id] = ArrayIntPtr ( new Array<int> ( (*_colLength)[local] ) );
     //_colPos[proc_id] = ArrayIntPtr ( new Array<int> ( (*_colLength)[local] ) );
      local++;
   }

   int indx = 0;
   local = 0;
   foreach ( const set<int>::value_type proc_id, _group ){
      for ( int i = 0; i < (*_rowLength)[local]; i++ )
         (*_row[proc_id])[i] = (*row)[indx++];
      local++;
   }

   indx  = 0;
   local = 0;
   foreach ( const set<int>::value_type proc_id, _group ){
      for ( int i = 0; i < (*_colLength)[local]; i++ )
         (*_col[proc_id])[i] = (*col)[indx++];
      local++;
   }

//   //create colPos data structure, you ask for cellID and returns position in col
//   indx  = 0;
//   local = 0;
//   foreach( const set<int>::value_type proc_id, _group ){
//      for ( int i = 0; i < (*_colLength)[local]; i++ ){
//         int cell_id = (*col)[indx++];
//         (*_colPos[proc_id])[cell_id] = i;
//      }
//      local++;
//   }



}


void 
LinearSystemMerger::get_local_to_global_map()
{
    //localToGlobal: assuming numbering start with lowest proc ID to the highest one
    //selfCounts is array from 0, 1, ... NPROCS-1 but groupID is not necessaryl holding that
    //future we might group processes so it will not be always from 0, 1, .. NPROCS-1
    int indx = 0;
    foreach ( const set<int>::value_type proc_id, _group){
       _localToGlobal[proc_id] = ArrayIntPtr( new Array<int> ( (*_selfCounts)[indx] ) );
       _totalCells +=  (*_selfCounts)[indx];
        indx++;
    }

   _globalToProc  = ArrayIntPtr( new Array<int> ( _totalCells ) );
   _globalToLocal = ArrayIntPtr( new Array<int> ( _totalCells ) );

    indx = 0;
    int global_id = 0;
    foreach ( const set<int>::value_type proc_id, _group ){
        for ( int n = 0; n < (*_selfCounts)[indx]; n++ ){
            (*_localToGlobal[proc_id])[n] = global_id;
            (*_globalToProc )[global_id] = proc_id;
            (*_globalToLocal)[global_id] = n;
            global_id++;
        }
        indx++;
    }


}


//merging connectivities
void  
LinearSystemMerger::set_merged_crconnectivity()
{
   //mergin sites
   const MultiField::ArrayIndexList& arrayIndices = _ls.getB().getArrayIndices();
   foreach( MultiField::ArrayIndex k, arrayIndices ){
       const StorageSite& cell_site = *k.second;
       _siteMerger = shared_ptr< StorageSiteMerger > ( new StorageSiteMerger( _targetID, _group, cell_site) );
    }

    //sev
    _site = _siteMerger->merge();
    update_gatherCells_from_scatterCells();

   _mergeCR = shared_ptr< CRConnectivity > ( new CRConnectivity(*_site, *_site) );
   //initCount
   _mergeCR->initCount();

   //addcount
   int indx = 0;
   foreach ( set<int>::value_type proc, _group ){
       const Array<int>& row = *_row[proc];
       const Array<int>& col = *(_col[proc]);
       const map<int,int>& gatherIDsLocalToGlobal  = _gatherIDsLocalToGlobalMap[proc];
       
       int ghost_endID   = row.getLength()-1;
       int ghost_startID = ghost_endID - int( gatherIDsLocalToGlobal.size() );
       int selfCount = (*_selfCounts)[proc];

       //loop over row
       for ( int i = 0; i < selfCount; i++ ){
            int count = row[i+1] - row[i];
//            _mergeCR->addCount(indx++,count);
            for ( int j = 0; j < count; j++){
              int localID  = col[ row[i]+j ]; //get surrounding local cellID
              //excluding boundary conditions
              if ( !(localID >= selfCount && localID < ghost_startID) )
                _mergeCR->addCount( indx, 1 );
            }
            indx++;
       }
   }

    //finish count
   _mergeCR->finishCount();

   //add
   indx = 0;
   foreach ( set<int>::value_type proc, _group ){
       const Array<int>& row = *(_row[proc]);
       const Array<int>& col = *(_col[proc]);
       const Array<int>& localToGlobal = *_localToGlobal[proc];
       const map<int,int>& gatherIDsLocalToGlobal  = _gatherIDsLocalToGlobalMap[proc];

       int ghost_endID   = row.getLength()-1;
       int ghost_startID = ghost_endID - int( gatherIDsLocalToGlobal.size() );
       int selfCount     = (*_selfCounts)[proc];

       //loop over row
       for ( int i = 0; i < selfCount; i++ ){
            int count = row[i+1] - row[i];
            //loop over surrounding cells
            for ( int j = 0; j < count; j++ ){
               int localID  = col[ row[i]+j ]; //get surrounding local cellID
               int globalID = -1;
               //ghost cell check if inner cells or ghost cells
               if ( localID < selfCount ){
                   globalID = localToGlobal[ localID ];
               } else {
                   globalID = gatherIDsLocalToGlobal.find(localID)->second;
               }

              //excluding boundary conditions
              if ( !(localID >= selfCount && localID < ghost_startID) )
                 _mergeCR->add( indx, globalID );
            }
            indx++;
       }
    }
    //finish add
    _mergeCR->finishAdd();

    //update ghost cells as globalID
    // update_col();

//    const Array<int>& col = _mergeCR->getCol();
//     cout << "llllllllllllllllllllllllllllllllllllllllcol =  " << col.getLength() << endl;
/*   _mergeColPos = ArrayIntPtr( new Array<int> ( col.getLength() ) );
   //store col position for given id
    for ( int i = 0; i < col.getLength(); i++ ){
       int id = col[i];
       (*_mergeColPos)[id] = i;
    }*/
}


void
LinearSystemMerger::update_gatherCells_from_scatterCells()
{
    //get scatter cells from neighbouring mesh and update the current mesh gather cells
    //loop over meshes

    foreach ( set<int>::value_type proc, _group ){
        const Array<int>& gatherInterfaceIDs   = *_gatherInterfaceIDs[proc];
        //const Array<int>& scatterInterfaceIDs  = *_scatterInterfaceIDs[proc];
        //loop over interfaces  on this mesh
        for ( int n = 0; n < gatherInterfaceIDs.getLength(); n++ ){
            int neighMeshID  = gatherInterfaceIDs[n];
            const Array<int>& gatherCells  = *_gatherCells [proc][neighMeshID];
            const Array<int>& scatterCells = *_scatterCells[neighMeshID][proc];
            //update gatherCells
            for ( int i = 0; i < gatherCells.getLength(); i++ ){
                int cellID = (*_localToGlobal[neighMeshID])[ scatterCells[i] ]; //get global numbering of scatter cells on other mesh
                int gatherID = gatherCells[i];
               _gatherIDsLocalToGlobalMap[proc][ gatherID ] = cellID;
                //gatherCells[i] = cellID;
            }
        }
    }




}

//update col
void
LinearSystemMerger::update_col()
{
   //update _col for gather Cells
   foreach( set<int>::value_type proc_id, _group ){
      Array<int>& col = *_col[proc_id];
      const Array<int>& row = *_row[proc_id];
      const map<int,int>& gatherIDsLocalToGLobal = _gatherIDsLocalToGlobalMap[proc_id];

      int ghost_endID   = row.getLength()-1;
      int ghost_startID = ghost_endID - int( gatherIDsLocalToGLobal.size() );
      int selfCount = (*_selfCounts)[proc_id];

      for ( int i = 0; i < col.getLength(); i++ ){
         //check if it is ghost cells
         if ( col[i] >= ghost_startID && col[i] <= ghost_endID )
            col[i] = gatherIDsLocalToGLobal.find( col[i] )->second;
         
        //mark boundary as -1 
         if( col[i] >= selfCount && col[i] < ghost_startID ){
            col[i] = -1;
         }
      }


   }

}

//merging x, b, residual  from coarse linear systems
void
LinearSystemMerger::set_ls_vectors()
{
    _mergeLS->_b = shared_ptr<MultiField> ( new MultiField() );

    const MultiField&  delta = _ls.getDelta();
    const MultiField::ArrayIndexList& arrayIndices = delta.getArrayIndices();
    foreach ( const MultiField::ArrayIndex& ai, arrayIndices){
         const StorageSite&  mergeSite = *_site;
         MultiField::ArrayIndex mergeIndex( ai.first, &mergeSite );
          //cout << " _procIDDDD=  " << _procID << " mergeSite.getCount() = " << mergeSite.getCount() << endl;
         _mergeLS->_b->addArray( mergeIndex, delta[ai].newSizedClone( mergeSite.getCount()) );
/*
        const Array<double> & xarray = dynamic_cast< const Array<double> & > ( delta[ai] );
        cout << " procid = " << _procID << " xarray.getLength() = " << xarray.getLength() << endl;
        for ( int i = 0; i < xarray.getLength(); i++ )
            cout << " xarray[" << i << "] = " << xarray[i] << endl;*/

    }

   _mergeLS->_b->zero();
   _mergeLS->_delta = dynamic_pointer_cast<MultiField> ( _mergeLS->_b->newClone() );
   _mergeLS->_delta->zero();
   _mergeLS->_residual = dynamic_pointer_cast<MultiField> ( _mergeLS->_b->newClone() );
   _mergeLS->_residual->zero();
 
   //_mergeLS->initAssembly();
   //_mergeLS->getMatrix().initAssembly();

}

//constructing mergin matrix
void
LinearSystemMerger::gatherMatrix()
{

    //getting  diag and off diag size (in BYTE )
    Array<int> sizeDiag   ( _totalProcs );
    Array<int> sizeOffDiag( _totalProcs );
    sizeDiag.zero();
    sizeOffDiag.zero();
    const MultiField::ArrayIndexList& arrayIndices = _ls.getB().getArrayIndices();
    foreach( MultiField::ArrayIndex k, arrayIndices ){
        int send_buffer = _ls.getMatrix().getMatrix(k,k).getDiagDataSize();
       _comm.Gather( &send_buffer , 1, MPI::INT, sizeDiag.getData(), 1, MPI::INT, _targetID );
        send_buffer = _ls.getMatrix().getMatrix(k,k).getOffDiagDataSize();
       _comm.Gather( &send_buffer , 1, MPI::INT, sizeOffDiag.getData(), 1, MPI::INT, _targetID );
    }

    //estimate size for diag and offdiag arrays for gathering
    int size_diag    = 0;
    int size_offdiag = 0;
    foreach ( set<int>::value_type proc, _group ){
        sizeDiag[proc]    = sizeDiag[proc]    / sizeof( double );
        sizeOffDiag[proc] = sizeOffDiag[proc] / sizeof( double );
        size_diag    += sizeDiag[proc];
        size_offdiag += sizeOffDiag[proc];
    }

    //
    int displs_diag[ _totalProcs ];
    int displs_offdiag[ _totalProcs ];
    displs_diag[0]    = 0;
    displs_offdiag[0] = 0;
    for ( int i = 1; i < _totalProcs; i++ ){
       displs_diag[i]    = displs_diag[i-1]    + sizeDiag[i-1];
       displs_offdiag[i] = displs_offdiag[i-1] + sizeOffDiag[i-1];
    }


    //allocating space for diag and offdiag
    ArrayDblePtr  diag   ( new Array<double> ( size_diag    ) );
    ArrayDblePtr  offdiag( new Array<double> ( size_offdiag ) );
//     cout << " procIDdDDDDDDDDDDDDD = " << _procID << " size_diag = " << size_diag << " size_offdiag = " << size_offdiag << endl;
//     if ( _procID == _targetID )
//        for ( int n = 0; n < _totalProcs; n++ )
//           cout << " sizeDiag[" << n << "] = " << sizeDiag[n] << " sizeOffDiag[" << n << "] = " << sizeOffDiag[n] << endl;
// 


    //gathering diag and offdiag 
    foreach( MultiField::ArrayIndex k, arrayIndices ){


//        //assigning fake value before shipping
//        double *diagPtr = reinterpret_cast < double* > (_ls.getMatrix().getMatrix(k,k).getDiagData() );
//        int diag_cnt = _ls.getMatrix().getMatrix(k,k).getDiagDataSize()    / sizeof( double );
//        int off_cnt  = _ls.getMatrix().getMatrix(k,k).getOffDiagDataSize() / sizeof( double );
//        for ( int i = 0; i < diag_cnt; i++)
//            *(diagPtr+i) = double(_procID * 100 + i);
// 
//        double *offDiagPtr = reinterpret_cast < double* > (_ls.getMatrix().getMatrix(k,k).getOffDiagData() );
//        for ( int i = 0; i < off_cnt; i++)
//            *(offDiagPtr+i) = double(_procID * 100 + i);


       int send_cnts_diag    = _ls.getMatrix().getMatrix(k,k).getDiagDataSize();
      _comm.Gatherv( _ls.getMatrix().getMatrix(k,k).getDiagData(), send_cnts_diag, MPI::BYTE, diag->getData(), &sizeDiag[0], displs_diag, MPI::DOUBLE, _targetID );
       int send_cnts_offdiag = _ls.getMatrix().getMatrix(k,k).getOffDiagDataSize();
      _comm.Gatherv( _ls.getMatrix().getMatrix(k,k).getOffDiagData(), send_cnts_offdiag, MPI::BYTE, offdiag->getData(), &sizeOffDiag[0], displs_offdiag, MPI::DOUBLE, _targetID );
    }

    //fill diag
    int indx = 0;
    foreach ( const set<int>::value_type proc_id, _group){
       _diag[proc_id]    = ArrayDblePtr( new Array<double> ( sizeDiag[proc_id]    ) );
        Array<double>& diagonal   = *_diag[proc_id];
       //fill data
        for ( int i = 0; i < sizeDiag[proc_id]; i++ ){
           diagonal[i]    = (*diag)[indx];
           indx++;
        }
    }
    //fill off diaogonal
    indx = 0;
    foreach ( const set<int>::value_type proc_id, _group){
       _offDiag[proc_id]    = ArrayDblePtr( new Array<double> ( sizeOffDiag[proc_id] ) );
        Array<double>& offDiagonal   = *_offDiag[proc_id];
       //fill data
        for ( int i = 0; i < sizeOffDiag[proc_id]; i++ ){
           offDiagonal[i]  = (*offdiag)[indx];
           indx++;
        }
    }



   int i = 0;
   MultiField& mergeMultiField = _mergeLS->getDelta();
   foreach( MultiField::ArrayIndex k, arrayIndices ){
       const MultiField::ArrayIndex& mergeIndex = mergeMultiField.getArrayIndex(i++);
       Matrix& mIJ = _ls.getMatrix().getMatrix( k, k );
      _mergeLS->getMatrix().addMatrix( mergeIndex, mergeIndex , mIJ.createMergeMatrix( *this) );
   }


}





void 
LinearSystemMerger::gatherB()
{

    int displs[ _totalProcs ];
    displs[0] = 0;
    for ( int i = 1; i < _totalProcs; i++ ){
         displs[i] = displs[i-1] + (*_selfCounts)[i-1];
     }

    //getting delta
   int i = 0;
   const MultiField::ArrayIndexList& arrayIndices = _ls.getB().getArrayIndices();
   const MultiField& multiField = _ls.getB();
   MultiField& mergeMultiField = _mergeLS->getB();
   _mergeLS->getDelta().zero();
   foreach( MultiField::ArrayIndex k, arrayIndices ){
       const StorageSite& site = *k.second;
       int send_cnts = site.getSelfCount();
       const ArrayBase&  B = multiField[k];
/*       Array<double>&  deltaDouble = dynamic_cast<  Array<double>&  > ( delta );
       for ( int indx = 0; indx < deltaDouble.getLength(); indx++ ){
           deltaDouble[indx] = double(_procID) * double(100) + double(indx) ;
           cout << "proc_id = " << _procID << " delta(" << indx << ") = " << deltaDouble[indx] << endl;
       }*/
       const MultiField::ArrayIndex& mergerIndex = mergeMultiField.getArrayIndex(i++);
       ArrayBase&  mergeB      = mergeMultiField[mergerIndex];
      _comm.Gatherv( B.getData(), send_cnts, MPI::DOUBLE, mergeB.getData(), &(*_selfCounts)[0], displs, MPI::DOUBLE, _targetID );

   }


}

void
LinearSystemMerger::scatterDelta()
{
   MPI::COMM_WORLD.Barrier();
   int displs[ _totalProcs ];
   displs[0] = 0;
   for ( int i = 1; i < _totalProcs; i++ )
       displs[i] = displs[i-1] + (*_selfCounts)[i-1];

   int i = 0;
   MultiField& multiField = _ls.getDelta();
   const MultiField::ArrayIndexList& arrayIndices = _ls.getDelta().getArrayIndices();
   foreach( MultiField::ArrayIndex k, arrayIndices ){
       const StorageSite& site = *k.second;
       int recv_cnts = site.getSelfCount();
       ArrayBase&  delta = multiField[k];

       const MultiField& mergeMultiField = _mergeLS->getDelta();
       const MultiField::ArrayIndex& mergerIndex = mergeMultiField.getArrayIndex(i++);
       const ArrayBase&  mergeDelta      = mergeMultiField[mergerIndex];
      _comm.Scatterv( mergeDelta.getData(), &(*_selfCounts)[0], displs, MPI::DOUBLE, delta.getData(), recv_cnts,  MPI::DOUBLE, _targetID );

   }


//     MPI::COMM_WORLD.Barrier();
//     cout << " _totalScatterCellss = " << _totalScatterCells << " procid = " << _procID << endl;
//     MPI::COMM_WORLD.Barrier();
// 
//     debug_print();

}


void
LinearSystemMerger::debug_print()
{

  if ( _procID == _targetID  ){
     stringstream ss;
     ss << "proc" << MPI::COMM_WORLD.Get_rank() << "_linear_system_merger.dat";
     ofstream  debug_file( (ss.str()).c_str() );

   
     //neighbourhd mesh counts
     for ( int n = 0; n < _totalProcs; n++ ) 
        debug_file << " _neighMeshCounts[" << n << "] = " << (*_neighMeshCounts)[n] << endl;
     debug_file << endl;

     //scatter size for each mesh
     for ( int n = 0; n < _totalProcs; n++ ) 
        debug_file << " _scatterSize[" << n << "] = " << (*_scatterSize)[n] << endl;
     debug_file << endl;

     //gather size for each mesh
     for ( int n = 0; n < _totalProcs; n++ ) 
        debug_file << " _gatherSize[" << n << "] = " << (*_gatherSize)[n] << endl;
     debug_file << endl;

      //total interfaces
      debug_file << " _totalInterfaces = " << _totalInterfaces  << endl;
      debug_file << endl;

      //total scatter cells
      debug_file << " _totalScatterCells = " << _totalScatterCells  << endl;
      debug_file << endl;

      //total gather cells
      debug_file << " _totalGatherCells = " << _totalGatherCells  << endl;
      debug_file << endl;


      //scatter and gather counts for each process
      foreach ( set<int>::value_type proc, _group ){
         debug_file << "procID : " << proc << endl;
         const Array<int>& scatter_counts = *_scatterInterfaceCounts[proc];
         const Array<int>&  gather_counts = *_gatherInterfaceCounts[proc];
         for ( int n = 0; n < scatter_counts.getLength(); n++ )
              debug_file << " scatter_counts[" << n << "] = " << scatter_counts[n] << endl;
         debug_file << endl;

         for ( int n = 0; n < gather_counts.getLength(); n++ )
              debug_file << " gather_counts[" << n << "] = " << gather_counts[n] << endl;
         debug_file << endl;
       }


     //scatter interface ids (neightbourhood mesh ids )
      int indx = 0;
      foreach ( set<int>::value_type proc, _group ){
         debug_file << "procID : " << proc << endl;
         const Array<int>& scatterInterfaceIDs  = *_scatterInterfaceIDs[proc];
         for ( int n = 0; n < scatterInterfaceIDs.getLength(); n++ )
             debug_file << " scatterInterfaceIDs[" << n << "] = " << scatterInterfaceIDs[n] << endl;
         debug_file  << endl;

      }

     //scatter cells
      foreach ( set<int>::value_type proc, _group ){
          debug_file << "procID : " << proc << endl;
          //loop over scatter segments
          const map<int, ArrayIntPtr>& cell_map = _scatterCells[proc];
          foreach ( const ArrayIntPtrMap::value_type& pos, cell_map ){
              debug_file << "    scatterID = " << pos.first <<  endl;
              for ( int n = 0; n < pos.second->getLength(); n++ )
                  debug_file << "         " << (*pos.second)[n] << endl;
              debug_file << endl;
          }
          debug_file << endl;
      } 

     //gather interface ids (neightbourhood mesh ids )
      indx = 0;
      foreach ( set<int>::value_type proc, _group ){
         debug_file << "procID : " << proc << endl;
         const Array<int>& gatherInterfaceIDs  = *_gatherInterfaceIDs[proc];
         for ( int n = 0; n < gatherInterfaceIDs.getLength(); n++ )
             debug_file << " gatherInterfaceIDs[" << n << "] = " << gatherInterfaceIDs[n] << endl;
         debug_file  << endl;
      }

      //gather cells
      foreach ( set<int>::value_type proc, _group ){
          debug_file << "procID : " << proc << endl;
          //loop over scatter segments
          const map<int, ArrayIntPtr>& cell_map = _gatherCells[proc];
          foreach ( const ArrayIntPtrMap::value_type& pos, cell_map ){
              debug_file << "    gatherID = " << pos.first <<  endl;
              for ( int n = 0; n < pos.second->getLength(); n++ )
                  debug_file << "         " << (*pos.second)[n] << endl;
              debug_file << endl;
          }
          debug_file << endl;
     } 

     //rowLength
     debug_file  << " _rowLength : " << endl;
     for ( int n = 0; n < _rowLength->getLength(); n++ )
          debug_file << " _rowLength[ " << n  << "] = " << (*_rowLength)[n] << endl;
     debug_file << endl;
     //colLength
     debug_file  << " _colLength : " << endl;
     for ( int n = 0; n < _colLength->getLength(); n++ )
          debug_file << " _colLength[ " << n  << "] = " << (*_colLength)[n] << endl;
     debug_file << endl;
     //selfCounts
     debug_file  << " _selfCounts : " << endl;
     for ( int n = 0; n < _selfCounts->getLength(); n++ )
          debug_file << " _selfCounts[ " << n  << "] = " << (*_selfCounts)[n] << endl;
     debug_file << endl;

     //local CRConnectivity
     foreach ( set<int>::value_type proc, _group ){
        debug_file << "procID : " << proc << endl;
        const Array<int>& row = *(_row[proc]);
        const Array<int>& col = *(_col[proc]);
        for ( int i = 0; i < row.getLength()-1; i++ ){
           int col_dim = row[i+1] - row[i];
           for ( int j = 0; j < col_dim; j++ ){
              debug_file << " coarseCR(" << i << "," << j << ") = " << col[ row[i]+j ] << endl;
           }
        }
     }
     debug_file << endl;

    //maps between local to Global, global to local
    indx = 0;
    foreach ( set<int>::value_type proc, _group ){
       debug_file << "procID: " << proc << endl;
       for ( int n = 0; n < (*_selfCounts)[indx]; n++ )
          debug_file << "localToGlobal[" << n << "] = " << (*_localToGlobal[proc])[n] << endl;
       indx++;
    }
    debug_file << endl;

    //_globalToProc, _globalToLocal
    for ( int n = 0; n < _totalCells; n++ )
       debug_file << "_globalToLocal[" << n << "] = " << (*_globalToLocal)[n] << " resides at procID = " << (*_globalToProc)[n] << endl;
    debug_file << endl;
 
      //gather cells
    foreach ( set<int>::value_type proc, _group ){
        debug_file << "procID (global numbering) : " << proc << endl;
        //loop over scatter segments
        const map<int, ArrayIntPtr>& cell_map = _gatherCells[proc];
        foreach ( const ArrayIntPtrMap::value_type& pos, cell_map ){
            debug_file << "    gatherID = " << pos.first <<  endl;
            for ( int n = 0; n < pos.second->getLength(); n++ )
                debug_file << "         " << (*pos.second)[n] << endl;
            debug_file << endl;
        }
        debug_file << endl;
    } 
    debug_file << endl;
    //merge CRconnectivity
    int row_dim = _mergeCR->getRowDim();
    for ( int i = 0; i < row_dim; i++ ){
       int col_dim = _mergeCR->getCount(i);
       for ( int j = 0; j < col_dim; j++ ){
          debug_file << " mergeCR(" << i << "," << j << ") = " << (*_mergeCR)(i,j) << endl;
       }
    }
    debug_file << endl;



   //mergeCol
    const Array<int>& mergeCol = _mergeCR->getCol();
    for ( int i = 0; i < mergeCol.getLength(); i++ )
       debug_file << " mergeCol(" << i << ") = " << mergeCol[i] << endl;
    debug_file << endl;

//    //mergeColPos
//     for ( int i = 0; i < _mergeColPos->getLength(); i++ )
//        debug_file << " mergeColPos(" << i << ") = " << (*_mergeColPos)[i] << endl;
//     debug_file << endl;

   //delta
   const MultiField::ArrayIndexList& arrayIndices = _mergeLS->getDelta().getArrayIndices();
    foreach( MultiField::ArrayIndex k, arrayIndices ){
       const MultiField& multiField = _mergeLS->getDelta();
       const Array<double>&  delta = dynamic_cast< const Array<double>&  > ( multiField[k] );
       for ( int i = 0; i < delta.getLength(); i++ )
           debug_file << " delta(" << i << ") = " << delta[i] << endl;

   }

    //matrix (diag)
    foreach ( set<int>::value_type proc, _group ){
         debug_file << "procID : " << proc << endl;
         const Array<double>& diag  = *_diag[proc];
         for ( int n = 0; n < diag.getLength(); n++ )
             debug_file << " diag[" << n << "] = " << diag[n] << endl;
         debug_file  << endl;
      }
     
     debug_file << endl;
    //matrix (diag)
    foreach ( set<int>::value_type proc, _group ){
         debug_file << "procID : " << proc << endl;
         const Array<double>& offDiag  = *_offDiag[proc];
         for ( int n = 0; n < offDiag.getLength(); n++ )
             debug_file << " offDiag[" << n << "] = " << offDiag[n] << endl;
         debug_file  << endl;
      }
     



     debug_file.close();
   }
   

}

#endif
