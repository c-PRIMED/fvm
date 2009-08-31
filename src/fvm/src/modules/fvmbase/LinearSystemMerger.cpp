
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "LinearSystemMerger.h"
#include "CRConnectivity.h"


using namespace std;

LinearSystemMerger::LinearSystemMerger(int target_proc_id, const set<int>& group, LinearSystem& ls )
:_targetID(target_proc_id), _groupID(target_proc_id), _group(group), _lsFine(ls)
{

  //cout << " matrix size = " << _lsFine.getMatrix().getSize()  << endl;


    init();
    merge();
}

LinearSystemMerger::~LinearSystemMerger()
{
  
}


void
LinearSystemMerger::init()
{
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
   //get number of all neightbourhood of a local mesh
   get_neigh_mesh_counts();
   get_scatter_cells();
   get_gather_cells();
   get_crconnectivity();
   get_local_to_global_map();
  
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
    const MultiFieldMatrix::MatrixMappersMap& scatterMap = _lsFine.getMatrix()._coarseScatterMaps;
    foreach( MultiFieldMatrix::MatrixMappersMap::value_type pos, scatterMap ){

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
   

   _comm.Gather( &nmesh, 1, MPI::INT,  _neighMeshCounts->getData(), 1, MPI::INT, _targetID );
   _comm.Gather( &_totalScatterCellsLocal, 1, MPI::INT, _scatterSize->getData(), 1, MPI::INT, _targetID );


    const MultiFieldMatrix::MatrixMappersMap& gatherMap = _lsFine.getMatrix()._coarseGatherMaps;
    foreach( MultiFieldMatrix::MatrixMappersMap::value_type pos, gatherMap ){
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
   foreach ( set<int>::value_type proc, _group ){
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
    foreach ( set<int>::value_type proc, _group ){
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
    const MultiFieldMatrix::MatrixMappersMap& scatterMap = _lsFine.getMatrix()._coarseScatterMaps;
    foreach( MultiFieldMatrix::MatrixMappersMap::value_type pos, scatterMap ){
         const StorageSite& site = *pos.first;
         int  neigh_proc_id  =  site.getGatherProcID();
        //if neigh proc belongs to the group, then ok
         if ( _group.count( neigh_proc_id ) > 0 ){
            const Array<int>& fromIndices = *pos.second;
            for ( int i = 0; i < fromIndices.getLength(); i++) 
                 (*scatterCellsLocal)[indx++] = fromIndices[i];
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
   foreach ( set<int>::value_type proc, _group ){
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
    const MultiFieldMatrix::MatrixMappersMap& gatherMap = _lsFine.getMatrix()._coarseGatherMaps;
    foreach( MultiFieldMatrix::MatrixMappersMap::value_type pos, gatherMap ){
         const StorageSite& site = *pos.first;
         int  neigh_proc_id  =  site.getGatherProcID();
        //if neigh proc belongs to the group, then ok
         if ( _group.count( neigh_proc_id ) > 0 ){
            const Array<int>& toIndices = *pos.second;
            for ( int i = 0; i < toIndices.getLength(); i++ ) 
                 (*gatherCellsLocal)[indx++] = toIndices[i];
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
   foreach ( set<int>::value_type proc, _group ){
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


   const MultiFieldMatrix::CoarseConnectivitiesMap&  coarseConnectivities = _lsFine.getMatrix()._coarseConnectivities;
   //first fill _rowLength and _colLength
   foreach( MultiFieldMatrix::CoarseConnectivitiesMap::value_type pos, coarseConnectivities ){
       //const MultiFieldMatrix::EntryIndex& entryIndx = pos.first;
       const CRConnectivity&  conn = *pos.second;
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
   foreach( MultiFieldMatrix::CoarseConnectivitiesMap::value_type pos, coarseConnectivities ){
       //const MultiFieldMatrix::EntryIndex& entryIndx = pos.first;
       const Array<int>&  row_local = pos.second->getRow();
       const Array<int>&  col_local = pos.second->getCol();
       int rowLength = row_local.getLength();
       int colLength = col_local.getLength();
      _comm.Gatherv( row_local.getData(), rowLength, MPI::INT, row->getData(), &(*_rowLength)[0], displs_row, MPI::INT, _targetID );
      _comm.Gatherv( col_local.getData(), colLength, MPI::INT, col->getData(), &(*_colLength)[0], displs_col, MPI::INT, _targetID );
   }

   //copying to map data structures, procid->row, or procid->col
   int local = 0;
   foreach ( set<int>::value_type proc_id, _group ){
     _row[proc_id] = ArrayIntPtr ( new Array<int> ( (*_rowLength)[local] ) );
     _col[proc_id] = ArrayIntPtr ( new Array<int> ( (*_colLength)[local] ) );
      local++;
   }

   int indx = 0;
   local = 0;
   foreach ( set<int>::value_type proc_id, _group ){
      for ( int i = 0; i < (*_rowLength)[local]; i++ )
         (*_row[proc_id])[i] = (*row)[indx++];
      local++;
   }

   indx  = 0;
   local = 0;
   foreach ( set<int>::value_type proc_id, _group ){
      for ( int i = 0; i < (*_colLength)[local]; i++ )
         (*_col[proc_id])[i] = (*col)[indx++];
      local++;
   }

}


void 
LinearSystemMerger::get_local_to_global_map()
{
    //localToGlobal: assuming numbering start with lowest proc ID to the highest one
    //selfCounts is array from 0, 1, ... NPROCS-1 but groupID is not necessaryl holding that
    //future we might group processes so it will not be always from 0, 1, .. NPROCS-1
    int indx = 0;
    foreach ( set<int>::value_type proc_id, _group){
       _localToGlobal[proc_id] = ArrayIntPtr( new Array<int> ( (*_selfCounts)[indx] ) );
       _totalCells +=  (*_selfCounts)[indx];
        indx++;
    }
    cout << " _totalCells = " << _totalCells << " procid = " << _procID << endl;

   _globalToProc  = ArrayIntPtr( new Array<int> ( _totalCells ) );
   _globalToLocal = ArrayIntPtr( new Array<int> ( _totalCells ) );


    indx = 0;
    int global_id = 0;
    foreach ( set<int>::value_type proc_id, _group ){
        for ( int n = 0; n < (*_selfCounts)[indx]; n++ ){
            cout << "proc_id = " << proc_id << " n = " << n << endl;
            (*_localToGlobal[proc_id])[n] = global_id;
            (*_globalToProc )[global_id] = proc_id;
            (*_globalToLocal)[global_id] = n;
            global_id++;
        }
        indx++;
    }

    MPI::COMM_WORLD.Barrier();
    cout << " _totalScatterCells = " << _totalScatterCells << " procid = " << _procID << endl;
    MPI::COMM_WORLD.Barrier();

   debug_print();
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

      //gather interface ids
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
 



     debug_file.close();
   }
   

}


void  
LinearSystemMerger::set_crconnectivity()
{

//     const MultiFieldMatrix::StorageSiteMap&  coarseSites = _lsFine.getMatrix()._coarseSites;
//     foreach( MultiFieldMatrix::StorageSiteMap::value_type pos, coarseSites ){
//        const StorageSite& cell_site = *pos.second;
//        _site = shared_ptr< StorageSite > ( new StorageSiteMerger( _targetID, _group, cell_site) );
//     }


   // _conn = share_ptr<CRConnectivity> ( new CRConnecit

}
