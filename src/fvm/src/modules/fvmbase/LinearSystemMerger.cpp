
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "LinearSystemMerger.h"
#include "CRConnectivity.h"

using namespace std;

LinearSystemMerger::LinearSystemMerger(int target_proc_id, const set<int>& group, LinearSystem& ls )
:_targetID(target_proc_id), _groupID(target_proc_id), _group(group), _lsCoarse(ls)
{

 // cout << " matrix size = " << _lsCoarse.getMatrix().getSize()  << endl;


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
   _neighMeshCounts = ArrayIntPtr( new Array<int>(_totalProcs) );
   _scatterSize     = ArrayIntPtr( new Array<int>(_totalProcs) );
   _gatherSize      = ArrayIntPtr( new Array<int>(_totalProcs) );
   _rowLength       = ArrayIntPtr( new Array<int>(_totalProcs) ); //CRConnecivity row_dimension,  size = _totalProcs
   _colLength       = ArrayIntPtr( new Array<int>(_totalProcs) );
   //get number of all neightbourhood of a local mesh
   get_neigh_mesh_counts();
   get_scatter_cells();
   get_scatter_cells();
   get_crconnectivity();
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

    cout << " matrix size ss= " << _lsCoarse.getMatrix().getSize()  << endl;
    cout << " get_neigh in " <<   " procID = " << _procID << endl;
    vector<int> siteScatterOffset;
    vector<int> siteGatherOffset;
    vector<int> scatterMapIDsLocal;
    vector<int> gatherMapIDsLocal;

    int nmesh = 0;
    const MultiFieldMatrix::MatrixMappersMap& scatterMap = _lsCoarse.getMatrix()._coarseScatterMaps;
    foreach( MultiFieldMatrix::MatrixMappersMap::value_type pos, scatterMap ){
         const StorageSite& site = *pos.first;
         int  neigh_proc_id  =  site.getScatterProcID();
        //if neigh proc belongs to the group, then ok
         if ( _group.count( neigh_proc_id ) > 0 ){
            const Array<int>& fromIndices = *pos.second;
            _totalScatterCellsLocal += fromIndices.getLength();
            siteScatterOffset.push_back( fromIndices.getLength() );
            scatterMapIDsLocal.push_back( neigh_proc_id );
            nmesh++;
         }
    }

   _comm.Gather( &nmesh, 1, MPI::INT,  _neighMeshCounts->getData(), 1, MPI::INT, _targetID );
   _comm.Gather( &_totalScatterCellsLocal, 1, MPI::INT, _scatterSize->getData(), 1, MPI::INT, _targetID );


    const MultiFieldMatrix::MatrixMappersMap& gatherMap = _lsCoarse.getMatrix()._coarseGatherMaps;
    foreach( MultiFieldMatrix::MatrixMappersMap::value_type pos, gatherMap ){
         const StorageSite& site = *pos.first;
         int  neigh_proc_id  =  site.getGatherProcID();
        //if neigh proc belongs to the group, then ok
         if ( _group.count( neigh_proc_id ) > 0 ){
            const Array<int>& toIndices = *pos.second;
            _totalGatherCellsLocal += toIndices.getLength();
            siteGatherOffset.push_back( toIndices.getLength() );
            gatherMapIDsLocal.push_back( neigh_proc_id );
         }
    }
   _comm.Gather( &_totalGatherCellsLocal, 1, MPI::INT, _gatherSize->getData(), 1, MPI::INT, _targetID );


     //total interfaces
     for ( int i = 0; i < _neighMeshCounts->getLength(); i++ )
          _totalInterfaces += (*_neighMeshCounts)[i];

     _scatterInterfaceCounts = ArrayIntPtr( new Array<int>( _totalInterfaces) );
     _gatherInterfaceCounts  = ArrayIntPtr( new Array<int>( _totalInterfaces) );

    int displs[ _totalInterfaces ];
    displs[0] = 0;
    for ( int i = 1; i < _totalInterfaces; i++ )
         displs[i] = displs[i-1] + (*_neighMeshCounts)[i-1];

    //recvcnts = _neighMeshCounts()->getData(),
    _comm.Gatherv( &siteScatterOffset[0], nmesh, MPI::INT, _scatterInterfaceCounts->getData(), &(*_neighMeshCounts)[0], displs, MPI::INT, _targetID );

    //recvcnts = _neighMeshCounts()->getData(),
    _comm.Gatherv( &siteGatherOffset[0],  nmesh, MPI::INT, _gatherInterfaceCounts->getData(), &(*_neighMeshCounts)[0], displs, MPI::INT, _targetID );


    _scatterMapIDs = ArrayIntPtr( new Array<int>( _totalInterfaces) );
    //scatter neightbour IDS
    _comm.Gatherv( &scatterMapIDsLocal[0], nmesh,  MPI::INT, _scatterMapIDs->getData(), &(*_neighMeshCounts)[0], displs, MPI::INT, _targetID );

    _gatherMapIDs = ArrayIntPtr( new Array<int>( _totalInterfaces) );
    //scatter neightbour IDS
    _comm.Gatherv( &gatherMapIDsLocal[0], nmesh,  MPI::INT, _gatherMapIDs->getData(), &(*_neighMeshCounts)[0], displs, MPI::INT, _targetID );





}

//merging scatterArray
void
LinearSystemMerger::get_scatter_cells()
{ 

     for ( int i = 0; i < _scatterInterfaceCounts->getLength(); i++ )
        _totalScatterCells += (*_scatterInterfaceCounts)[i];

   //fill scatter cells in local for shipping
    ArrayIntPtr  scatterCellsLocal( new Array<int> ( _totalScatterCells ) );
    int indx = 0;
    const MultiFieldMatrix::MatrixMappersMap& scatterMap = _lsCoarse.getMatrix()._coarseScatterMaps;
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
    _scatterCells = ArrayIntPtr( new Array<int> ( _totalScatterCells) );


    int displs[ _totalProcs ];
    displs[0] = 0;
    for ( int i = 1; i < _totalProcs; i++ )
         displs[i] = displs[i-1] + (*_scatterSize)[i-1];

     //scatterSize is receive count array, 
    _comm.Gatherv( scatterCellsLocal->getData(), _totalScatterCellsLocal, MPI::INT,
                  _scatterCells->getData(), &(*_scatterSize)[0], displs, MPI::INT, _targetID );

 
}


//merging scatterArray
void
LinearSystemMerger::get_gather_cells()
{ 

     for ( int i = 0; i < _gatherInterfaceCounts->getLength(); i++ )
        _totalGatherCells += (*_gatherInterfaceCounts)[i];

   //fill gather cells in local for shipping
    ArrayIntPtr  gatherCellsLocal( new Array<int> ( _totalGatherCells ) );
    int indx = 0;
    const MultiFieldMatrix::MatrixMappersMap& gatherMap = _lsCoarse.getMatrix()._coarseGatherMaps;
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
    _gatherCells = ArrayIntPtr( new Array<int> ( _totalGatherCells) );

    int displs[ _totalProcs ];
    displs[0] = 0;
    for ( int i = 1; i < _totalProcs; i++ )
         displs[i] = displs[i-1] + (*_gatherSize)[i-1];

     //scatterSize is receive count array, 
    _comm.Gatherv( gatherCellsLocal->getData(), _totalGatherCellsLocal, MPI::INT,
                  _gatherCells->getData(), &(*_gatherSize)[0], displs, MPI::INT, _targetID );

 
}


void
LinearSystemMerger::get_crconnectivity()
{


   const MultiFieldMatrix::CoarseConnectivitiesMap&  coarseConnectivities = _lsCoarse.getMatrix()._coarseConnectivities;
   //first fill _rowLength and _colLength
   foreach( MultiFieldMatrix::CoarseConnectivitiesMap::value_type pos, coarseConnectivities ){
       //const MultiFieldMatrix::EntryIndex& entryIndx = pos.first;
       const CRConnectivity&  conn = *pos.second;
       int rowLength = conn.getRow().getLength();
       int colLength = conn.getCol().getLength();
       _comm.Gather( &rowLength, 1, MPI::INT, _rowLength->getData(), 1, MPI::INT, _targetID );
       _comm.Gather( &colLength, 1, MPI::INT, _colLength->getData(), 1, MPI::INT, _targetID );
   }

  //compute aggloremation size
   int row_size = 0;
   int col_size = 0;
   for ( int i = 0; i < _totalProcs; i++ ){
       row_size  += (*_rowLength)[i];
       col_size  += (*_colLength)[i];
   }

   //forming arrays
   _row = ArrayIntPtr( new Array<int> ( row_size ) );
   _col = ArrayIntPtr( new Array<int> ( col_size ) );

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
       const Array<int>&  row = pos.second->getRow();
       const Array<int>&  col = pos.second->getCol();
       int rowLength = row.getLength();
       int colLength = col.getLength();
      _comm.Gatherv( row.getData(), rowLength, MPI::INT, _row->getData(), &(*_rowLength)[0], displs_row, MPI::INT, _targetID );
      _comm.Gatherv( col.getData(), colLength, MPI::INT, _col->getData(), &(*_colLength)[0], displs_col, MPI::INT, _targetID );
   }






}


void
LinearSystemMerger::debug_print()
{

   stringstream ss;
   ss << "proc" << MPI::COMM_WORLD.Get_rank() << "_linear_system_merger.dat";
   ofstream  debug_file( (ss.str()).c_str() );
  
   debug_file << " neigh Mesh Counts   = " << (*_neighMeshCounts)[_procID] << endl;
   //debug_file << " GhostCount  = " << _mergeSiteGhostSize << endl;
   //debug_file << " count       = " << _mergeSiteSize + _mergeSiteGhostSize << endl;

  debug_file.close();

}

