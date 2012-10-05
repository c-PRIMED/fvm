// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include "CRConnectivity.h"
#include "SpikeStorage.h"
#include "Field.h"
#include <iostream>

using namespace std;


SpikeStorage::SpikeStorage(const CRConnectivity& conn, int semi_bandwidth):
_conn(conn), _bandwidth(semi_bandwidth) 
{
	//logCtor();
  
  init();
}

SpikeStorage::~SpikeStorage()
{
	//logDtor();
}

                    //PRIVATE
void
SpikeStorage::init()
{
   //reserve minumum space for data structures
   _LSPK_INTERIOR.reserve(_bandwidth);
   _LSPK_GHOST.reserve   (_bandwidth);
   _LSPK_OFFD_PTR.reserve(_bandwidth);
   _RSPK_INTERIOR.reserve(_bandwidth);
   _RSPK_GHOST.reserve   (_bandwidth);
   _RSPK_OFFD_PTR.reserve(_bandwidth);
   _LSPK_I.reserve(_bandwidth);
   _LSPK_J.reserve(_bandwidth);
   _RSPK_I.reserve(_bandwidth);
   _RSPK_J.reserve(_bandwidth);
   //keeping left and right spike matrices counter, we will know if there is more than one ghost cell
    int ncells = _conn.getRowSite().getSelfCount();
   _LSPKCountGhost.resize(ncells,0);
   _RSPKCountGhost.resize(ncells,0);


#ifdef FVM_PARALLEL
   _cellSelfCounts.resize( MPI::COMM_WORLD.Get_size() );
   _procID = MPI::COMM_WORLD.Get_rank();
#endif
   //get size of cells from other process
   gatherCellSizes();
   //syn operation to get neighbour local indices to this ghost cells
   syncCellIDs();
   //create local array to store global indices 
   setGlobalIndices();
   //spk_offd_ptr
   setOffDiagPtr();
}

//gather cells size
void
SpikeStorage::gatherCellSizes()
{
   const StorageSite& cellSite = _conn.getRowSite();
   _localCellSelfCount = cellSite.getSelfCount(); 
#ifdef FVM_PARALLEL
   MPI::COMM_WORLD.Allgather(&_localCellSelfCount, 1, MPI::INT, &_cellSelfCounts[0], 1, MPI::INT);
#endif
}
//sync operation to get cellIDs 
void 
SpikeStorage::syncCellIDs()
{
   shared_ptr<Field> cellIndicesField( new Field("cellID") );

   const StorageSite& cells = _conn.getRowSite();	
   const int cellCount = cells.getCount();
   shared_ptr< Array<int> > indPtr( new Array<int>(cellCount) );
   Array<int>& indices = *indPtr;
   Array<int>  indicesOld( cellCount ); //keep original indices
   //zeroing
   indices.zero(); 
   //filling indices
   for ( int n = 0; n < cells.getCount(); n++ ){
       indices[n] = n;
       indicesOld[n] = n;
   }
   //volume add
   cellIndicesField->addArray(cells, indPtr);
   //synLocal to get neighbourhood
   cellIndicesField->syncLocal();
   //create mapping between old and new ghost ids
   int ibeg = cells.getSelfCount();
   int iend = cells.getCount();
   for ( int i = ibeg; i < iend; i++ ){
       _ghostMap[ indicesOld[i] ] = indices[i];
   }


}

//create local array to store global indices
void
SpikeStorage::setGlobalIndices()
{
  int indx_base = 0; 
  for ( int i = 1; i <= _procID; i++)
     indx_base += _cellSelfCounts[i-1];
  
  const StorageSite& cellSite = _conn.getRowSite();
  const StorageSite::GatherMap& gatherMap = cellSite.getGatherMap();
  //loop over gather maps
  foreach( const StorageSite::GatherMap::value_type& mpos, gatherMap ){
     const StorageSite& oSite = *mpos.first;
     const int oRank = oSite.getGatherProcID();
     const Array<int>& ghostIndices = *(mpos.second);
     //get inner indices 
     Array<int> innerIndices( ghostIndices.getLength() );
     for ( int n = 0; n < ghostIndices.getLength(); n++ ){
         innerIndices[n] = _conn( ghostIndices[n], 0 ); 
     }
     //finding base for neighbour
     int indx_base_other = 0; 
     for (  int i = 1; i <= oRank; i++)
         indx_base_other += _cellSelfCounts[i-1];

     //check
     for ( int n = 0; n < ghostIndices.getLength(); n++ ){
         const int iGlbIndx = indx_base + innerIndices[n];
	 const int jGlbIndx = indx_base_other + _ghostMap[ ghostIndices[n] ];
         //check left or right
         if ( (iGlbIndx > jGlbIndx) && ((iGlbIndx - jGlbIndx) <= _bandwidth) ){ 
            _LSPK_INTERIOR.push_back( innerIndices[n] );
	    _LSPK_GHOST.push_back   ( ghostIndices[n] );
            _LSPK_I.push_back( innerIndices[n] );
            _LSPK_J.push_back( _bandwidth + jGlbIndx - indx_base );
            _LSPKCountGhost[innerIndices[n]] += 1;
	 } else if ( (iGlbIndx < jGlbIndx) && (jGlbIndx - iGlbIndx) <= _bandwidth ) {
            _RSPK_INTERIOR.push_back( innerIndices[n] );
	    _RSPK_GHOST.push_back   ( ghostIndices[n] );
            _RSPK_I.push_back( innerIndices[n] );
            _RSPK_J.push_back( jGlbIndx - indx_base_other );
            _RSPKCountGhost[innerIndices[n]] += 1;
         }
     }
  }
  
}
//compute RSP_OFFD_PTR and LSPK_OFFD_PTR
void
SpikeStorage::setOffDiagPtr()
{
    //loop over LSPK
    const Array<int>& row = _conn.getRow();
    const Array<int>& col = _conn.getCol();
    for( unsigned int i = 0; i < _LSPK_INTERIOR.size(); i++ ){
        const int rowIndx = _LSPK_INTERIOR[i]; 
	for ( int j = row[rowIndx]; j < row[rowIndx+1]; j++ ){
           const int neighCellID = col[j];
	   if ( neighCellID == _LSPK_GHOST[i] )
              _LSPK_OFFD_PTR.push_back( j ); 
	}
    }
    //loop over RSPK
    for( unsigned int i = 0; i < _RSPK_INTERIOR.size(); i++ ){
        const int rowIndx = _RSPK_INTERIOR[i]; 
	for ( int j = row[rowIndx]; j < row[rowIndx+1]; j++ ){
           const int neighCellID = col[j];
	   if ( neighCellID == _RSPK_GHOST[i] )
              _RSPK_OFFD_PTR.push_back( j ); 
	}
    }


}
